from contextlib import contextmanager
from dataclasses import fields
from dataclasses import FrozenInstanceError
from decimal import Decimal
import inspect
import os
import stat
import subprocess

import pytest
import pysam

import viral_ngs.core.misc
from viral_ngs.classify import lyra


def _segment(query_name, flag=0x4):
    record = pysam.AlignedSegment()
    record.query_name = query_name
    record.flag = flag
    record.query_sequence = "A" * 50
    record.query_qualities = pysam.qualitystring_to_array("I" * 50)
    record.reference_id = -1
    record.reference_start = -1
    return record


def _write_bam(tmp_path, name, query_name="read"):
    path = tmp_path / name
    header = {
        "HD": {"VN": "1.6", "SO": "unsorted"},
        "PG": [{"ID": "source", "PN": "fixture"}],
        "SQ": [],
    }
    with pysam.AlignmentFile(str(path), "wb", header=header) as bam:
        bam.write(_segment(query_name))
    return path, header


def _write_scores(tmp_path, name="scores.tsv", query_name="read"):
    path = tmp_path / name
    path.write_text(
        "read_id\tscore\tcall\n{}\t0.9\t1\n".format(query_name),
        encoding="utf-8",
    )
    return path


def _write_score_rows(tmp_path, name, rows):
    path = tmp_path / name
    path.write_text(
        "read_id\tscore\tcall\n"
        + "".join("{}\t{}\t{}\n".format(*row) for row in rows),
        encoding="utf-8",
    )
    return path


def _write_bam_records(tmp_path, name, records, header=None):
    path = tmp_path / name
    header = header or {
        "HD": {"VN": "1.6", "SO": "unsorted"},
        "PG": [{"ID": "source", "PN": "fixture"}],
        "SQ": [],
    }
    with pysam.AlignmentFile(str(path), "wb", header=header) as bam:
        for record in records:
            bam.write(record)
    return path, header


@contextmanager
def _publication_store(tmp_path, name="publication", rows=None, records=None):
    rows = rows if rows is not None else [
        ("alpha", "0.9", "0"),
        ("zeta", "0.1", "1"),
    ]
    records = records if records is not None else [
        _segment("zeta"),
        _segment("alpha"),
    ]
    score_path = _write_score_rows(tmp_path, name + "-scores.tsv", rows)
    bam_path, header = _write_bam_records(
        tmp_path,
        name + "-source.bam",
        records,
    )
    with lyra.reconcile_lyra_fragments(
        score_path,
        bam_path,
        "sample exact",
        Decimal("0.8000"),
        work_dir=tmp_path,
    ) as store:
        yield store, score_path, bam_path, header


def _staged_plan(tmp_path, score_path, bam_path, normalized_suffix=".tsv"):
    parents = [tmp_path / role for role in ("normalized", "summary", "bam")]
    for parent in parents:
        parent.mkdir(exist_ok=True)
    outputs = (
        parents[0] / ("result" + normalized_suffix),
        parents[1] / "result.tsv",
        parents[2] / "result.bam",
    )
    return lyra._build_artifact_path_plan(score_path, bam_path, *outputs), outputs


def _generate_staged_artifacts(store, path_plan, work_dir):
    transaction = lyra.LyraArtifactTransaction(store, path_plan, work_dir=work_dir)
    transaction._create_stages()
    stages = transaction.stages
    normalized_rows = lyra._write_normalized(store, stages[0])
    bam_records = lyra._write_viral_bam(
        store,
        stages[2],
        work_dir=work_dir,
    )
    lyra._validate_artifact_counts(store.counts, normalized_rows, bam_records)
    lyra._write_summary(store, stages[1], bam_records)
    return stages, normalized_rows, bam_records


def _remove_stages(stages):
    for stage in stages:
        if os.path.lexists(stage.stage_path):
            os.unlink(stage.stage_path)
    for stage in reversed(stages):
        os.close(stage.descriptor)
    closed_parents = set()
    for stage in reversed(stages):
        if stage.parent.descriptor not in closed_parents:
            os.close(stage.parent.descriptor)
            closed_parents.add(stage.parent.descriptor)


def _artifact_paths(tmp_path, prefix):
    return (
        tmp_path / (prefix + "-normalized.tsv"),
        tmp_path / (prefix + "-summary.tsv"),
        tmp_path / (prefix + "-viral.bam"),
    )


def _entry_snapshot(path):
    path = os.fspath(path)
    try:
        entry_stat = os.lstat(path)
    except (FileNotFoundError, NotADirectoryError):
        return ("absent",)
    if stat.S_ISLNK(entry_stat.st_mode):
        return ("symlink", os.readlink(path))
    if stat.S_ISDIR(entry_stat.st_mode):
        return ("directory", entry_stat.st_dev, entry_stat.st_ino)
    with open(path, "rb") as stream:
        contents = stream.read()
    return ("file", entry_stat.st_dev, entry_stat.st_ino, contents)


def _invoke_invalid_postprocess(monkeypatch, paths, expected_category):
    reconcile_calls = []

    @contextmanager
    def fail_if_called(*args, **kwargs):
        reconcile_calls.append((args, kwargs))
        pytest.fail("reconciliation must not run after invalid path preflight")
        yield

    monkeypatch.setattr(lyra, "reconcile_lyra_fragments", fail_if_called)
    snapshots = {os.fspath(path): _entry_snapshot(path) for path in paths}

    with pytest.raises(lyra.LyraPathError) as exc_info:
        lyra.postprocess_lyra(
            paths[0],
            paths[1],
            "sample",
            "0.8",
            paths[2],
            paths[3],
            paths[4],
        )

    assert exc_info.value.category == expected_category
    assert reconcile_calls == []
    assert {
        os.fspath(path): _entry_snapshot(path)
        for path in paths
    } == snapshots
    return exc_info.value


@pytest.mark.parametrize("parent_kind", ["missing", "regular_file"])
def test_preflight_rejects_invalid_output_parent_before_reconciliation(
    tmp_path,
    monkeypatch,
    parent_kind,
):
    score_path = _write_scores(tmp_path, "parent-scores.tsv")
    bam_path, _ = _write_bam(tmp_path, "parent-source.bam")
    parent = tmp_path / "invalid-parent"
    if parent_kind == "regular_file":
        parent.write_bytes(b"caller parent")
    paths = (
        score_path,
        bam_path,
        parent / "normalized.tsv",
        tmp_path / "summary.tsv",
        tmp_path / "viral.bam",
    )

    error = _invoke_invalid_postprocess(monkeypatch, paths, "output_parent")

    assert error.role == "normalized"
    assert error.path == os.fspath(paths[2])
    assert error.conflicting_role is None
    assert error.conflicting_path is None
    if parent_kind == "regular_file":
        assert parent.read_bytes() == b"caller parent"


@pytest.mark.parametrize(
    "alias_kind",
    ["direct", "relative", "dotdot", "symlink_parent"],
)
def test_preflight_rejects_output_alias_spellings_before_reconciliation(
    tmp_path,
    monkeypatch,
    alias_kind,
):
    score_path = _write_scores(tmp_path, "output-alias-scores.tsv")
    bam_path, _ = _write_bam(tmp_path, "output-alias-source.bam")
    parent = tmp_path / "outputs"
    parent.mkdir()
    normalized = parent / "result.tsv"
    if alias_kind == "direct":
        summary = normalized
    elif alias_kind == "relative":
        summary = os.path.relpath(normalized, os.getcwd())
    elif alias_kind == "dotdot":
        nested = parent / "nested"
        nested.mkdir()
        summary = nested / ".." / normalized.name
    else:
        parent_alias = tmp_path / "outputs-link"
        parent_alias.symlink_to(parent, target_is_directory=True)
        summary = parent_alias / normalized.name
    paths = (
        score_path,
        bam_path,
        normalized,
        summary,
        tmp_path / "viral.bam",
    )

    error = _invoke_invalid_postprocess(monkeypatch, paths, "output_output_alias")

    assert error.role == "summary"
    assert error.conflicting_role == "normalized"
    assert error.path == os.fspath(summary)
    assert error.conflicting_path == os.fspath(normalized)


@pytest.mark.parametrize("input_role", ["score", "bam"])
@pytest.mark.parametrize("alias_kind", ["exact", "symlink", "hardlink"])
def test_preflight_rejects_input_output_aliases_before_reconciliation(
    tmp_path,
    monkeypatch,
    input_role,
    alias_kind,
):
    score_path = _write_scores(tmp_path, "input-alias-scores.tsv")
    bam_path, _ = _write_bam(tmp_path, "input-alias-source.bam")
    source = score_path if input_role == "score" else bam_path
    suffix = ".tsv" if input_role == "score" else ".bam"
    if alias_kind == "exact":
        alias_output = source
    elif alias_kind == "symlink":
        alias_output = tmp_path / ("input-alias" + suffix)
        alias_output.symlink_to(source.name)
    else:
        alias_output = tmp_path / ("input-hardlink" + suffix)
        os.link(source, alias_output)
    normalized = (
        alias_output
        if input_role == "score"
        else tmp_path / "input-alias-normalized.tsv"
    )
    viral_bam = (
        alias_output
        if input_role == "bam"
        else tmp_path / "input-alias-viral.bam"
    )
    paths = (
        score_path,
        bam_path,
        normalized,
        tmp_path / "summary.tsv",
        viral_bam,
    )

    error = _invoke_invalid_postprocess(monkeypatch, paths, "input_output_alias")

    assert error.role == ("normalized" if input_role == "score" else "viral_bam")
    assert error.conflicting_role == input_role
    assert error.path == os.fspath(alias_output)
    assert error.conflicting_path == os.fspath(source)


@pytest.mark.parametrize(
    "entry_kind",
    ["regular_file", "directory", "symlink", "dangling_symlink"],
)
def test_preflight_rejects_every_existing_final_entry_before_reconciliation(
    tmp_path,
    monkeypatch,
    entry_kind,
):
    score_path = _write_scores(tmp_path, "existing-scores.tsv")
    bam_path, _ = _write_bam(tmp_path, "existing-source.bam")
    normalized = tmp_path / "existing-normalized.tsv"
    if entry_kind == "regular_file":
        normalized.write_bytes(b"caller artifact")
    elif entry_kind == "directory":
        normalized.mkdir()
    elif entry_kind == "symlink":
        target = tmp_path / "existing-target"
        target.write_bytes(b"caller target")
        normalized.symlink_to(target.name)
    else:
        normalized.symlink_to("missing-target")
    paths = (
        score_path,
        bam_path,
        normalized,
        tmp_path / "summary.tsv",
        tmp_path / "viral.bam",
    )

    error = _invoke_invalid_postprocess(monkeypatch, paths, "output_exists")

    assert error.role == "normalized"
    assert error.path == os.fspath(normalized)


def test_preflight_existing_input_alias_precedes_generic_existing_error(
    tmp_path,
    monkeypatch,
):
    score_path = _write_scores(tmp_path, "precedence-scores.tsv")
    bam_path, _ = _write_bam(tmp_path, "precedence-source.bam")
    normalized = tmp_path / "precedence-normalized.tsv"
    normalized.symlink_to(score_path.name)
    paths = (
        score_path,
        bam_path,
        normalized,
        tmp_path / "summary.tsv",
        tmp_path / "viral.bam",
    )

    error = _invoke_invalid_postprocess(monkeypatch, paths, "input_output_alias")

    assert error.conflicting_role == "score"


def test_preflight_existing_output_alias_precedes_generic_existing_error(
    tmp_path,
    monkeypatch,
):
    score_path = _write_scores(tmp_path, "output-precedence-scores.tsv")
    bam_path, _ = _write_bam(tmp_path, "output-precedence-source.bam")
    normalized = tmp_path / "output-precedence-normalized.tsv"
    summary = tmp_path / "output-precedence-summary.tsv"
    normalized.write_bytes(b"caller output")
    os.link(normalized, summary)
    paths = (
        score_path,
        bam_path,
        normalized,
        summary,
        tmp_path / "output-precedence-viral.bam",
    )

    error = _invoke_invalid_postprocess(monkeypatch, paths, "output_output_alias")

    assert error.role == "summary"
    assert error.conflicting_role == "normalized"


@pytest.mark.parametrize(
    ("normalized_name", "expected_suffix"),
    [
        ("normalized.tsv", ".tsv"),
        ("normalized.tsv.gz", ".tsv.gz"),
        ("normalized.tsv.zst", ".tsv.zst"),
    ],
)
def test_preflight_valid_distinct_parents_produce_immutable_resolved_plan(
    tmp_path,
    normalized_name,
    expected_suffix,
):
    score_path = _write_scores(tmp_path, "valid-scores.tsv")
    bam_path, _ = _write_bam(tmp_path, "valid-source.bam")
    parents = [tmp_path / role for role in ("normalized", "summary", "bam")]
    for parent in parents:
        parent.mkdir()
    outputs = (
        parents[0] / normalized_name,
        parents[1] / "summary.tsv",
        parents[2] / "viral.bam",
    )

    path_plan = lyra._build_artifact_path_plan(
        score_path,
        bam_path,
        *outputs,
    )

    assert [field.name for field in fields(lyra.ArtifactDestination)] == [
        "role",
        "caller_path",
        "final_path",
        "parent_path",
        "destination_key",
        "suffix",
    ]
    assert [field.name for field in fields(lyra.LyraArtifactPathPlan)] == [
        "score_path",
        "bam_path",
        "normalized",
        "summary",
        "viral_bam",
    ]
    assert path_plan.score_path == os.fspath(score_path)
    assert path_plan.bam_path == os.fspath(bam_path)
    assert path_plan.normalized.suffix == expected_suffix
    assert path_plan.summary.suffix == ".tsv"
    assert path_plan.viral_bam.suffix == ".bam"
    for destination, output, parent in zip(
        (path_plan.normalized, path_plan.summary, path_plan.viral_bam),
        outputs,
        parents,
    ):
        parent_stat = os.stat(parent)
        assert destination.caller_path == os.fspath(output)
        assert destination.parent_path == os.path.realpath(parent, strict=True)
        assert destination.final_path == os.path.join(
            destination.parent_path,
            output.name,
        )
        assert destination.destination_key == (
            parent_stat.st_dev,
            parent_stat.st_ino,
            output.name,
        )
        assert not os.path.lexists(destination.final_path)
        with pytest.raises(FrozenInstanceError):
            destination.final_path = "redirected"
    with pytest.raises(FrozenInstanceError):
        path_plan.score_path = "redirected"


def test_coordinator_converts_each_pathlike_once_and_routes_one_plan(
    tmp_path,
    monkeypatch,
):
    class StatefulPath(os.PathLike):
        def __init__(self, valid_path, bypass_path):
            self.valid_path = os.fspath(valid_path)
            self.bypass_path = os.fspath(bypass_path)
            self.fspath_calls = 0

        def __fspath__(self):
            self.fspath_calls += 1
            if self.fspath_calls == 1:
                return self.valid_path
            return self.bypass_path

    score_file = _write_scores(tmp_path, "coordinator-scores.tsv")
    bam_file, _ = _write_bam(tmp_path, "coordinator-source.bam")
    pathlikes = (
        StatefulPath(score_file, tmp_path / "bypass-scores.tsv"),
        StatefulPath(bam_file, tmp_path / "bypass-source.bam"),
        StatefulPath(tmp_path / "normalized.tsv", tmp_path / "bypass-normalized.txt"),
        StatefulPath(tmp_path / "summary.tsv", tmp_path / "bypass-summary.txt"),
        StatefulPath(tmp_path / "viral.bam", tmp_path / "bypass-viral.sam"),
    )
    calls = []
    store = object()
    received = {}
    original_sample_validator = lyra.validate_sample_id
    original_threshold_validator = lyra.validate_lyra_threshold
    original_path_builder = lyra._build_artifact_path_plan

    def validate_sample(value):
        calls.append("sample")
        return original_sample_validator(value)

    def validate_threshold(value):
        calls.append("threshold")
        return original_threshold_validator(value)

    def build_path_plan(*args):
        calls.append("preflight")
        plan = original_path_builder(*args)
        received["built_plan"] = plan
        return plan

    @contextmanager
    def reconcile(*args, **kwargs):
        calls.append("reconcile")
        received["reconcile"] = (args, kwargs)
        yield store

    def write(current_store, path_plan, work_dir=None):
        calls.append("writer")
        received["writer"] = (current_store, path_plan, work_dir)

    monkeypatch.setattr(lyra, "validate_sample_id", validate_sample)
    monkeypatch.setattr(lyra, "validate_lyra_threshold", validate_threshold)
    monkeypatch.setattr(lyra, "_build_artifact_path_plan", build_path_plan)
    monkeypatch.setattr(lyra, "reconcile_lyra_fragments", reconcile)
    monkeypatch.setattr(lyra, "write_lyra_artifacts", write)

    lyra.postprocess_lyra(
        pathlikes[0],
        pathlikes[1],
        "sample",
        "0.8",
        pathlikes[2],
        pathlikes[3],
        pathlikes[4],
        work_dir=tmp_path,
    )

    path_plan = received["built_plan"]
    assert calls == ["sample", "threshold", "preflight", "reconcile", "writer"]
    assert received["reconcile"] == (
        (
            path_plan.score_path,
            path_plan.bam_path,
            "sample",
            Decimal("0.8"),
        ),
        {"work_dir": tmp_path},
    )
    assert received["writer"] == (store, path_plan, tmp_path)
    assert all(path.fspath_calls == 1 for path in pathlikes)
    assert not any(os.path.lexists(path.bypass_path) for path in pathlikes)


def test_writer_boundary_race_rejects_final_before_any_producer(
    tmp_path,
    monkeypatch,
):
    score_path = _write_scores(tmp_path, "writer-race-scores.tsv")
    bam_path, _ = _write_bam(tmp_path, "writer-race-source.bam")
    outputs = _artifact_paths(tmp_path, "writer-race")
    path_plan = lyra._build_artifact_path_plan(
        score_path,
        bam_path,
        *outputs,
    )
    outputs[0].write_bytes(b"raced caller artifact")
    collaborator_calls = []
    monkeypatch.setattr(
        lyra,
        "_assert_source_bam_identity",
        lambda *args: collaborator_calls.append("source"),
    )
    monkeypatch.setattr(
        lyra,
        "_write_normalized",
        lambda *args: collaborator_calls.append("normalized"),
    )

    with pytest.raises(lyra.LyraPublicationError) as exc_info:
        lyra.write_lyra_artifacts(object(), path_plan, work_dir=tmp_path)

    assert exc_info.value.category == "publication"
    assert exc_info.value.stage == "pre_generation"
    assert isinstance(exc_info.value.__cause__, lyra.LyraPathError)
    assert exc_info.value.__cause__.category == "output_exists"
    assert exc_info.value.__cause__.role == "normalized"
    assert collaborator_calls == []
    assert outputs[0].read_bytes() == b"raced caller artifact"
    assert not outputs[1].exists()
    assert not outputs[2].exists()


def test_link_no_clobber_links_absent_final_and_preserves_raced_final(tmp_path):
    stage = tmp_path / ".lyra-stage.tsv"
    stage.write_bytes(b"staged artifact")
    final = tmp_path / "final.tsv"

    lyra._link_no_clobber(stage, final)

    assert final.read_bytes() == b"staged artifact"
    assert os.stat(stage).st_dev == os.stat(final).st_dev
    assert os.stat(stage).st_ino == os.stat(final).st_ino

    second_stage = tmp_path / ".lyra-second-stage.tsv"
    second_stage.write_bytes(b"second staged artifact")
    raced_final = tmp_path / "raced-final.tsv"
    raced_final.write_bytes(b"racer bytes")
    raced_identity = (os.stat(raced_final).st_dev, os.stat(raced_final).st_ino)

    with pytest.raises(FileExistsError):
        lyra._link_no_clobber(second_stage, raced_final)

    assert raced_final.read_bytes() == b"racer bytes"
    assert (os.stat(raced_final).st_dev, os.stat(raced_final).st_ino) == raced_identity


def _replace_with_valid_bam(tmp_path, source_path):
    replacement, _ = _write_bam(tmp_path, "replacement.bam")
    os.replace(replacement, source_path)


def _append_source(_tmp_path, source_path):
    with open(source_path, "ab") as stream:
        stream.write(b"mutation")


def _truncate_source(_tmp_path, source_path):
    with open(source_path, "r+b") as stream:
        stream.truncate(max(1, source_path.stat().st_size - 8))


def _touch_source(_tmp_path, source_path):
    current = source_path.stat()
    os.utime(
        source_path,
        ns=(current.st_atime_ns, current.st_mtime_ns + 1_000_000_000),
    )


def _finalize_source_store(store, source_path, identity, header):
    store._counts = lyra.LyraReconciliationCounts(0, 0, 0, 0, 0, 0, 0, 0, 0)
    store._threshold = Decimal("0.8")
    source_bam_fd = os.open(source_path, os.O_RDONLY)
    try:
        store._install_source_metadata(
            source_bam_fd=source_bam_fd,
            source_bam_path=str(source_path.absolute()),
            source_bam_display_path=str(source_path),
            source_bam_identity=identity,
            source_bam_header=header,
        )
    except BaseException:
        os.close(source_bam_fd)
        raise


def test_source_file_identity_is_frozen_with_exact_descriptor_fields(tmp_path):
    source_path = tmp_path / "source.bam"
    source_path.write_bytes(b"bam identity")
    stat_result = source_path.stat()

    identity = lyra._file_identity(stat_result)

    assert [field.name for field in fields(lyra.FileIdentity)] == [
        "device",
        "inode",
        "size",
        "mtime_ns",
    ]
    assert identity == lyra.FileIdentity(
        device=stat_result.st_dev,
        inode=stat_result.st_ino,
        size=stat_result.st_size,
        mtime_ns=stat_result.st_mtime_ns,
    )
    assert all(type(getattr(identity, field.name)) is int for field in fields(identity))
    with pytest.raises(FrozenInstanceError):
        identity.size = 0


def test_source_properties_require_finalized_live_store_and_copy_header(tmp_path):
    source_path = tmp_path / "reads.bam"
    source_path.write_bytes(b"bam metadata fixture")
    identity = lyra.FileIdentity(1, 2, 3, 4)
    retained_header = {
        "HD": {"VN": "1.6"},
        "SQ": [{"SN": "ref", "LN": 100}],
    }
    store = lyra.LyraFragmentStore("sample", work_dir=tmp_path)

    for property_name in (
        "source_bam_fd",
        "source_bam_path",
        "source_bam_display_path",
        "source_bam_identity",
        "source_bam_header",
    ):
        with pytest.raises(RuntimeError, match="not finalized"):
            getattr(store, property_name)

    _finalize_source_store(store, source_path, identity, retained_header)

    assert type(store.source_bam_fd) is int
    os.fstat(store.source_bam_fd)
    assert store.source_bam_path == str(source_path.absolute())
    assert store.source_bam_display_path == str(source_path)
    assert store.source_bam_identity is identity
    returned_header = store.source_bam_header
    returned_header["HD"]["VN"] = "mutated"
    returned_header["SQ"].append({"SN": "other", "LN": 1})
    assert store.source_bam_header == retained_header

    with pytest.raises(RuntimeError, match="already installed"):
        _finalize_source_store(store, source_path, identity, retained_header)

    store.close()
    for property_name in (
        "source_bam_fd",
        "source_bam_path",
        "source_bam_display_path",
        "source_bam_identity",
        "source_bam_header",
    ):
        with pytest.raises(RuntimeError, match="closed"):
            getattr(store, property_name)


def test_source_descriptor_matches_identity_and_closes_once_with_store(
    tmp_path,
    monkeypatch,
):
    source_path, _ = _write_bam(tmp_path, "descriptor-lifecycle.bam")
    score_path = _write_scores(tmp_path, "descriptor-lifecycle.tsv")
    close_calls = []
    real_close = os.close

    with lyra.reconcile_lyra_fragments(
        score_path,
        source_path,
        "sample",
        "0.8",
        work_dir=tmp_path,
    ) as store:
        source_bam_fd = store.source_bam_fd
        assert lyra._file_identity(os.fstat(source_bam_fd)) == (
            store.source_bam_identity
        )

        def record_close(descriptor):
            if descriptor == source_bam_fd:
                close_calls.append(descriptor)
            return real_close(descriptor)

        monkeypatch.setattr(lyra.os, "close", record_close)

    assert close_calls == [source_bam_fd]
    with pytest.raises(OSError):
        os.fstat(source_bam_fd)


@pytest.mark.parametrize(
    "failing_helper",
    ["_validate_reconciliation", "_finalize_fragments"],
)
def test_source_descriptor_created_before_failure_is_closed(
    tmp_path,
    monkeypatch,
    failing_helper,
):
    source_path, _ = _write_bam(
        tmp_path,
        "descriptor-failure-{}.bam".format(failing_helper),
    )
    score_path = _write_scores(
        tmp_path,
        "descriptor-failure-{}.tsv".format(failing_helper),
    )
    duplicated_descriptors = []
    real_dup = os.dup

    def record_duplicate(descriptor):
        duplicate = real_dup(descriptor)
        duplicated_descriptors.append(duplicate)
        return duplicate

    primary = RuntimeError("injected {} failure".format(failing_helper))
    monkeypatch.setattr(lyra.os, "dup", record_duplicate)
    monkeypatch.setattr(lyra, failing_helper, lambda *args: _raise(primary))

    with pytest.raises(RuntimeError) as exc_info:
        with lyra.reconcile_lyra_fragments(
            score_path,
            source_path,
            "sample",
            "0.8",
            work_dir=tmp_path,
        ):
            pytest.fail("failed reconciliation must not expose a store")

    assert exc_info.value is primary
    assert len(duplicated_descriptors) == 1
    with pytest.raises(OSError):
        os.fstat(duplicated_descriptors[0])


def test_source_descriptor_close_failure_still_cleans_store_once(
    tmp_path,
    monkeypatch,
):
    source_path, header = _write_bam(tmp_path, "descriptor-close-failure.bam")
    identity = lyra._file_identity(os.stat(source_path))
    store = lyra.LyraFragmentStore("sample", work_dir=tmp_path)
    _finalize_source_store(store, source_path, identity, header)
    source_bam_fd = store.source_bam_fd
    real_descriptor_close = os.close
    real_connection = store._connection
    real_temporary_cleanup = store._temporary_directory.cleanup
    operations = []

    class RecordingConnection:
        def close(self):
            operations.append("database")
            real_connection.close()

    def fail_source_close(descriptor):
        if descriptor == source_bam_fd:
            operations.append("source")
            raise OSError(5, "injected source descriptor close failure")
        return real_descriptor_close(descriptor)

    def cleanup_temporary_directory():
        operations.append("temporary_directory")
        real_temporary_cleanup()

    store._connection = RecordingConnection()
    store._temporary_directory.cleanup = cleanup_temporary_directory
    with monkeypatch.context() as patch:
        patch.setattr(lyra.os, "close", fail_source_close)
        with pytest.raises(OSError, match="source descriptor close failure"):
            store.close()
        store.close()

    assert operations == ["source", "database", "temporary_directory"]
    real_descriptor_close(source_bam_fd)


def test_source_identity_rejects_retargeted_retained_descriptor(tmp_path):
    source_path, _ = _write_bam(tmp_path, "descriptor-identity-source.bam")
    replacement_path, _ = _write_bam(
        tmp_path,
        "descriptor-identity-replacement.bam",
        query_name="replacement",
    )
    score_path = _write_scores(tmp_path, "descriptor-identity.tsv")

    with lyra.reconcile_lyra_fragments(
        score_path,
        source_path,
        "sample",
        "0.8",
        work_dir=tmp_path,
    ) as store:
        replacement_fd = os.open(replacement_path, os.O_RDONLY)
        try:
            os.dup2(replacement_fd, store.source_bam_fd)
        finally:
            os.close(replacement_fd)

        with pytest.raises(lyra.LyraSourceIdentityError) as exc_info:
            lyra._assert_source_bam_identity(store, "descriptor_check")

    assert exc_info.value.stage == "descriptor_check"
    assert exc_info.value.path == os.fspath(source_path)
    assert exc_info.value.expected == lyra._file_identity(os.stat(source_path))
    assert exc_info.value.actual == lyra._file_identity(os.stat(replacement_path))


def test_source_transient_path_swap_cannot_change_filtered_bam(tmp_path, monkeypatch):
    header = {
        "HD": {"VN": "1.6", "SO": "unsorted"},
        "PG": [{"ID": "source", "PN": "fixture"}],
        "SQ": [],
    }
    reconciled_record = _segment("read")
    reconciled_record.set_tag("NM", 1)
    replacement_record = _segment("read")
    replacement_record.set_tag("NM", 99)
    source_path, _ = _write_bam_records(
        tmp_path,
        "transient-source.bam",
        [reconciled_record],
        header,
    )
    replacement_path, _ = _write_bam_records(
        tmp_path,
        "transient-replacement.bam",
        [replacement_record],
        header,
    )
    displaced_path = tmp_path / "transient-displaced.bam"
    score_path = _write_scores(tmp_path, "transient-scores.tsv")
    output_path = tmp_path / "transient-output.bam"
    original_filter = lyra.util_misc.ReadIdStore.filter_bam_by_ids
    filter_calls = []

    def swap_during_filter(read_ids, in_bam, out_bam, *args, **kwargs):
        filter_calls.append((in_bam, kwargs.get("in_bam_fd")))
        os.replace(source_path, displaced_path)
        os.replace(replacement_path, source_path)
        try:
            return original_filter(read_ids, in_bam, out_bam, *args, **kwargs)
        finally:
            os.replace(source_path, replacement_path)
            os.replace(displaced_path, source_path)

    monkeypatch.setattr(
        lyra.util_misc.ReadIdStore,
        "filter_bam_by_ids",
        swap_during_filter,
    )

    with lyra.reconcile_lyra_fragments(
        score_path,
        source_path,
        "sample",
        "0.8",
        work_dir=tmp_path,
    ) as store:
        source_bam_fd = store.source_bam_fd
        assert lyra._write_viral_bam(
            store,
            output_path,
            work_dir=tmp_path,
        ) == 1
        assert filter_calls == [
            (store.source_bam_display_path, source_bam_fd)
        ]

    with pysam.AlignmentFile(output_path, "rb", check_sq=False) as output:
        output_records = list(output.fetch(until_eof=True))
    assert len(output_records) == 1
    assert output_records[0].get_tag("NM") == 1


def test_source_identity_error_exposes_stable_bounded_mismatch_facts():
    expected = lyra.FileIdentity(1, 2, 3, 4)
    actual = lyra.FileIdentity(5, 6, 7, 8)
    long_path = "/caller/" + "x" * 300 + ".bam"

    error = lyra.LyraSourceIdentityError(
        stage="pre_generation",
        path=long_path,
        expected=expected,
        actual=actual,
    )

    assert error.category == "source_bam_identity"
    assert error.stage == "pre_generation"
    assert error.path == long_path
    assert error.expected is expected
    assert error.actual is actual
    assert error.actual_status is None
    assert long_path not in str(error)
    assert "x" * 160 not in str(error)
    assert len(str(error)) < 500


def test_source_identity_error_exposes_bounded_missing_status():
    expected = lyra.FileIdentity(1, 2, 3, 4)

    error = lyra.LyraSourceIdentityError(
        stage="reconciliation",
        path="missing.bam",
        expected=expected,
        actual_status="missing",
    )

    assert error.category == "source_bam_identity"
    assert error.stage == "reconciliation"
    assert error.path == "missing.bam"
    assert error.expected is expected
    assert error.actual is None
    assert error.actual_status == "missing"
    assert "actual_status='missing'" in str(error)


def test_source_owned_artifact_apis_accept_no_second_bam_argument():
    assert list(inspect.signature(lyra._write_viral_bam).parameters) == [
        "store",
        "viral_bam_output",
        "work_dir",
    ]
    assert list(inspect.signature(lyra.write_lyra_artifacts).parameters) == [
        "store",
        "path_plan",
        "work_dir",
    ]


@pytest.mark.parametrize("alias_kind", ["direct", "symlink", "hardlink"])
def test_source_direct_symlink_and_hardlink_aliases_produce_artifacts(
    tmp_path,
    alias_kind,
):
    source_path, header = _write_bam(tmp_path, "target-{}.bam".format(alias_kind))
    access_path = source_path
    if alias_kind == "symlink":
        access_path = tmp_path / "source-symlink.bam"
        access_path.symlink_to(source_path.name)
    elif alias_kind == "hardlink":
        access_path = tmp_path / "source-hardlink.bam"
        os.link(source_path, access_path)
    score_path = _write_scores(tmp_path, "scores-{}.tsv".format(alias_kind))
    outputs = _artifact_paths(tmp_path, alias_kind)
    path_plan = lyra._build_artifact_path_plan(
        score_path,
        access_path,
        *outputs,
    )

    with lyra.reconcile_lyra_fragments(
        score_path,
        access_path,
        "sample",
        "0.8",
        work_dir=tmp_path,
    ) as store:
        assert store.source_bam_path == os.path.abspath(os.fspath(access_path))
        assert store.source_bam_display_path == os.fspath(access_path)
        assert store.source_bam_identity == lyra._file_identity(os.stat(access_path))
        assert store.source_bam_header == header
        lyra.write_lyra_artifacts(store, path_plan, work_dir=tmp_path)

    assert all(path.exists() for path in outputs)


@pytest.mark.parametrize(
    "mutation",
    [
        _replace_with_valid_bam,
        _append_source,
        _truncate_source,
        _touch_source,
    ],
    ids=["replacement", "append", "truncate", "mtime"],
)
def test_source_mutation_during_reconciliation_is_rejected(
    tmp_path,
    monkeypatch,
    mutation,
):
    source_path, _ = _write_bam(tmp_path, "mutated-source.bam")
    score_path = _write_scores(tmp_path, "mutated-scores.tsv")
    original_collect = lyra._collect_bam_evidence

    def mutate_after_traversal(store, bam, display_path):
        original_collect(store, bam, display_path)
        mutation(tmp_path, source_path)

    monkeypatch.setattr(lyra, "_collect_bam_evidence", mutate_after_traversal)

    with pytest.raises(lyra.LyraSourceIdentityError) as exc_info:
        with lyra.reconcile_lyra_fragments(
            score_path,
            source_path,
            "sample",
            "0.8",
            work_dir=tmp_path,
        ):
            pytest.fail("a mutated source must not yield a store")

    assert exc_info.value.stage == "reconciliation"
    assert exc_info.value.path == os.fspath(source_path)


def test_source_symlink_retarget_during_reconciliation_is_rejected(
    tmp_path,
    monkeypatch,
):
    first_target, _ = _write_bam(tmp_path, "first-target.bam")
    second_target, _ = _write_bam(tmp_path, "second-target.bam")
    source_link = tmp_path / "retargeted-source.bam"
    source_link.symlink_to(first_target.name)
    score_path = _write_scores(tmp_path, "retargeted-scores.tsv")
    original_collect = lyra._collect_bam_evidence

    def retarget_after_traversal(store, bam, display_path):
        original_collect(store, bam, display_path)
        source_link.unlink()
        source_link.symlink_to(second_target.name)

    monkeypatch.setattr(lyra, "_collect_bam_evidence", retarget_after_traversal)

    with pytest.raises(lyra.LyraSourceIdentityError) as exc_info:
        with lyra.reconcile_lyra_fragments(
            score_path,
            source_link,
            "sample",
            "0.8",
            work_dir=tmp_path,
        ):
            pytest.fail("a retargeted source must not yield a store")

    assert exc_info.value.stage == "reconciliation"
    assert exc_info.value.path == os.fspath(source_link)


@pytest.mark.parametrize("mutation", [_replace_with_valid_bam, _touch_source])
def test_source_mutation_before_generation_opens_no_output_or_producer(
    tmp_path,
    monkeypatch,
    mutation,
):
    source_path, _ = _write_bam(tmp_path, "pregeneration-source.bam")
    score_path = _write_scores(tmp_path, "pregeneration-scores.tsv")
    outputs = _artifact_paths(tmp_path, "pregeneration")
    path_plan = lyra._build_artifact_path_plan(
        score_path,
        source_path,
        *outputs,
    )
    producer_calls = []

    with lyra.reconcile_lyra_fragments(
        score_path,
        source_path,
        "sample",
        "0.8",
        work_dir=tmp_path,
    ) as store:
        mutation(tmp_path, source_path)
        monkeypatch.setattr(
            lyra,
            "_write_normalized",
            lambda *args, **kwargs: producer_calls.append("normalized"),
        )
        monkeypatch.setattr(
            lyra,
            "_write_viral_bam",
            lambda *args, **kwargs: producer_calls.append("bam"),
        )

        with pytest.raises(lyra.LyraPublicationError) as exc_info:
            lyra.write_lyra_artifacts(store, path_plan, work_dir=tmp_path)

    assert exc_info.value.stage == "pre_generation"
    assert isinstance(exc_info.value.__cause__, lyra.LyraSourceIdentityError)
    assert producer_calls == []
    assert not any(path.exists() for path in outputs)


def test_source_missing_before_generation_preserves_filesystem_cause(tmp_path):
    source_path, _ = _write_bam(tmp_path, "missing-source.bam")
    score_path = _write_scores(tmp_path, "missing-scores.tsv")
    outputs = _artifact_paths(tmp_path, "missing")
    path_plan = lyra._build_artifact_path_plan(
        score_path,
        source_path,
        *outputs,
    )

    with lyra.reconcile_lyra_fragments(
        score_path,
        source_path,
        "sample",
        "0.8",
        work_dir=tmp_path,
    ) as store:
        source_path.unlink()
        with pytest.raises(lyra.LyraPublicationError) as exc_info:
            lyra.write_lyra_artifacts(store, path_plan, work_dir=tmp_path)

    assert exc_info.value.stage == "pre_generation"
    source_error = exc_info.value.__cause__
    assert isinstance(source_error, lyra.LyraSourceIdentityError)
    assert source_error.actual is None
    assert source_error.actual_status == "missing"
    assert isinstance(source_error.__cause__, FileNotFoundError)
    assert not any(path.exists() for path in outputs)


def test_source_mutation_after_bam_filter_prevents_counts_and_summary(
    tmp_path,
    monkeypatch,
):
    source_path, _ = _write_bam(tmp_path, "postfilter-source.bam")
    score_path = _write_scores(tmp_path, "postfilter-scores.tsv")
    outputs = _artifact_paths(tmp_path, "postfilter")
    path_plan = lyra._build_artifact_path_plan(
        score_path,
        source_path,
        *outputs,
    )
    original_filter = lyra.util_misc.ReadIdStore.filter_bam_by_ids
    accepted_calls = []

    def mutate_after_filter(read_ids, *args, **kwargs):
        count = original_filter(read_ids, *args, **kwargs)
        _touch_source(tmp_path, source_path)
        return count

    monkeypatch.setattr(
        lyra.util_misc.ReadIdStore,
        "filter_bam_by_ids",
        mutate_after_filter,
    )
    monkeypatch.setattr(
        lyra,
        "_validate_artifact_counts",
        lambda *args: accepted_calls.append("counts"),
    )
    monkeypatch.setattr(
        lyra,
        "_write_summary",
        lambda *args: accepted_calls.append("summary"),
    )

    with lyra.reconcile_lyra_fragments(
        score_path,
        source_path,
        "sample",
        "0.8",
        work_dir=tmp_path,
    ) as store:
        with pytest.raises(lyra.LyraPublicationError) as exc_info:
            lyra.write_lyra_artifacts(store, path_plan, work_dir=tmp_path)

    assert exc_info.value.stage == "after_bam_filter"
    assert isinstance(exc_info.value.__cause__, lyra.LyraSourceIdentityError)
    assert accepted_calls == []
    assert not outputs[1].exists()


@pytest.mark.parametrize(
    "normalized_suffix",
    [".tsv", ".tsv.gz", ".tsv.zst"],
)
def test_staged_paths_are_hidden_same_parent_and_preserve_exact_suffix(
    tmp_path,
    normalized_suffix,
):
    score_path = _write_scores(tmp_path, "stage-scores.tsv")
    bam_path, _ = _write_bam(tmp_path, "stage-source.bam")
    path_plan, _ = _staged_plan(
        tmp_path,
        score_path,
        bam_path,
        normalized_suffix,
    )
    destinations = (
        path_plan.normalized,
        path_plan.summary,
        path_plan.viral_bam,
    )
    transaction = lyra.LyraArtifactTransaction(object(), path_plan)
    transaction._create_stages()
    stages = transaction.stages
    try:
        assert [field.name for field in fields(lyra.ArtifactStage)] == [
            "role",
            "basename",
            "display_path",
            "descriptor",
            "object_identity",
            "destination",
            "parent",
        ]
        for stage, destination in zip(stages, destinations):
            stage_parent = os.stat(os.path.dirname(stage.stage_path))
            final_parent = os.stat(destination.parent_path)
            assert (stage_parent.st_dev, stage_parent.st_ino) == (
                final_parent.st_dev,
                final_parent.st_ino,
            )
            assert os.path.basename(stage.stage_path).startswith(
                ".lyra-{}-".format(destination.role)
            )
            assert stage.stage_path.endswith(destination.suffix)
            assert stage.destination is destination
    finally:
        _remove_stages(stages)


def test_parent_anchor_and_stage_descriptors_match_exact_contract(
    tmp_path,
    monkeypatch,
):
    score_path = _write_scores(tmp_path, "descriptor-stage-scores.tsv")
    bam_path, _ = _write_bam(tmp_path, "descriptor-stage-source.bam")
    outputs = _artifact_paths(tmp_path, "descriptor-stage")
    path_plan = lyra._build_artifact_path_plan(
        score_path,
        bam_path,
        *outputs,
    )
    opened = []
    real_open = os.open

    def record_open(path, flags, mode=0o777, *, dir_fd=None):
        descriptor = real_open(path, flags, mode, dir_fd=dir_fd)
        opened.append((path, flags, mode, dir_fd, descriptor))
        return descriptor

    monkeypatch.setattr(lyra.os, "open", record_open)
    transaction = lyra.LyraArtifactTransaction(object(), path_plan)
    transaction._create_stages()
    stages = transaction.stages

    try:
        assert [field.name for field in fields(lyra.ParentDirectoryAnchor)] == [
            "parent_path",
            "identity",
            "descriptor",
        ]
        assert [field.name for field in fields(lyra.ArtifactStage)] == [
            "role",
            "basename",
            "display_path",
            "descriptor",
            "object_identity",
            "destination",
            "parent",
        ]
        assert len(stages) == 3
        assert all(stage.parent is stages[0].parent for stage in stages)

        parent = stages[0].parent
        parent_stat = os.stat(path_plan.normalized.parent_path)
        assert parent.parent_path == path_plan.normalized.parent_path
        assert parent.identity == (parent_stat.st_dev, parent_stat.st_ino)
        retained_parent_stat = os.fstat(parent.descriptor)
        assert parent.identity == (
            retained_parent_stat.st_dev,
            retained_parent_stat.st_ino,
        )

        parent_opens = [
            item
            for item in opened
            if item[0] == parent.parent_path and item[3] is None
        ]
        assert len(parent_opens) == 1
        parent_flags = parent_opens[0][1]
        assert parent_flags & os.O_DIRECTORY
        assert parent_flags & os.O_NOFOLLOW
        assert parent_flags & os.O_CLOEXEC

        stage_opens = [item for item in opened if item[3] == parent.descriptor]
        assert len(stage_opens) == 3
        for stage, opened_stage in zip(stages, stage_opens):
            path, flags, mode, dir_fd, descriptor = opened_stage
            stage_stat = os.fstat(stage.descriptor)
            entry_stat = os.stat(
                stage.basename,
                dir_fd=parent.descriptor,
                follow_symlinks=False,
            )
            assert path == stage.basename
            assert flags & os.O_CREAT
            assert flags & os.O_EXCL
            assert flags & os.O_NOFOLLOW
            assert flags & os.O_CLOEXEC
            assert mode == 0o600
            assert dir_fd == parent.descriptor
            assert descriptor == stage.descriptor
            assert stage.object_identity == (stage_stat.st_dev, stage_stat.st_ino)
            assert stage.object_identity == (entry_stat.st_dev, entry_stat.st_ino)
            assert stat.S_IMODE(stage_stat.st_mode) == 0o600
            assert stage.basename.startswith(".lyra-{}-".format(stage.role))
            assert stage.basename.endswith(stage.destination.suffix)
            assert stage.display_path == os.path.join(
                parent.parent_path,
                stage.basename,
            )
    finally:
        transaction.rollback_and_cleanup()

    for _, _, _, _, descriptor in opened:
        with pytest.raises(OSError):
            os.fstat(descriptor)
    assert not list(tmp_path.rglob(".lyra-*"))


def test_distinct_parent_anchors_match_each_planned_parent(tmp_path):
    score_path = _write_scores(tmp_path, "distinct-parent-scores.tsv")
    bam_path, _ = _write_bam(tmp_path, "distinct-parent-source.bam")
    path_plan, _ = _staged_plan(tmp_path, score_path, bam_path)
    transaction = lyra.LyraArtifactTransaction(object(), path_plan)
    transaction._create_stages()

    try:
        anchors = [stage.parent for stage in transaction.stages]
        assert len({anchor.descriptor for anchor in anchors}) == 3
        assert len({anchor.identity for anchor in anchors}) == 3
        for stage in transaction.stages:
            expected = stage.destination.destination_key[:2]
            retained = os.fstat(stage.parent.descriptor)
            assert stage.parent.identity == expected
            assert expected == (retained.st_dev, retained.st_ino)
    finally:
        transaction.rollback_and_cleanup()


def test_parent_replacement_is_rejected_while_descriptor_anchors_original(
    tmp_path,
):
    score_path = _write_scores(tmp_path, "replaced-parent-scores.tsv")
    bam_path, _ = _write_bam(tmp_path, "replaced-parent-source.bam")
    path_plan, _ = _staged_plan(tmp_path, score_path, bam_path)
    transaction = lyra.LyraArtifactTransaction(object(), path_plan)
    transaction._create_stages()
    normalized_stage = next(
        stage for stage in transaction.stages if stage.role == "normalized"
    )
    parent = normalized_stage.parent
    moved_parent = tmp_path / "moved-normalized-parent"
    os.rename(parent.parent_path, moved_parent)
    os.mkdir(parent.parent_path)

    try:
        with pytest.raises(lyra.LyraPathError) as exc_info:
            lyra._assert_parent_anchor(parent, "parent_recheck")

        error = exc_info.value
        assert error.category == "output_parent_identity"
        assert error.stage == "parent_recheck"
        assert error.path == parent.parent_path
        assert normalized_stage.basename in os.listdir(parent.descriptor)
        assert not os.path.exists(
            os.path.join(parent.parent_path, normalized_stage.basename)
        )
        assert len(str(error)) < 700
    finally:
        transaction.rollback_and_cleanup()


def test_stage_symlink_substitution_never_follows_or_changes_target(tmp_path):
    score_path = _write_scores(tmp_path, "stage-symlink-scores.tsv")
    bam_path, _ = _write_bam(tmp_path, "stage-symlink-source.bam")
    path_plan, _ = _staged_plan(tmp_path, score_path, bam_path)
    transaction = lyra.LyraArtifactTransaction(object(), path_plan)
    transaction._create_stages()
    stage = transaction.stages[0]
    victim = tmp_path / "stage-symlink-victim"
    victim.write_bytes(b"caller bytes")
    victim_identity = lyra._file_identity(os.stat(victim))
    os.unlink(stage.basename, dir_fd=stage.parent.descriptor)
    os.symlink(
        victim,
        stage.basename,
        dir_fd=stage.parent.descriptor,
    )

    try:
        with pytest.raises(lyra.LyraArtifactConsistencyError) as exc_info:
            lyra._assert_stage_handle(stage)

        error = exc_info.value
        assert error.category == "stage_identity"
        assert error.field == stage.role
        assert error.expected == stage.object_identity
        assert victim.read_bytes() == b"caller bytes"
        assert lyra._file_identity(os.stat(victim)) == victim_identity
        assert stat.S_ISLNK(
            os.stat(
                stage.basename,
                dir_fd=stage.parent.descriptor,
                follow_symlinks=False,
            ).st_mode
        )
    finally:
        os.unlink(stage.basename, dir_fd=stage.parent.descriptor)
        transaction.rollback_and_cleanup()


def test_partial_stage_acquisition_closes_every_owned_descriptor_once(
    tmp_path,
    monkeypatch,
):
    score_path = _write_scores(tmp_path, "partial-stage-scores.tsv")
    bam_path, _ = _write_bam(tmp_path, "partial-stage-source.bam")
    outputs = _artifact_paths(tmp_path, "partial-stage")
    path_plan = lyra._build_artifact_path_plan(
        score_path,
        bam_path,
        *outputs,
    )
    real_open = os.open
    real_close = os.close
    owned_descriptors = []
    close_calls = []
    stage_open_count = 0

    def fail_second_stage_open(path, flags, mode=0o777, *, dir_fd=None):
        nonlocal stage_open_count
        if dir_fd is not None and flags & os.O_CREAT:
            stage_open_count += 1
            if stage_open_count == 2:
                raise OSError(5, "injected stage acquisition failure")
        descriptor = real_open(path, flags, mode, dir_fd=dir_fd)
        if flags & os.O_DIRECTORY or dir_fd is not None:
            owned_descriptors.append(descriptor)
        return descriptor

    def record_close(descriptor):
        if descriptor in owned_descriptors:
            close_calls.append(descriptor)
        return real_close(descriptor)

    monkeypatch.setattr(lyra.os, "open", fail_second_stage_open)
    monkeypatch.setattr(lyra.os, "close", record_close)
    transaction = lyra.LyraArtifactTransaction(object(), path_plan)

    with pytest.raises(OSError, match="stage acquisition failure"):
        transaction._create_stages()

    assert stage_open_count == 2
    assert close_calls == list(reversed(owned_descriptors))
    assert len(close_calls) == len(set(close_calls))
    for descriptor in owned_descriptors:
        with pytest.raises(OSError):
            os.fstat(descriptor)
    assert not list(tmp_path.rglob(".lyra-*"))


@pytest.mark.parametrize(
    (
        "case",
        "normalized_suffix",
        "rows",
        "records",
        "expected_rows",
        "expected_bam",
    ),
    [
        ("empty", ".tsv", [], [], 0, 0),
        (
            "no-hit",
            ".tsv.gz",
            [("alpha", "0.1", "1")],
            [_segment("alpha")],
            1,
            0,
        ),
        ("rich-zstd", ".tsv.zst", None, None, 2, 1),
    ],
)
def test_staged_readback_accepts_empty_no_hit_and_compressed_outputs(
    tmp_path,
    case,
    normalized_suffix,
    rows,
    records,
    expected_rows,
    expected_bam,
):
    with _publication_store(tmp_path, case, rows=rows, records=records) as (
        store,
        score_path,
        bam_path,
        _,
    ):
        path_plan, outputs = _staged_plan(
            tmp_path,
            score_path,
            bam_path,
            normalized_suffix,
        )
        stages, producer_rows, producer_bam = _generate_staged_artifacts(
            store,
            path_plan,
            tmp_path,
        )
        try:
            assert producer_rows == expected_rows
            assert producer_bam == expected_bam
            assert lyra._validate_staged_normalized(
                store,
                stages[0],
            ) == expected_rows
            os.fstat(stages[0].descriptor)
            assert lyra._validate_staged_bam(
                store,
                stages[2],
            ) == expected_bam
            os.fstat(stages[2].descriptor)
            lyra._validate_staged_summary(
                store,
                stages[1],
                expected_rows,
                expected_bam,
            )
            os.fstat(stages[1].descriptor)
            assert not any(path.exists() for path in outputs)
        finally:
            _remove_stages(stages)


@pytest.mark.parametrize(
    ("role", "normalized_suffix"),
    [
        ("normalized", ".tsv"),
        ("normalized", ".tsv.gz"),
        ("normalized", ".tsv.zst"),
        ("summary", ".tsv"),
        ("viral_bam", ".tsv"),
    ],
)
def test_stage_symlink_before_producer_cannot_change_victim(
    tmp_path,
    monkeypatch,
    role,
    normalized_suffix,
):
    with _publication_store(tmp_path, "producer-symlink-" + role) as (
        store,
        score_path,
        bam_path,
        _,
    ):
        path_plan, _ = _staged_plan(
            tmp_path,
            score_path,
            bam_path,
            normalized_suffix,
        )
        transaction = lyra.LyraArtifactTransaction(
            store,
            path_plan,
            work_dir=tmp_path,
        )
        transaction._create_stages()
        stage = next(item for item in transaction.stages if item.role == role)
        victim = tmp_path / ("producer-victim-" + role)
        victim.write_bytes(b"caller bytes")
        victim_identity = lyra._file_identity(os.stat(victim))
        os.unlink(stage.basename, dir_fd=stage.parent.descriptor)
        os.symlink(victim, stage.basename, dir_fd=stage.parent.descriptor)
        filter_descriptors = []
        original_filter = lyra.util_misc.ReadIdStore.filter_bam_by_ids

        def record_filter(read_ids, in_bam, out_bam, *args, **kwargs):
            filter_descriptors.append(
                (kwargs.get("in_bam_fd"), kwargs.get("out_bam_fd"))
            )
            return original_filter(read_ids, in_bam, out_bam, *args, **kwargs)

        monkeypatch.setattr(
            lyra.util_misc.ReadIdStore,
            "filter_bam_by_ids",
            record_filter,
        )

        try:
            if role == "normalized":
                assert lyra._write_normalized(store, stage) == 2
            elif role == "summary":
                lyra._write_summary(store, stage, 1)
            else:
                assert lyra._write_viral_bam(
                    store,
                    stage,
                    work_dir=tmp_path,
                ) == 1
                assert filter_descriptors == [
                    (store.source_bam_fd, stage.descriptor)
                ]

            os.fstat(stage.descriptor)
            assert victim.read_bytes() == b"caller bytes"
            assert lyra._file_identity(os.stat(victim)) == victim_identity
            with pytest.raises(lyra.LyraArtifactConsistencyError) as exc_info:
                lyra._assert_stage_handle(stage)
            assert exc_info.value.category == "stage_identity"
        finally:
            os.unlink(stage.basename, dir_fd=stage.parent.descriptor)
            transaction.rollback_and_cleanup()


def test_transaction_records_full_stage_identities_after_handle_readback(tmp_path):
    with _publication_store(tmp_path, "validated-identities") as (
        store,
        score_path,
        bam_path,
        _,
    ):
        path_plan, _ = _staged_plan(
            tmp_path,
            score_path,
            bam_path,
            ".tsv.zst",
        )
        transaction = lyra.LyraArtifactTransaction(
            store,
            path_plan,
            work_dir=tmp_path,
        )
        transaction._create_stages()
        try:
            transaction._generate_and_validate()

            assert transaction._validated_stage_identities == {
                stage.role: lyra._file_identity(os.fstat(stage.descriptor))
                for stage in transaction.stages
            }
            assert set(transaction._validated_stage_identities) == {
                "normalized",
                "summary",
                "viral_bam",
            }
        finally:
            transaction.rollback_and_cleanup()


@pytest.mark.parametrize(
    ("corruption", "mutate"),
    [
        ("header", lambda value: value.replace(b"SAMPLE_ID", b"sample_id", 1)),
        ("width", lambda value: value.replace(b"\tViral\n", b"\n", 1)),
        ("sample", lambda value: value.replace(b"sample exact", b"other", 1)),
        ("read-id", lambda value: value.replace(b"\talpha\t", b"\tbeta\t", 1)),
        (
            "order",
            lambda value: b"\n".join(
                value.split(b"\n")[:1]
                + list(reversed(value.split(b"\n")[1:-1]))
            )
            + b"\n",
        ),
        (
            "score-count",
            lambda value: value.replace(b"\talpha\t1\t", b"\talpha\t2\t", 1),
        ),
        ("pairing", lambda value: value.replace(b"Single-end", b"single-end", 1)),
        ("call", lambda value: value.replace(b"Viral\n", b"viral\n", 1)),
        ("decimal", lambda value: value.replace(b"\t0.9\t0.9\t", b"\t0.90\t0.9\t", 1)),
        ("row-count", lambda value: b"\n".join(value.split(b"\n")[:-2]) + b"\n"),
        ("missing-lf", lambda value: value.rstrip(b"\n")),
    ],
)
def test_staged_normalized_readback_rejects_exact_contract_corruption(
    tmp_path,
    corruption,
    mutate,
):
    with _publication_store(tmp_path, "normalized-corrupt-" + corruption) as (
        store,
        score_path,
        bam_path,
        _,
    ):
        path_plan, _ = _staged_plan(tmp_path, score_path, bam_path)
        stages, _, _ = _generate_staged_artifacts(store, path_plan, tmp_path)
        try:
            normalized_stage = stages[0].stage_path
            with open(normalized_stage, "rb") as stream:
                contents = stream.read()
            with open(normalized_stage, "wb") as stream:
                stream.write(mutate(contents))

            with pytest.raises(lyra.LyraArtifactConsistencyError) as exc_info:
                lyra._validate_staged_normalized(store, normalized_stage)

            assert exc_info.value.category == "normalized_readback"
        finally:
            _remove_stages(stages)


@pytest.mark.parametrize("normalized_suffix", [".tsv.gz", ".tsv.zst"])
def test_staged_normalized_readback_rejects_truncated_compression(
    tmp_path,
    normalized_suffix,
):
    with _publication_store(tmp_path, "truncated" + normalized_suffix) as (
        store,
        score_path,
        bam_path,
        _,
    ):
        path_plan, _ = _staged_plan(
            tmp_path,
            score_path,
            bam_path,
            normalized_suffix,
        )
        stages, _, _ = _generate_staged_artifacts(store, path_plan, tmp_path)
        try:
            normalized_stage = stages[0].stage_path
            with open(normalized_stage, "r+b") as stream:
                stream.truncate(os.path.getsize(normalized_stage) - 2)

            with pytest.raises(lyra.LyraArtifactConsistencyError) as exc_info:
                lyra._validate_staged_normalized(store, normalized_stage)

            assert exc_info.value.category == "normalized_readback"
            assert isinstance(exc_info.value.__cause__, lyra._COMPRESSION_EXCEPTIONS)
        finally:
            _remove_stages(stages)


@pytest.mark.parametrize(
    ("corruption", "mutate"),
    [
        ("header", lambda value: value.replace(b"SAMPLE_ID", b"sample_id", 1)),
        ("width", lambda value: value.replace(b"\t1\n", b"\n", 1)),
        ("value", lambda value: value.replace(b"sample exact", b"other", 1)),
        ("equation", lambda value: value.rsplit(b"\t1\n", 1)[0] + b"\t0\n"),
        ("extra-row", lambda value: value + b"extra\n"),
        ("trailing-bytes", lambda value: value + b"trailing"),
    ],
)
def test_staged_summary_readback_rejects_corruption(
    tmp_path,
    corruption,
    mutate,
):
    with _publication_store(tmp_path, "summary-corrupt-" + corruption) as (
        store,
        score_path,
        bam_path,
        _,
    ):
        path_plan, _ = _staged_plan(tmp_path, score_path, bam_path)
        stages, normalized_rows, bam_records = _generate_staged_artifacts(
            store,
            path_plan,
            tmp_path,
        )
        try:
            summary_stage = stages[1].stage_path
            with open(summary_stage, "rb") as stream:
                contents = stream.read()
            with open(summary_stage, "wb") as stream:
                stream.write(mutate(contents))

            with pytest.raises(lyra.LyraArtifactConsistencyError) as exc_info:
                lyra._validate_staged_summary(
                    store,
                    summary_stage,
                    normalized_rows,
                    bam_records,
                )

            assert exc_info.value.category == "summary_readback"
        finally:
            _remove_stages(stages)


@pytest.mark.parametrize("corruption", ["header", "record-count", "truncated"])
def test_staged_bam_readback_rejects_header_count_and_truncation(
    tmp_path,
    corruption,
):
    with _publication_store(tmp_path, "bam-corrupt-" + corruption) as (
        store,
        score_path,
        bam_path,
        header,
    ):
        path_plan, _ = _staged_plan(tmp_path, score_path, bam_path)
        stages, normalized_rows, expected_bam_records = _generate_staged_artifacts(
            store,
            path_plan,
            tmp_path,
        )
        try:
            bam_stage = stages[2].stage_path
            if corruption == "header":
                wrong_header = dict(header)
                wrong_header["CO"] = ["changed"]
                _write_bam_records(
                    tmp_path / "bam",
                    os.path.basename(bam_stage),
                    [_segment("alpha")],
                    header=wrong_header,
                )
            elif corruption == "record-count":
                _write_bam_records(
                    tmp_path / "bam",
                    os.path.basename(bam_stage),
                    [],
                    header=header,
                )
            else:
                with open(bam_stage, "r+b") as stream:
                    stream.truncate(max(1, os.path.getsize(bam_stage) - 16))

            if corruption == "record-count":
                with pytest.raises(lyra.LyraArtifactConsistencyError) as exc_info:
                    lyra._validate_staged_artifacts(
                        store,
                        stages[0],
                        stages[1],
                        stages[2],
                        normalized_rows,
                        expected_bam_records,
                    )
                assert exc_info.value.category == "bam_readback"
            else:
                with pytest.raises(lyra.LyraArtifactConsistencyError) as exc_info:
                    lyra._validate_staged_bam(store, bam_stage)
                assert exc_info.value.category == "bam_readback"
        finally:
            _remove_stages(stages)


def test_all_staged_readback_finishes_before_first_final_link(
    tmp_path,
    monkeypatch,
):
    with _publication_store(tmp_path, "readback-before-link") as (
        store,
        score_path,
        bam_path,
        _,
    ):
        path_plan, outputs = _staged_plan(
            tmp_path,
            score_path,
            bam_path,
            ".tsv.zst",
        )
        readbacks = []
        original_normalized = lyra._validate_staged_normalized
        original_bam = lyra._validate_staged_bam
        original_summary = lyra._validate_staged_summary
        original_link = lyra._link_no_clobber

        def validate_normalized(*args):
            assert not any(path.exists() for path in outputs)
            result = original_normalized(*args)
            readbacks.append("normalized")
            return result

        def validate_bam(*args):
            assert not any(path.exists() for path in outputs)
            result = original_bam(*args)
            readbacks.append("bam")
            return result

        def validate_summary(*args):
            assert not any(path.exists() for path in outputs)
            result = original_summary(*args)
            readbacks.append("summary")
            return result

        def link(stage_path, final_path):
            assert readbacks == ["normalized", "bam", "summary"]
            return original_link(stage_path, final_path)

        monkeypatch.setattr(lyra, "_validate_staged_normalized", validate_normalized)
        monkeypatch.setattr(lyra, "_validate_staged_bam", validate_bam)
        monkeypatch.setattr(lyra, "_validate_staged_summary", validate_summary)
        monkeypatch.setattr(lyra, "_link_no_clobber", link)

        lyra.write_lyra_artifacts(store, path_plan, work_dir=tmp_path)

    assert all(path.exists() for path in outputs)
    assert not list(tmp_path.rglob(".lyra-*"))


@pytest.mark.parametrize(
    ("helper_name", "expected_flags"),
    [
        ("_fsync_file", os.O_RDONLY),
        ("_fsync_directory", os.O_RDONLY | os.O_DIRECTORY),
    ],
)
def test_fsync_helpers_open_read_only_close_and_propagate(
    monkeypatch,
    helper_name,
    expected_flags,
):
    operations = []

    def open_path(path, flags):
        operations.append(("open", path, flags))
        return 17

    def fail_fsync(descriptor):
        operations.append(("fsync", descriptor))
        raise OSError(5, "injected fsync failure")

    def close_descriptor(descriptor):
        operations.append(("close", descriptor))

    with monkeypatch.context() as patch:
        patch.setattr(lyra.os, "open", open_path)
        patch.setattr(lyra.os, "fsync", fail_fsync)
        patch.setattr(lyra.os, "close", close_descriptor)

        with pytest.raises(OSError, match="injected fsync failure"):
            getattr(lyra, helper_name)("/artifact/path")

    assert operations == [
        ("open", "/artifact/path", expected_flags),
        ("fsync", 17),
        ("close", 17),
    ]


def _operation_role(path, path_plan):
    for destination in (
        path_plan.normalized,
        path_plan.viral_bam,
        path_plan.summary,
    ):
        if path == destination.final_path:
            return destination.role + "_final"
        if os.path.basename(path).startswith(
            ".lyra-{}-".format(destination.role)
        ):
            return destination.role + "_stage"
        if path == destination.parent_path:
            return destination.role + "_parent"
    raise AssertionError("unknown transaction path: {}".format(path))


def _raise(error):
    raise error


@pytest.mark.parametrize(
    ("boundary", "expected_stage", "normalized_suffix"),
    [
        ("normalized_producer", "generate_normalized", ".tsv"),
        ("normalized_producer", "generate_normalized", ".tsv.gz"),
        ("normalized_producer", "generate_normalized", ".tsv.zst"),
        ("bam_filter", "generate_viral_bam", ".tsv"),
        ("source_recheck", "after_bam_filter", ".tsv"),
        ("generation_equations", "validate_generation_counts", ".tsv"),
        ("normalized_readback", "validate_normalized_readback", ".tsv"),
        ("bam_readback", "validate_bam_readback", ".tsv"),
        ("summary_readback", "validate_summary_readback", ".tsv"),
        ("fsync_normalized", "fsync_normalized_stage", ".tsv"),
        ("fsync_viral_bam", "fsync_viral_bam_stage", ".tsv"),
        ("fsync_summary", "fsync_summary_stage", ".tsv"),
        ("fsync_stage_parent", "fsync_normalized_stage_parent", ".tsv"),
        ("fsync_final_parent", "sync_normalized_final", ".tsv"),
        ("link_normalized", "publish_normalized", ".tsv"),
        ("link_viral_bam", "publish_viral_bam", ".tsv"),
        ("link_summary", "publish_summary", ".tsv"),
        ("stage_cleanup", "remove_normalized_stage", ".tsv"),
        ("final_summary_sync", "sync_summary_final", ".tsv"),
    ],
)
def test_publication_failure_reports_exact_stage_and_original_cause(
    tmp_path,
    monkeypatch,
    boundary,
    expected_stage,
    normalized_suffix,
):
    with _publication_store(tmp_path, "failure-" + boundary) as (
        store,
        score_path,
        bam_path,
        _,
    ):
        path_plan, outputs = _staged_plan(
            tmp_path,
            score_path,
            bam_path,
            normalized_suffix,
        )
        primary = OSError(5, "injected publication boundary failure")

        if boundary == "normalized_producer":
            monkeypatch.setattr(lyra, "_write_normalized", lambda *args: _raise(primary))
        elif boundary == "bam_filter":
            monkeypatch.setattr(
                lyra.util_misc.ReadIdStore,
                "filter_bam_by_ids",
                lambda *args, **kwargs: _raise(primary),
            )
        elif boundary == "source_recheck":
            original = lyra._assert_source_bam_identity
            calls = []

            def fail_second_source_check(current_store, stage):
                calls.append(stage)
                if len(calls) == 2:
                    raise lyra.LyraSourceIdentityError(
                        stage="after_bam_filter",
                        path=current_store.source_bam_display_path,
                        expected=current_store.source_bam_identity,
                        actual_status="unreadable",
                    ) from primary
                return original(current_store, stage)

            monkeypatch.setattr(
                lyra,
                "_assert_source_bam_identity",
                fail_second_source_check,
            )
        elif boundary == "generation_equations":
            monkeypatch.setattr(
                lyra,
                "_validate_artifact_counts",
                lambda *args: _raise(primary),
            )
        elif boundary.endswith("_readback"):
            monkeypatch.setattr(
                lyra,
                "_validate_staged_" + boundary.removesuffix("_readback"),
                lambda *args: _raise(primary),
            )
        elif boundary.startswith("fsync_") and boundary not in {
            "fsync_stage_parent",
            "fsync_final_parent",
        }:
            role = boundary[len("fsync_"):]
            original = lyra._fsync_file

            def fail_selected_stage(path):
                if os.path.basename(path).startswith(".lyra-{}-".format(role)):
                    raise primary
                return original(path)

            monkeypatch.setattr(lyra, "_fsync_file", fail_selected_stage)
        elif boundary in {"fsync_stage_parent", "fsync_final_parent"}:
            original = lyra._fsync_directory
            calls = []

            def fail_selected_directory_sync(path):
                if path == path_plan.normalized.parent_path:
                    calls.append(path)
                    target_call = 1 if boundary == "fsync_stage_parent" else 2
                    if len(calls) == target_call:
                        raise primary
                return original(path)

            monkeypatch.setattr(
                lyra,
                "_fsync_directory",
                fail_selected_directory_sync,
            )
        elif boundary.startswith("link_"):
            role = boundary[len("link_"):]
            destination = getattr(path_plan, role).final_path
            original = lyra._link_no_clobber

            def fail_selected_link(stage_path, final_path):
                if final_path == destination:
                    raise primary
                return original(stage_path, final_path)

            monkeypatch.setattr(lyra, "_link_no_clobber", fail_selected_link)
        elif boundary == "stage_cleanup":
            original = lyra._unlink_path
            failed = []

            def fail_first_normalized_stage_unlink(path):
                if (
                    not failed
                    and os.path.basename(path).startswith(".lyra-normalized-")
                ):
                    failed.append(path)
                    raise primary
                return original(path)

            monkeypatch.setattr(
                lyra,
                "_unlink_path",
                fail_first_normalized_stage_unlink,
            )
        else:
            assert boundary == "final_summary_sync"
            original = lyra._fsync_directory

            def fail_final_summary_sync(path):
                if (
                    path == path_plan.summary.parent_path
                    and os.path.lexists(path_plan.summary.final_path)
                ):
                    raise primary
                return original(path)

            monkeypatch.setattr(
                lyra,
                "_fsync_directory",
                fail_final_summary_sync,
            )

        with pytest.raises(lyra.LyraPublicationError) as exc_info:
            lyra.write_lyra_artifacts(store, path_plan, work_dir=tmp_path)

        error = exc_info.value
        expected_primary = error.__cause__
        if boundary == "source_recheck":
            assert isinstance(expected_primary, lyra.LyraSourceIdentityError)
            assert expected_primary.__cause__ is primary
            assert error.primary_category == "source_bam_identity"
        else:
            assert expected_primary is primary
            assert error.primary_category is None
        assert error.category == "publication"
        assert error.stage == expected_stage
        assert error.primary_type == type(expected_primary).__name__
        assert error.primary_errno == getattr(expected_primary, "errno", None)
        assert error.primary_returncode is None
        assert not any(path.exists() for path in outputs)
        assert not list(tmp_path.rglob(".lyra-*"))


@pytest.mark.parametrize("role", ["normalized", "viral_bam", "summary"])
def test_raced_final_survives_real_no_clobber_link_and_owned_finals_roll_back(
    tmp_path,
    monkeypatch,
    role,
):
    with _publication_store(tmp_path, "race-" + role) as (
        store,
        score_path,
        bam_path,
        _,
    ):
        path_plan, outputs = _staged_plan(tmp_path, score_path, bam_path)
        destination = getattr(path_plan, role).final_path
        original_link = lyra._link_no_clobber
        racer_bytes = ("racer-" + role).encode("ascii")
        raced_identity = []

        def race_selected_link(stage_path, final_path):
            if final_path == destination:
                with open(final_path, "wb") as stream:
                    stream.write(racer_bytes)
                raced_identity.append(lyra._file_identity(os.lstat(final_path)))
            return original_link(stage_path, final_path)

        monkeypatch.setattr(lyra, "_link_no_clobber", race_selected_link)

        with pytest.raises(lyra.LyraPublicationError) as exc_info:
            lyra.write_lyra_artifacts(store, path_plan, work_dir=tmp_path)

        error = exc_info.value

    assert isinstance(error.__cause__, FileExistsError)
    assert error.stage == "publish_{}".format(role)
    assert os.path.exists(destination)
    assert open(destination, "rb").read() == racer_bytes
    assert lyra._file_identity(os.lstat(destination)) == raced_identity[0]
    assert all(
        not path.exists()
        for path in outputs
        if os.fspath(path) != destination
    )
    assert not list(tmp_path.rglob(".lyra-*"))


def test_ambiguous_link_error_claims_exact_stage_identity_for_rollback(
    tmp_path,
    monkeypatch,
):
    with _publication_store(tmp_path, "ambiguous-link") as (
        store,
        score_path,
        bam_path,
        _,
    ):
        path_plan, outputs = _staged_plan(tmp_path, score_path, bam_path)
        original_link = lyra._link_no_clobber
        primary = OSError(5, "ambiguous remote link status")

        def ambiguous_normalized_link(stage_path, final_path):
            original_link(stage_path, final_path)
            raise primary

        monkeypatch.setattr(
            lyra,
            "_link_no_clobber",
            ambiguous_normalized_link,
        )

        with pytest.raises(lyra.LyraPublicationError) as exc_info:
            lyra.write_lyra_artifacts(store, path_plan, work_dir=tmp_path)

        error = exc_info.value

    assert error.__cause__ is primary
    assert error.stage == "publish_normalized"
    assert not any(path.exists() for path in outputs)
    assert not list(tmp_path.rglob(".lyra-*"))


@pytest.mark.parametrize("role", ["normalized", "viral_bam", "summary"])
def test_replacement_before_rollback_identity_observation_survives_and_is_reported(
    tmp_path,
    monkeypatch,
    role,
):
    with _publication_store(tmp_path, "replacement-" + role) as (
        store,
        score_path,
        bam_path,
        _,
    ):
        path_plan, outputs = _staged_plan(tmp_path, score_path, bam_path)
        transaction = lyra.LyraArtifactTransaction(
            store,
            path_plan,
            work_dir=tmp_path,
        )
        target = getattr(path_plan, role).final_path
        replacement_bytes = ("replacement-" + role).encode("ascii")
        real_lstat = os.lstat
        real_unlink = os.unlink
        replacement_identity = []
        publication_failed = []
        original_fsync_directory = lyra._fsync_directory

        def fail_final_success_barrier(path):
            if (
                transaction.stage == "sync_summary_final"
                and path == path_plan.summary.parent_path
            ):
                publication_failed.append(True)
                raise OSError(5, "injected final success barrier failure")
            return original_fsync_directory(path)

        def replace_at_rollback_observation(path, *args, **kwargs):
            if (
                not args
                and not kwargs
                and publication_failed
                and not replacement_identity
                and os.fspath(path) == target
            ):
                real_unlink(target)
                with open(target, "wb") as stream:
                    stream.write(replacement_bytes)
                replacement_identity.append(
                    lyra._file_identity(real_lstat(target))
                )
            return real_lstat(path, *args, **kwargs)

        monkeypatch.setattr(
            lyra,
            "_fsync_directory",
            fail_final_success_barrier,
        )
        monkeypatch.setattr(lyra.os, "lstat", replace_at_rollback_observation)

        with pytest.raises(lyra.LyraPublicationError) as exc_info:
            transaction.generate_validate_and_publish()

        error = exc_info.value

    assert open(target, "rb").read() == replacement_bytes
    assert lyra._file_identity(real_lstat(target)) == replacement_identity[0]
    mismatch = [
        failure
        for failure in error.cleanup_failures
        if failure.category == "rollback_identity_mismatch"
    ]
    assert [(failure.role, failure.path) for failure in mismatch] == [
        (role, target)
    ]
    assert all(
        not path.exists()
        for path in outputs
        if os.fspath(path) != target
    )
    assert not list(tmp_path.rglob(".lyra-*"))


def test_cleanup_continues_after_failure_and_ignores_tmpkeep(
    tmp_path,
    monkeypatch,
):
    monkeypatch.setenv("VIRAL_NGS_TMP_DIRKEEP", "1")
    with _publication_store(tmp_path, "cleanup-continues") as (
        store,
        score_path,
        bam_path,
        _,
    ):
        path_plan, outputs = _staged_plan(tmp_path, score_path, bam_path)
        primary = RuntimeError("primary validation failure")
        original_unlink = lyra._unlink_path
        attempts = []

        monkeypatch.setattr(
            lyra,
            "_validate_staged_artifacts",
            lambda *args, **kwargs: _raise(primary),
        )

        def fail_one_stage_unlink(path):
            attempts.append(path)
            if os.path.basename(path).startswith(".lyra-normalized-"):
                raise OSError(5, "injected cleanup unlink failure")
            return original_unlink(path)

        monkeypatch.setattr(lyra, "_unlink_path", fail_one_stage_unlink)

        with pytest.raises(lyra.LyraPublicationError) as exc_info:
            lyra.write_lyra_artifacts(store, path_plan, work_dir=tmp_path)

        error = exc_info.value

    assert error.__cause__ is primary
    attempted_roles = [
        role
        for role in ("normalized", "summary", "viral_bam")
        if any(
            os.path.basename(path).startswith(".lyra-{}-".format(role))
            for path in attempts
        )
    ]
    assert attempted_roles == ["normalized", "summary", "viral_bam"]
    assert [failure.category for failure in error.cleanup_failures] == [
        "stage_cleanup_unlink_failed"
    ]
    remaining_stages = list(tmp_path.rglob(".lyra-*"))
    assert len(remaining_stages) == 1
    assert remaining_stages[0].name.startswith(".lyra-normalized-")
    assert not any(path.exists() for path in outputs)
    remaining_stages[0].unlink()


def test_bam_filter_and_transaction_cleanup_diagnostics_are_combined(
    tmp_path,
    monkeypatch,
):
    core_primary = OSError(5, "injected BAM-filter wait failure")
    core_primary.cleanup_failures = (
        viral_ngs.core.misc.ReadIdCleanupFailure(
            operation="writer_stdin_close",
            error_type="BrokenPipeError",
            errno=None,
        ),
        viral_ngs.core.misc.ReadIdCleanupFailure(
            operation="reader_wait",
            error_type="OSError",
            errno=5,
        ),
    )
    original_unlink = lyra._unlink_path

    def fail_bam_filter(*args, **kwargs):
        raise core_primary

    def fail_one_stage_unlink(path):
        if os.path.basename(path).startswith(".lyra-normalized-"):
            raise OSError(6, "injected stage cleanup failure")
        return original_unlink(path)

    monkeypatch.setattr(
        viral_ngs.core.misc.ReadIdStore,
        "filter_bam_by_ids",
        fail_bam_filter,
    )
    monkeypatch.setattr(lyra, "_unlink_path", fail_one_stage_unlink)

    with _publication_store(tmp_path, "bam-diagnostic") as (
        store,
        score_path,
        bam_path,
        _,
    ):
        path_plan, outputs = _staged_plan(tmp_path, score_path, bam_path)
        with pytest.raises(lyra.LyraPublicationError) as exc_info:
            lyra.write_lyra_artifacts(store, path_plan, work_dir=tmp_path)

    error = exc_info.value
    assert error.__cause__ is core_primary
    assert [
        (failure.category, failure.operation, failure.error_type, failure.errno)
        for failure in error.cleanup_failures
    ] == [
        (
            "stage_cleanup_unlink_failed",
            "cleanup",
            "OSError",
            6,
        ),
        (
            "bam_filter_cleanup",
            "writer_stdin_close",
            "BrokenPipeError",
            None,
        ),
        (
            "bam_filter_cleanup",
            "reader_wait",
            "OSError",
            5,
        ),
    ]
    assert error.cleanup_failures_truncated is False
    assert not any(path.exists() for path in outputs)
    remaining_stages = list(tmp_path.rglob(".lyra-*"))
    assert len(remaining_stages) == 1
    remaining_stages[0].unlink()


def test_primary_and_multiple_rollback_failures_share_one_ordered_error(
    tmp_path,
    monkeypatch,
):
    with _publication_store(tmp_path, "multiple-rollback-failures") as (
        store,
        score_path,
        bam_path,
        _,
    ):
        path_plan, _ = _staged_plan(tmp_path, score_path, bam_path)
        transaction = lyra.LyraArtifactTransaction(
            store,
            path_plan,
            work_dir=tmp_path,
        )
        primary = OSError(5, "injected completion-barrier failure")
        original_fsync_directory = lyra._fsync_directory
        original_unlink = lyra._unlink_path
        failed_roles = {"summary", "viral_bam"}

        def fail_completion_barrier(path):
            if (
                transaction.stage == "sync_summary_final"
                and path == path_plan.summary.parent_path
            ):
                raise primary
            return original_fsync_directory(path)

        def fail_selected_rollbacks(path):
            for role in failed_roles:
                if path == getattr(path_plan, role).final_path:
                    raise OSError(5, "injected rollback unlink failure")
            return original_unlink(path)

        monkeypatch.setattr(
            lyra,
            "_fsync_directory",
            fail_completion_barrier,
        )
        monkeypatch.setattr(lyra, "_unlink_path", fail_selected_rollbacks)

        with pytest.raises(lyra.LyraPublicationError) as exc_info:
            transaction.generate_validate_and_publish()

        error = exc_info.value

    assert error.__cause__ is primary
    assert [
        (failure.category, failure.role)
        for failure in error.cleanup_failures
    ] == [
        ("rollback_unlink_failed", "summary"),
        ("rollback_unlink_failed", "viral_bam"),
    ]
    for role in failed_roles:
        path = getattr(path_plan, role).final_path
        assert os.path.exists(path)
        os.unlink(path)
    assert not os.path.exists(path_plan.normalized.final_path)
    assert not list(tmp_path.rglob(".lyra-*"))


def test_system_exiting_base_exception_is_preserved_after_cleanup(
    tmp_path,
    monkeypatch,
):
    with _publication_store(tmp_path, "system-exit") as (
        store,
        score_path,
        bam_path,
        _,
    ):
        path_plan, outputs = _staged_plan(tmp_path, score_path, bam_path)
        primary = SystemExit(17)
        monkeypatch.setattr(
            lyra,
            "_validate_staged_artifacts",
            lambda *args, **kwargs: _raise(primary),
        )

        with pytest.raises(SystemExit) as exc_info:
            lyra.write_lyra_artifacts(store, path_plan, work_dir=tmp_path)

    assert exc_info.value is primary
    assert not any(path.exists() for path in outputs)
    assert not list(tmp_path.rglob(".lyra-*"))


def test_publication_diagnostic_is_bounded_and_caps_ordered_cleanup_facts():
    long_value = "x" * 10000
    primary = subprocess.CalledProcessError(
        23,
        ["samtools", long_value],
    )
    cleanup_failures = tuple(
        lyra.CleanupFailure(
            operation="cleanup",
            role="role-{}".format(index),
            path="/output/{}".format(long_value),
            error_type="OSError",
            errno=5,
            category="stage_cleanup_unlink_failed",
        )
        for index in range(20)
    )

    error = lyra.LyraPublicationError(
        stage="publish_summary",
        primary=primary,
        cleanup_failures=cleanup_failures,
    )

    assert error.category == "publication"
    assert error.stage == "publish_summary"
    assert error.primary_type == "CalledProcessError"
    assert error.primary_category is None
    assert error.primary_errno is None
    assert error.primary_returncode == 23
    assert len(error.cleanup_failures) == 16
    assert [failure.role for failure in error.cleanup_failures] == [
        "role-{}".format(index) for index in range(16)
    ]
    assert error.cleanup_failures_truncated is True
    assert long_value not in str(error)
    assert len(str(error)) < 6000
    with pytest.raises(FrozenInstanceError):
        error.cleanup_failures[0].role = "changed"


def test_publication_diagnostic_never_stringifies_arbitrary_primary():
    class HostileError(RuntimeError):
        def __str__(self):
            raise AssertionError("primary __str__ must not run")

        def __repr__(self):
            raise AssertionError("primary __repr__ must not run")

    primary = HostileError()

    error = lyra.LyraPublicationError(
        stage="generate_normalized",
        primary=primary,
    )

    assert error.primary_type == "HostileError"
    assert "HostileError" in str(error)


def test_publish_order_flushes_every_stage_and_syncs_summary_parent_last(
    tmp_path,
    monkeypatch,
):
    with _publication_store(tmp_path, "publish-order") as (
        store,
        score_path,
        bam_path,
        _,
    ):
        path_plan, outputs = _staged_plan(tmp_path, score_path, bam_path)
        operations = []
        original_validate = lyra._validate_staged_artifacts
        original_fsync_file = lyra._fsync_file
        original_fsync_directory = lyra._fsync_directory
        original_link = lyra._link_no_clobber
        original_unlink = lyra._unlink_path

        def validate(*args, **kwargs):
            result = original_validate(*args, **kwargs)
            operations.append(("validated", "all"))
            return result

        def fsync_file(path):
            operations.append(("fsync_file", _operation_role(path, path_plan)))
            return original_fsync_file(path)

        def fsync_directory(path):
            operations.append(("fsync_dir", _operation_role(path, path_plan)))
            return original_fsync_directory(path)

        def link(stage_path, final_path):
            operations.append(("link", _operation_role(final_path, path_plan)))
            return original_link(stage_path, final_path)

        def unlink(path):
            operations.append(("unlink", _operation_role(path, path_plan)))
            return original_unlink(path)

        monkeypatch.setattr(lyra, "_validate_staged_artifacts", validate)
        monkeypatch.setattr(lyra, "_fsync_file", fsync_file)
        monkeypatch.setattr(lyra, "_fsync_directory", fsync_directory)
        monkeypatch.setattr(lyra, "_link_no_clobber", link)
        monkeypatch.setattr(lyra, "_unlink_path", unlink)

        lyra.write_lyra_artifacts(store, path_plan, work_dir=tmp_path)

    assert operations == [
        ("validated", "all"),
        ("fsync_file", "normalized_stage"),
        ("fsync_file", "viral_bam_stage"),
        ("fsync_file", "summary_stage"),
        ("fsync_dir", "normalized_parent"),
        ("fsync_dir", "viral_bam_parent"),
        ("fsync_dir", "summary_parent"),
        ("link", "normalized_final"),
        ("fsync_dir", "normalized_parent"),
        ("link", "viral_bam_final"),
        ("fsync_dir", "viral_bam_parent"),
        ("unlink", "normalized_stage"),
        ("fsync_dir", "normalized_parent"),
        ("unlink", "viral_bam_stage"),
        ("fsync_dir", "viral_bam_parent"),
        ("link", "summary_final"),
        ("unlink", "summary_stage"),
        ("fsync_dir", "summary_parent"),
    ]
    assert all(path.exists() for path in outputs)
    assert not list(tmp_path.rglob(".lyra-*"))


def test_same_parent_publication_keeps_every_directory_transition_barrier(
    tmp_path,
    monkeypatch,
):
    with _publication_store(tmp_path, "same-parent") as (
        store,
        score_path,
        bam_path,
        _,
    ):
        outputs = _artifact_paths(tmp_path, "same-parent")
        path_plan = lyra._build_artifact_path_plan(
            score_path,
            bam_path,
            *outputs,
        )
        synced_directories = []
        original_fsync_directory = lyra._fsync_directory

        def record_directory_sync(path):
            synced_directories.append(path)
            return original_fsync_directory(path)

        monkeypatch.setattr(lyra, "_fsync_directory", record_directory_sync)

        lyra.write_lyra_artifacts(store, path_plan, work_dir=tmp_path)

    assert synced_directories == [os.path.realpath(tmp_path, strict=True)] * 8
    assert all(path.exists() for path in outputs)


@pytest.mark.parametrize("failure_kind", ["file", "directory"])
def test_fsync_failure_prevents_summary_and_cleans_all_stages(
    tmp_path,
    monkeypatch,
    failure_kind,
):
    with _publication_store(tmp_path, "fsync-failure-" + failure_kind) as (
        store,
        score_path,
        bam_path,
        _,
    ):
        path_plan, outputs = _staged_plan(tmp_path, score_path, bam_path)
        helper_name = "_fsync_file" if failure_kind == "file" else "_fsync_directory"

        def fail(_path):
            raise OSError(5, "injected {} fsync failure".format(failure_kind))

        monkeypatch.setattr(lyra, helper_name, fail)

        with pytest.raises(lyra.LyraPublicationError) as exc_info:
            lyra.write_lyra_artifacts(store, path_plan, work_dir=tmp_path)

    assert isinstance(exc_info.value.__cause__, OSError)
    assert not any(path.exists() for path in outputs)
    assert not list(tmp_path.rglob(".lyra-*"))


def test_final_summary_parent_fsync_failure_rolls_back_every_final(
    tmp_path,
    monkeypatch,
):
    with _publication_store(tmp_path, "summary-fsync-failure") as (
        store,
        score_path,
        bam_path,
        _,
    ):
        path_plan, outputs = _staged_plan(tmp_path, score_path, bam_path)
        original_fsync_directory = lyra._fsync_directory

        def fail_final_summary_sync(path):
            if (
                path == path_plan.summary.parent_path
                and os.path.lexists(path_plan.summary.final_path)
            ):
                raise OSError(5, "injected final summary fsync failure")
            return original_fsync_directory(path)

        monkeypatch.setattr(
            lyra,
            "_fsync_directory",
            fail_final_summary_sync,
        )

        with pytest.raises(lyra.LyraPublicationError) as exc_info:
            lyra.write_lyra_artifacts(store, path_plan, work_dir=tmp_path)

    assert exc_info.value.stage == "sync_summary_final"
    assert isinstance(exc_info.value.__cause__, OSError)
    assert not any(path.exists() for path in outputs)
    assert not list(tmp_path.rglob(".lyra-*"))


def test_bam_link_failure_rolls_back_normalized_and_leaves_no_summary(
    tmp_path,
    monkeypatch,
):
    with _publication_store(tmp_path, "bam-link-failure") as (
        store,
        score_path,
        bam_path,
        _,
    ):
        path_plan, outputs = _staged_plan(tmp_path, score_path, bam_path)
        original_link = lyra._link_no_clobber

        def fail_bam_link(stage_path, final_path):
            if final_path == path_plan.viral_bam.final_path:
                raise OSError(5, "injected BAM link failure")
            return original_link(stage_path, final_path)

        monkeypatch.setattr(lyra, "_link_no_clobber", fail_bam_link)

        with pytest.raises(lyra.LyraPublicationError) as exc_info:
            lyra.write_lyra_artifacts(store, path_plan, work_dir=tmp_path)

    assert exc_info.value.stage == "publish_viral_bam"
    assert isinstance(exc_info.value.__cause__, OSError)
    assert not any(path.exists() for path in outputs)
    assert not list(tmp_path.rglob(".lyra-*"))


def test_stage_cleanup_ignores_tmpkeep_after_validation_failure(
    tmp_path,
    monkeypatch,
):
    monkeypatch.setenv("VIRAL_NGS_TMP_DIRKEEP", "1")
    with _publication_store(tmp_path, "tmpkeep-failure") as (
        store,
        score_path,
        bam_path,
        _,
    ):
        path_plan, outputs = _staged_plan(tmp_path, score_path, bam_path)

        def fail_validation(*args, **kwargs):
            raise RuntimeError("injected validation failure")

        monkeypatch.setattr(lyra, "_validate_staged_artifacts", fail_validation)

        with pytest.raises(lyra.LyraPublicationError) as exc_info:
            lyra.write_lyra_artifacts(store, path_plan, work_dir=tmp_path)

    assert isinstance(exc_info.value.__cause__, RuntimeError)
    assert not any(path.exists() for path in outputs)
    assert not list(tmp_path.rglob(".lyra-*"))


def test_rollback_preserves_observed_replacement_and_reports_identity_mismatch(
    tmp_path,
    monkeypatch,
):
    with _publication_store(tmp_path, "rollback-replacement") as (
        store,
        score_path,
        bam_path,
        _,
    ):
        path_plan, outputs = _staged_plan(tmp_path, score_path, bam_path)
        transaction = lyra.LyraArtifactTransaction(
            store,
            path_plan,
            work_dir=tmp_path,
        )
        original_link = lyra._link_no_clobber

        def replace_before_failure(stage_path, final_path):
            if final_path == path_plan.viral_bam.final_path:
                os.unlink(path_plan.normalized.final_path)
                with open(path_plan.normalized.final_path, "wb") as stream:
                    stream.write(b"caller replacement")
                raise OSError(5, "injected BAM link failure")
            return original_link(stage_path, final_path)

        monkeypatch.setattr(lyra, "_link_no_clobber", replace_before_failure)

        with pytest.raises(lyra.LyraPublicationError) as exc_info:
            transaction.generate_validate_and_publish()

        outcomes = transaction.cleanup_outcomes

    assert exc_info.value.stage == "publish_viral_bam"
    assert isinstance(exc_info.value.__cause__, OSError)
    assert outputs[0].read_bytes() == b"caller replacement"
    assert not outputs[1].exists()
    assert not outputs[2].exists()
    assert any(
        outcome.operation == "rollback"
        and outcome.role == "normalized"
        and outcome.status == "identity_mismatch"
        for outcome in outcomes
    )
    assert not list(tmp_path.rglob(".lyra-*"))


def test_publication_records_are_frozen_and_transaction_state_is_observable(
    tmp_path,
):
    score_path = _write_scores(tmp_path, "transaction-scores.tsv")
    bam_path, _ = _write_bam(tmp_path, "transaction-source.bam")
    path_plan, _ = _staged_plan(tmp_path, score_path, bam_path)
    transaction = lyra.LyraArtifactTransaction(object(), path_plan)

    assert [field.name for field in fields(lyra.PublishedArtifact)] == [
        "role",
        "final_path",
        "identity",
    ]
    assert transaction.stage == "initialized"
    assert transaction.stages == ()
    assert transaction.published == ()
    assert transaction.cleanup_outcomes == ()
