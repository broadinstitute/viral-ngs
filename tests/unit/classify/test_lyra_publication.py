from contextlib import contextmanager
from dataclasses import fields
from dataclasses import FrozenInstanceError
from decimal import Decimal
import inspect
import os
import stat

import pytest
import pysam

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
    stages = tuple(
        lyra._create_artifact_stage(destination)
        for destination in (
            path_plan.normalized,
            path_plan.summary,
            path_plan.viral_bam,
        )
    )
    normalized_rows = lyra._write_normalized(store, stages[0].stage_path)
    bam_records = lyra._write_viral_bam(
        store,
        stages[2].stage_path,
        work_dir=work_dir,
    )
    lyra._validate_artifact_counts(store.counts, normalized_rows, bam_records)
    lyra._write_summary(store, stages[1].stage_path, bam_records)
    return stages, normalized_rows, bam_records


def _remove_stages(stages):
    for stage in stages:
        if os.path.lexists(stage.stage_path):
            os.unlink(stage.stage_path)


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

    with pytest.raises(lyra.LyraPathError) as exc_info:
        lyra.write_lyra_artifacts(object(), path_plan, work_dir=tmp_path)

    assert exc_info.value.category == "output_exists"
    assert exc_info.value.role == "normalized"
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
    store._install_source_metadata(
        source_bam_path=str(source_path.absolute()),
        source_bam_display_path=str(source_path),
        source_bam_identity=identity,
        source_bam_header=header,
    )


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
    identity = lyra.FileIdentity(1, 2, 3, 4)
    retained_header = {
        "HD": {"VN": "1.6"},
        "SQ": [{"SN": "ref", "LN": 100}],
    }
    store = lyra.LyraFragmentStore("sample", work_dir=tmp_path)

    for property_name in (
        "source_bam_path",
        "source_bam_display_path",
        "source_bam_identity",
        "source_bam_header",
    ):
        with pytest.raises(RuntimeError, match="not finalized"):
            getattr(store, property_name)

    _finalize_source_store(store, source_path, identity, retained_header)

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
        "source_bam_path",
        "source_bam_display_path",
        "source_bam_identity",
        "source_bam_header",
    ):
        with pytest.raises(RuntimeError, match="closed"):
            getattr(store, property_name)


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

        with pytest.raises(lyra.LyraSourceIdentityError) as exc_info:
            lyra.write_lyra_artifacts(store, path_plan, work_dir=tmp_path)

    assert exc_info.value.stage == "pre_generation"
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
        with pytest.raises(lyra.LyraSourceIdentityError) as exc_info:
            lyra.write_lyra_artifacts(store, path_plan, work_dir=tmp_path)

    assert exc_info.value.stage == "pre_generation"
    assert exc_info.value.actual is None
    assert exc_info.value.actual_status == "missing"
    assert isinstance(exc_info.value.__cause__, FileNotFoundError)
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
        with pytest.raises(lyra.LyraSourceIdentityError) as exc_info:
            lyra.write_lyra_artifacts(store, path_plan, work_dir=tmp_path)

    assert exc_info.value.stage == "after_bam_filter"
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
    stages = tuple(lyra._create_artifact_stage(item) for item in destinations)
    try:
        assert [field.name for field in fields(lyra.ArtifactStage)] == [
            "role",
            "stage_path",
            "destination",
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
                stages[0].stage_path,
            ) == expected_rows
            assert lyra._validate_staged_bam(
                store,
                stages[2].stage_path,
            ) == expected_bam
            lyra._validate_staged_summary(
                store,
                stages[1].stage_path,
                expected_rows,
                expected_bam,
            )
            assert not any(path.exists() for path in outputs)
        finally:
            _remove_stages(stages)


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
        stages, _, expected_bam_records = _generate_staged_artifacts(
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
                assert lyra._validate_staged_bam(store, bam_stage) == 0
                assert expected_bam_records == 1
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
