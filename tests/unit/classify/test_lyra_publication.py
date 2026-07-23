from dataclasses import fields
from dataclasses import FrozenInstanceError
from decimal import Decimal
import inspect
import os

import pytest
import pysam

from viral_ngs.classify import lyra


def _segment(query_name):
    record = pysam.AlignedSegment()
    record.query_name = query_name
    record.flag = 0x4
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


def _artifact_paths(tmp_path, prefix):
    return (
        tmp_path / (prefix + "-normalized.tsv"),
        tmp_path / (prefix + "-summary.tsv"),
        tmp_path / (prefix + "-viral.bam"),
    )


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
        "normalized_output",
        "summary_output",
        "viral_bam_output",
        "work_dir",
    ]


@pytest.mark.parametrize("alias_kind", ["direct", "symlink", "hardlink"])
def test_source_direct_symlink_and_hardlink_aliases_produce_artifacts(
    tmp_path,
    alias_kind,
):
    source_path, header = _write_bam(tmp_path, "source-{}.bam".format(alias_kind))
    access_path = source_path
    if alias_kind == "symlink":
        access_path = tmp_path / "source-symlink.bam"
        access_path.symlink_to(source_path.name)
    elif alias_kind == "hardlink":
        access_path = tmp_path / "source-hardlink.bam"
        os.link(source_path, access_path)
    score_path = _write_scores(tmp_path, "scores-{}.tsv".format(alias_kind))
    outputs = _artifact_paths(tmp_path, alias_kind)

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
        lyra.write_lyra_artifacts(store, *outputs, work_dir=tmp_path)

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
            lyra.write_lyra_artifacts(store, *outputs, work_dir=tmp_path)

    assert exc_info.value.stage == "pre_generation"
    assert producer_calls == []
    assert not any(path.exists() for path in outputs)


def test_source_missing_before_generation_preserves_filesystem_cause(tmp_path):
    source_path, _ = _write_bam(tmp_path, "missing-source.bam")
    score_path = _write_scores(tmp_path, "missing-scores.tsv")
    outputs = _artifact_paths(tmp_path, "missing")

    with lyra.reconcile_lyra_fragments(
        score_path,
        source_path,
        "sample",
        "0.8",
        work_dir=tmp_path,
    ) as store:
        source_path.unlink()
        with pytest.raises(lyra.LyraSourceIdentityError) as exc_info:
            lyra.write_lyra_artifacts(store, *outputs, work_dir=tmp_path)

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
            lyra.write_lyra_artifacts(store, *outputs, work_dir=tmp_path)

    assert exc_info.value.stage == "after_bam_filter"
    assert accepted_calls == []
    assert not outputs[1].exists()
