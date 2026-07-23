from dataclasses import fields
from dataclasses import FrozenInstanceError
from decimal import Decimal

import pytest

from viral_ngs.classify import lyra


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
