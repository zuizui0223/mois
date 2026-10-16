"""REAL original FCP ROI-v4 adapter interface verification, no photo download.

The earlier Mois ROI tests injected a fake backend: only this separate
integration test actually imports the byte-pinned ORIGINAL FCP repository
wrapper and invokes its actual terminal/schema/image-URL security checks.
It does NOT fabricate flowering-photo pixel measurements.
"""
from __future__ import annotations
import hashlib
import os
import inspect
from pathlib import Path

import pandas as pd
import pytest

from scripts.acquisition.mois_execute_all_dated_original_roi_workunit_20261016 import (
    FCP_BACKEND_GIT_BLOB,BIO as MOIS_BIO
)
from scripts.acquisition.mois_execute_sample_roi_v4_20261016 import (
    backend_from_fcp
)

SOURCE_COMMIT="9fae6ccdf684a46026f72ba12e98de2c5c54bf2a"
FCP_WRAPPER_BLOB="45de671444c0fac4fb4611375d19ed7e5a1608f3"


def fcp_original():
    value=os.environ["MOIS_ORIGINAL_FCP_REPOSITORY"]
    root=Path(value).resolve()
    assert root.is_dir()
    assert (root/"fcp_pipeline/fcp_fresh_study.py").is_file()
    return backend_from_fcp(root)


def test_actual_original_FCP_backend_byte_identity_and_full_real_function_abi():
    backend=fcp_original()
    assert FCP_BACKEND_GIT_BLOB==FCP_WRAPPER_BLOB
    assert backend._git_blob(Path(backend.__file__))==FCP_WRAPPER_BLOB
    assert tuple(backend.BIO)==MOIS_BIO
    assert tuple(backend.WORKER)==("measurement_id","image_filename","photo_license")
    assert tuple(backend.MORPHS)==(
        "white","yellow_orange","red_pink","blue_purple")
    assert backend.SOURCE_METHOD_COMMIT==SOURCE_COMMIT
    assert "verify_frozen_method" in dir(backend)
    assert "image_loader" in dir(backend)
    assert "evaluate_with_frozen_source" in dir(backend)
    assert "validate_terminal" in dir(backend)
    assert "PermanentOriginalPhotoError" in dir(backend)
    assert "ProviderTransientAccessError" in dir(backend)
    assert "ProviderRateGateError" in dir(backend)
    assert list(inspect.signature(backend.verify_frozen_method).parameters)==[
        "source","detector","sam"
    ]
    assert list(inspect.signature(backend.evaluate_with_frozen_source).parameters)==[
        "worker","images","source","detector","sam","out"
    ]
    assert list(inspect.signature(backend.validate_terminal).parameters)==[
        "rows","expected"
    ]
    assert list(inspect.signature(backend.image_loader).parameters)==["url"]


def test_original_FCP_terminal_accepts_true_colour_and_failure_without_fake_nonwhite():
    backend=fcp_original()
    rows=[
        {
            "measurement_id":"FCPR-"+"A"*24,
            "morph":"white",
            "measurement_status":"classified_four_state_morph",
            **{"flower_fraction_"+x:float(x=="white") for x in backend.BIO}
        },
        {
            "measurement_id":"FCPR-"+"B"*24,
            "morph":"mixed_uncertain",
            "measurement_status":"image_acquisition_failed",
            **{"flower_fraction_"+x:0.0 for x in backend.BIO}
        }
    ]
    terminal=pd.DataFrame(rows)
    expected=[r["measurement_id"] for r in rows]
    assert backend.validate_terminal(terminal,expected) is None
    changed=terminal.copy()
    changed.loc[1,"morph"]="blue_purple"
    with pytest.raises(ValueError,match="Original 4-state classification gate"):
        backend.validate_terminal(changed,expected)
    changed=terminal.copy()
    changed.loc[0,"flower_fraction_white"]=.7
    with pytest.raises(ValueError,match="nine-colour composition invalid"):
        backend.validate_terminal(changed,expected)
    with pytest.raises(ValueError,match="every fresh source image exactly once"):
        backend.validate_terminal(terminal,expected[:1])


def test_original_FCP_provider_uri_rejects_non_original_images_before_any_network():
    backend=fcp_original()
    for url in (
        "https://example.com/photos/200/large.jpg",
        "http://inaturalist-open-data.s3.amazonaws.com/photos/200/large.jpg",
        "https://inaturalist-open-data.s3.amazonaws.com/photos/200/large.jpg?unsafe=1",
        "https://inaturalist-open-data.s3.amazonaws.com/photos/hello/large.jpg",
    ):
        with pytest.raises(ValueError,match="not a frozen original provider"):
            backend.image_loader(url)


def test_original_ROI_method_hash_gate_fails_with_nonoriginal_source_before_pixels():
    backend=fcp_original()
    absent=Path(os.environ["MOIS_ORIGINAL_FCP_REPOSITORY"])/"README.md"
    with pytest.raises(ValueError,match="Original frozen ROI-v4 measurement source changed"):
        backend.verify_frozen_method(absent,absent,absent)
