"""No invented flower labels or cross-reviewer leakage in the 400-photo handoff."""
import csv
from pathlib import Path
import pytest

from scripts.analysis.mois_prepare_400_botanical_review_handoff import (
    FIELDS, SOURCE_MAP, STAGE, TISSUE, make
)


def _csv(tmp_path:Path,duplicate:bool=False)->Path:
    file=tmp_path/"sources.csv"
    with file.open("w",newline="") as f:
        fields=[k for k in SOURCE_MAP if k!="review_token"]+["attribution"]
        w=csv.DictWriter(f,fieldnames=fields,lineterminator="\n")
        w.writeheader()
        for i in (1,2):
            w.writerow({
              "inat_taxon_id":str(i),"species":f"Plant species {i}",
              "photo_id":str(i),"observation_id":f"o{i}",
              "observer_id":"same" if duplicate else f"u{i}",
              "observed_on":"2025-06-01","photo_license":"cc-by",
              "photo_url_large":f"https://example.invalid/only-for-synthetic-test/{i}",
              "attribution":"synthetic CC-BY test attribution",
              "source_rank":"abc","site_latcell":"510","site_loncell":"300",
              "season_year":"2025","selected_season":"summer",
            })
    return file


def test_reviewer_forms_both_unfilled_and_independent(tmp_path:Path):
    f=_csv(tmp_path)
    d=make(f,tmp_path/"review",fixture=True)
    assert d["n_original_photos"]==2 and d["distinct_observers"]==2
    assert d["all_review_response_fields_blank"] is True
    assert d["seasonal_white_fraction_estimated"] is False
    for reviewer in ("A","B"):
        with (tmp_path/"review"/f"mois_400_botanical_open_flower_reviewer_{reviewer}_UNFILLED.csv").open() as fin:
            r=list(csv.DictReader(fin))
        assert len(r)==2
        assert set(r[0])==set(FIELDS)
        assert all(not row["open_flower_stage"] and not row["roi_organ"]
                   and not row["reviewer_id"] for row in r)
        assert all("observed_on" not in row and "source_rank" not in row
                   and "selected_season" not in row for row in r)
    assert "open_flower" in STAGE and "petal_tepal" in TISSUE


def test_reused_original_observer_is_blocked(tmp_path:Path):
    with pytest.raises(ValueError,match="observer reused"):
        make(_csv(tmp_path,duplicate=True),tmp_path/"no",fixture=True)


def test_immutable_source_sha_blocks_fake_full_review(tmp_path:Path):
    with pytest.raises(ValueError,match="not the exact"):
        make(_csv(tmp_path),tmp_path/"no",fixture=False)
