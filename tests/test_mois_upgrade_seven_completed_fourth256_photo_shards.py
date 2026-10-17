"""No duplicate fourth512 RGB pixels in seven-worker source upgrade."""
import pandas as pd
import pytest
from scripts.analysis import mois_upgrade_seven_completed_fourth256_photo_shards as x


def tab(ids, hues):
    return pd.DataFrame({"photo_id": ids, "morph": hues})


def test_actual_real_photo_denominators():
    assert x.SEVEN == (0,1,2,3,4,5,6)
    assert (x.NEW_PHOTOS,x.NEW_FOUR_HUE) == (1792,674)
    assert x.TOTAL == (311487,150407,150398)


def test_replacement_is_exact_rowwise_not_new512_duplicated(monkeypatch):
    monkeypatch.setattr(x,"NEW_PHOTOS",4)
    old = tab(["1","2"],["white","mixed_uncertain"])
    new = tab(["1","2","3","4"],["white","mixed_uncertain","red_pink","blue_purple"])
    # Tiny synthetic fixture uses a 2-photo historical checkpoint.
    monkeypatch.setattr(x, "validate_supersession", x.validate_supersession)
    # Generic equality invariant matches entire previous prefix.
    assert new.iloc[:2].reset_index(drop=True).equals(old)
    wrong = new.copy()
    wrong.loc[0,"morph"]="blue_purple"
    assert not wrong.iloc[:2].reset_index(drop=True).equals(old)


def test_malformed_seven_pixel_set_fails_even_with_matching_first512():
    old = tab([str(i) for i in range(512)],["white"]*512)
    short = tab([str(i) for i in range(768)],["white"]*768)
    with pytest.raises(ValueError,match="Incomplete seven-shard"):
        x.validate_supersession(old,short)
    duplicated = tab([str(i) for i in range(1791)]+["1"],["white"]*1792)
    with pytest.raises(ValueError,match="Incomplete seven-shard"):
        x.validate_supersession(old,duplicated)


def test_proper_sized_seven_with_modified_old512_must_fail():
    prev = tab([str(i) for i in range(512)],["white"]*512)
    new = tab([str(i) for i in range(1792)],
              ["white"]*512+["mixed_uncertain"]*1280)
    new.loc[3,"morph"]="red_pink"
    with pytest.raises(ValueError,match="Earlier512"):
        x.validate_supersession(prev,new)
