"""Global MOIS source selection rewards breadth rather than FCP taxon depth."""
from datetime import date
from scripts.acquisition import mois_global_angiosperm_breadth_source as m


def o(pid,tid,name,observer,lat=34.5,lon=140.3,grid=""):
    return {
        "photo_id":str(pid),"taxon_id":str(tid),
        "species":name,"observer_id":str(observer),
        "observation_uuid":"uuid"+str(pid),
        "latitude":str(lat),"longitude":str(lon),
        "observed_on":"2025-05-04","positional_accuracy":"8",
        "photo_license":"cc-by",
        "photo_url_large":
            "https://inaturalist-open-data.s3.amazonaws.com/photos/"+
            str(pid)+"/large.jpg",
        "grid10":grid or m.grid10(lat,lon),
        "quality_grade":"research",
    }


def test_flowering_plant_taxonomy_must_be_active_species_not_generic_plant():
    assert m.ANGIOSPERMAE=="47125"
    assert m.new_eligible_taxon("123","1/47125/111","species",
                               "Cirsium japonicum","t")
    assert not m.new_eligible_taxon("123","1/47126/111","species",
                                   "Picea jezoensis","t")
    assert not m.new_eligible_taxon("123","1/47125/111","genus",
                                   "Cirsium","t")
    assert not m.new_eligible_taxon("123","1/47125/111","species",
                                   "Cirsium japonicum","f")


def test_research_observation_validation_keeps_original_day_and_geo():
    species={"22":("Campanula punctata","Campanula")}
    one=m.eligible_observation("22","research","57","unique-uuid",
        "2025-07-18","43.1","141.3","500",species)
    assert one["taxon_id"]=="22"
    assert one["quality_grade"]=="research"
    assert one["observed_on"]=="2025-07-18"
    assert m.eligible_observation("22","needs_id","57","uuid","2025-07-18",
                                  "43.1","141.3","500",species) is None
    assert m.eligible_observation("22","research","57","uuid","2025-07-18",
                                  "43.1","141.3","50000",species) is None
    assert m.eligible_observation("22","research","57","uuid","2025-07-18",
                                  "93.1","141.3","",species) is None


def test_new_genus_gets_priority_over_many_photos_in_same_old_species():
    taxa={}
    rows=[]
    for k in range(15):
        tid=str(100+k)
        gen=f"Novel{k}"
        taxa[tid]=(gen+" newflora",gen)
        rows.append(o(1000+k,tid,gen+" newflora",10000+k))
    taxa["30"]=("Campanula existing","Campanula")
    for k in range(30):
        rows.append(o(2000+k,"30","Campanula existing",30000+k))
    selected,receipt=m.rank_by_breadth(
        rows,taxa,old_taxa={"30"},old_canon={"campanula existing"},
        old_genera={"Campanula"},old_cells={},limit=16)
    assert len(selected)==16
    assert len({x["taxon_id"] for x in selected})==16
    assert set(str(i) for i in range(100,115)).issubset(
        {x["taxon_id"] for x in selected})
    assert receipt.get("new_genus_new_species")==15
    assert receipt.get("existing_species_new_region")==1


def test_no_hue_or_climate_used_in_hash_rank():
    items=[o(i,i,f"Flowering{i} newflora",i) for i in range(10,24)]
    taxa={str(i):(f"Flowering{i} newflora",f"Flowering{i}") for i in range(10,24)}
    a,_=m.rank_by_breadth(items,taxa,set(),set(),set(),{},limit=9)
    for x in items:
        x["morph"]="white" if int(x["photo_id"])%2 else "blue_purple"
        x["temperature"]=55
    b,_=m.rank_by_breadth(items,taxa,set(),set(),set(),{},limit=9)
    assert [x["photo_id"] for x in a]==[x["photo_id"] for x in b]
    assert len({x["photo_id"] for x in b})==9


def test_provider_photo_licence_and_minimum_size_constants():
    assert "cc-by" in m.PHOTO_LICENSES
    assert "cc0" in m.PHOTO_LICENSES
    assert m.PHOTO_EXTENSIONS=={"jpg","jpeg","png"}
    assert (m.POOL_GOAL,m.FIRST_WAVE,m.PER_SHARD,m.N_SHARDS)==(
        32768,4096,256,16)
