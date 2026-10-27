"""Synthetic 128 source-shard selection: zero fake photos, authenticated priority and no P5 targeted."""
from __future__ import annotations
import collections,hashlib,unittest
from scripts.acquisition.mois_bulk128_tier10_source import (
    choose_part,take_priority,rank,add_original_photo,photo_worker,
    NO_COLOUR_PRIORITY,THIN_PRIORITY,P5,GAP_ARM,REFERENCE_ARM,
    TARGET_PER_SHARD,REFERENCE_PER_SHARD,SHARDS,TIER,
    source_five_degree_cell,
)
from scripts.acquisition.mois_onepass_full93m_target_tiers import obs_draw

class Bulk128PhotoSourceContractTests(unittest.TestCase):
    def make_population(self,*,p2=1100,p3=300,p4=150,p5=400):
        observations={};priorities={}
        categories=[(NO_COLOUR_PRIORITY[1],p2),(THIN_PRIORITY[0],p3),
                    (THIN_PRIORITY[1],p4),(P5,p5)]
        num=0
        z=("N_temperate","N_tropical","S_tropical","S_temperate")
        t=("summer","winter","autumn","spring")
        for label,n in categories:
            for _ in range(n):
                num+=1;taxon=str(10000+num)
                uuid="bulk-obs-"+str(num)
                region=z[num%len(z)];season=t[num%len(t)]
                observations[uuid]={
                    "uuid":uuid,"taxon":taxon,"photo_id":str(9000000+num),
                    "ref_rank":rank("uniform_reference_observation",uuid),
                    "target_rank":rank("independent_target_observation",uuid),
                    "zone":region,"season":season,
                }
                priorities[taxon]=label
        return observations,priorities

    def test_shard_selects_1024_independent_photos_without_P5_targeted(self):
        obs,prior=self.make_population()
        selected,receipt=choose_part(obs,prior,3)
        self.assertEqual(len(selected),TARGET_PER_SHARD)
        self.assertEqual(receipt["selected_reference_photo_attempts"],256)
        self.assertEqual(receipt["selected_targeted_photo_attempts"],768)
        self.assertEqual(receipt["targeted_unfilled_from_unrepresented_source_support"],0)
        self.assertEqual(len({x["uuid"] for _,_,x in selected}),1024)
        self.assertEqual(len({x["photo_id"] for _,_,x in selected}),1024)
        self.assertEqual(sum(role=="targeted" and label in NO_COLOUR_PRIORITY for role,label,_ in selected),
                         500)  # ceil(0.65 * 768)
        self.assertEqual(sum(role=="targeted" and label==P5 for role,label,_ in selected),0)

    def test_unavailable_source_photos_do_not_become_1024_fake_targets(self):
        obs,prior=self.make_population(p2=70,p3=80,p4=10,p5=500)
        selected,receipt=choose_part(obs,prior,17)
        self.assertLess(len(selected),1024)
        self.assertGreater(receipt["targeted_unfilled_from_unrepresented_source_support"],0)
        self.assertEqual(receipt["selected_reference_photo_attempts"],256)
        self.assertEqual(sum(role=="targeted" and label==P5 for role,label,_ in selected),0)

    def test_actual_source_taxon_modulo_and_original_observation_photo_multiplicity(self):
        table={}
        taxon="257"  # part 1 because 257 mod 128 == 1
        uuid="bulk128-one-original-observation"
        for photo in ["7001","7002","7003"]:
            row={
               "photo_id":photo,"observation_uuid":uuid,"inat_taxon_id":taxon,
               "tier":str(TIER),"observation_inclusion_score":format(obs_draw(uuid),".17g"),
               "observed_on":"2025-08-20","latitude":"-18","longitude":"130",
               "observer_id":"o1","positional_accuracy":"9",
               "photo_license":"cc-by","image_extension":"jpg",
            }
            add_original_photo(table,row,part=1)
        self.assertEqual(len(table),1)
        item=table[uuid]
        self.assertEqual(item["nphotos"],3)
        self.assertEqual(item["zone"],"S_tropical")
        self.assertEqual(item["season"],"winter")
        item["species"]="Source species"
        worker=photo_worker(item,"targeted",1,4,100)
        self.assertEqual(worker["source_arm"],GAP_ARM)
        self.assertTrue(0<float(worker["selection_probability"])<1)
        self.assertEqual(worker["inat_taxon_id"],taxon)
        ref=photo_worker(item,"reference",1,4,100)
        self.assertEqual(ref["source_arm"],REFERENCE_ARM)
        self.assertTrue(0<float(ref["selection_probability"])<=1)
        row["inat_taxon_id"]="258"
        with self.assertRaisesRegex(ValueError,"modulo128"):
            add_original_photo({},row,part=1)

    def test_unverified_geographic_coordinate_does_not_invent_5degree_cell(self):
        self.assertEqual(source_five_degree_cell("-18.3","130.2"),"14:62")
        self.assertEqual(source_five_degree_cell("", "120"),"unknown")
        self.assertEqual(source_five_degree_cell("44", "not-coordinate"),"unknown")
        self.assertEqual(source_five_degree_cell("95", "140"),"unknown")
        self.assertEqual(source_five_degree_cell("10", "inf"),"unknown")

    def test_priority_caps_do_not_reset_across_multiple_passes(self):
        cand=[]
        for i in range(35):
            cand.append({"uuid":f"obs{i}","taxon":"1001","zone":"S_tropical",
                         "season":"summer","target_rank":rank("t",i)})
        used=set();counts=collections.Counter();served=collections.Counter()
        chosen1=take_priority(cand,used,4,served,{("S_tropical","summer"):2.0},8,counts)
        chosen2=take_priority(cand,used,25,served,{("S_tropical","summer"):2.0},8,counts)
        self.assertEqual(len(chosen1),4)
        self.assertEqual(len(chosen2),4)
        self.assertEqual(counts["1001"],8)
        self.assertEqual(len(used),8)

if __name__=="__main__":
    unittest.main()
