#!/usr/bin/env python3
"""Pre-outcome validity gates for MOIS five-hypothesis *future* photo study.

DO NOT call historical original pixels confirmatory on the basis of a new
protocol. Require a signed complete 93m source photo/observation index,
outcome-blind source queue, complete all-repository prior photo/observation
exposure roster and a disjoint whole-observation holdout BEFORE colours exist.
This file never reads or predicts flower hues and never fetches photos.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

FRAME_PHOTOS=93209351
FRAME_TAXA=42111
ALREADY_VIEWED_MOIS_PHOTOS=350399
SOURCE_SCHEMA="mois_full_42111_93209351_uncapped_original_photo_observation_index_v1"
QUEUE_SCHEMA="mois_full_source_observation_sampling_cohort_v1"
SPEC_SCHEMA="mois_five_photo_hypotheses_holdout_protocol_v1"
REQUIRED_HYPOTHESES=(
    "H1_ALTITUDE","H2_DRY_SEASON","H3_NATIVE_ISLAND","H4_ARCTIC",
    "H5_ATMOSPHERIC_MOISTURE",
)


def require(value, message):
    if not value:
        raise ValueError(message)


def sha256(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as file:
        for block in iter(lambda:file.read(1<<20),b""):
            h.update(block)
    return h.hexdigest()


def load(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def validate_protocol(protocol):
    p=load(protocol) if isinstance(protocol,(str,Path)) else protocol
    require(p.get("schema")==SPEC_SCHEMA,"Wrong protocol schema")
    require("NOT_retrospective_preregistration" in p.get("status",""),
            "This cannot be called a retrospective pre-registration")
    prior=p["original_images_currently_seen"]
    require(prior["source_photos_attempted"]==ALREADY_VIEWED_MOIS_PHOTOS
            and prior["prior_outcomes_exposed"] is True,
            "Existing photographed outcomes must be acknowledged as exposed")
    ids=[h["id"] for h in p["directional_hypotheses"]]
    require(ids==list(REQUIRED_HYPOTHESES),"Unexpected switched hypothesis IDs")
    require(all(h["previously_explored"] is True and
                h["white_direction"] in ("negative","positive")
                for h in p["directional_hypotheses"]),
            "No pre-existing hypothesis may be labelled independently prospective")
    holdout=p["holdout"]
    require(holdout["previously_seen_photo_ids_forbidden"] and
            holdout["previously_seen_observation_ids_forbidden"] and
            holdout["all_photos_from_selected_observation_retained"] and
            holdout["no_per_taxon_or_per_genus_photo_cap"] and
            holdout["no_result_disclosure_until_frozen_holdout_manifest_and_covariate_source_digests_exist"],
            "Prior-photo and source observation disclosure guard disabled")
    require(p["photo_universe"]["photo_candidates"]==FRAME_PHOTOS and
            p["photo_universe"]["source_taxa"]==FRAME_TAXA and
            holdout["proposed_initial_probability"]==0.0002,
            "Source frame/observation-level baseline probability was modified")
    require(p["external_source_gates"]["WCVP"]["same_island_native_assertion_permitted"] is False
            and p["external_source_gates"]["bee_richness"]["visitation_causality_claim_permitted"] is False
            and p["inference"]["output_if_unseen_holdout_unavailable"].startswith("EXPLORATORY_ONLY"),
            "Missing native or bee information cannot be promoted to causal inference")
    return p


def csv_ids(path, required_fields, *, unique_photo=True):
    fields=None
    rows=[]
    with Path(path).open(newline="",encoding="utf-8") as f:
        reader=csv.DictReader(f)
        fields=set(reader.fieldnames or [])
        require(set(required_fields)<=fields,"Missing mandatory observation/photo ID fields")
        forbidden={"white","morph","hue","classified","measurement_status",
                   "flower_fraction_white","predicted_hue"}
        # Blind queue must contain identities only, never outcome labels.
        if unique_photo: require(not (fields & forbidden),
                                 "Holdout roster has photographic colour or ROI outcomes")
        seen=set()
        for row in reader:
            pid=str(row["photo_id"]).strip()
            oid=str(row["observation_id"]).strip()
            require(pid.isdecimal() and oid,"Invalid original photo/observation ID")
            if unique_photo:
                require(pid not in seen,"Holdout photo ID repeated")
            seen.add(pid)
            rows.append((pid,oid,row))
    require(rows,"Empty signed source observation/photo roster")
    return rows


def gate(protocol, full_index_manifest, cohort_manifest,
         holdout_roster_csv, previously_viewed_csv, prior_receipt,
         *, expected_photos=FRAME_PHOTOS, expected_taxa=FRAME_TAXA,
         expected_shards=128, min_prior_attempts=ALREADY_VIEWED_MOIS_PHOTOS):
    """Return a *pre-outcome eligibility audit*, not a statistical confirmation."""
    p=validate_protocol(protocol)
    index=load(full_index_manifest); queue=load(cohort_manifest)
    prior=load(prior_receipt)
    require(index.get("schema")==SOURCE_SCHEMA and
            index.get("photos_emitted")==expected_photos and
            index.get("taxa_in_frame")==expected_taxa and
            index.get("photo_caps") is None and
            index.get("per_observation_photo_cap") is None and
            index.get("photo_colours_measured") is False and
            len(index.get("shards",[]))==expected_shards,
            "Full verified original photo-and-observation index has not finished")
    require(sum(s["photo_count"] for s in index["shards"])==expected_photos
            and all("photos_sha256" in s and "observations_sha256" in s
                    and "observation_count" in s for s in index["shards"]),
            "Source photo/observation shard audit incomplete")
    require(queue.get("schema")==QUEUE_SCHEMA and
            queue.get("full_frame_photo_ids")==expected_photos and
            queue.get("full_source_manifest_sha256")==sha256(full_index_manifest) and
            queue.get("selection_unit")=="observation_uuid" and
            queue.get("selection_probability")==p["holdout"]["proposed_initial_probability"] and
            queue.get("sampled_before_any_hue_measurement") is True and
            queue.get("colour_measured") is False and
            queue.get("source_original_images_measured") is False and
            queue.get("all_unselected_source_photos_eligible_at_higher_probability") is True and
            queue.get("no_per_taxon_or_observation_photo_count_cap") is True,
            "Not a frozen complete-source pre-colour observation cohort")
    require(prior.get("schema")=="mois_prior_exposed_photo_and_observation_roster_v1"
            and prior.get("includes_all_mois_original_attempts") is True
            and prior.get("includes_other_repositories_and_failed_original_attempts") is True
            and prior.get("photo_id_observation_id_csv_sha256")==sha256(previously_viewed_csv),
            "Prior exposed pixels in MOIS and other repositories are not source-attested")
    old=csv_ids(previously_viewed_csv,("photo_id","observation_id"),unique_photo=False)
    require(len(old)>=min_prior_attempts and
            len({a for a,_,_ in old})==len(old),
            "Prior full original-photo and observation exclusion roster incomplete/duplicated")
    previously_photos={a for a,_,_ in old}
    previously_obs={b for _,b,_ in old}
    future=csv_ids(holdout_roster_csv,("photo_id","observation_id"),unique_photo=True)
    future_ids={a for a,_,_ in future}
    future_obs={b for _,b,_ in future}
    require(not (future_ids & previously_photos),
            "Previously seen original photo IDs leak into the confirmatory holdout")
    require(not (future_obs & previously_obs),
            "A different PHOTO of a previously seen OBSERVATION leaked into the holdout")
    # This is exactly the queue's frozen source image roster, no early stopping.
    require(len(future)==queue["new_photo_attempts_planned"]
            and len(future_obs)==queue["selected_observations_with_new_photos"]
            and not queue.get("ecological_analysis_run",False),
            "Future source roster/number of observation units is not complete")
    result={
        "schema":"mois_five_hypotheses_future_photo_holdout_preoutcome_gate_v1",
        "status":"OUTCOME_BLIND_PHOTO_ROSTER_ELIGIBLE_NOT_ECOLOGICAL_CONFIRMATION",
        "protocol_sha256":sha256(protocol),
        "full_index_manifest_sha256":sha256(full_index_manifest),
        "queue_manifest_sha256":sha256(cohort_manifest),
        "previously_exposed_source_csv_sha256":sha256(previously_viewed_csv),
        "future_unseen_source_csv_sha256":sha256(holdout_roster_csv),
        "future_photos":len(future),
        "future_distinct_observations":len(future_obs),
        "previously_exposed_photo_ids":len(old),
        "previously_exposed_observation_ids":len(previously_obs),
        "independent_source_observation_and_photo_IDs_verified":True,
        "all_photo_pixels_in_the_future_roster_have_colours_measured":False,
        "native_UVB_bee_and_botanical_review_external_data_ready":False,
        "hypothesis_confirmatory_inference_permitted":False,
        "reason":"Independent original photo outcomes, vetted open flower organs, source-native ranges, UV and independent bee predictors have not yet been obtained and verified."
    }
    return result


def main():
    x=argparse.ArgumentParser(description=__doc__)
    for key in ("protocol","full-index-manifest","cohort-manifest",
                "holdout-roster-csv","previously-viewed-csv","prior-receipt","output"):
        x.add_argument("--"+key,type=Path,required=True)
    a=x.parse_args()
    result=gate(a.protocol,a.full_index_manifest,a.cohort_manifest,
                a.holdout_roster_csv,a.previously_viewed_csv,a.prior_receipt)
    require(not a.output.exists(),"Do not overwrite outcome-blind protocol audit")
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    print("MOIS_FIVE_HYPOTHESIS_HOLDOUT_PRECOLOUR_SOURCE_VALID",
          json.dumps(result,sort_keys=True),flush=True)


if __name__=="__main__":
    main()
