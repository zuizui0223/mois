#!/usr/bin/env python3
"""MOIS epoch9 offset1408: new observation-created interval after complete independently audited epoch8 (2,197 taxa).

This is a new temporal source interval, not a replay of any epoch6 source SHA.
The source selects the next 128 taxon IDs from the 2,197 still-zero-strict9 taxa in a
frozen epoch7 hash order; previously queried taxa are eligible only for
observations created after the epoch6 final source receipt. All older photo
IDs/observation UUIDs are excluded. Raw image SHA disjointness against the
full historical ~1.23 million images is NOT certified.
"""
import argparse
import hashlib
import json
from pathlib import Path
import pandas as pd

from scripts.acquisition import (
    mois_epoch6_post_complete2199_new_created_original_photo9_offset2176 as predecessor
)
from scripts.acquisition import mois_20261115_created_epoch2_disjoint_original_photo9_after2048 as core

from scripts.acquisition import mois_epoch7_after_epoch6_new_created_photo9_offset2048 as prior_2048
m = prior_2048.m
m.SCHEMA = "mois_incremental_20261119_epoch9_offset1408_post_complete_epoch8_1231356_original9_v1"
m.SOURCE_ARM = "MOIS_EPOCH9_OFFSET1408_NEW_CREATED_AFTER_EPOCH8_FINAL2197_TERMINAL"
# Epoch6 final source artifact was frozen at 2026-11-17T11:47:05Z.
m.MIN_CREATED = "2026-11-18T19:07:09+00:00"
m.CREATED_D1 = "2026-11-18"
m.MAX_SOURCE_TAXA = 128
m.ATLAS_SHA = "bc646ee29c02fa02111e073f096b2c7a63594d7294d7e658c2c7d6e8424bc154"
m.OLD_PHOTOS = 1231356
m.OLD_NUMERIC_PHOTOS = 467568
m.OLD_NUMERIC_TAXA = 39914
m.OLD_CLEAR_TAXA = 39174

def select_taxa_epoch8(atlas, n=128, offset=0):
    m.require(len(atlas) == 42111 and atlas.inat_taxon_id.is_unique,
              "nominal atlas must contain exactly 42,111 distinct taxa")
    z = atlas.loc[atlas.after_strict9.eq(0)].copy()
    m.require(len(z) == 2197 and z.inat_taxon_id.str.fullmatch(r"\d+").all(),
              "epoch7 must select only the 2,197 presently unmeasured taxa")
    z["priority"] = z.inat_taxon_id.map(
        lambda tid: hashlib.sha256(
            ("MOIS-20261118-EPOCH8-AFTER-EPOCH7-2197-TEMPORAL-NEW-CREATED|" + tid).encode()
        ).hexdigest())
    ordered = z.sort_values("priority", kind="stable").inat_taxon_id.tolist()
    selected = ordered[offset:offset+n]
    m.require(offset == 1408 and n == 128
              and len(selected) == len(set(selected)) == n
              and not set(selected).intersection(ordered[:offset]),
              "epoch9 offset1408 frozen taxon selection changed")
    return selected

# Frozen source-only 2048 remains immutable: its only independent audit is the
# read-only recovery run 44234354022 (byte verified, 0 candidates).
# Exactly the same 1,231,355 IDs and 111 known recent raw hashes carry forward.
# Prior epoch8 offsets0–1664 are sealed zero-source; its photo/UUID IDs do not change history.
def history_after_epoch8_offset1920():
    # Start from the independently signed and immutable old 1,231,355 ID and 111
    # recent-raw-SHA ledger. Add only the audit-proven new technical-unknown image.
    path = Path("old/epoch8_offset1920/new_real_original_nine_photo9_delta.parquet")
    expected_sha = "5666543e87beee8ab4e63df5b3107fbea5b01b49b4ee5982c71fd33210bc5745"
    m.require(path.is_file() and m.sha256(path) == expected_sha,
              "previous epoch8 offset1920 independently audited delta absent/drifted")
    current=m.OLD_PHOTOS
    try:
        m.OLD_PHOTOS=1231355
        ids, observations, raw = prior_2048.history_after_epoch7_offset1920()
    finally:
        m.OLD_PHOTOS=current
    m.require(len(ids)==1231355 and len(raw)==111,
              "old 1,231,355 original-ID and 111 raw-SHA history drift")
    d=pd.read_parquet(path)
    m.require(len(d)==1 and {"photo_id","observation_id","image_sha256",
                               "strict9","inat_taxon_id"}.issubset(d.columns),
              "old independently audited delta missing required data")
    row=d.iloc[0]
    pid,oid,h,tid=(str(row.photo_id),str(row.observation_id),
                   str(row.image_sha256),str(row.inat_taxon_id))
    m.require(pid=="764874941"
              and oid=="a5b58834-4657-4f90-95c0-9f4d0292ceb6"
              and h=="7cbe598b5e4b1e2cbc2ad2a4f653cecd225ee27c171fb5b42b69632df557d8bf"
              and tid=="774116" and not bool(row.strict9)
              and pid not in ids and oid not in observations and h not in raw,
              "new raw SHA, photo ID, observation UUID, strict9 or taxon drift")
    m.require(all(pd.isna(row["strict9_"+p]) for p in m.PALETTE),
              "previous technical unknown must have all nine palette components missing")
    ids.add(pid);observations.add(oid);raw.add(h)
    m.require(len(ids)==1231356 and len(raw)==112,
              "complete 1,231,356 original-ID / 112 raw SHA history mismatch")
    return ids,observations,raw

m.history = history_after_epoch8_offset1920

m.select_taxa = select_taxa_epoch8
m.source = core.source

def main():
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest="stage", required=True)
    for stage, keys in (
        ("source", ("atlas", "output")),
        ("audit", ("source", "worker", "output")),
        ("append", ("atlas", "source", "worker", "audit", "output")),
    ):
        q = sub.add_parser(stage)
        for key in keys:
            q.add_argument("--" + key, required=True)
        if stage == "source":
            q.add_argument("--offset", type=int, default=1408)
    args = parser.parse_args()
    if args.stage == "source":
        m.source(args.atlas, args.output, quota=128, offset=args.offset)
        receipt_path = Path(args.output) / "receipt.json"
        r = json.loads(receipt_path.read_bytes())
        m.require(
            (r["hash_disjoint_taxon_offset"], r["searched_zero_numeric_taxa"],
             r["source_API_error_count"], r["old_original_photo_attempt_ID_count"],
             r["old_original_photo_obs_ID_overlap"]) == (1408, 128, 0, 1231356, 0)
            and r["old_all_image_raw_SHA_disjointness_proven"] is False
            and r["created_utc_lower_bound"] == m.MIN_CREATED
            and r["original_pixels_fetched"] == r["strict9_numeric_photos"] == 0,
            "epoch9 offset1408 new original source-only receipt gate failed")
        r["recent_known_original_image_raw_SHA_count"] = 112
        r["preceding_epoch8_final_source_run"] = 44375274312
        r["preceding_epoch8_final_source_artifact"] = 14460963430
        r["preceding_epoch8_final_source_zip_SHA256"] = "ee3348942a5f1c4455783c097dd28a666201347fdf15a4e6751742f6afa88931"
        r["preceding_epoch8_complete2197_recovery_run"] = 44376974882
        r["preceding_epoch8_complete2197_recovery_artifact"] = 14462440335
        r["preceding_epoch8_complete2197_recovery_zip_SHA256"] = "efe3668865c363026fc48c330f3a999e50c6ea6d80ba92eddfdbdb5a774c5dc2"
        r["preceding_epoch9_offset0_run"] = 44414467074
        r["preceding_epoch9_offset0_source_artifact"] = 14479865043
        r["preceding_epoch9_offset0_source_ZIP_SHA256"] = "8667ba3638dd7212d963e72e04eb73119ea4f426f29aa506668582882368bb5c"
        r["preceding_epoch9_offset0_new_candidates"] = 0
        r["preceding_epoch9_offset128_run"] = 44418086333
        r["preceding_epoch9_offset128_source_artifact"] = 14480986460
        r["preceding_epoch9_offset128_source_ZIP_SHA256"] = "517c7b4e26bfe099397029b94459efa11bf39b318d78453b4f7eb06c54da24bf"
        r["preceding_epoch9_offset128_source_receipt_SHA256"] = "415e3ee10b640cd8963f5218912d60ef6103b9df248e51a45904d2ac94ca9ad9"
        r["preceding_epoch9_offset128_new_candidates"] = 0
        r["preceding_epoch9_offset256_run"] = 44418896209
        r["preceding_epoch9_offset256_source_artifact"] = 14480333609
        r["preceding_epoch9_offset256_source_ZIP_SHA256"] = "a323884f12d322f5073c1f7b6846e4294f62143ee4a1bcd526f0354b196b2762"
        r["preceding_epoch9_offset256_source_receipt_SHA256"] = "01306b6aaa1b028e8a5a8a0b3c704d4c1277a72831dd1e170c6167ffee2c2e2d"
        r["preceding_epoch9_offset256_new_candidates"] = 0
        r["preceding_epoch9_offset384_run"] = 44427759101
        r["preceding_epoch9_offset384_source_artifact"] = 14484529215
        r["preceding_epoch9_offset384_source_ZIP_SHA256"] = "e6799199dc5d156740563d7c97d6d0e6bf4456856760957f7b4b379ef2b557d2"
        r["preceding_epoch9_offset384_source_receipt_SHA256"] = "4ff7e5d7bcc2731461408c3e0a9d5da52415467eabfd95aef7dbf519781bb1ed"
        r["preceding_epoch9_offset384_new_candidates"] = 0
        r["preceding_epoch9_offset512_run"] = 44428331270
        r["preceding_epoch9_offset512_source_artifact"] = 14484559749
        r["preceding_epoch9_offset512_source_ZIP_SHA256"] = "f4f0e94fb4e729f27ab851eebfd04d76642395161003d067b2e4cac776b2bb98"
        r["preceding_epoch9_offset512_source_receipt_SHA256"] = "b02a8387c5e470a609aec2d34a252b447467d8a3e5d3e382f6460520a67f6dd3"
        r["preceding_epoch9_offset512_new_candidates"] = 0
        r["preceding_epoch9_offset640_run"] = 44438815760
        r["preceding_epoch9_offset640_source_artifact"] = 14490710603
        r["preceding_epoch9_offset640_source_ZIP_SHA256"] = "2bc8a0730c232cd308043dfd3baa00338ef0890148e13207ea643be226f4fa1c"
        r["preceding_epoch9_offset640_source_receipt_SHA256"] = "81cbdceedfe83a63b41b03d3ae97faa71bf6078583ddba91ee19103b8ac6c1a0"
        r["preceding_epoch9_offset640_source_status_SHA256"] = "5dcd8c99ce227b1bfedb5d0c137e16738772ecc142c826cdecb9f7ecc64b2ea9"
        r["preceding_epoch9_offset640_new_candidates"] = 0
        r["preceding_epoch9_offset768_run"] = 44444368945
        r["preceding_epoch9_offset768_source_artifact"] = 14492723162
        r["preceding_epoch9_offset768_source_ZIP_SHA256"] = "191e89608ec0f6331405860ad64f63961fea642998f2f7b913628a58f66c8092"
        r["preceding_epoch9_offset768_source_receipt_SHA256"] = "5c458002edb6bde2e16e8fe7e12426c7a03f56ec8e089234bd8a742bdd483624"
        r["preceding_epoch9_offset768_source_status_SHA256"] = "8e87be07cadf709b7ca556c1a689d1264b2190897811bd086d8a7004f30e93ec"
        r["preceding_epoch9_offset768_new_candidates"] = 0
        r["preceding_epoch9_offset896_run"] = 44453284367
        r["preceding_epoch9_offset896_source_artifact"] = 14496856997
        r["preceding_epoch9_offset896_source_ZIP_SHA256"] = "7fcb92a02e7ab26486634a8e29ee6efac0b1a95423631b2e991cf68875093bba"
        r["preceding_epoch9_offset896_source_receipt_SHA256"] = "da72be008899ea5a6acd8672c13bb6049e054e0e8a0f830944ca80bffa8011b0"
        r["preceding_epoch9_offset896_source_status_SHA256"] = "fb64a5c35410ba3e63dc42b4455ea830e16ef28442d379c1d5d92b0a17164e08"
        r["preceding_epoch9_offset896_new_candidates"] = 0
        r["preceding_epoch9_offset1024_run"] = 44453805964
        r["preceding_epoch9_offset1024_source_artifact"] = 14496442990
        r["preceding_epoch9_offset1024_source_ZIP_SHA256"] = "10df3eb3b92988d466bbc43e3d39c9adeaa27b2f93381ce81ce653842299a4f7"
        r["preceding_epoch9_offset1024_source_receipt_SHA256"] = "89c8e4a501e5f144eb089a916e4372990dcffd0075bb935e4911f74598d81a1f"
        r["preceding_epoch9_offset1024_source_status_SHA256"] = "07040ef9fa4291753ae856ab3451a8065d423e20e27359296bf51209b18c91f5"
        r["preceding_epoch9_offset1024_new_candidates"] = 0
        r["preceding_epoch9_offset1152_run"] = 44459892412
        r["preceding_epoch9_offset1152_source_artifact"] = 14499662068
        r["preceding_epoch9_offset1152_source_ZIP_SHA256"] = "9dac2ab05fa42a0f10f709a4e260dd148555e6aec89b879d74fe65de4d3f7bac"
        r["preceding_epoch9_offset1152_source_receipt_SHA256"] = "0943917a4eba2ba745cea889049afb29e3ec474818632f6714cae7b08b09ebd4"
        r["preceding_epoch9_offset1152_source_status_SHA256"] = "289ec1f182885e5c30a311dbacc2c375a26aa9005ba179dbb92a2cd9f6bdd7f6"
        r["preceding_epoch9_offset1152_new_candidates"] = 0
        r["preceding_epoch9_offset1280_run"] = 44460557777
        r["preceding_epoch9_offset1280_source_artifact"] = 14499827914
        r["preceding_epoch9_offset1280_source_ZIP_SHA256"] = "48a3471f04a81cf6ac357cda6a42187741cd4c80349f2f3ea467b0b4181c9121"
        r["preceding_epoch9_offset1280_source_receipt_SHA256"] = "31f793f754f73972dee4abf8983bca9c8b526a388d1c649dc396d4ddfa7a539d"
        r["preceding_epoch9_offset1280_source_status_SHA256"] = "f9356dbf61aa9d26c574222f8eca1f4e302adfe6a517cce9e26a131b64b4c69d"
        r["preceding_epoch9_offset1280_new_candidates"] = 0
        r["new_temporal_interval_after_epoch8_final_terminal"] = True
        r["prior_signed_epoch6_taxon_frame_complete"] = 2199
        r["currently_zero_numeric_taxa"] = 2197
        r["previous_epoch6_terminal_run"] = 44076781021
        r["previous_epoch7_offset0_run"] = 44084420377
        r["previous_epoch7_offset0_source_ZIP_SHA256"] = "8eaee465a4d4b31c318a0abd26303cebcedaa23345fa32d57c16cbdd84ac28d2"
        r["previous_epoch7_offset0_source_candidate_count"] = 0
        r["previous_epoch7_offset128_run"] = 44088047232
        r["previous_epoch7_offset128_source_ZIP_SHA256"] = "a1fdd2b04e03a7d8455c7bd09c0959016a75427bd1c6d2a6c2c35dc63d5ddf46"
        r["previous_epoch7_offset128_source_candidate_count"] = 0
        r["previous_epoch7_offset256_run"] = 44105932217
        r["previous_epoch7_offset256_source_ZIP_SHA256"] = "7b9a159cc772dbd1b859139ff498cc900947971602a855617ace3d6322b0bf37"
        r["previous_epoch7_offset256_source_candidate_count"] = 0
        r["previous_epoch7_offset384_run"] = 44108050010
        r["previous_epoch7_offset384_source_ZIP_SHA256"] = "9b181cbed5ba5d33e5f5accc1036751c9aa120060319331601d7d96650b4788c"
        r["previous_epoch7_offset384_source_candidate_count"] = 0
        r["previous_epoch7_offset512_run"] = 44125684318
        r["previous_epoch7_offset512_source_ZIP_SHA256"] = "dc1f2e0171661170ed8b8d393b279169d8843f06f8def56778f436327dbaa9aa"
        r["previous_epoch7_offset512_source_candidate_count"] = 0
        r["previous_epoch7_offset640_run"] = 44138261989
        r["previous_epoch7_offset640_source_ZIP_SHA256"] = "02c9c1b60c375f1450fe4f299ce9c793f20924f939d3dd2d27b1d5334a3c5153"
        r["previous_epoch7_offset640_source_candidate_count"] = 0
        r["previous_epoch7_offset768_run"] = 44140649923
        r["previous_epoch7_offset768_source_ZIP_SHA256"] = "d797a36f97b6f26db6ebfeeb33ca4e8e9c7444ec5157580730af3a8be238ccd7"
        r["previous_epoch7_offset768_original_photo_ID"] = "764609465"
        r["previous_epoch7_offset768_audited_delta_SHA256"] = "2f5d7d4b905e071caf1198ecdbac9835b3793c0e0d37f2fc73aa5d41ab8229b8"
        r["previous_epoch7_offset768_metadata_bug_corrected_this_receipt_only"] = True
        r["previous_epoch7_offset896_run"] = 44159159412
        r["previous_epoch7_offset896_source_artifact"] = 14363370860
        r["previous_epoch7_offset896_source_ZIP_SHA256"] = "c2b2d53d678ff6590abaa7cc5d5bdabbb72d5a3829b2212a6439d02c6df1647a"
        r["previous_epoch7_offset896_source_candidate_count"] = 0
        r["previous_epoch7_offset1024_run"] = 44170221242
        r["previous_epoch7_offset1024_source_artifact"] = 14367587413
        r["previous_epoch7_offset1024_source_ZIP_SHA256"] = "0fcc9363f7598ec278e4eb8c5cbc621abac9d617b79591fd0f1f46135539d7a4"
        r["previous_epoch7_offset1024_source_candidate_count"] = 0
        r["previous_epoch7_offset1152_run"] = 44173253159
        r["previous_epoch7_offset1152_source_artifact"] = 14368388310
        r["previous_epoch7_offset1152_source_ZIP_SHA256"] = "e00fc9d0215f0fdb172b6fa7731db2a68307b7d4274b7e667c72a84bc2605520"
        r["previous_epoch7_offset1152_source_candidate_count"] = 0
        r["previous_epoch7_offset1280_run"] = 44192320978
        r["previous_epoch7_offset1280_source_artifact"] = 14376667948
        r["previous_epoch7_offset1280_source_ZIP_SHA256"] = "d28b743dbff4832fd85d6c2cb4b950ddbd88a23f962615643fecd56c46913664"
        r["previous_epoch7_offset1280_source_candidate_count"] = 0
        r["previous_epoch7_offset1408_run"] = 44201443264
        r["previous_epoch7_offset1408_source_artifact"] = 14381214117
        r["previous_epoch7_offset1408_source_ZIP_SHA256"] = "cdcba19d017db4420bddf874b7ebd5ee6e668474af0c5609f5b938013bb5f89b"
        r["previous_epoch7_offset1408_source_candidate_count"] = 0
        r["previous_epoch7_offset1536_run"] = 44208024637
        r["previous_epoch7_offset1536_source_artifact"] = 14385435693
        r["previous_epoch7_offset1536_source_ZIP_SHA256"] = "53e6f7e83be2b7f0ff68fd800c9d6675d695efad3ec4a4a8667cefd3ff1207fe"
        r["previous_epoch7_offset1536_audited_delta_SHA256"] = "c653d80677f10947a361ca61584450f68806392322cd9b3fcd13cc7dd00eada1"
        r["previous_epoch7_offset1536_photo_ID"] = "764630362"
        r["previous_epoch7_offset1664_run"] = 44215477957
        r["previous_epoch7_offset1664_source_artifact"] = 14387762866
        r["previous_epoch7_offset1664_source_ZIP_SHA256"] = "8675bd7989a38132c214d66980424d361cbb73646a15a9c0378a01b3688c0b5e"
        r["previous_epoch7_offset1664_receipt_SHA256"] = "684b182b0f4e5aec4a031ba9e165cb84d0bc3e713ac3cedbf29c2977120dd7f5"
        r["previous_epoch7_offset1664_source_candidate_count"] = 0
        r["previous_epoch7_offset1792_run"] = 44219184371
        r["previous_epoch7_offset1792_source_artifact"] = 14390790013
        r["previous_epoch7_offset1792_source_ZIP_SHA256"] = "6042d82ea0bca4894bd6b32ef9d48334c277ef1e73b8b06a527d42536181c6bb"
        r["previous_epoch7_offset1792_receipt_SHA256"] = "83614ad3a1c308fec0c13ab98aa4ae27437e8cc6c9aaff53977635dc1ceac76c"
        r["previous_epoch7_offset1792_source_candidate_count"] = 0
        r["previous_epoch7_offset1920_run"] = 44231043954
        r["previous_epoch7_offset1920_source_artifact"] = 14395566128
        r["previous_epoch7_offset1920_source_ZIP_SHA256"] = "3930f94aa5e4ab8c3578f3f6708053dddd932f5954a782d7f2251f825a49f511"
        r["previous_epoch7_offset1920_worker_artifact"] = 14394733765
        r["previous_epoch7_offset1920_worker_ZIP_SHA256"] = "4759f24bf28bee6210ce3f3d6717ae4f2118106ca68ac3448b93e8670b99a3b9"
        r["previous_epoch7_offset1920_independent_audit_artifact"] = 14395023981
        r["previous_epoch7_offset1920_independent_audit_ZIP_SHA256"] = "fb5d5c10d20b28c808cc172a92d687eb5c2d0dd7023d2286d9bf1cadfd7037e9"
        r["previous_epoch7_offset1920_formal_atlas_artifact"] = 14394934674
        r["previous_epoch7_offset1920_formal_atlas_ZIP_SHA256"] = "bafee36929746dba5b00f380e21a10ffa7785e1ff4bf4d882a91b8e6242967db"
        r["previous_epoch7_offset1920_audited_delta_SHA256"] = "278cf722d023644025c4a4b7745753e05739a5c1d3e93f225c9441161025012c"
        r["previous_epoch7_offset1920_photo_ID"] = "764700473"
        r["previous_epoch7_offset2048_source_run"] = 44233454465
        r["previous_epoch7_offset2048_source_artifact"] = 14395854062
        r["previous_epoch7_offset2048_source_ZIP_SHA256"] = "e88a282b4c05906399cbd59ffc02a3c7babb2b8091d79130e21ffe04ae537698"
        r["previous_epoch7_offset2048_source_receipt_SHA256"] = "e167b51ac82dda189fc846b5f4ff5412414a1dc5a6f0a35357cf60b3b1c3776d"
        r["previous_epoch7_offset2048_readonly_recovery_audit_run"] = 44234354022
        r["previous_epoch7_offset2048_readonly_recovery_audit_artifact"] = 14396826157
        r["previous_epoch7_offset2048_readonly_recovery_audit_ZIP_SHA256"] = "c3588a93f1d34acd7180b863a7e09adcbcc91c2a0ac1f059d9c3dca106fbd11c"
        r["previous_epoch7_offset2048_source_candidate_count"] = 0
        r["previous_epoch7_final_offset2176_source_run"] = 44235063258
        r["previous_epoch7_final_offset2176_source_artifact"] = 14396233558
        r["previous_epoch7_final_offset2176_source_ZIP_SHA256"] = "5fcd21233648d4a5983c709002548d5f34c861a3b6cb25c9bbabb6ad3d61dca9"
        r["previous_epoch7_final_offset2176_receipt_SHA256"] = "afd81c5c1225259c8ee0e98d1153e45d92f3dfcb2b0fd1acd1a0a15ce82eac03"
        r["previous_epoch7_final_offset2176_source_candidate_count"] = 0
        r["epoch8_temporal_after_epoch7_final_source"] = True
        r["previous_epoch8_offset0_terminal_run"] = 44241913276
        r["previous_epoch8_offset0_source_artifact"] = 14400002460
        r["previous_epoch8_offset0_source_ZIP_SHA256"] = "2c88f44781e6f53b0ef33a72bd353c6abbd8cb223e15fdd09c4a3fadc744fd9e"
        r["previous_epoch8_offset0_new_source_candidates"] = 0
        r["previous_epoch8_offset128_terminal_run"] = 44251161824
        r["previous_epoch8_offset128_source_artifact"] = 14403982782
        r["previous_epoch8_offset128_source_ZIP_SHA256"] = "9871cc7677400205f0c2e79ccf9fb49ba8cb38e8ae1b0f1979d10238ca239470"
        r["previous_epoch8_offset128_source_receipt_SHA256"] = "a79a5600ea02302f8622b41260affb81affeb54872c92b517595932b4bffacbc"
        r["previous_epoch8_offset128_source_status_SHA256"] = "bd738e1f468595d34b6e30783107227ccb4c1bba8098677d1ed897df6b1395a3"
        r["previous_epoch8_offset128_new_source_candidates"] = 0
        r["previous_epoch8_offset256_terminal_run"] = 44252212051
        r["previous_epoch8_offset256_source_artifact"] = 14404457828
        r["previous_epoch8_offset256_source_ZIP_SHA256"] = "4caa58596d9bd9be8e676801279af01c694d23d0cca25a0a652b6f560395d958"
        r["previous_epoch8_offset256_source_receipt_SHA256"] = "bc4f46b9e9abef2fbde711f861597a9d20f6e4ca19cd5d0912c237957350fd34"
        r["previous_epoch8_offset256_source_status_SHA256"] = "779e4b081f92928f75ff49b9794dc3c1101ece1bfa7a742b93609d435857fd41"
        r["previous_epoch8_offset256_new_source_candidates"] = 0
        r["previous_epoch8_offset384_terminal_run"] = 44259582461
        r["previous_epoch8_offset384_source_artifact"] = 14407887198
        r["previous_epoch8_offset384_source_ZIP_SHA256"] = "21b59d080bd98de4975dccbb3c54f12832594d3c329c6470eb04cc88d6433566"
        r["previous_epoch8_offset384_source_receipt_SHA256"] = "8d837ab0a5241498d4448c702698b04c8fe462c8a70202fda7bd8ca2a6511be9"
        r["previous_epoch8_offset384_source_status_SHA256"] = "100e62ec3fa0196c9007ea522a2bbfc7bef899cd26d8120caba61d999129f16c"
        r["previous_epoch8_offset384_new_source_candidates"] = 0
        r["previous_epoch8_offset512_terminal_run"] = 44265122905
        r["previous_epoch8_offset512_source_artifact"] = 14409618707
        r["previous_epoch8_offset512_source_ZIP_SHA256"] = "6421bc22193167f97b383091bf5550a879308224c2369738cbef94007fd2f53f"
        r["previous_epoch8_offset512_source_receipt_SHA256"] = "825565ab7cd2465150453604ca0562352898f7aa3ac09a3d64d9482fab608ad8"
        r["previous_epoch8_offset512_source_status_SHA256"] = "261b2742e00dc2b1e671a3b934cabf0a72800a35916865648dd9c99d0f3f88f0"
        r["previous_epoch8_offset512_new_source_candidates"] = 0
        r["previous_epoch8_offset640_terminal_run"] = 44278444204
        r["previous_epoch8_offset640_source_artifact"] = 14415992331
        r["previous_epoch8_offset640_source_ZIP_SHA256"] = "d0bc47d2b05caaf78f25717844bd9ea3763e94d201b923dd2511d1c3195b647b"
        r["previous_epoch8_offset640_source_receipt_SHA256"] = "18d89bf27645926a0f221ce05607c8088db82ab5d2eacb234d24a6308540bfc7"
        r["previous_epoch8_offset640_source_status_SHA256"] = "5dcd8c99ce227b1bfedb5d0c137e16738772ecc142c826cdecb9f7ecc64b2ea9"
        r["previous_epoch8_offset640_new_source_candidates"] = 0
        r["previous_epoch8_offset768_terminal_run"] = 44282241887
        r["previous_epoch8_offset768_source_artifact"] = 14418046010
        r["previous_epoch8_offset768_source_ZIP_SHA256"] = "f348f3dbc4504efce7d5208b69dfba058d6d0ed4aabe6f8a0b752809d7b56919"
        r["previous_epoch8_offset768_source_receipt_SHA256"] = "aaef87450b4cf035959c054fd994b6ac75cdea2f75f315b680afb7a0d68b58b3"
        r["previous_epoch8_offset768_source_status_SHA256"] = "8e87be07cadf709b7ca556c1a689d1264b2190897811bd086d8a7004f30e93ec"
        r["previous_epoch8_offset768_new_source_candidates"] = 0
        r["previous_epoch8_offset896_terminal_run"] = 44291835760
        r["previous_epoch8_offset896_source_artifact"] = 14420949834
        r["previous_epoch8_offset896_source_ZIP_SHA256"] = "eb15cac0bb89436d900f3f2cb9cc69d25cca8232509cfe6dba5224ca3df5ec1e"
        r["previous_epoch8_offset896_source_receipt_SHA256"] = "eeae4cccd72670fdb2d695e5197b61bf32244070604945de8aff5b72270ce726"
        r["previous_epoch8_offset896_source_status_SHA256"] = "fb64a5c35410ba3e63dc42b4455ea830e16ef28442d379c1d5d92b0a17164e08"
        r["previous_epoch8_offset896_new_source_candidates"] = 0
        r["previous_epoch8_offset1024_terminal_run"] = 44306614585
        r["previous_epoch8_offset1024_source_artifact"] = 14428593193
        r["previous_epoch8_offset1024_source_ZIP_SHA256"] = "86bab8dc138421be03470e2e55e80b8832d2f0c16b9f37003be86cc21763efd6"
        r["previous_epoch8_offset1024_source_receipt_SHA256"] = "77dc55c5fa198c551f4b7bac5ecbbc3912130612b4d44d8c126badb370483697"
        r["previous_epoch8_offset1024_source_status_SHA256"] = "07040ef9fa4291753ae856ab3451a8065d423e20e27359296bf51209b18c91f5"
        r["previous_epoch8_offset1024_new_source_candidates"] = 0
        r["previous_epoch8_offset1152_terminal_run"] = 44309830390
        r["previous_epoch8_offset1152_source_artifact"] = 14431030548
        r["previous_epoch8_offset1152_source_ZIP_SHA256"] = "26718dd109a02faccd23e7294485594a1610411dc43a8a5c6a2accc2311c4169"
        r["previous_epoch8_offset1152_source_receipt_SHA256"] = "9076c7269b517c98bf60015230ddc4806540a90065e0cf6139bee1f1f8160961"
        r["previous_epoch8_offset1152_source_status_SHA256"] = "289ec1f182885e5c30a311dbacc2c375a26aa9005ba179dbb92a2cd9f6bdd7f6"
        r["previous_epoch8_offset1152_new_source_candidates"] = 0
        r["previous_epoch8_offset1280_terminal_run"] = 44315766630
        r["previous_epoch8_offset1280_source_artifact"] = 14433646347
        r["previous_epoch8_offset1280_source_ZIP_SHA256"] = "c8a0feb237b41480eb56f7668895f17f1d4cf4c4b146e29914faef0497225e49"
        r["previous_epoch8_offset1280_source_receipt_SHA256"] = "e6e408c401bceddf184daa68b8628c59459b1a0b2445aec9d30747189dae7fb4"
        r["previous_epoch8_offset1280_source_status_SHA256"] = "a9b3e3b699bea06bc8f0ebdc29f864dc62dcde3cb3edd1bc94de5c562460e8fb"
        r["previous_epoch8_offset1280_new_source_candidates"] = 0
        r["previous_epoch8_offset1408_terminal_run"] = 44327640093
        r["previous_epoch8_offset1408_source_artifact"] = 14437289063
        r["previous_epoch8_offset1408_source_ZIP_SHA256"] = "828096ecb857e48878718f2a0b2843ca5e4727021bb7e7ebbebaf3d34f3d7059"
        r["previous_epoch8_offset1408_source_receipt_SHA256"] = "d51db154ede3c0d7b9a2f7df2c2c737adc71aa91af9a823bdd619a0f2e70972c"
        r["previous_epoch8_offset1408_source_status_SHA256"] = "3ad512539806f4cf8360bf578925302a64300bc674e0250f99efc05ecc095596"
        r["previous_epoch8_offset1408_new_source_candidates"] = 0
        r["previous_epoch8_offset1536_terminal_run"] = 44329037825
        r["previous_epoch8_offset1536_source_artifact"] = 14439091393
        r["previous_epoch8_offset1536_source_ZIP_SHA256"] = "415c6f1d2075d41612cc59ff2d2215c06975cc703e30b9713224e31f77a9a68a"
        r["previous_epoch8_offset1536_source_receipt_SHA256"] = "f2c919ff75abb00fd678dbb91cd65d6e3082ae720a026c7c90bc8c2072a4ad9a"
        r["previous_epoch8_offset1536_source_status_SHA256"] = "4a8914d379ebb7dc7431945b18db8855abf45ca129e387d6e8d97c98dd6db939"
        r["previous_epoch8_offset1536_new_source_candidates"] = 0
        r["previous_epoch8_offset1664_terminal_run"] = 44330471831
        r["previous_epoch8_offset1664_source_artifact"] = 14440065636
        r["previous_epoch8_offset1664_source_ZIP_SHA256"] = "cf6ef77f1e87264eef94eaf76345e7fef425c18a4c73c98e5161e0c8f959c0aa"
        r["previous_epoch8_offset1664_source_receipt_SHA256"] = "56acfa0cec965c04f05d562665a8533a6167d5c56743f83b71ff8de0aed3de4b"
        r["previous_epoch8_offset1664_source_status_SHA256"] = "be7de4701df8f17144f3c000052754f076a6b757690aca76a402c35c2f0e27b5"
        r["previous_epoch8_offset1664_new_source_candidates"] = 0
        r["previous_epoch8_offset1792_terminal_run"] = 44340566644
        r["previous_epoch8_offset1792_source_artifact"] = 14443567563
        r["previous_epoch8_offset1792_source_ZIP_SHA256"] = "d8392cf143b17eacc5a05724f64ba4ab05aa5f83f6eb00ed8e70def151e2dfcf"
        r["previous_epoch8_offset1792_source_receipt_SHA256"] = "b85eae13eca430b4994e2d24334759ddf8c05846629ab1ae1cb455e9016761d9"
        r["previous_epoch8_offset1792_source_status_SHA256"] = "77ef2681228159e3caf5ca11dc8f72dce0dc6abda51f301f72be0204c94e9cdf"
        r["previous_epoch8_offset1792_new_source_candidates"] = 0
        r["previous_epoch8_offset1920_terminal_run"] = 44364417366
        r["previous_epoch8_offset1920_source_artifact"] = 14454878851
        r["previous_epoch8_offset1920_source_ZIP_SHA256"] = "ae8d26c6ddad4d06390d02787d664a34c179e3b03a1402e1281d19cbb82a9548"
        r["previous_epoch8_offset1920_worker_artifact"] = 14457020969
        r["previous_epoch8_offset1920_worker_ZIP_SHA256"] = "fd7090e02cf657568d94357d55571baf33ec7f77865ad106c4a91388f25e98fe"
        r["previous_epoch8_offset1920_audit_artifact"] = 14456249618
        r["previous_epoch8_offset1920_audit_ZIP_SHA256"] = "02123c6600d094509ba59b6c622ae7feba07ff76296022d3be135406410e9704"
        r["previous_epoch8_offset1920_formal_artifact"] = 14456399686
        r["previous_epoch8_offset1920_formal_ZIP_SHA256"] = "74c5200467ed65ecfd858e093f38c1a38f92a81ea0909e4f328e3aeea49fceed"
        r["previous_epoch8_offset1920_audited_delta_SHA256"] = "5666543e87beee8ab4e63df5b3107fbea5b01b49b4ee5982c71fd33210bc5745"
        r["previous_epoch8_offset1920_original_image_raw_SHA256"] = "7cbe598b5e4b1e2cbc2ad2a4f653cecd225ee27c171fb5b42b69632df557d8bf"
        r["previous_epoch8_offset1920_technical_unknown"] = 1
        r["previous_epoch8_offset1920_new_strict9"] = 0
        r["previous_epoch6_terminal_source_ZIP_SHA256"] = (
            "3ce764b6d837348d02009173513539772d71de569088832a2f768d09c22fb87e"
        )
        r.pop("recent_97_nonempty_raw_image_SHA_count", None)
        receipt_path.write_text(json.dumps(r, sort_keys=True, indent=2) + "\n")
        print("MOIS_EPOCH9_OFFSET1408_NEW_TEMPORAL_SOURCE_SIGNED",
              json.dumps(r, sort_keys=True), flush=True)
    elif args.stage == "audit":
        m.audit(args.source, args.worker, args.output)
    else:
        m.append(args.atlas, args.source, args.worker, args.audit, args.output)

if __name__ == "__main__":
    main()
