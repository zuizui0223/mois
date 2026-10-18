# FCP H9 photo acquisition applied to MOIS without the FCP within-species design

**18 October 2026.** The broad MOIS first256 original RGB-image pilot had
69 successful visible colour ROIs and 187 classification-unknown source
photos (173 ROI/flip-gate failures and six original image access failures).
All 69 successful source taxa were new angiosperm genera/species; keep
those colours and do NOT acquire redundant photos for them.

The original frozen FCP acquisition code at commit
9fae6ccdf684a46026f72ba12e98de2c5c54bf2a (the code file
fcp_pipeline/random_photo_pool.py has Git blob SHA
a46341c98aa7ada422ed0e2799ccb863f01a3626)
implements a measured and considerate iNaturalist API client
with at least 1.05 seconds between requests, up to four bounded retries,
per_page=200/page=1, research-grade licensed photos, geo=true,
<=5000m positional uncertainty, unobscured coordinates and
plant phenology Flowers and Fruits/Flowers annotation 12:13.

The new MOIS collector uses EXACTLY this method for no more than
the 187 species whose first source photograph did not yield four-hue
colour data; it does not need 100 images for every species as FCP does.
It queries one FCP-style randomly ordered flowering-annotated API
page per failed source taxon, locally rechecks named taxon ID/species,
licence, original image URL and observer/observation identity, and
retains at most one alternative photo selected by outcome-blind hash.
The original 313535 measured photo IDs, all 32768 future source
photo IDs, existing first4096 source, independent fifth4096
and separate seasonal/SI4096 photo reservations are excluded.
No phenotype hue label is looked up for the candidate photo.

Sources:
https://help.inaturalist.org/en/support/solutions/articles/151000170738
https://github.com/zuizui0223/fcp/blob/9fae6ccdf684a46026f72ba12e98de2c5c54bf2a/fcp_pipeline/random_photo_pool.py

**Crucial limitation.** A flowering annotation is attached to an
OBSERVATION, not necessarily to each individual PHOTO in that observation.
The client metadata is source-only. It DOES NOT claim that a photo has
verified petals or that new photos have classifiable white/pigmented hue.
A follow-up real RGB original YOLO11/EfficientSAM ROI-v4 step must report
known four-hue photo labels, unknowns, source-photo IDs and
measured new species/genera per photographed attempt.

Flowering-annotated records are not an unbiased random subsample of
the worldwide photo population. Their later visible-colour fractions
must be kept separate from the original global 313535-source ratio
and compared with a non-annotated source arm at matched species/regions
before inferring ecological pigment composition. Pigment production cost,
UV colours, selection and reproductive fitness are not measured.

The FCP method's greatest transferable benefit is using the API
to preselect likely flowers BEFORE decoding and evaluating expensive
RGB pixels. This follow-up requires up to 187 deliberately spaced
API calls, versus rescanning ~500 million monthly provider photo
metadata records and rerunning the original flower ROI estimator on
all 256 old source photos.
