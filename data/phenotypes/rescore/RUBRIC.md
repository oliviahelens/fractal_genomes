# Leaf form rescoring rubric (blind to earlier scores)

Score the **typical mature foliage leaf of an adult wild-type plant**: for rosette plants the
basal/rosette leaves; for trees and shrubs ordinary adult leaves (not juvenile or bract leaves);
for ferns with dimorphic fronds the sterile (trophophyll) fronds. Score the species as usually
described in floras; if the genus is uniform and the species is poorly known, score from the genus
and set confidence to low.

Columns (tab-separated):
- `sp`, `species`: copied from the input.
- `leaf_form`, one of:
  - `simple_entire`: undivided blade; entire, toothed, crenate or with shallow lobes (< 1/2 way to midrib/base)
  - `dissected`: continuous blade cut >= 1/2 way to the midrib/base (pinnatifid, pinnatisect, palmatisect, bipinnatifid, finely dissected); segments are NOT separate stalked/articulated leaflets
  - `compound`: blade divided into separate leaflets (pinnate, palmate, trifoliolate, bipinnate, decompound; includes fern fronds with discrete pinnae)
  - `reduced`: needles, scales, microphylls, tiny bilobed or awl-shaped leaves
  - `not_comparable`: no ordinary leaves (thalloid, leafless, phyllodes or cladodes as main photosynthetic organs, duckweed fronds, fully parasitic, cactus spines)
- `typical_order`: nested levels of division of the typical leaf. A level counts only if cuts reach >= 1/2 way to the axis they divide, or produce separate leaflets/pinnae.
  0 = simple (entire/toothed/shallowly lobed), 1 = once divided (pinnatifid; or once compound: pinnate, palmate, trifoliolate),
  2 = twice (bipinnate, bipinnatifid, pinnate with pinnatifid/deeply lobed leaflets), 3 = three or more (tripinnate and beyond, finely dissected).
- `max_order`: same scale, for the most divided mature foliage leaf the plant regularly produces (e.g. basal vs cauline leaves, fertile fronds). >= typical_order.
- `include`: 1 if leaves are comparable (all forms except `not_comparable`), else 0.
- `confidence`: high / medium / low.
- `note`: short botanical description supporting the score (e.g. "2-pinnate, pinnules pinnatifid").
