# fractal_genomes

Can gene families (not just known leaf-shape genes) predict whether a plant has "fractal"
(repeatedly subdivided) leaves, across independent origins of that trait?

## Pipeline
| step | script | output |
|---|---|---|
| 1 | `scripts/01_fetch_plaza.sh` | PLAZA gene families (dicots 5, monocots 5, basal 1 [ferns], gymno 1) |
| 2 | `scripts/02_build_family_matrix.py` | `data/processed/family_copy_number.tsv.gz`: 146 genomes x 22k families, cross-mapped onto the dicots-5 backbone via shared species |
| 3 | `scripts/03_species_tree.py` | `data/processed/species_tree.nwk` from Open Tree of Life |
| - | `data/phenotypes/leaf_dissection_curated.tsv` | leaf dissection order 0-3 per species, hand-curated, with confidence flags |
| 4 | `scripts/04_associate_and_predict.py [K] [n_nulls]` | `results/`: PGLS association per family, positive controls, leave-one-clade-out classifier vs null traits |

Run: `pip install -r requirements.txt && scripts/01_fetch_plaza.sh && python3 scripts/02_build_family_matrix.py && python3 scripts/03_species_tree.py && python3 scripts/04_associate_and_predict.py 50 20`

## First-pass result (`results/summary.txt`)
130 species with comparable leaves; 9 "fractal" (order >= 2) in 8 independent clades; 10,067 families present in >= 10 species.
- **No predictive signal.** Leave-one-clade-out AUROC 0.38 (null traits: 0.44-0.50).
- **No credible association.** 250 families pass q<0.05, but random traits simulated on the same tree give a median of ~150 (up to ~800), so the test is inflated and the count is within noise.
- **Positive controls don't stand out.** KNOX1, CUC/NAM, LFY, UFO, ARP, RCO-family copy numbers are not enriched. Consistent with the literature: leaf dissection evolves mostly through changes in *where/when* these genes are expressed, which copy number cannot see.

## Caveats
- Phenotype scores are curated from botanical descriptions, not measured; image-based fractal dimension is the planned check.
- Open Tree topology has no branch lengths (Grafen lengths used).
- Gymnosperm families map less completely onto the backbone (63% of genes).
