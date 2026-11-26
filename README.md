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
| 5 | `scripts/05_fetch_proteomes.sh` | PLAZA proteomes for all species (`data/raw/proteomes/`) |
| 6 | `scripts/06_gene_trees.py [n\|all] [workers]` | per family: MAFFT, one copy per species, trimAl, IQ-TREE branch lengths on the fixed species tree (3,989 near-single-copy families; ~4 h on 4 CPUs) |
| 7 | `scripts/07_rer.py [n_null] [n_null_cv] [K]` | `results/rer/`: relative evolutionary rates (RERconverge-style), per-gene tests vs permulated traits, leave-one-clade-out classifier |

Run: `pip install -r requirements.txt && scripts/01_fetch_plaza.sh && python3 scripts/02_build_family_matrix.py && python3 scripts/03_species_tree.py && python3 scripts/04_associate_and_predict.py 50 20`

## First-pass result (`results/summary.txt`)
130 species with comparable leaves; 9 "fractal" (order >= 2) in 8 independent clades; 10,067 families present in >= 10 species.
- **No predictive signal.** Leave-one-clade-out AUROC 0.38 (null traits: 0.44-0.50).
- **No credible association.** 250 families pass q<0.05, but random traits simulated on the same tree give a median of ~150 (up to ~800), so the test is inflated and the count is within noise.
- **Positive controls don't stand out.** KNOX1, CUC/NAM, LFY, UFO, ARP, RCO-family copy numbers are not enriched. Consistent with the literature: leaf dissection evolves mostly through changes in *where/when* these genes are expressed, which copy number cannot see.

## Second pass: evolutionary rates (`results/rer/summary.txt`)
Does any gene evolve faster or slower in fractal lineages (a common signature of convergent change in function or constraint)?
3,989 families present in >= 80% of species, mostly single-copy; terminal-branch rates relative to the genome-wide expectation.
- **No gene passes FDR** (0 at q<0.1, for fractal-vs-not and for correlation with dissection order). Small p-values are rarer than expected by chance (17 at p<0.01 vs 40 expected), so there is no hidden excess either.
- **No predictive signal.** Leave-one-clade-out AUROC 0.43; random traits on the same tree give 0.40 on average (95th percentile 0.68).

## Caveats
- Phenotype scores are curated from botanical descriptions, not measured; image-based fractal dimension is the planned check.
- Open Tree topology has no branch lengths (Grafen lengths used).
- Only 9 fractal species in 8 clades: power to detect a gene shared by most origins is limited, and a gene used in only 2-3 origins is undetectable here.
- Gymnosperm families map less completely onto the backbone (63% of genes).
