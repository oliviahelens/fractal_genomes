# fractal_genomes

Effort towards discovering gene constellations that predict whether a plant has "fractal"
or repeatedly subdivided leaves, across independent origins.

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
| 8 | `scripts/08_fetch_onekp.py` | 1KP transcriptome protein sets (1,350 samples) via the CyVerse Data Commons (Terrain API, anonymous) |
| 9 | `scripts/09_fetch_ncbi.py` | NCBI annotated vascular-plant genomes not in PLAZA (423), via NCBI Datasets |
| 10 | `scripts/10_species_table.py` | `data/processed/species_table.tsv`: one entry per species (PLAZA > RefSeq > GenBank > 1KP), 1KP consortium-flagged samples removed, NCBI taxonomy |
| - | `data/phenotypes/curation_parts/scored_*.tsv` | leaf dissection scores for the new species (same 0-3 scale, confidence flags) |
| 11 | `scripts/11_assign_families.py [threads]` | DIAMOND assignment of new proteins to the 3,989 families (reference: 30 PLAZA genomes); best hit per family; `data/processed/completeness.tsv` |
| 12 | `scripts/12_species_tree_expanded.py` | Open Tree topology for all species (1,225 tips; missing species grafted next to congeners) |
| 13 | `scripts/13_gene_trees_expanded.py [n\|all] [workers]` | new sequences added to the PLAZA alignments (MAFFT --addfragments --keeplength), lineage-mismatch and fragment filters, FastTree branch lengths on the fixed topology (~7 h on 4 CPUs) |
| 14 | `scripts/14_rer_expanded.py [n_null] [n_null_cv] [K]` | `results/rer_expanded/`: step 7 on the expanded set, leave-one-order-out, genome/transcriptome artifact checks |
| 15 | `scripts/15_within_order_test.py [n_perm]` | within-order comparison (labels shuffled within orders) |

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

## Third pass: 1,133 species with genomes and transcriptomes (`results/rer_expanded/`)
Added 417 NCBI genomes and 655 1KP transcriptomes to the same 3,989 families, to get power: 1,133 species pass all filters
(527 genomes, 606 transcriptomes), **90 fractal species in 24 orders** (vs 9 in 8 clades before).

Keeping the sources comparable:
- Same features for every species: evolutionary rates only, never copy number (transcriptomes miss unexpressed genes and split genes into isoforms).
- One sequence per species and family: the best DIAMOND hit, which must cover >= 50% of the reference protein. Transcript fragments are added to the fixed PLAZA alignments, so they cannot reshape them.
- Filters: species need >= 50% of families (1,054 of 1,072 new species pass; 1KP median 89%); sequences need >= 50% of alignment columns filled; sequences whose closest PLAZA relative is from another major lineage are dropped (contamination or a different paralog; ~10%); samples the 1KP consortium flagged are excluded.
- Checks: every test is repeated on genomes only. The same features are also used to try to predict genome vs transcriptome.

Results:
- **No gene passes FDR** (q<0.1: 0; also 0 for dissection order, genomes only, and without low-confidence scores). Small p-values are rarer than chance under the phylogenetic null (6 at p<0.01 vs 40 expected).
- **No predictive signal.** Leave-one-order-out AUROC 0.64; random traits on the same tree give 0.59 on average (95th percentile 0.87; P = 0.33).
- **Within-order test** (`within_order_summary.txt`; 19 orders that contain both fractal and non-fractal species, labels shuffled within each order): 0 genes at q<0.1. There is a mild excess of small p-values (77 at p<0.01 vs 40 expected; 8 at p<0.001 vs 4). That is expected from clustering of fractal species within families, which this test does not remove. The most consistent gene, HOM05D004396 (*DAD1/DAD2*, slower in fractal species, p = 0.001 in genomes and 0.01 in transcriptomes), is a weak lead, not a finding.
- **Data-type warning:** the rate features distinguish genomes from transcriptomes with AUROC 0.98. Source leaves a strong mark, so any future hit must replicate within genomes.

Interpretation: three independent passes (copy number; rates on 130 species; rates on 1,133 species) find no genome-wide gene pattern shared by fractal leaves. The power is now adequate for genes used by many origins. What remains open:
- Origin-specific changes: genes used in only a few origins.
- Regulatory and expression changes, which are known to drive leaf dissection (KNOX, CUC, RCO).
- Changes on the ancestral branches where dissection arose. The tests here use terminal branches, and fern dissection is ancient.

## Caveats
- Phenotype scores are curated from botanical descriptions, not measured; image-based fractal dimension is the planned check.
- Open Tree topology has no branch lengths (Grafen lengths used).
- First two passes: only 9 fractal species in 8 clades (low power). Third pass: 90 in 24 orders, but ferns and Ranunculales contribute many of them.
- Terminal-branch rates only; an ancestral-state test on internal branches (where dissection arose) is the natural next design.
- Gymnosperm families map less completely onto the backbone (63% of genes).
