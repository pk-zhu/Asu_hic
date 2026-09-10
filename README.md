# Code release: subgenome-resolved 3D genome and expression analysis in *Arabidopsis suecica*

Custom analysis scripts underlying the manuscript. These cover the steps that are
not reproducible from a single standard-tool command: syntenic coordinate mapping
between the two *A. suecica* subgenomes and their parental genomes, gene-neighborhood
contact similarity, genome-wide cross-genome contact correlation (SCC-like),
compartment/TAD conservation, 3D remodeling trajectories, and the joint
expression–3D analyses.

Standard upstream processing (read trimming, alignment, contact-matrix building,
TAD/compartment calling, RNA-seq quantification) is **not** included; it uses
published tools with the parameters given in the Methods and is listed under
"Dependencies" below.

## Repository layout

The scripts expect the original project directory tree under a single root
(`ROOT`). Each script locates `ROOT` as follows:

1. the `ARABIDOPSIS_HIC_ROOT` environment variable, if set; otherwise
2. three levels up from the script's own location.

This resolves correctly both in the original layout
(`<ROOT>/<module>/scripts/<script>.py`) and in this release
(`<ROOT>/code_release/<module>/<script>.py`).

Point `ARABIDOPSIS_HIC_ROOT` at a directory with the following structure:

```
<ROOT>/
├── 01.single/<Species>/                 # per-species Hi-C processing
│   ├── Matrix/Final/<Species>.mcool     # balanced multi-resolution matrices
│   ├── TADs/HiCExplorer/<res>/<Species>_domains.bed
│   └── samples.txt
├── 02.collinearity/
│   ├── 02_genome_bins/                  # fixed-resolution syntenic bins
│   ├── 03_mummer/                       # nucmer/delta-filter/show-coords/show-snps output
│   ├── fastANI/                         # ANI between corresponding chromosomes
│   └── 04_collinear_pairs/              # syntenic bin-pair tables (orth_bins_*.tsv)
├── 03.pairwise_subgenome_conservation/results/
├── 04.cds_synteny_3d_conservation/results/
└── 09.expression_3d_integration/results/
```

Species abbreviations: `Asu` (*A. suecica*, natural allopolyploid),
`syn_Asu` (synthetic Allo738), `Ath` (*A. thaliana*), `Aar` (*A. arenosa*).
Subgenomes: `sT1`–`sT5` (Ath-derived), `sA6`–`sA13` (Aar-derived).

## Dependencies

Python 3 with `numpy`, `pandas`, `scipy`, `cooler` (v0.10.4 used for reading
matrices). External tools that must be on `PATH` for the relevant steps:

- HiCExplorer v3.7.3 and cooler v0.9.3 (matrix construction / TAD calling; upstream)
- MUMmer v4.0.1 (`nucmer`, `delta-filter`, `show-coords`, `show-snps`); the
  boundary-InDel script locates binaries via `PATH`, then `$MUMMER_HOME`
- MMseqs2 v15-6f452 (`easy-rbh`, one-to-one sT–sA homologs; upstream)
- R package HiCRep (within-sample SCC; upstream)
- Coolsecture (https://github.com/pk-zhu/Coolsecture); the SCC-like logic used
  here is reimplemented in `03_structure_conservation/01_bin_matrix_conservation.py`

## Modules and correspondence to Methods

### `02_synteny_mapping/` — Subgenome-aware syntenic coordinate mapping
- `02_generate_genome_bins.py` — build fixed-resolution bins for each subgenome/parent.
- `04_coords_to_collinear_synteny.py` — convert MUMmer `show-coords` alignments
  into syntenic bin pairs at 5/10/25/100 kb. Output bin pairs drive every later
  cross-genome coordinate projection.

### `03_structure_conservation/` — compartment / TAD / matrix conservation and 3D trajectories
- `01_bin_matrix_conservation.py` — genome-wide cross-genome contact similarity:
  per-distance-stratum Pearson correlations on synteny-aligned balanced matrices,
  weighted by valid-pixel count (the SCC-like statistic), at 10/25/100 kb.
- `02_compartment_conservation.py` — 100-kb bin classification
  (conserved / A_to_B / B_to_A / unmapped) after CALDER2 A/B reorientation.
- `03_tad_conservation.py` — TAD conservation score
  `min(1, max overlap / total overlapping bins)`; ≥0.8 = highly conserved.
- `06_syn_asu_conservation.py` — same conservation pipeline for the synthetic stage.
- `07_trajectory_classification.py` — PR/IR/TP/GR remodeling trajectories from
  gene-neighborhood `matrix_corr`, 33rd-percentile threshold (primary).
- `07b_trajectory_sensitivity.py` — threshold sensitivity at 25/33/40/50th percentiles.
- `07_ab_features.py` — gene density / TE coverage versus conservation (Spearman).
- `09_tad_boundary_indel.py` — InDel density at TAD boundaries (±250 bp) vs interiors.
- Sensitivity / enrichment / plotting: `03b_tad_40k.py`, `03c_tad_threshold_sensitivity.py`,
  `07c_tp_gene_enrichment.py`, `07d_trajectory_enrichment.py`, `05_plot_all.py`.
- `chr_mapping.py` — shared sT↔Ath / sA↔Aar chromosome correspondence.

### `04_gene_contact_similarity/` — gene-neighborhood contact similarity
- `03_gene_region_3d_conservation.py` — for each gene, ±50-kb balanced submatrices
  around the Asu interval and its syntenically projected parental interval
  (linear interpolation between neighboring syntenic bins); truncate to common
  minimum dimension, flatten in full (diagonal and symmetric entries included),
  masked values set to zero, Pearson correlation via NumPy. No log / O-E /
  distance normalization; genes with common dimension <4 bins or zero variance
  are excluded. Runs at 5/10/25 kb.
- `04_syn_asu_gene_3d.py` — identical with the Asu matrix replaced by the syn_Asu
  matrix (same coordinates, mappings, parental matrices, parameters).
- Sensitivity / plotting: `04b_2kb_sensitivity.py`, `04c_multires_sensitivity.py`,
  `04d_window_sensitivity_runner.py`, `04e_cds_tad_sensitivity.py`,
  `04_cds_density_vs_conservation.py`, `05_plot_all.py`.
- `chr_mapping.py` — shared chromosome correspondence.

### `09_expression/` — RNA-seq integration and expression trajectories
- `05_dominance_direct.py`, `06_dominance_allo738.py` — sT–sA expression bias
  on MMseqs2 reciprocal-best-hit one-to-one homologs; `log2FC` and biased-gene
  classification (natural Asu and synthetic Allo738).
- `07b_expression_trajectory_parental.py` — parental expression conservation:
  per-gene Δ of mean log1p-TPM versus the mapped parental gene, standardized
  `z = (Δ − median)/SD` within each comparison; |z|>1 = divergent; four
  expression trajectories.
- `02_gene_3d_vs_expr.py`, `01_bin_3d_vs_expr.py` — joint expression × 3D
  (compartment × TAD region) classification.
- `07d_expr_vs_3d_trajectory.py` — cross-tabulation of expression and 3D trajectories.
- `08_sample_spearman.py` — between-replicate RNA-seq Spearman correlations.
- `09_gene_region_expression_joint.py` — gene-level joint expression/structure table.
- `common.py` — shared paths and helpers.

## Notes on reproducibility

- Set `ARABIDOPSIS_HIC_ROOT` (or run inside the original tree) before executing.
- Scripts are organized by analysis step; within a module the numeric prefixes
  indicate run order (mapping → conservation → trajectories → expression integration).
- Figures are generated by the `05_plot_all.py` scripts in each module.
- Intermediate result tables are written under each module's `results/` directory;
  create these directories or let the scripts create them as needed.
