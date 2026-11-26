# Scientific Reports revision reproducibility

This directory supplies the completed revision experiment runtime, frozen design and deterministic seeds, the immutable completed scalar results, and the scientific derivations and result tables supporting the revised TSRClust paper. It supplements the original submitted study; it does not replace the package or change its algorithm.

## Versions and scope

- The original R package 1.0.0 and original submitted reproduction code/results remain in the repository root, `reproducibility/` and `results/`. Their exact historical identity is commit `61812b21a7abc04e518133c4593d5f9c15dae99d` (26 August 2026). The root README is extended only to locate this bundle.
- `runtime/R/` contains the separate revision adapters, including `legacy_rule`, `pvalue_only`, conditional-uniform `recalibrated_rule`, macro-reference alternatives and the hierarchical all-extra evaluator. These are experimental/diagnostic variants, not silent replacements of root `R/`.
- `runtime/config/design.csv` and `tasks_tier1.csv` preserve the frozen 282 configurations and 9,880 ten-fit task blocks. Completed Tier 1 contains 98,800 successful scalar method records. Configuration-level repetitions and reused datasets must not be counted as independent evidence across parameter settings.
- `results/TSRClust_Tier1_Analysis.tar.gz` is the exact completed archive (SHA256 `c1852811d7dd94cf8f7b2dbae66d6c17f4284a6cdcf6abbfd8721d73076f94f6`). It contains replications, summaries, paired differences, accepted-count diagnostics and aggregation tables/figures. Full fitted arrays, original observation memberships and final worker checkpoints were not exported and are not claimed available.
- `analysis_sources/` contains the existing read-only scientific analyses. Reviewer/comment directory names are stable analysis identifiers only; no correspondence, response documents or private status reports are included. `tables/` supplies gzip-compressed numerical CSV outputs; `figures/` supplies their retained scientific plots.
- Exact typed original RDS result exports in `inputs/original_rds/` support the historical sufficient-count audit. They correspond to the public original CSV exports; they are not new fits. The small original automatic-K diagnostic is separately identified in `inputs/legacy_diagnostics/`.

## Verify completed results without rerunning simulations

From the repository root, with R 4.5.2 and Python 3:

```sh
D=reproducibility/scientific_reports_revision
python3 "$D/verify_manifest.py"
Rscript --vanilla "$D/verify_exports.R" "$D"
python3 "$D/prepare_workspace.py" /tmp/tsr-revision-reproduction
export TSR_TIER1_ARCHIVE="$PWD/$D/results/TSRClust_Tier1_Analysis.tar.gz"
cd /tmp/tsr-revision-reproduction
Rscript --vanilla Revision1/reviewer1/comment2/code/audit_novelty.R
Rscript --vanilla Revision1/reviewer1/comment4/code/reproduce_comment4.R
Rscript --vanilla Revision1/reviewer1/comment5/code/reproduce_comment5.R
Rscript --vanilla Revision1/reviewer2/comment1/code/reproduce_comment1.R
Rscript --vanilla Revision1/reviewer2/comment2/code/reproduce_comment2.R
Rscript --vanilla Revision1/reviewer2/comment3/code/reproduce_comment3.R
Rscript --vanilla Revision1/reviewer2/comment4/code/reproduce_comment4.R
Rscript --vanilla Revision1/reviewer1/comment3/finalization/reproduce_final_bounds.R
```

Choose a new empty workspace, not an existing results directory. The workspace adapter restores the historical relative directory layout and supplies public-only input guards. It does not require the author's scratch filesystem. Source hashes and every operational path adaptation are recorded in `source_manifest.csv`. Original private guards covered manuscripts and internal evidence as well; they are deliberately not republished. Historical local-object inventories and cross-installation Shuttle comparisons are outside public reproduction scope and are not evidence of access to missing memberships.

For a portable Shuttle count-only reproduction (no raw dataset, refit or private historical check), run from the repository root:

```sh
Rscript --vanilla reproducibility/scientific_reports_revision/reproduce_shuttle_counts.R . /tmp/tsr-shuttle-counts
```

This entry point preserves the original derivation's count/ARI expressions and removes only checks against unavailable local raw-data installations and the excluded old fit. It writes into the chosen separate output directory. Use the class trimming and microcluster summaries in root `results/real_data_results/` and the published seven-by-twelve table in `tables/reviewer2/comment5/`: the unchanged count/ARI calculations are in `analysis_sources/reviewer2/comment5/code/reproduce_comment5.R`. That archival source also records excluded checks against a historical fit and two local raw-data installations; these checks are not portable and are not required for its published sufficient-count calculations. Obtain raw observations separately through `mlbench::Shuttle`/UCI and `pdfCluster::oliveoil`; raw real datasets are not redistributed here.

## Manuscript result mapping

| Item | Producing scientific source | Archived inputs and retained outputs |
|---|---|---|
| Main Table 1 original design | root `reproducibility/summaries/generate_assets.R` | root `results/simulation_design/`; `publication_inputs/original/tables/simulation_scale.csv` |
| Main Table 2 competitor information | `analysis_sources/reviewer1/comment4/code/reproduce_comment4.R` and `analysis_sources/reviewer2/comment4/code/reproduce_comment4.R` | original and Tier 1 scalar settings; `tables/reviewer1/comment4/` and `tables/reviewer2/comment4/` |
| Main Table 3 evaluation definitions | root `reproducibility/evaluation/evaluation_metrics.R`, root `R/metrics.R`, `runtime/R/metrics.R` | `tables/reviewer1/comment1/metric_definitions.csv.gz` and metric audit exports; mathematical definitions are retained separately by version |
| Main Table 4 Tier 1 scoring audit | `analysis_sources/reviewer1/comment1/code/audit_metrics.R` | completed scalar archive; `tables/reviewer1/comment1/` |
| Main Table 5 parameter sensitivity | `analysis_sources/reviewer2/comment4/code/reproduce_comment4.R` | Tier 1 B/C; `tables/reviewer2/comment4/` |
| Main Table 6 selected K | `analysis_sources/reviewer1/comment4/code/reproduce_comment4.R` | Tier 1 D1/D2; `tables/reviewer1/comment4/` |
| Main Table 7 original headline results | root `reproducibility/summaries/20_final_summary_v2.R`, `generate_assets.R` | original replication exports; `publication_inputs/original/tables/` |
| Main Table 8 expanded null controls | `analysis_sources/reviewer2/comment1/code/reproduce_comment1.R` | Tier 1 A-END/A-COND; `tables/reviewer2/comment1/` |
| Main Table 9 Shuttle microclusters | `analysis_sources/reviewer2/comment5/code/reproduce_comment5.R` | original class/microcluster summaries; `tables/reviewer2/comment5/` |
| Main Tables 10–11 original real-data settings/agreement | root `reproducibility/real_data/02_realdata_v2.R`, `summaries/generate_assets.R` | root `results/real_data_results/` |
| Main Figures 1–4 tradeoff, stress, null controls, runtime | root `reproducibility/summaries/generate_assets.R` | original summary exports; exact retained vectors in `publication_inputs/figures/` |
| SI Table S1 combined pipeline ablation | root `reproducibility/summaries/generate_assets.R` | `publication_inputs/original/supplementary/validation_ablation_summary.csv` |
| SI Table S2 Shuttle class flow | `analysis_sources/reviewer2/comment5/code/reproduce_comment5.R` | `tables/reviewer2/comment5/shuttle_class_flow.csv.gz` |
| SI Figure S1 real-data sizes | root `reproducibility/summaries/generate_assets.R` | `publication_inputs/figures/fig6_realdata_microclusters.pdf` |
| SI Figure S2 Bypass coverage | `analysis_sources/reviewer2/comment5/code/reproduce_comment5.R` | `figures/reviewer2/comment5/shuttle_class_coverage.pdf` |
| Supplementary Data 1 | original summaries and `publication_inputs/original/supplementary/` | all original cell, paired and real-data CSV sources |
| Supplementary Data 2 | parameter/noise/Shuttle sources above | `tables/reviewer2/comment4/{cell_summaries,paired_sensitivity,density_minpts_interactions,hyperparameter_sensitivity_summary}.csv.gz`, reviewer2/comment1 null controls and reviewer2/comment5 sufficient counts |
| Supplementary Data 3 | `analysis_sources/reviewer1/comment3/finalization/reproduce_final_bounds.R` | original RDS, completed Tier 1 D1-E, `inputs/reviewer1/comment3/identified_aggregate_counts.rds`; `tables/reviewer1/comment3/finalization/` |
| Additional novelty/reference/support/MC evidence | reviewer1/comment2, reviewer1/comment5, reviewer2/comment2 and reviewer2/comment3 scientific R sources | completed Tier 1 archive and corresponding `tables/` directories; native and external events remain distinct |

The four main figures correspond to source names fig2_macro_minor_tradeoff, fig3_stress_factors, fig4_negative_control and fig5_runtime; numbering in the journal is 1–4. Publication CSV sources and the exact vectors are supplied rather than internal manuscript build scripts or reviewer text. The supplementary workbooks are editorial renderings of these numerical sources, not additional raw data. Unavailable per-observation evidence is not reconstructed.

## Dependencies and HPC reproduction

The completed revision used R 4.5.2 with the 17 versions, source URLs and SHA256 pins in `runtime/config/dependencies/source_lock.csv`. No installed binaries or private R startup files are redistributed. `fetch_dependencies.py` downloads public sources; the installer uses an isolated new library and verifies all pins. On Linux it must execute inside a SLURM compute-node allocation. The original package remains installable with `remotes::install_github("hsnbulut/TSRClust")` or `R CMD INSTALL .` with its declared dependencies.

For independent reruns, inspect `runtime/scripts/run_task.R`, the frozen task manifest, `runtime/scripts/aggregate.R` and the SLURM scripts. Work from `runtime/`; set `TSR_SOURCE_ROOT` to the repository root and choose an isolated `TSR_R_LIB` and a new output directory. Public runtime pins were regenerated only for this portable layout; `historical_runtime_md5.csv` retains the executed layout's pin bytes. Scientific generation, fitting, validation and metric function bodies are unchanged. `R/load.R` resolves the original repository's config_v2.R location. Never compare local runtime elapsed times as identical to TRUBA's recorded times.

The measured TRUBA environment was GCC `comp/gcc/14.1.0`, R `apps/R/4.5.2-gcc-14.1.0`, account `hbulut`, QoS `normal`; debug was used for acceptance. Readers must use their own allocated account/partition and scheduler limits. The same R scripts can run serially without TRUBA. Set BLAS/OpenMP threads to one for each independent fit and use unique `task_%06d.rds` paths; the worker checks source/design identity, deterministic seeds and successful-record restartability. Full reruns are opt-in; merely extracting or verifying this bundle launches no job. No new simulations were run for publication of this bundle.

## Data and limits

Synthetic observations are reproducible from frozen generators and seeds; generated observation matrices are not archived here. Original and revision scalar replication exports, summaries and diagnostic tables are public. Public real datasets retain their original attribution and must be obtained from their original providers. Row-level fitted memberships and unexported historical arrays remain unavailable; high microcluster purity, local class coverage and empirical null frequencies retain the scientific qualifications in the manuscript. These materials support the reported analyses without claiming a universal error guarantee or redistributing third-party datasets.
