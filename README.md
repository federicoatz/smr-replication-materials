# Replication materials

This package reproduces the SCORE sample-size analyses and the controlled simulation reported in the manuscript.

The included public SCORE input is `g5sny-osfstorage-archive/extracted/analyst data.RData`, from https://doi.org/10.17605/OSF.IO/G5SNY. The analysis uses the non-COVID, version-of-record claim-level data. No individual participant records are used. The simulation uses seed 20260916; full scenario definitions are in `analysis/smr_revision/README.md`.

Requirements: R and packages dplyr, pwr, sandwich, lmtest, splines, ggplot2, and grid. Run the SCORE scripts from `analysis/` in this order:

```sh
Rscript manuscript_tables.R
Rscript compute_safeguard_inflation.R
Rscript regression_and_robustness.R
Rscript make_figures.R
Rscript make_protocol_figure.R
```

Run the simulation from the package root:

```sh
Rscript analysis/smr_revision/simulation.R
```

CSV files retain the analysis output names used by the scripts. The two essential empirical tables correspond to `table1_gradient.csv` and `table2_elasticity.csv`. Observation accounting is `table3_resources.csv` and `reallocation_benchmark.csv`; remaining `appendix_*.csv` files contain supplementary analyses. The simulation table is based on `analysis/smr_revision/simulation_summary.csv` and reports Monte Carlo standard errors, not empirical model standard errors. Existing numerical figure outputs are in `manuscript/`; the code overwrites only those outputs within this extracted package.

The archive contains code and derived outputs for anonymous review. Permanent archival of this article's replication package remains to be arranged before acceptance.
