# Species diversity reduces risk in tropical forest restoration: a portfolio effect across heterogeneous sites

*Title TBD*

Code for Quaglino, Cole, Quirós Cedeño, Rodriguez, Tingle, Manyard, Weissenhofer, Crowther & Werden — targeting Global Change Biology.

## Data

- `FDiv_PlantedTreeData_final_excel.xlsx` — raw tree and canopy cover data (4 sites, 144 plots, year 3)

## Analysis scripts

Run in order. All scripts source `01_preprocess.R` for cleaned data.

| Script | Purpose | Outputs |
|---|---|---|
| `01_preprocess.R` | Load raw data, fix errors, back-calculate DBH, attach wood density from 5 databases, standardize treatments | `data_cleaned.rds`, `canopy_cleaned.rds` |
| `02_main_figures.R` | Fig 2 (height, canopy, carbon, survival) and Fig 3 (Gini coefficients) — pooled bar plots with CLD letters and eta-squared | `fig2_height_canopy_carbon.png`, `fig3_gini.png` |
| `03_synthesis_figures.R` | Fig 4: radar plot, PCA, betadisper multivariate dispersion analysis | `fig4_synthesis.png`, `supp_betadisper_pairwise.csv` |
| `04_sem.R` | Fig 5: piecewise SEM path diagrams (main text: log richness; supplement: Inga proportion) | `fig5_sem.png`, `supp_sem_inga.png`, `supp_sem_coefficients.csv` |
| `05_supplement_figures.R` | By-site bar plots (7 metrics) + combined Inga performance figure | `supp_*_by_site.jpeg`, `supp_inga_performance.png` |
| `06_tables.R` | Table 3 (ANOVA results) and supplement tables (site means, wood density, DHARMa diagnostics) | `table2_anova.csv`, `supp_*.csv` |
| `07_sensitivity_analyses.R` | LOO, permutation test, risk-return, temporal trajectory, replanting sensitivity | `supp_loo_sensitivity.csv`, `supp_permutation_test.csv`, `supp_risk_return.png`, `supp_temporal_dispersion.png`, `supp_replanting_sensitivity.csv` |

## exploratory/

Early analysis scripts and their outputs, moved out of the main pipeline. Includes Rachele's original Rmd reports (canopy cover, heterogeneity, ground cover, species performance), an early biomass script (`04_biomass_carbon.R`), the site map script (`map.R`), and intermediate cached .rds files. Not needed for reproduction.

## Model specification

All plot-level mixed models: `response ~ Treatment * Site + (1|Site:Plot_number)`

- Site is a fixed effect (4 levels)
- Block (= Plot_number) nested within Site as random intercept
- Survival: binomial GLMM via glmmTMB
- Carbon and basal area: log-transformed
- Pairwise comparisons: Tukey-adjusted via emmeans

## Key R packages

lme4, lmerTest, glmmTMB, emmeans, piecewiseSEM, vegan (betadisper), BIOMASS, ineq (Gini), ggplot2, patchwork, DHARMa

## Reproduction

```r
# Install dependencies
install.packages(c("readxl", "dplyr", "ggplot2", "lme4", "lmerTest",
  "glmmTMB", "emmeans", "multcomp", "ineq", "BIOMASS", "vegan",
  "piecewiseSEM", "DHARMa", "patchwork", "ggrepel"))

# Run full pipeline
source("01_preprocess.R")
source("02_main_figures.R")
source("03_synthesis_figures.R")
source("04_sem.R")
source("05_supplement_figures.R")
source("06_tables.R")
source("07_sensitivity_analyses.R")
```
