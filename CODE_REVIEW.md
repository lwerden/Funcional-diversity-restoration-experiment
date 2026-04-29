# Code review and changes

## Summary

All analysis scripts were rewritten with corrected model specifications. The original scripts used `Treatment + (1|Site)` (Site as random effect) but with only 4 sites, Site must be a fixed effect. The correct model is `Treatment * Site + (1 | Site:Plot_number)` where Plot_number represents the block in the RCBD design.

All conclusions hold in the same direction and are generally **stronger** with the correct models, but exact F-values, p-values, and some descriptive statistics change. Several significant Treatment:Site interactions were uncovered that need to be acknowledged in the manuscript.

## Script changes

### All scripts
- **Model specification**: Changed from `Treatment + (1|Site)` or `Treatment * Site + (1|Site/Plot)` to `Treatment * Site + (1 | Site:Plot_number)` — Site as fixed effect, block (Plot_number) nested within Site as random
- **Data paths**: Changed from hardcoded `~/Desktop/Master/...` to relative paths (`FDiv_PlantedTreeData_final_excel.xlsx`)
- **Output paths**: Changed from `~/Desktop/...` to relative paths
- **ggplot2**: Replaced deprecated `size` with `linewidth` in line geoms
- **Removed unused package loads** (glmmTMB, performance where not used)

### Canopy Cover 2026.Rmd
- Fixed bug: Gini canopy diagnostics ran on wrong model object (`canopy_model` instead of `gini_model_canopy`)
- Fixed bug: Gini plot y-axis label said "Mean Canopy Cover (%)" instead of "Mean Gini Coefficient of Canopy Cover"
- Removed duplicate library() calls and redundant code
- Refactored time series site plots into a loop

### Heterogeneity 2026.Rmd
- Fixed bug: `summary(mod3)` and `Anova(mod3)` referenced undefined object — fixed to `gini_model_3`
- Added `library(car)` (was missing, needed for `Anova()`)
- Applied `cld()` to `$emmeans` component consistently

### Species performance across treatments 2026.Rmd
- **Complete rewrite**: original script read from a CSV not in the repo; now uses the xlsx `SurvivalGrowth` sheet
- Fixed bug: RGR response variables were swapped (`lm_RGR_ba` used Height, `lm_RGR_height` used BA)
- Fixed bug: referenced undefined `SubsetSummary` data frame
- Changed `lm()` to `lmer()` with proper random effects (was pseudoreplication)
- Removed incorrect SD averaging (`mean(sd_...)`)
- Negative RGR values are now **retained** (were silently dropped; 49 of 1163 = 4.2%)
- RGR analysis limited to 3 sites (FAB, FCEA, LLFS) since GAM has no year 2 data

### Ground Cover year 3.Rmd
- Fixed data path (was missing `~/` prefix, inconsistent with other scripts)
- Cleaned up code, removed duplicate library() calls
- Added TODO note: no statistical model exists — currently descriptive only

### map.R
- Removed ~150 lines of dead code (abandoned Google Maps, MODIS, ggmap, tmap attempts)

---

## Manuscript changes required

### Statistical methods section (needs rewriting)

The methods should state:
> "We built linear mixed-effects models using the lme4 package in R with treatment, site, and their interaction as fixed effects, and block nested within site as a random intercept: `Response ~ Treatment * Site + (1 | Site:Block)`."

### Results — updated F-values

All F-values change because the model specification changed. Conclusions are the same or stronger.

| Metric | Old F | New F | Old p | New p | Notes |
|--------|-------|-------|-------|-------|-------|
| Canopy cover (Treatment) | 29.04 | **32.28** | < 0.001 | < 2.2e-16 | |
| Height Gini (Treatment) | 85.39 | **135.82** | < 0.001 | < 2.2e-16 | **Treatment:Site now significant (F=5.23, p=2.2e-06)** |
| Canopy Gini (Treatment) | 21.98 | **24.59** | < 0.001 | 8.1e-16 | |
| Total BA (Treatment) | 16.88 | **20.01** | < 0.001 | 1.9e-11 | |
| Mean height (Treatment) | — | **46.34** | — | < 2.2e-16 | Not reported as F-value in ms |
| Inga BA (Treatment) | 3.81 | **11.54** | 0.02 | 1.7e-07 | **Treatment:Site now significant (F=8.78, p=5.3e-13)** |
| Inga height (Treatment) | 3.05 | **65.94** | 0.04 | < 2.2e-16 | **Treatment:Site now significant (F=10.96, p<2.2e-16)** |
| Inga RGR height (Treatment) | 0.27 | **2.05** | > 0.05 | 0.11 | Still non-significant, consistent |

### Results — new significant interactions to report

1. **Height Gini: Treatment:Site interaction** (F = 5.23, p = 2.2e-06). The manuscript says height heterogeneity "increased sequentially along the functional diversity gradient across all sites" — this needs to be qualified. The pattern holds overall but varies in magnitude by site.

2. **Inga BA and Height: Treatment:Site interactions** are highly significant. The manuscript describes a simple monotonic decrease from M1 to diverse treatments, but with 4 sites this pattern is site-dependent.

### Results — Site effects now visible

Site is significant in all models (was hidden when Site was a random effect):
- Canopy cover: F = 9.50, p = 0.0004
- Height Gini: F = 8.03, p = 0.001
- Canopy Gini: F = 20.06, p = 4.8e-06
- Total BA: F = 8.05, p = 0.001
- Mean height: F = 12.28, p = 9.4e-05

### Results — descriptive values confirmed

These values match the manuscript and do not need changing:
- M1 mean heights: FCEA = 323 cm, FAB = 290 cm, LLFS = 463 cm, GAM = 402 cm, pooled = 370 cm
- M1 has highest canopy cover at all sites
- M1 has highest total BA

### Results — Inga individual-level pattern may have changed

With all 4 sites, Inga mean individual BA by treatment:
- M1: 32.2 cm²
- 2SP: 27.3 cm² (15% lower than M1)
- 6SP: 36.7 cm² (14% **higher** than M1)
- 12SP: 46.0 cm² (43% **higher** than M1)

The manuscript states "mean BA of Inga edulis decreased" along the diversity gradient. This no longer holds for 6SP and 12SP. **This needs careful re-examination** — the significant Treatment:Site interaction means the pattern differs by site, so the pooled means may be misleading. Site-specific emmeans should be reported.

### Discussion — new limitation text needed

Add to "Importance and limitations" section:
> "Our design manipulated functional diversity by increasing species richness along a trait-informed gradient, but practical constraints in seedling availability prevented inclusion of high-richness/low-functional-diversity treatments that would have allowed us to fully disentangle species richness from functional diversity effects."

### Ground cover — no model

The ground cover section is currently descriptive only (stacked barplots). If the manuscript makes causal claims about treatment effects on ground cover composition, a statistical model (e.g., beta regression via glmmTMB) should be added.

### Reproducibility

- `FDiv_GroundCover_2025.xlsx` needs to be added to the repo
- The old `Height and BA all sites with mortality.csv` is no longer needed (Species performance now reads from the xlsx)
