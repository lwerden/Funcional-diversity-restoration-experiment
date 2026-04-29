# Code review and changes

## Changes made in this PR

### Bug fixes (would produce wrong results)

1. **Heterogeneity 2026.Rmd** (line 55-56): `summary(mod3)` and `Anova(mod3, type=3)` referenced an undefined object `mod3`. Fixed to `gini_model_3`.

2. **Canopy Cover 2026.Rmd** (lines 186-189): After fitting `gini_model_canopy`, the diagnostics block ran `summary(canopy_model)` and `simulateResiduals(fittedModel = canopy_model)` — checking the *wrong model*. Fixed to reference `gini_model_canopy`.

3. **Species performance across treatments 2026.Rmd** (lines 410-416): Two issues in one block:
   - Referenced `SubsetSummary` which does not exist. Fixed to `SubsetSummaryHeight` or `SubsetSummaryBA` as appropriate for each model.
   - RGR response variables were swapped: `lm_RGR_ba` used `Mean_RGR_Height` and `lm_RGR_height` used `Mean_RGR_BA`. Fixed to match variable names to model names.

### Label errors

4. **Canopy Cover 2026.Rmd** (line 225): Y-axis label on the Gini coefficient of canopy cover plot read "Mean Canopy Cover (%)". Fixed to "Mean Gini Coefficient of Canopy Cover".

### Dead code removal

5. **map.R**: Removed ~150 lines of unused code (lines 56-200) including abandoned Google Maps API calls, MODIS queries, and duplicate mapping attempts. Only the final working version (rnaturalearth + ggspatial) was retained.

---

## Issues flagged for discussion (not changed)

These require decisions from the analysis team before fixing.

### Statistical methodology — high priority

6. **Missing block random effect in all models**: The experiment uses a randomized complete block design (6 blocks per site), but no model in any script includes block as a random effect. This is a significant omission for a blocked design — failing to account for blocking inflates residual variance and can produce anticonservative tests. All models should include block, e.g., `(1 | Site:Block/Plot)` or similar.

7. **Site as both fixed and random effect** (Canopy Cover 2026.Rmd line 61, Heterogeneity 2026.Rmd line 53): Site appears as a fixed effect in `Treatment * Site` and simultaneously as a random effect `(1 | Site/...)`. With only 4 sites, this is statistically problematic. Recommendation: keep Site as fixed (4 levels is too few for reliable random effect variance estimation) and use `(1 | Block:Plot)` as the random structure.

8. **`lm()` used instead of `lmer()`** (Species performance 2026.Rmd lines 410-416): Six models use ordinary `lm()` on hierarchically structured data (individuals in plots in sites). This is pseudoreplication. The manuscript claims mixed-effects models were used throughout. These should use `lmer()` with appropriate random effects.

9. **Negative RGR values silently removed** (Species performance 2026.Rmd lines 255-258): All negative relative growth rates are filtered out without reporting how many observations this removes or justifying why shrinkage is biologically implausible. This introduces upward bias in growth estimates. Recommendation: report the number/percentage removed, justify the decision, and consider a sensitivity analysis with and without negative values.

10. **Standard deviations averaged incorrectly** (Species performance 2026.Rmd lines 303-313): `mean(sd_...)` is used to compute treatment-level SD, which is mathematically invalid — SDs do not average linearly. This underestimates true variability by ignoring between-plot variance. Should be recalculated from individual-level data or pooled correctly.

11. **No statistical model for ground cover** (Ground Cover year 3.Rmd): Only descriptive stacked barplots are produced. If the manuscript makes claims about treatment effects on ground cover, they are currently unsupported by any hypothesis test.

12. **Canopy cover modeled as Gaussian** (Canopy Cover 2026.Rmd): `glmmTMB` is loaded but never used. Canopy cover (bounded 0-100%) would be better modeled with beta regression, especially if values cluster near boundaries.

### Statistical methodology — moderate priority

13. **Multiple comparison correction inconsistency**: Canopy Cover and Heterogeneity use Tukey adjustment; Species performance uses Benjamini-Hochberg. These control different error rates (familywise vs. FDR). The manuscript should use one consistent approach.

14. **Error bars show raw SD, not model-based uncertainty**: All bar plots use raw standard deviations from `group_by() %>% summarise()`. Model-based confidence intervals from `emmeans()` would be more appropriate and consistent with the mixed-effects framework.

15. **Time series model lacks random slopes** (Canopy Cover 2026.Rmd line 273): Repeated measures on the same plots across years, but only a random intercept is included. A random slope for Year `(1 + Year | Plot)` should be considered.

16. **`library(car)` not loaded in Heterogeneity 2026.Rmd**: `Anova()` (capital A, Type III) requires the car package, which is not loaded in this script.

17. **Replacing zeros/NAs with 1e-6** (Species performance 2026.Rmd lines 211-228): Creates extreme artificial RGR values. The assumption that NA "meant 0" (line 223) should be verified against field protocols.

### Reproducibility

18. **Hardcoded file paths in all scripts**: All scripts use absolute paths like `~/Desktop/Master/3. master thesis/Data/...`. These should be made relative to the repo root.

19. **Missing data files**: `Height and BA all sites with mortality.csv` (used in Species performance) and `FDiv_GroundCover_2025.xlsx` (used in Ground Cover) are not in the repository.

20. **Treatment factor level inconsistency**: Canopy Cover uses lowercase (`"2sp"`, `"6sp"`, `"12sp"`) while all other scripts use uppercase (`"2SP"`, `"6SP"`, `"12SP"`).
