# Sensitivity analyses: LOO, permutation test, risk-return, temporal trajectory, replanting, mixed model
library(dplyr)
library(ggplot2)
library(ineq)
library(BIOMASS)
library(vegan)
library(lme4)
library(lmerTest)
library(patchwork)

source("01_preprocess.R")
CC3 <- CanopyCover %>% filter(Year == "3")

# Rebuild the same 6 plot-level metrics used in 05_synthesis_figures.R
ht <- TreeData %>% filter(!is.na(Height_year_3)) %>%
  group_by(Site, Treatment, Plot_number) %>%
  summarise(mean_height = mean(Height_year_3 / 100), .groups = "drop")

bio <- TreeData %>% filter(!is.na(DBH_year_3), !is.na(Height_year_3), Height_year_3 > 0, !is.na(WD))
bio$AGB_Mg <- computeAGB(D = bio$DBH_year_3, WD = bio$WD, H = bio$Height_year_3 / 100)
bio$C_kg <- bio$AGB_Mg * 1000 * 0.47

c_plot <- bio %>% group_by(Site, Treatment, Plot_number) %>%
  summarise(C_Mg_ha = sum(AGB_Mg) / 0.0225 * 0.47, .groups = "drop")

gini_ht <- TreeData %>% filter(!is.na(Height_year_3)) %>%
  group_by(Site, Treatment, Plot_number) %>%
  summarize(height_gini = ineq(Height_year_3, type = "Gini"), .groups = "drop")

carbon_gini <- bio %>% group_by(Site, Treatment, Plot_number) %>%
  summarize(carbon_gini = ineq(C_kg, type = "Gini"), .groups = "drop") %>% filter(!is.nan(carbon_gini))

cc_plot <- CC3 %>% group_by(Site, Treatment, Plot_number) %>%
  summarise(canopy_cover = mean(Percent_canopy_cover, na.rm = TRUE), .groups = "drop")

gini_cc <- CC3 %>% group_by(Site, Treatment, Plot_number) %>%
  summarize(canopy_gini = ineq(Percent_canopy_cover, type = "Gini"), .groups = "drop") %>% filter(!is.nan(canopy_gini))

pm <- ht %>%
  inner_join(c_plot, by = c("Site", "Treatment", "Plot_number")) %>%
  inner_join(gini_ht, by = c("Site", "Treatment", "Plot_number")) %>%
  inner_join(carbon_gini, by = c("Site", "Treatment", "Plot_number")) %>%
  left_join(cc_plot, by = c("Site", "Treatment", "Plot_number")) %>%
  left_join(gini_cc, by = c("Site", "Treatment", "Plot_number")) %>%
  filter(!is.na(canopy_gini))

pca_mat <- pm %>%
  select(mean_height, C_Mg_ha, canopy_cover, height_gini, canopy_gini, carbon_gini) %>%
  scale()

no_m2 <- pm$Treatment != "M2"
pm_no_m2 <- pm[no_m2, ]
mat_no_m2 <- pca_mat[no_m2, ]

trt_colors <- c("M1" = "#D55E00", "2SP" = "darkseagreen2", "6SP" = "darkseagreen", "12SP" = "darkolivegreen")
richness_map <- c("M1" = 1, "2SP" = 2, "6SP" = 6, "12SP" = 12)
richness <- richness_map[as.character(pm_no_m2$Treatment)]

bd <- betadisper(dist(mat_no_m2), droplevels(pm_no_m2$Treatment))

# ============================================================
# A) LOO SENSITIVITY — drop each site, re-run richness-dispersion correlation
# ============================================================
sites <- unique(pm_no_m2$Site)
loo_results <- data.frame()
for (drop_site in c("none", sites)) {
  if (drop_site == "none") {
    keep <- rep(TRUE, nrow(pm_no_m2))
  } else {
    keep <- pm_no_m2$Site != drop_site
  }
  bd_loo <- betadisper(dist(mat_no_m2[keep, ]), droplevels(pm_no_m2$Treatment[keep]))
  rich_loo <- richness[keep]
  ct <- suppressWarnings(cor.test(rich_loo, bd_loo$distances, method = "spearman"))
  loo_results <- rbind(loo_results, data.frame(
    Dropped = ifelse(drop_site == "none", "None (full)", drop_site),
    n_plots = sum(keep), rho = round(ct$estimate, 3),
    p_value = round(ct$p.value, 4),
    Significant = ifelse(ct$p.value < 0.05, "*", "")))
}
write.csv(loo_results, "supp_loo_sensitivity.csv", row.names = FALSE)
cat("Saved: supp_loo_sensitivity.csv\n")
print(loo_results)

# ============================================================
# B) PERMUTATION TEST — null model for richness-dispersion correlation
# Permute richness labels 9999x, compute SES against null distribution
# ============================================================
set.seed(42)
n_perm <- 9999
obs_rho <- cor(richness, bd$distances, method = "spearman")

perm_rhos <- replicate(n_perm, {
  cor(sample(richness), bd$distances, method = "spearman")
})
p_perm <- (sum(perm_rhos <= obs_rho) + 1) / (n_perm + 1)
ses <- (obs_rho - mean(perm_rhos)) / sd(perm_rhos)

cat(sprintf("\nPermutation test: rho = %.4f, p = %.4f, SES = %.2f\n", obs_rho, p_perm, ses))

perm_results <- data.frame(
  Test = c("Spearman rho", "Permutation p (9999 perms, one-tailed)", "SES"),
  Value = c(round(obs_rho, 4), round(p_perm, 4), round(ses, 2)))
write.csv(perm_results, "supp_permutation_test.csv", row.names = FALSE)
cat("Saved: supp_permutation_test.csv\n")

# ============================================================
# C) RISK-RETURN FRONTIER — mean vs CV for each metric by treatment
# ============================================================
rr <- pm_no_m2 %>%
  group_by(Treatment) %>%
  summarise(
    Height_mean = mean(mean_height), Height_cv = sd(mean_height) / mean(mean_height) * 100,
    Carbon_mean = mean(C_Mg_ha), Carbon_cv = sd(C_Mg_ha) / mean(C_Mg_ha) * 100,
    Canopy_mean = mean(canopy_cover), Canopy_cv = sd(canopy_cover) / mean(canopy_cover) * 100,
    .groups = "drop")

mv_rr <- data.frame(Treatment = names(bd$group.distances),
                     Dispersion = bd$group.distances,
                     Richness = richness_map[names(bd$group.distances)])

p_rr <- ggplot(mv_rr, aes(x = Dispersion, y = Richness, color = Treatment)) +
  geom_point(size = 5) +
  geom_text(aes(label = Treatment), hjust = -0.3, fontface = "bold", size = 4, show.legend = FALSE) +
  scale_color_manual(values = trt_colors, guide = "none") +
  labs(x = "Multivariate dispersion\n(distance to centroid)", y = "Planted species richness") +
  theme_minimal()

ggsave("supp_risk_return.png", plot = p_rr, width = 6, height = 5, dpi = 300, bg = "white")
cat("Saved: supp_risk_return.png\n")

# ============================================================
# D) TEMPORAL TRAJECTORY — does dispersion-richness trend strengthen over years 1-3?
# Only 3 sites have all 3 years (GAM planted last, only has year 3)
# ============================================================
temporal_sites <- c("FAB", "FCEA", "LLFS")
temporal_results <- data.frame()

for (yr in 1:3) {
  ht_col <- paste0("Height_year_", yr)

  ht_yr <- TreeData %>%
    filter(!is.na(.data[[ht_col]]), Site %in% temporal_sites) %>%
    group_by(Site, Treatment, Plot_number) %>%
    summarise(mean_height = mean(.data[[ht_col]] / 100), .groups = "drop")

  gini_yr <- TreeData %>%
    filter(!is.na(.data[[ht_col]]), Site %in% temporal_sites) %>%
    group_by(Site, Treatment, Plot_number) %>%
    summarize(height_gini = ineq(.data[[ht_col]], type = "Gini"), .groups = "drop") %>%
    filter(!is.nan(height_gini))

  yr_data <- ht_yr %>%
    inner_join(gini_yr, by = c("Site", "Treatment", "Plot_number")) %>%
    filter(Treatment != "M2")

  if (nrow(yr_data) < 10) next

  yr_mat <- yr_data %>% select(mean_height, height_gini) %>% scale()
  yr_rich <- richness_map[as.character(yr_data$Treatment)]

  bd_yr <- betadisper(dist(yr_mat), droplevels(yr_data$Treatment))
  ct_yr <- suppressWarnings(cor.test(yr_rich, bd_yr$distances, method = "spearman"))

  dists_yr <- bd_yr$group.distances
  temporal_results <- rbind(temporal_results, data.frame(
    Year = yr, n_plots = nrow(yr_data),
    rho = round(ct_yr$estimate, 3), p_value = round(ct_yr$p.value, 4),
    M1_disp = round(dists_yr["M1"], 3),
    SP2_disp = round(dists_yr["2SP"], 3),
    SP6_disp = round(dists_yr["6SP"], 3),
    SP12_disp = round(dists_yr["12SP"], 3)))
}

write.csv(temporal_results, "supp_temporal_dispersion.csv", row.names = FALSE)
cat("Saved: supp_temporal_dispersion.csv\n")
print(temporal_results)

# Temporal dispersion figure dropped from supplement — results in Table S3E only

# ============================================================
# E) REPLANTING SENSITIVITY
# ============================================================
# Identify replanted trees: died year 1 or 2 but alive at year 3
TreeData <- TreeData %>%
  mutate(is_replant = (Mortality_year_1 == 1 | Mortality_year_2 == 1) & !is.na(Height_year_3))
TreeData$is_replant[is.na(TreeData$is_replant)] <- FALSE
TD_norep <- TreeData %>% filter(!is_replant)

cat(sprintf("\n=== REPLANTING SENSITIVITY ===\n"))
cat(sprintf("Replants: %d of %d trees (%.1f%%)\n",
    sum(TreeData$is_replant), nrow(TreeData),
    sum(TreeData$is_replant) / nrow(TreeData) * 100))

# Height Gini WITH replants (full dataset baseline)
ht_wr <- TreeData %>% filter(!is.na(Height_year_3), Treatment != "M2") %>%
  group_by(Site, Treatment, Plot_number) %>%
  summarise(height_gini = ineq(Height_year_3, type = "Gini"), .groups = "drop") %>%
  filter(!is.nan(height_gini))

m_gini_wr <- suppressWarnings(lmer(height_gini ~ Treatment * Site + (1 | Site:Plot_number), data = ht_wr))
a_wr <- anova(m_gini_wr)
eta_wr <- a_wr["Treatment", "F value"] * a_wr["Treatment", "NumDF"] /
  (a_wr["Treatment", "F value"] * a_wr["Treatment", "NumDF"] + a_wr["Treatment", "DenDF"])

# Height Gini WITHOUT replants
ht_nr <- TD_norep %>% filter(!is.na(Height_year_3)) %>%
  group_by(Site, Treatment, Plot_number) %>%
  summarise(height_gini = ineq(Height_year_3, type = "Gini"), .groups = "drop") %>%
  filter(!is.nan(height_gini), Treatment != "M2")

m_gini_nr <- suppressWarnings(lmer(height_gini ~ Treatment * Site + (1 | Site:Plot_number), data = ht_nr))
a_nr <- anova(m_gini_nr)
eta_nr <- a_nr["Treatment", "F value"] * a_nr["Treatment", "NumDF"] /
  (a_nr["Treatment", "F value"] * a_nr["Treatment", "NumDF"] + a_nr["Treatment", "DenDF"])

# Betadisper without replants
bio_nr <- TD_norep %>% filter(!is.na(DBH_year_3), !is.na(Height_year_3), Height_year_3 > 0, !is.na(WD))
bio_nr$AGB_Mg <- computeAGB(D = bio_nr$DBH_year_3, WD = bio_nr$WD, H = bio_nr$Height_year_3 / 100)
bio_nr$C_kg <- bio_nr$AGB_Mg * 1000 * 0.47

ht_nr2 <- TD_norep %>% filter(!is.na(Height_year_3)) %>%
  group_by(Site, Treatment, Plot_number) %>%
  summarise(mean_height = mean(Height_year_3 / 100), .groups = "drop")
gini_ht_nr <- TD_norep %>% filter(!is.na(Height_year_3)) %>%
  group_by(Site, Treatment, Plot_number) %>%
  summarize(height_gini = ineq(Height_year_3, type = "Gini"), .groups = "drop") %>%
  filter(!is.nan(height_gini))
c_nr <- bio_nr %>% group_by(Site, Treatment, Plot_number) %>%
  summarise(C_Mg_ha = sum(AGB_Mg) / 0.0225 * 0.47, .groups = "drop")
cg_nr <- bio_nr %>% group_by(Site, Treatment, Plot_number) %>%
  summarize(carbon_gini = ineq(C_kg, type = "Gini"), .groups = "drop") %>% filter(!is.nan(carbon_gini))

pm_nr <- ht_nr2 %>%
  inner_join(c_nr, by = c("Site", "Treatment", "Plot_number")) %>%
  inner_join(gini_ht_nr, by = c("Site", "Treatment", "Plot_number")) %>%
  inner_join(cg_nr, by = c("Site", "Treatment", "Plot_number")) %>%
  left_join(cc_plot, by = c("Site", "Treatment", "Plot_number")) %>%
  left_join(gini_cc, by = c("Site", "Treatment", "Plot_number")) %>%
  filter(!is.na(canopy_gini), Treatment != "M2")

mat_nr <- pm_nr %>% select(mean_height, C_Mg_ha, canopy_cover, height_gini, canopy_gini, carbon_gini) %>% scale()
rich_nr <- richness_map[as.character(pm_nr$Treatment)]
bd_nr <- betadisper(dist(mat_nr), droplevels(pm_nr$Treatment))
ct_nr <- suppressWarnings(cor.test(rich_nr, bd_nr$distances, method = "spearman"))

# Compare
ct_wr <- suppressWarnings(cor.test(richness, bd$distances, method = "spearman"))
rep_compare <- data.frame(
  Analysis = c("Height Gini", "Height Gini", "Betadisper", "Betadisper"),
  Dataset = c("With replants", "Without replants", "With replants", "Without replants"),
  n_trees = c(nrow(TreeData), nrow(TD_norep), nrow(TreeData), nrow(TD_norep)),
  n_replants_excluded = c(0, sum(TreeData$is_replant), 0, sum(TreeData$is_replant)),
  Statistic = c(
    round(a_wr["Treatment", "F value"], 1),
    round(a_nr["Treatment", "F value"], 1),
    round(ct_wr$estimate, 3),
    round(ct_nr$estimate, 3)),
  p_value = c(a_wr["Treatment", "Pr(>F)"],
              a_nr["Treatment", "Pr(>F)"],
              ct_wr$p.value,
              ct_nr$p.value),
  Effect_size = c(round(eta_wr, 3), round(eta_nr, 3),
                  round(bd$group.distances["M1"], 2), round(bd_nr$group.distances["M1"], 2))
)

write.csv(rep_compare, "supp_replanting_sensitivity.csv", row.names = FALSE)
cat("Saved: supp_replanting_sensitivity.csv\n")
cat(sprintf("  Height Gini eta2p: %.3f (without replants)\n", eta_nr))
cat(sprintf("  Betadisper rho: %.3f, p = %.4f (without replants)\n", ct_nr$estimate, ct_nr$p.value))
cat(sprintf("  M1 disp: %.2f, 12SP disp: %.2f (without replants)\n",
    bd_nr$group.distances["M1"], bd_nr$group.distances["12SP"]))

# ============================================================
# F) MIXED-MODEL ROBUSTNESS — accounts for site-level clustering
# Addresses potential pseudoreplication in the plot-level Spearman test
# ============================================================
bd_df <- data.frame(
  distance = bd$distances,
  richness = richness,
  Site = pm_no_m2$Site,
  Block = pm_no_m2$Plot_number,
  Treatment = droplevels(pm_no_m2$Treatment)
)

mm_site <- lmer(distance ~ richness + (1 | Site), data = bd_df)
mm_summary <- summary(mm_site)
mm_coef <- coef(mm_summary)["richness", ]

mm_full <- lmer(distance ~ richness + (1 | Site / Block), data = bd_df)
mm_full_summary <- summary(mm_full)
mm_full_coef <- coef(mm_full_summary)["richness", ]

cat("\n=== MIXED-MODEL ROBUSTNESS (betadisper distances) ===\n")
cat(sprintf("  Model: distance ~ richness + (1|Site)\n"))
cat(sprintf("  Richness β = %.4f, t = %.2f, p = %.4f\n",
    mm_coef["Estimate"], mm_coef["t value"], mm_coef["Pr(>|t|)"]))
cat(sprintf("  Site variance: %.4f (SD = %.4f)\n",
    as.numeric(VarCorr(mm_site)$Site), sqrt(as.numeric(VarCorr(mm_site)$Site))))
cat(sprintf("\n  Model: distance ~ richness + (1|Site/Block)\n"))
cat(sprintf("  Richness β = %.4f, t = %.2f, p = %.4f\n",
    mm_full_coef["Estimate"], mm_full_coef["t value"], mm_full_coef["Pr(>|t|)"]))

mm_results <- data.frame(
  Model = c("distance ~ richness + (1|Site)",
            "distance ~ richness + (1|Site/Block)",
            "Spearman (no hierarchy)"),
  Beta_or_rho = c(round(mm_coef["Estimate"], 4),
                  round(mm_full_coef["Estimate"], 4),
                  round(obs_rho, 4)),
  t_or_S = c(round(mm_coef["t value"], 2),
             round(mm_full_coef["t value"], 2),
             NA),
  p_value = c(round(mm_coef["Pr(>|t|)"], 4),
              round(mm_full_coef["Pr(>|t|)"], 4),
              round(suppressWarnings(cor.test(richness, bd$distances, method = "spearman"))$p.value, 4)),
  n_plots = rep(nrow(bd_df), 3),
  n_sites = rep(length(unique(bd_df$Site)), 3)
)
write.csv(mm_results, "supp_mixed_model_robustness.csv", row.names = FALSE)
cat("Saved: supp_mixed_model_robustness.csv\n")

cat("\nAll sensitivity analyses done.\n")
