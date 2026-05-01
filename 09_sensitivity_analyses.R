library(dplyr)
library(ggplot2)
library(ineq)
library(BIOMASS)
library(vegan)
library(patchwork)

source("00_preprocess.R")
CC3 <- CanopyCover %>% filter(Year == "3")

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

bd <- betadisper(dist(mat_no_m2), pm_no_m2$Treatment)

# ============================================================
# A) LOO SENSITIVITY
# ============================================================
sites <- unique(pm_no_m2$Site)
loo_results <- data.frame()
for (drop_site in c("none", sites)) {
  if (drop_site == "none") {
    keep <- rep(TRUE, nrow(pm_no_m2))
  } else {
    keep <- pm_no_m2$Site != drop_site
  }
  bd_loo <- betadisper(dist(mat_no_m2[keep, ]), pm_no_m2$Treatment[keep])
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
# B) PERMUTATION LINEAR CONTRAST + NULL MODEL
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
# C) RISK-RETURN FRONTIER
# ============================================================
rr <- pm_no_m2 %>%
  group_by(Treatment) %>%
  summarise(
    Height_mean = mean(mean_height), Height_cv = sd(mean_height) / mean(mean_height) * 100,
    Carbon_mean = mean(C_Mg_ha), Carbon_cv = sd(C_Mg_ha) / mean(C_Mg_ha) * 100,
    Canopy_mean = mean(canopy_cover), Canopy_cv = sd(canopy_cover) / mean(canopy_cover) * 100,
    .groups = "drop")

p_rr_carbon <- ggplot(rr, aes(x = Carbon_cv, y = Carbon_mean, color = Treatment)) +
  geom_point(size = 5) +
  geom_text(aes(label = Treatment), vjust = -1.2, fontface = "bold", size = 3.5, show.legend = FALSE) +
  scale_color_manual(values = trt_colors, guide = "none") +
  labs(x = "Cross-site variability (CV %)", y = expression("Mean carbon (Mg C " * ha^-1 * ")"),
       subtitle = "Carbon") +
  theme_minimal() + theme(plot.subtitle = element_text(face = "bold"))

p_rr_height <- ggplot(rr, aes(x = Height_cv, y = Height_mean, color = Treatment)) +
  geom_point(size = 5) +
  geom_text(aes(label = Treatment), vjust = -1.2, fontface = "bold", size = 3.5, show.legend = FALSE) +
  scale_color_manual(values = trt_colors, guide = "none") +
  labs(x = "Cross-site variability (CV %)", y = "Mean height (m)",
       subtitle = "Height") +
  theme_minimal() + theme(plot.subtitle = element_text(face = "bold"))

p_rr_canopy <- ggplot(rr, aes(x = Canopy_cv, y = Canopy_mean, color = Treatment)) +
  geom_point(size = 5) +
  geom_text(aes(label = Treatment), vjust = -1.2, fontface = "bold", size = 3.5, show.legend = FALSE) +
  scale_color_manual(values = trt_colors, guide = "none") +
  labs(x = "Cross-site variability (CV %)", y = "Mean canopy cover (%)",
       subtitle = "Canopy cover") +
  theme_minimal() + theme(plot.subtitle = element_text(face = "bold"))

mv_rr <- data.frame(Treatment = names(bd$group.distances),
                     Dispersion = bd$group.distances,
                     Richness = richness_map[names(bd$group.distances)])

p_rr_mv <- ggplot(mv_rr, aes(x = Dispersion, y = Richness, color = Treatment)) +
  geom_point(size = 5) +
  geom_text(aes(label = Treatment), hjust = -0.3, fontface = "bold", size = 3.5, show.legend = FALSE) +
  scale_color_manual(values = trt_colors, guide = "none") +
  labs(x = "Multivariate dispersion\n(distance to centroid)", y = "Planted species richness",
       subtitle = "Portfolio effect") +
  theme_minimal() + theme(plot.subtitle = element_text(face = "bold"))

fig_rr <- (p_rr_carbon + p_rr_height) / (p_rr_canopy + p_rr_mv) +
  plot_annotation(tag_levels = "a")
ggsave("supp_risk_return.png", plot = fig_rr, width = 10, height = 8, dpi = 300, bg = "white")
cat("Saved: supp_risk_return.png\n")

# ============================================================
# D) TEMPORAL DISPERSION TRAJECTORY
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

  bd_yr <- betadisper(dist(yr_mat), yr_data$Treatment)
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

# Reshape for plotting
temp_long <- temporal_results %>%
  tidyr::pivot_longer(cols = c(M1_disp, SP2_disp, SP6_disp, SP12_disp),
                      names_to = "Treatment", values_to = "Dispersion") %>%
  mutate(Treatment = case_when(
    Treatment == "M1_disp" ~ "M1", Treatment == "SP2_disp" ~ "2SP",
    Treatment == "SP6_disp" ~ "6SP", Treatment == "SP12_disp" ~ "12SP"))

p_temporal <- ggplot(temp_long, aes(x = Year, y = Dispersion, color = Treatment)) +
  geom_line(linewidth = 1.2) + geom_point(size = 4) +
  scale_color_manual(values = trt_colors) +
  scale_x_continuous(breaks = 1:3, labels = paste("Year", 1:3)) +
  labs(x = "Year after planting", y = "Multivariate dispersion\n(distance to centroid)",
       subtitle = "Dispersion trajectory (FAB + FCEA + LLFS, height + height Gini)") +
  theme_minimal() +
  theme(plot.subtitle = element_text(face = "bold", size = 10))

ggsave("supp_temporal_dispersion.png", plot = p_temporal, width = 7, height = 5, dpi = 300, bg = "white")
cat("Saved: supp_temporal_dispersion.png\n")

cat("\nAll sensitivity analyses done.\n")
