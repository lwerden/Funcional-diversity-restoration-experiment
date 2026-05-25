# Fig 4: radar plot, PCA, and betadisper multivariate dispersion analysis
library(dplyr)
library(ggplot2)
library(ggrepel)
library(ggtext)
library(ineq)
library(BIOMASS)
library(patchwork)
library(fmsb)
library(png)
library(grid)
library(vegan)

source("01_preprocess.R")

trt_colors_all <- c("M1" = "#D55E00", "M2" = "#E69F00",
                     "2SP" = "darkseagreen2", "6SP" = "darkseagreen", "12SP" = "darkolivegreen")
trt_colors_no_m2 <- c("M1" = "#D55E00",
                       "2SP" = "darkseagreen2", "6SP" = "darkseagreen", "12SP" = "darkolivegreen")

# ============================================================
# Compute plot-level metrics (7 variables for PCA and betadisper)
# ============================================================
# Mean height per plot (cm -> m)
ht_plot <- TreeData %>%
  filter(!is.na(Height_year_3)) %>%
  group_by(Site, Treatment, Plot_number) %>%
  summarise(mean_height = mean(Height_year_3 / 100), .groups = "drop")

# Biomass via Chave 2014 allometry, then carbon = AGB * 0.47
bio <- TreeData %>%
  filter(!is.na(DBH_year_3), !is.na(Height_year_3), Height_year_3 > 0, !is.na(WD))
bio$AGB_Mg <- computeAGB(D = bio$DBH_year_3, WD = bio$WD, H = bio$Height_year_3 / 100)
bio$C_kg <- bio$AGB_Mg * 1000 * 0.47

# Plot-level carbon scaled to Mg C per hectare (plot = 0.0225 ha)
c_plot <- bio %>%
  group_by(Site, Treatment, Plot_number) %>%
  summarise(C_Mg_ha = sum(AGB_Mg) / 0.0225 * 0.47, .groups = "drop")

# Within-plot heterogeneity: Gini coefficients (0 = uniform, 1 = max inequality)
gini_ht <- TreeData %>%
  filter(!is.na(Height_year_3)) %>%
  group_by(Site, Treatment, Plot_number) %>%
  summarize(height_gini = ineq(Height_year_3, type = "Gini"), .groups = "drop")

carbon_gini <- bio %>%
  group_by(Site, Treatment, Plot_number) %>%
  summarize(carbon_gini = ineq(C_kg, type = "Gini"), .groups = "drop") %>%
  filter(!is.nan(carbon_gini))

# Canopy cover: mean of 6 densiometer readings per plot
cc_plot <- CanopyCover %>%
  filter(Year == "3") %>%
  group_by(Site, Treatment, Plot_number) %>%
  summarise(canopy_cover = mean(Percent_canopy_cover, na.rm = TRUE), .groups = "drop")

gini_cc <- CanopyCover %>%
  filter(Year == "3") %>%
  group_by(Site, Treatment, Plot_number) %>%
  summarize(canopy_gini = ineq(Percent_canopy_cover, type = "Gini"), .groups = "drop") %>%
  filter(!is.nan(canopy_gini))

# Plot-level survival (proportion alive)
surv_plot <- TreeData %>%
  group_by(Site, Treatment, Plot_number) %>%
  summarise(survival = sum(!is.na(Height_year_3)) / n(), .groups = "drop")

# Join all 7 metrics into one plot-level dataframe
plot_metrics <- ht_plot %>%
  inner_join(c_plot, by = c("Site", "Treatment", "Plot_number")) %>%
  inner_join(gini_ht, by = c("Site", "Treatment", "Plot_number")) %>%
  inner_join(carbon_gini, by = c("Site", "Treatment", "Plot_number")) %>%
  left_join(cc_plot, by = c("Site", "Treatment", "Plot_number")) %>%
  left_join(gini_cc, by = c("Site", "Treatment", "Plot_number")) %>%
  left_join(surv_plot, by = c("Site", "Treatment", "Plot_number")) %>%
  filter(!is.na(canopy_gini), !is.na(survival))

# ============================================================
# RADAR PLOT helper
# ============================================================
make_radar <- function(data, colors, filename) {
  radar_data <- data %>%
    group_by(Treatment) %>%
    summarise(
      `Carbon\n(Mg C ha⁻¹)` = mean(C_Mg_ha),
      `Canopy\ncover (%)` = mean(canopy_cover),
      `Survival` = mean(survival),
      `Height\nGini` = mean(height_gini),
      `Carbon\nGini` = mean(carbon_gini),
      `Canopy\nGini` = mean(canopy_gini),
      `Mean\nheight (m)` = mean(mean_height),
      .groups = "drop")

  treatments <- as.character(radar_data$Treatment)
  radar_norm <- radar_data %>%
    select(-Treatment) %>%
    mutate(across(everything(), ~ . / max(.)))
  radar_fmsb <- rbind(rep(1, 7), rep(0, 7), radar_norm)

  png(filename, width = 6, height = 6, units = "in", res = 300)
  par(mar = c(1, 1.5, 2, 1.5), family = "Helvetica")
  radarchart(radar_fmsb, axistype = 1,
    pcol = unname(colors[treatments]),
    pfcol = NA,
    plwd = 2.5, plty = 1, cglcol = "grey75", cglty = 1, cglwd = 0.8,
    axislabcol = "grey30", vlcex = 0.95, calcex = 1.05,
    caxislabels = c("0", "0.25", "0.5", "0.75", "1"))
  legend("topright", legend = treatments, col = unname(colors[treatments]),
         lwd = 2.5, bty = "n", cex = 0.9)
  dev.off()
  cat("Saved:", filename, "\n")
}

# Main text radar (all treatments including M2)
make_radar(plot_metrics, trt_colors_all, "fig4a_radar.png")

# ============================================================
# PCA (all treatments including M2)
# ============================================================
pca_data <- plot_metrics

pca_mat <- pca_data %>%
  select(mean_height, C_Mg_ha, canopy_cover, height_gini, canopy_gini, carbon_gini, survival) %>%
  scale()
colnames(pca_mat) <- c("Mean height (m)", "Carbon\n(Mg C ha⁻¹)", "Canopy\ncover (%)",
                        "Height Gini", "Canopy Gini", "Carbon Gini", "Survival")
pca_result <- prcomp(pca_mat)
ve <- summary(pca_result)$importance[2, ]

scores <- as.data.frame(pca_result$x[, 1:2])
scores$Treatment <- pca_data$Treatment
scores$Site <- pca_data$Site

loadings <- as.data.frame(pca_result$rotation[, 1:2])
loadings$Variable <- rownames(loadings)

# Multivariate dispersion (betadisper) on full 7-variable space
mv_dist_all <- dist(pca_mat)
bd_all <- betadisper(mv_dist_all, pca_data$Treatment)
bd_perm_all <- permutest(bd_all, pairwise = TRUE)
pw_all <- bd_perm_all$pairwise$permuted
dists_all <- bd_all$group.distances
pair_names_all <- names(pw_all)

pw_table <- data.frame(
  Pair = pair_names_all,
  Distance_1 = sapply(pair_names_all, function(nm) dists_all[strsplit(nm, "-")[[1]][1]]),
  Distance_2 = sapply(pair_names_all, function(nm) dists_all[strsplit(nm, "-")[[1]][2]]),
  Difference = sapply(pair_names_all, function(nm) {
    parts <- strsplit(nm, "-")[[1]]
    dists_all[parts[1]] - dists_all[parts[2]]
  }),
  p_value = as.numeric(pw_all),
  Significant = ifelse(pw_all < 0.05, "*", ""),
  row.names = NULL
)
pw_table$Distance_1 <- round(pw_table$Distance_1, 3)
pw_table$Distance_2 <- round(pw_table$Distance_2, 3)
pw_table$Difference <- round(pw_table$Difference, 3)
write.csv(pw_table, "supp_betadisper_pairwise.csv", row.names = FALSE)
cat("Saved: supp_betadisper_pairwise.csv\n")

# Primary analysis 1: monoculture vs polyculture (all 5 treatments, no exclusions)
mono_poly <- ifelse(pca_data$Treatment %in% c("M1", "M2"), "monoculture", "polyculture")
bd_mp <- betadisper(mv_dist_all, mono_poly)
bd_mp_perm <- permutest(bd_mp)
mp_F <- round(bd_mp_perm$tab[1, "F"], 2)
mp_p <- bd_mp_perm$tab[1, "Pr(>F)"]
cat(sprintf("Mono vs poly: F = %s, p = %.3f\n", mp_F, mp_p))
cat(sprintf("  Mono dispersion = %.3f, Poly dispersion = %.3f\n",
    bd_mp$group.distances["monoculture"], bd_mp$group.distances["polyculture"]))

# Primary analysis 2: richness-dispersion trend along the designed gradient (M1→2SP→6SP→12SP)
# M2 is a parallel monoculture testing species identity, not part of the richness gradient
no_m2 <- pca_data$Treatment != "M2"
bd_no_m2 <- betadisper(dist(pca_mat[no_m2, ]), droplevels(pca_data$Treatment[no_m2]))
richness <- case_when(
  pca_data$Treatment[no_m2] == "M1" ~ 1,
  pca_data$Treatment[no_m2] == "2SP" ~ 2,
  pca_data$Treatment[no_m2] == "6SP" ~ 6,
  pca_data$Treatment[no_m2] == "12SP" ~ 12
)
trend_test <- cor.test(richness, bd_no_m2$distances, method = "spearman")
rho <- round(trend_test$estimate, 2)
trend_p <- trend_test$p.value

# Pairwise contrasts along the richness gradient
bd_perm_no_m2 <- permutest(bd_no_m2, pairwise = TRUE)
pw_no_m2 <- bd_perm_no_m2$pairwise$permuted
m1_12sp_p <- pw_no_m2["M1-12SP"]
dists_no_m2 <- bd_no_m2$group.distances

trend_p_txt <- ifelse(trend_p < 0.001, "< 0.001", formatC(trend_p, format = "f", digits = 3))
m1_12sp_p_txt <- formatC(m1_12sp_p, format = "f", digits = 3)
disp_pct <- round((dists_no_m2["M1"] / dists_no_m2["12SP"] - 1) * 100)

disp_html <- paste0(
  "<span style='font-family:Helvetica;font-size:11pt'>",
  "<b>Multivariate dispersion</b><br>",
  "M1 & M2 vs 2SP, 6SP & 12SP: <i>F</i> = ", mp_F,
  ", <i>p</i> = ", formatC(mp_p, format = "f", digits = 3), "<br>",
  "Dispersion decreases with richness: ρ = ", rho,
  ", <i>p</i> ", trend_p_txt, "<br>",
  "M1 ", disp_pct, "% more dispersed than 12SP (<i>p</i> = ", m1_12sp_p_txt, ")",
  "</span>"
)

pca_plot <- ggplot(scores, aes(x = PC1, y = PC2)) +
  # Crosshairs at 0,0
  geom_hline(yintercept = 0, color = "grey80", linewidth = 0.3) +
  geom_vline(xintercept = 0, color = "grey80", linewidth = 0.3) +
  # Ellipses
  stat_ellipse(aes(color = Treatment, fill = Treatment),
               geom = "polygon", alpha = 0.05, level = 0.95,
               linewidth = 0.5, linetype = "dashed") +
  # Points
  geom_point(aes(color = Treatment, shape = Site), size = 2, alpha = 0.6) +
  # Loading arrows
  geom_segment(data = loadings,
               aes(x = 0, y = 0, xend = PC1 * 3, yend = PC2 * 3),
               arrow = arrow(length = unit(0.2, "cm")),
               color = "grey30", linewidth = 0.5) +
  # Loading labels — not italic, black
  geom_text_repel(data = loadings,
                  aes(x = PC1 * 3.3, y = PC2 * 3.3, label = Variable),
                  size = 3.7, color = "black", fontface = "bold",
                  max.overlaps = 15, family = "Helvetica") +
  scale_color_manual(values = trt_colors_all, name = "Treatment") +
  scale_fill_manual(values = trt_colors_all, guide = "none") +
  scale_shape_manual(values = c("FAB" = 16, "FCEA" = 17, "GAM" = 15, "LLFS" = 3)) +
  geom_richtext(data = data.frame(x = -Inf, y = Inf, label = disp_html),
                aes(x = x, y = y, label = label),
                hjust = 0, vjust = 1, size = 3.5, color = "grey20",
                fill = alpha("white", 0.9), label.color = "grey50",
                label.padding = unit(c(5, 6, 5, 6), "pt"),
                label.r = unit(0, "pt"),
                family = "Helvetica", inherit.aes = FALSE) +
  labs(x = paste0("Higher productivity ←              PC1 (", round(ve[1] * 100, 1), "%)              → Canopy heterogeneity"),
       y = paste0("→ Height & carbon heterogeneity              PC2 (", round(ve[2] * 100, 1), "%)              Structural uniformity ←"),
       shape = "Site", tag = "b") +
  guides(
    color = guide_legend(override.aes = list(shape = 15, size = 5, alpha = 1,
                                              color = unname(trt_colors_all)),
                         nrow = 1, title = "Treatment"),
    shape = guide_legend(override.aes = list(size = 3), nrow = 1, title = "Site")
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank()
  )

# ============================================================
# COMBINED FIG 4
# ============================================================
radar_gg <- ggplot() +
  annotation_custom(rasterGrob(readPNG("fig4a_radar.png"), interpolate = TRUE)) +
  theme_void() + labs(tag = "a")

# PCA with site legend inside upper right, treatment legend inside upper left
pca_final <- pca_plot +
  guides(fill = "none",
         color = guide_legend(override.aes = list(shape = 15, size = 4, alpha = 1),
                              title = "Treatment"),
         shape = guide_legend(override.aes = list(size = 3), title = "Site")) +
  theme(
    legend.position = c(0.98, 0.98),
    legend.justification = c(1, 1),
    legend.box = "vertical",
    legend.background = element_rect(fill = alpha("white", 0.8), color = NA),
    legend.key.size = unit(0.35, "cm"),
    legend.key = element_rect(fill = NA, color = NA),
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 10, face = "bold"),
    legend.spacing.y = unit(0.1, "cm")
  )

fig4 <- radar_gg + pca_final + plot_layout(widths = c(1, 1.2))
ggsave("Fig4_multivariate_synthesis.png", plot = fig4, width = 13, height = 7, dpi = 300, bg = "white")
file.remove("fig4a_radar.png")
cat("Saved: Fig4_multivariate_synthesis.png\n")
