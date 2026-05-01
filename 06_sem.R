library(dplyr)
library(ggplot2)
library(lme4)
library(lmerTest)
library(ineq)
library(BIOMASS)
library(piecewiseSEM)
library(readxl)

source("00_preprocess.R")

# ============================================================
# Build plot-level metrics for SEM
# ============================================================
ht <- TreeData %>% filter(!is.na(Height_year_3)) %>%
  group_by(Site, Treatment, Plot_number) %>%
  summarise(mean_height = mean(Height_year_3), .groups = "drop")

bio <- TreeData %>%
  filter(!is.na(DBH_year_3), !is.na(Height_year_3), Height_year_3 > 0, !is.na(WD))
bio$AGB_Mg <- computeAGB(D = bio$DBH_year_3, WD = bio$WD, H = bio$Height_year_3 / 100)
c_plot <- bio %>%
  group_by(Site, Treatment, Plot_number) %>%
  summarise(C_Mg_ha = sum(AGB_Mg) / 0.0225 * 0.47, .groups = "drop")

gini_ht <- TreeData %>% filter(!is.na(Height_year_3)) %>%
  group_by(Site, Treatment, Plot_number) %>%
  summarize(height_gini = ineq(Height_year_3, type = "Gini"), .groups = "drop")

cc <- CanopyCover %>% filter(Year == "3") %>%
  group_by(Site, Treatment, Plot_number) %>%
  summarise(canopy_cover = mean(Percent_canopy_cover, na.rm = TRUE), .groups = "drop")

gini_cc <- CanopyCover %>% filter(Year == "3") %>%
  group_by(Site, Treatment, Plot_number) %>%
  summarize(canopy_gini = ineq(Percent_canopy_cover, type = "Gini"), .groups = "drop")

sem_data <- ht %>%
  inner_join(c_plot, by = c("Site", "Treatment", "Plot_number")) %>%
  inner_join(gini_ht, by = c("Site", "Treatment", "Plot_number")) %>%
  left_join(cc, by = c("Site", "Treatment", "Plot_number")) %>%
  left_join(gini_cc, by = c("Site", "Treatment", "Plot_number")) %>%
  filter(!is.na(Treatment), !is.na(canopy_gini), !is.na(C_Mg_ha)) %>%
  mutate(
    sp_richness = case_when(
      Treatment == "M1" ~ 1, Treatment == "M2" ~ 1,
      Treatment == "2SP" ~ 2, Treatment == "6SP" ~ 6, Treatment == "12SP" ~ 12
    ),
    log_richness = log(sp_richness)
  )

sem_df <- sem_data %>%
  select(Site, log_richness, mean_height, C_Mg_ha, canopy_cover,
         height_gini, canopy_gini) %>%
  mutate(Site = as.factor(Site)) %>%
  mutate(across(c(log_richness:canopy_gini), ~ as.numeric(scale(.)))) %>%
  as.data.frame()

# ============================================================
# Fit piecewise SEM — using log(richness) as predictor
# ============================================================
m1 <- lmer(height_gini ~ log_richness + (1 | Site), data = sem_df)
m2 <- lmer(C_Mg_ha ~ log_richness + (1 | Site), data = sem_df)
m3 <- lmer(canopy_cover ~ log_richness + height_gini + C_Mg_ha + (1 | Site), data = sem_df)
m4 <- lmer(canopy_gini ~ log_richness + canopy_cover + (1 | Site), data = sem_df)

sem_model <- psem(m1, m2, m3, m4)
sem_summary <- summary(sem_model, .progressBar = FALSE)
print(sem_summary)

# ============================================================
# SEM path diagram
# ============================================================
nodes <- data.frame(
  name = c("Planted\nspecies\nrichness", "Height\nGini", "Carbon\n(Mg C ha⁻¹)",
           "Canopy\ncover (%)", "Canopy\nGini"),
  x = c(0, 1.5, 1.5, 3, 4.5),
  y = c(0, 1, -1, 0, 0),
  r2 = c(NA, 0.53, 0.02, 0.66, 0.75))

coefs <- sem_summary$coefficients
get_beta <- function(resp, pred) {
  row <- coefs[coefs$Response == resp & coefs$Predictor == pred, ]
  list(beta = row$Std.Estimate, pval = row$P.Value)
}

b1 <- get_beta("height_gini", "log_richness")
b2 <- get_beta("C_Mg_ha", "log_richness")
b3 <- get_beta("canopy_cover", "log_richness")
b4 <- get_beta("canopy_cover", "height_gini")
b5 <- get_beta("canopy_cover", "C_Mg_ha")
b6 <- get_beta("canopy_gini", "canopy_cover")

edges <- data.frame(
  from_x = c(0, 0, 0, 1.5, 1.5, 3),
  from_y = c(0, 0, 0, 1, -1, 0),
  to_x = c(1.5, 1.5, 3, 3, 3, 4.5),
  to_y = c(1, -1, 0, 0, 0, 0),
  beta = c(b1$beta, b2$beta, b3$beta, b4$beta, b5$beta, b6$beta),
  pval = c(b1$pval, b2$pval, b3$pval, b4$pval, b5$pval, b6$pval))

shorten <- 0.38
edges <- edges %>%
  mutate(
    dx = to_x - from_x, dy = to_y - from_y,
    len = sqrt(dx^2 + dy^2),
    x1 = from_x + shorten * dx / len, y1 = from_y + shorten * dy / len,
    x2 = to_x - shorten * dx / len, y2 = to_y - shorten * dy / len,
    sig = ifelse(pval < 0.001, "***", ifelse(pval < 0.05, "*", "ns")),
    color = ifelse(pval >= 0.05, "black", ifelse(beta > 0, "#009E73", "#D55E00")),
    lwd = ifelse(pval >= 0.05, 0.4, abs(beta) * 2.5 + 0.5),
    ltype = ifelse(pval >= 0.05, "dashed", "solid"),
    label = paste0(round(beta, 2), sig),
    # Place labels at midpoint of the shortened arrow, offset perpendicular
    mid_x = (x1 + x2) / 2, mid_y = (y1 + y2) / 2,
    # Perpendicular offset for diagonal arrows
    perp_x = -dy / len, perp_y = dx / len,
    label_x = mid_x + perp_x * 0.18,
    label_y = mid_y + perp_y * 0.18)

fig5 <- ggplot() +
  geom_segment(data = edges,
    aes(x = x1, y = y1, xend = x2, yend = y2,
        linewidth = lwd, linetype = ltype),
    arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
    color = edges$color, lineend = "round", show.legend = FALSE) +
  scale_linewidth_identity() + scale_linetype_identity() +
  geom_label(data = edges,
    aes(x = label_x, y = label_y, label = label),
    size = 3.5, fill = "white", label.size = 0,
    label.padding = unit(0.15, "lines"),
    color = edges$color, fontface = "bold", family = "Helvetica") +
  geom_point(data = nodes, aes(x = x, y = y),
    size = 24, shape = 21, fill = "white", color = "grey30", stroke = 1.2) +
  geom_text(data = nodes, aes(x = x, y = y, label = name),
    size = 3.2, lineheight = 0.85, fontface = "bold", family = "Helvetica") +
  geom_text(data = nodes %>% filter(!is.na(r2)),
    aes(x = x, y = y - 0.45,
        label = paste0("R² = ", sprintf("%.2f", r2))),
    size = 3.1, color = "black", family = "Helvetica") +
  coord_cartesian(xlim = c(-0.9, 5.4), ylim = c(-1.8, 1.8)) +
  theme_void() +
  theme(plot.margin = margin(10, 10, 10, 10),
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA))

ggsave("fig5_sem.png", plot = fig5, width = 9, height = 5, dpi = 300, bg = "white")
cat("Saved: fig5_sem.png\n")
