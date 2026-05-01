library(dplyr)
library(ggplot2)
library(lme4)
library(lmerTest)
library(glmmTMB)
library(emmeans)
library(multcomp)
library(ineq)
library(BIOMASS)
library(patchwork)

source("00_preprocess.R")

trt_colors <- c("M1" = "#D55E00", "M2" = "#E69F00",
                "2SP" = "darkseagreen2", "6SP" = "darkseagreen", "12SP" = "darkolivegreen")

CC3 <- CanopyCover %>% filter(Year == "3")

# ============================================================
# Helper: single pooled bar panel with CLD + effect size
# ============================================================
pooled_panel <- function(data, val_col, ylab, log_transform = FALSE) {
  data <- data %>% rename(val = !!val_col)

  if (log_transform) {
    m <- suppressWarnings(lmer(log(val) ~ Treatment * Site + (1 | Site:Plot_number), data = data))
  } else {
    m <- suppressWarnings(lmer(val ~ Treatment * Site + (1 | Site:Plot_number), data = data))
  }

  em <- emmeans(m, pairwise ~ Treatment)
  cld_df <- as.data.frame(cld(em$emmeans, Letters = letters, adjust = "tukey",
                               alpha = 0.05, decreasing = TRUE))

  aov_tbl <- anova(m)
  f_val <- aov_tbl["Treatment", "F value"]
  num_df <- aov_tbl["Treatment", "NumDF"]
  den_df <- aov_tbl["Treatment", "DenDF"]
  p_val <- aov_tbl["Treatment", "Pr(>F)"]
  eta2p <- f_val * num_df / (f_val * num_df + den_df)
  p_txt <- ifelse(p_val < 0.001, "< 0.001", sprintf("= %.3f", p_val))
  subtitle <- paste0("η²p = ", formatC(eta2p, format = "f", digits = 2), ", p ", p_txt)

  sm <- data %>%
    group_by(Treatment) %>%
    summarise(m = mean(val, na.rm = TRUE), s = sd(val, na.rm = TRUE), .groups = "drop")

  y_top <- max(c(max(data$val, na.rm = TRUE), max(sm$m + sm$s, na.rm = TRUE))) * 1.08
  cld_df$y_top <- y_top

  ggplot(sm, aes(x = Treatment, y = m, fill = Treatment)) +
    geom_bar(stat = "identity") +
    geom_jitter(data = data, aes(y = val),
                width = 0.2, size = 1.2, alpha = 0.4, color = "black") +
    geom_errorbar(aes(ymin = pmax(m - s, 0), ymax = m + s),
                  width = 0.2, linetype = "dashed", linewidth = 0.4) +
    geom_text(data = cld_df, aes(y = y_top, label = .group),
              size = 4.5, fontface = "bold", hjust = 0.5) +
    scale_fill_manual(values = trt_colors, guide = "none") +
    labs(x = NULL, y = ylab, subtitle = subtitle) +
    theme_classic() +
    theme(
      plot.subtitle = element_text(size = 10, color = "black"),
      axis.text.x = element_text(size = 11, face = "bold"),
      axis.text.y = element_text(size = 10),
      axis.title = element_text(size = 12)
    )
}

# ============================================================
# Fig 2: Height + Canopy cover + Carbon + Survival (pooled)
# ============================================================

# a) Mean height
ht <- TreeData %>% filter(!is.na(Height_year_3)) %>%
  group_by(Site, Treatment, Plot_number) %>%
  summarise(mean_height = mean(Height_year_3 / 100), .groups = "drop")
p_ht <- pooled_panel(ht, "mean_height", "Mean height (m)")

# b) Canopy cover
cc <- CC3 %>%
  group_by(Site, Treatment, Plot_number) %>%
  summarise(canopy_cover = mean(Percent_canopy_cover, na.rm = TRUE), .groups = "drop")
p_cc <- pooled_panel(cc, "canopy_cover", "Canopy cover (%)")

# c) Carbon (log-transformed model)
bio <- TreeData %>%
  filter(!is.na(DBH_year_3), !is.na(Height_year_3), Height_year_3 > 0, !is.na(WD))
bio$AGB_Mg <- computeAGB(D = bio$DBH_year_3, WD = bio$WD, H = bio$Height_year_3 / 100)
c_dat <- bio %>%
  group_by(Site, Treatment, Plot_number) %>%
  summarise(carbon = sum(AGB_Mg) / 0.0225 * 0.47, .groups = "drop")
p_carbon <- pooled_panel(c_dat, "carbon", expression("Carbon (Mg C " * ha^-1 * ")"),
                         log_transform = TRUE)

# d) Survival (binomial GLMM — custom panel)
surv_counts <- TreeData %>%
  group_by(Site, Treatment, Plot_number) %>%
  summarise(n_alive = sum(!is.na(Height_year_3)),
            n_total = n(), .groups = "drop") %>%
  mutate(n_dead = n_total - n_alive,
         survival = n_alive / n_total * 100)

surv_model <- glmmTMB(cbind(n_alive, n_dead) ~ Treatment * Site + (1 | Site:Plot_number),
                      data = surv_counts, family = binomial)

em_surv <- emmeans(surv_model, pairwise ~ Treatment, type = "response")
cld_surv_pooled <- as.data.frame(cld(em_surv$emmeans, Letters = letters, adjust = "tukey",
                                      alpha = 0.05, decreasing = TRUE))

aov_surv <- car::Anova(surv_model, type = 3)
chisq_trt <- aov_surv["Treatment", "Chisq"]
df_trt <- aov_surv["Treatment", "Df"]
p_surv_val <- aov_surv["Treatment", "Pr(>Chisq)"]
p_surv_txt <- ifelse(p_surv_val < 0.001, "< 0.001", sprintf("= %.3f", p_surv_val))
subtitle_surv <- paste0("χ²(", df_trt, ") = ",
                        formatC(chisq_trt, format = "f", digits = 1),
                        ", p ", p_surv_txt)

sm_surv <- surv_counts %>%
  group_by(Treatment) %>%
  summarise(m = mean(survival), s = sd(survival), .groups = "drop")

y_top_surv <- max(c(max(surv_counts$survival),
                    max(sm_surv$m + sm_surv$s, na.rm = TRUE))) * 1.08
cld_surv_pooled$y_top <- y_top_surv

p_surv <- ggplot(sm_surv, aes(x = Treatment, y = m, fill = Treatment)) +
  geom_bar(stat = "identity") +
  geom_jitter(data = surv_counts, aes(y = survival),
              width = 0.2, size = 1.2, alpha = 0.4, color = "black") +
  geom_errorbar(aes(ymin = pmax(m - s, 0), ymax = m + s),
                width = 0.2, linetype = "dashed", linewidth = 0.4) +
  geom_text(data = cld_surv_pooled, aes(y = y_top_surv, label = .group),
            size = 4.5, fontface = "bold", hjust = 0.5) +
  scale_fill_manual(values = trt_colors, guide = "none") +
  labs(x = NULL, y = "Survival (%)", subtitle = subtitle_surv) +
  theme_classic() +
  theme(
    plot.subtitle = element_text(size = 9, color = "grey30"),
    axis.text.x = element_text(size = 11, face = "bold"),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12)
  )

fig2 <- (p_ht + p_cc) / (p_carbon + p_surv) +
  plot_annotation(tag_levels = "a")

ggsave("fig2_height_canopy_carbon.png", plot = fig2, width = 10, height = 9, dpi = 300, bg = "white")
cat("Saved: fig2_height_canopy_carbon.png\n")

# ============================================================
# Fig 3: Height Gini + Canopy Gini + Carbon Gini (pooled)
# ============================================================

gini_ht <- TreeData %>% filter(!is.na(Height_year_3)) %>%
  group_by(Site, Treatment, Plot_number) %>%
  summarize(height_gini = ineq(Height_year_3, type = "Gini"), .groups = "drop") %>%
  filter(!is.nan(height_gini))
p_gini_ht <- pooled_panel(gini_ht, "height_gini", "Gini coefficient of height")

gini_cc <- CC3 %>%
  group_by(Site, Treatment, Plot_number) %>%
  summarize(canopy_gini = ineq(Percent_canopy_cover, type = "Gini"), .groups = "drop") %>%
  filter(!is.nan(canopy_gini))
p_gini_cc <- pooled_panel(gini_cc, "canopy_gini", "Gini coefficient of canopy cover")

bio$C_kg <- bio$AGB_Mg * 1000 * 0.47
carbon_gini <- bio %>%
  group_by(Site, Treatment, Plot_number) %>%
  summarize(carbon_gini = ineq(C_kg, type = "Gini"), .groups = "drop") %>%
  filter(!is.nan(carbon_gini))
p_gini_c <- pooled_panel(carbon_gini, "carbon_gini", "Gini coefficient of carbon")

fig3 <- p_gini_ht + p_gini_cc + p_gini_c +
  plot_layout(ncol = 3) +
  plot_annotation(tag_levels = "a")

ggsave("fig3_gini.png", plot = fig3, width = 14, height = 5, dpi = 300, bg = "white")
cat("Saved: fig3_gini.png\n")
