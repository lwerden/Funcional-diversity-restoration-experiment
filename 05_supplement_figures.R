# Supplement figures: by-site bar plots with CLD and eta-squared
library(dplyr)
library(ggplot2)
library(lme4)
library(lmerTest)
library(glmmTMB)
library(emmeans)
library(multcomp)
library(ineq)
library(BIOMASS)
library(DHARMa)
library(patchwork)

source("01_preprocess.R")

trt_colors <- c("M1" = "#D55E00", "M2" = "#E69F00",
                "2SP" = "darkseagreen2", "6SP" = "darkseagreen", "12SP" = "darkolivegreen")

CC3 <- CanopyCover %>% filter(Year == "3")

# ============================================================
# Helper: 4-panel bar plot by site with CLD letters + stats
# ============================================================
site_bar <- function(data, val_col, ylab, filename, log_transform = FALSE) {
  data <- data %>% rename(val = !!val_col)

  if (log_transform) {
    m <- suppressWarnings(lmer(log(val) ~ Treatment * Site + (1 | Site:Plot_number), data = data))
  } else {
    m <- suppressWarnings(lmer(val ~ Treatment * Site + (1 | Site:Plot_number), data = data))
  }

  em <- emmeans(m, pairwise ~ Treatment | Site)
  cld_df <- as.data.frame(cld(em$emmeans,
    Letters = letters, adjust = "tukey", alpha = 0.05, decreasing = TRUE))

  sm <- data %>%
    group_by(Site, Treatment) %>%
    summarise(m = mean(val, na.rm = TRUE), s = sd(val, na.rm = TRUE), .groups = "drop")

  site_max <- data %>%
    group_by(Site) %>%
    summarise(y_max = max(val, na.rm = TRUE), .groups = "drop")
  site_err_max <- sm %>%
    group_by(Site) %>%
    summarise(y_err = max(m + s, na.rm = TRUE), .groups = "drop")
  site_tops <- site_max %>%
    left_join(site_err_max, by = "Site") %>%
    mutate(y_top = pmax(y_max, y_err) * 1.08)
  cld_df <- cld_df %>% left_join(site_tops %>% dplyr::select(Site, y_top), by = "Site")

  jt <- as.data.frame(joint_tests(m, by = "Site"))
  jt <- jt %>%
    mutate(
      eta2p = F.ratio * df1 / (F.ratio * df1 + df2),
      p_txt = ifelse(p.value < 0.001, "< 0.001", sprintf("= %.3f", p.value))
    )
  strip_labs <- setNames(
    paste0(jt$Site, "  (η²p = ",
           formatC(jt$eta2p, format = "f", digits = 2),
           ", p ", jt$p_txt, ")"),
    jt$Site)

  p <- ggplot(sm, aes(x = Treatment, y = m, fill = Treatment)) +
    geom_bar(stat = "identity") +
    geom_jitter(data = data, aes(y = val),
                width = 0.2, size = 1.2, alpha = 0.4, color = "black") +
    geom_errorbar(aes(ymin = pmax(m - s, 0), ymax = m + s),
                  width = 0.2, linetype = "dashed", linewidth = 0.4) +
    geom_text(data = cld_df, aes(y = y_top, label = .group),
              size = 4, fontface = "bold", hjust = 0.5) +
    scale_fill_manual(values = trt_colors, guide = "none") +
    facet_wrap(~ Site, ncol = 2, scales = "fixed",
               labeller = labeller(Site = strip_labs)) +
    labs(x = "Treatment", y = ylab) +
    theme_classic() +
    theme(
      strip.text = element_text(face = "bold", size = 9),
      axis.text.x = element_text(size = 11, face = "bold"),
      axis.text.y = element_text(size = 10),
      axis.title = element_text(size = 12)
    )

  ggsave(filename, plot = p, width = 8, height = 7, dpi = 300)
  cat("Saved:", filename, "\n")
}

# ============================================================
# Build all metrics
# ============================================================

# Mean height
ht <- TreeData %>% filter(!is.na(Height_year_3)) %>%
  group_by(Site, Treatment, Plot_number) %>%
  summarise(mean_height = mean(Height_year_3 / 100), .groups = "drop")
site_bar(ht, "mean_height", "Mean height (m)", "FigS3_height_by_site.jpeg")

# Canopy cover
cc <- CC3 %>%
  group_by(Site, Treatment, Plot_number) %>%
  summarise(canopy_cover = mean(Percent_canopy_cover, na.rm = TRUE), .groups = "drop")
site_bar(cc, "canopy_cover", "Canopy cover (%)", "FigS4_canopy_cover_by_site.jpeg")

# Carbon
bio <- TreeData %>%
  filter(!is.na(DBH_year_3), !is.na(Height_year_3), Height_year_3 > 0, !is.na(WD))
bio$AGB_Mg <- computeAGB(D = bio$DBH_year_3, WD = bio$WD, H = bio$Height_year_3 / 100)
c_dat <- bio %>%
  group_by(Site, Treatment, Plot_number) %>%
  summarise(carbon = sum(AGB_Mg) / 0.0225 * 0.47, .groups = "drop")
site_bar(c_dat, "carbon", expression("Carbon (Mg C " * ha^-1 * ")"), "FigS5_carbon_by_site.jpeg", log_transform = TRUE)

# Height Gini
gini_ht <- TreeData %>% filter(!is.na(Height_year_3)) %>%
  group_by(Site, Treatment, Plot_number) %>%
  summarize(height_gini = ineq(Height_year_3, type = "Gini"), .groups = "drop") %>%
  filter(!is.nan(height_gini))
site_bar(gini_ht, "height_gini", "Gini coefficient of height", "FigS7_height_gini_by_site.jpeg")

# Canopy Gini
gini_cc <- CC3 %>%
  group_by(Site, Treatment, Plot_number) %>%
  summarize(canopy_gini = ineq(Percent_canopy_cover, type = "Gini"), .groups = "drop") %>%
  filter(!is.nan(canopy_gini))
site_bar(gini_cc, "canopy_gini", "Gini coefficient of canopy cover", "FigS8_canopy_gini_by_site.jpeg")

# Carbon Gini
bio$C_kg <- bio$AGB_Mg * 1000 * 0.47
carbon_gini <- bio %>%
  group_by(Site, Treatment, Plot_number) %>%
  summarize(carbon_gini = ineq(C_kg, type = "Gini"), .groups = "drop") %>%
  filter(!is.nan(carbon_gini))
site_bar(carbon_gini, "carbon_gini", "Gini coefficient of carbon", "FigS9_carbon_gini_by_site.jpeg")

# Survival — binomial GLMM on plot-level alive/dead counts
surv_counts <- TreeData %>%
  group_by(Site, Treatment, Plot_number) %>%
  summarise(n_alive = sum(!is.na(Height_year_3)),
            n_total = n(), .groups = "drop") %>%
  mutate(n_dead = n_total - n_alive,
         survival = n_alive / n_total * 100)

surv_model <- glmmTMB(cbind(n_alive, n_dead) ~ Treatment * Site + (1 | Site:Plot_number),
                      data = surv_counts, family = binomial)

em_surv <- emmeans(surv_model, pairwise ~ Treatment | Site, type = "response")
cld_surv <- as.data.frame(cld(em_surv$emmeans, Letters = letters,
                               adjust = "tukey", alpha = 0.05, decreasing = TRUE))

sm_surv <- surv_counts %>%
  group_by(Site, Treatment) %>%
  summarise(m = mean(survival), s = sd(survival), .groups = "drop")

site_max_surv <- surv_counts %>%
  group_by(Site) %>%
  summarise(y_max = max(survival), .groups = "drop")
site_err_surv <- sm_surv %>%
  group_by(Site) %>%
  summarise(y_err = max(m + s, na.rm = TRUE), .groups = "drop")
site_tops_surv <- site_max_surv %>%
  left_join(site_err_surv, by = "Site") %>%
  mutate(y_top = pmax(y_max, y_err) * 1.08)
cld_surv <- cld_surv %>%
  left_join(site_tops_surv %>% dplyr::select(Site, y_top), by = "Site")

jt_surv <- as.data.frame(joint_tests(surv_model, by = "Site"))
jt_surv <- jt_surv %>%
  mutate(
    eta2p = F.ratio * df1 / (F.ratio * df1 + df2),
    p_txt = ifelse(p.value < 0.001, "< 0.001", sprintf("= %.3f", p.value))
  )
strip_labs_surv <- setNames(
  paste0(jt_surv$Site, "  (η²p = ",
         formatC(jt_surv$eta2p, format = "f", digits = 2),
         ", p ", jt_surv$p_txt, ")"),
  jt_surv$Site)

p_surv <- ggplot(sm_surv, aes(x = Treatment, y = m, fill = Treatment)) +
  geom_bar(stat = "identity") +
  geom_jitter(data = surv_counts, aes(y = survival),
              width = 0.2, size = 1.2, alpha = 0.4, color = "black") +
  geom_errorbar(aes(ymin = pmax(m - s, 0), ymax = m + s),
                width = 0.2, linetype = "dashed", linewidth = 0.4) +
  geom_text(data = cld_surv, aes(y = y_top, label = .group),
            size = 4, fontface = "bold", hjust = 0.5) +
  scale_fill_manual(values = trt_colors, guide = "none") +
  facet_wrap(~ Site, ncol = 2, scales = "fixed",
             labeller = labeller(Site = strip_labs_surv)) +
  labs(x = "Treatment", y = "Survival (%)") +
  theme_classic() +
  theme(
    strip.text = element_text(face = "bold", size = 9),
    axis.text.x = element_text(size = 11, face = "bold"),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12)
  )

ggsave("FigS6_survival_by_site.jpeg", plot = p_surv, width = 8, height = 7, dpi = 300)
cat("Saved: FigS6_survival_by_site.jpeg\n")

# Inga individual performance — combined BA + height
inga_ba <- TreeData %>%
  filter(Species == "Inga edulis", !is.na(BA_year_3)) %>%
  group_by(Site, Treatment, Plot_number) %>%
  summarise(inga_BA = mean(BA_year_3, na.rm = TRUE), .groups = "drop")

inga_ht <- TreeData %>%
  filter(Species == "Inga edulis", !is.na(Height_year_3)) %>%
  group_by(Site, Treatment, Plot_number) %>%
  summarise(inga_height = mean(Height_year_3 / 100, na.rm = TRUE), .groups = "drop")

make_inga_panel <- function(data, val_col, ylab, log_transform = FALSE) {
  data <- data %>% rename(val = !!val_col)
  if (log_transform) {
    m <- suppressWarnings(lmer(log(val) ~ Treatment * Site + (1 | Site:Plot_number), data = data))
  } else {
    m <- suppressWarnings(lmer(val ~ Treatment * Site + (1 | Site:Plot_number), data = data))
  }
  em <- emmeans(m, pairwise ~ Treatment | Site)
  cld_df <- as.data.frame(cld(em$emmeans, Letters = letters, adjust = "tukey", alpha = 0.05, decreasing = TRUE))
  sm <- data %>% group_by(Site, Treatment) %>%
    summarise(m = mean(val, na.rm = TRUE), s = sd(val, na.rm = TRUE), .groups = "drop")
  site_tops <- sm %>% group_by(Site) %>%
    summarise(y_top = max(m + s, na.rm = TRUE) * 1.15, .groups = "drop")
  cld_df <- cld_df %>% left_join(site_tops, by = "Site")

  ggplot(sm, aes(x = Treatment, y = m, fill = Treatment)) +
    geom_bar(stat = "identity") +
    geom_jitter(data = data, aes(y = val), width = 0.2, size = 1, alpha = 0.4, color = "black") +
    geom_errorbar(aes(ymin = pmax(m - s, 0), ymax = m + s), width = 0.2, linetype = "dashed", linewidth = 0.4) +
    geom_text(data = cld_df, aes(y = y_top, label = .group), size = 3.5, fontface = "bold") +
    scale_fill_manual(values = trt_colors, guide = "none") +
    facet_wrap(~ Site, ncol = 2, scales = "fixed") +
    labs(x = "Treatment", y = ylab) +
    theme_classic() +
    theme(strip.text = element_text(face = "bold", size = 9),
          axis.text.x = element_text(size = 10, face = "bold"),
          axis.text.y = element_text(size = 9),
          axis.title = element_text(size = 11))
}

p_ba <- make_inga_panel(inga_ba, "inga_BA", expression("Mean BA (" * cm^2 * ")"), log_transform = TRUE)
p_ht <- make_inga_panel(inga_ht, "inga_height", "Mean height (m)")
p_inga <- p_ba / p_ht + plot_annotation(tag_levels = "a")

ggsave("FigS10_inga_performance.png", plot = p_inga, width = 8, height = 12, dpi = 300, bg = "white")
cat("Saved: FigS10_inga_performance.png\n")

cat("\nAll supplement figures done.\n")
