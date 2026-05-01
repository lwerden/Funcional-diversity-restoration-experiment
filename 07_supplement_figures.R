library(dplyr)
library(ggplot2)
library(lme4)
library(lmerTest)
library(emmeans)
library(multcomp)
library(ineq)
library(BIOMASS)

source("00_preprocess.R")

trt_colors <- c("M1" = "#D55E00", "M2" = "#E69F00",
                "2SP" = "darkseagreen2", "6SP" = "darkseagreen", "12SP" = "darkolivegreen")

CC3 <- CanopyCover %>% filter(Year == "3")

# ============================================================
# Helper: 4-panel bar plot by site with CLD letters
# ============================================================
site_bar <- function(data, val_col, ylab, filename) {
  data <- data %>% rename(val = !!val_col)

  m <- suppressWarnings(lmer(val ~ Treatment * Site + (1 | Site:Plot_number), data = data))
  em <- emmeans(m, pairwise ~ Treatment | Site)
  cld_df <- as.data.frame(cld(em$emmeans,
    Letters = letters, adjust = "tukey", alpha = 0.05, decreasing = TRUE))

  sm <- data %>%
    group_by(Site, Treatment) %>%
    summarise(m = mean(val, na.rm = TRUE), s = sd(val, na.rm = TRUE), .groups = "drop")

  p <- ggplot(sm, aes(x = Treatment, y = m, fill = Treatment)) +
    geom_bar(stat = "identity") +
    geom_jitter(data = data, aes(y = val),
                width = 0.2, size = 1.2, alpha = 0.4, color = "black") +
    geom_errorbar(aes(ymin = pmax(m - s, 0), ymax = m + s),
                  width = 0.2, linetype = "dashed", linewidth = 0.4) +
    geom_text(data = cld_df, aes(y = Inf, label = .group),
              vjust = 1, size = 2.5) +
    scale_fill_manual(values = trt_colors, guide = "none") +
    facet_wrap(~ Site, ncol = 2) +
    labs(x = "Treatment", y = ylab) +
    theme_classic() +
    theme(strip.text = element_text(face = "bold"))

  ggsave(filename, plot = p, width = 8, height = 7, dpi = 300)
  cat("Saved:", filename, "\n")
}

# ============================================================
# Build all metrics
# ============================================================

# Mean height
ht <- TreeData %>% filter(!is.na(Height_year_3)) %>%
  group_by(Site, Treatment, Plot_number) %>%
  summarise(mean_height = mean(Height_year_3), .groups = "drop")
site_bar(ht, "mean_height", "Mean height (cm)", "supp_height_by_site.jpeg")

# Canopy cover
cc <- CC3 %>%
  group_by(Site, Treatment, Plot_number) %>%
  summarise(canopy_cover = mean(Percent_canopy_cover, na.rm = TRUE), .groups = "drop")
site_bar(cc, "canopy_cover", "Canopy cover (%)", "supp_canopy_cover_by_site.jpeg")

# Carbon
bio <- TreeData %>%
  filter(!is.na(DBH_year_3), !is.na(Height_year_3), Height_year_3 > 0, !is.na(WD))
bio$AGB_Mg <- computeAGB(D = bio$DBH_year_3, WD = bio$WD, H = bio$Height_year_3 / 100)
c_dat <- bio %>%
  group_by(Site, Treatment, Plot_number) %>%
  summarise(carbon = sum(AGB_Mg) / 0.0225 * 0.47, .groups = "drop")
site_bar(c_dat, "carbon", expression("Carbon (Mg C " * ha^-1 * ")"), "supp_carbon_by_site.jpeg")

# Total BA
ba <- TreeData %>% filter(!is.na(BA_year_3)) %>%
  group_by(Site, Treatment, Plot_number) %>%
  summarise(total_BA = sum(BA_year_3, na.rm = TRUE), .groups = "drop")
site_bar(ba, "total_BA", expression("Total basal area (" * cm^2 * ")"), "supp_BA_by_site.jpeg")

# Height Gini
gini_ht <- TreeData %>% filter(!is.na(Height_year_3)) %>%
  group_by(Site, Treatment, Plot_number) %>%
  summarize(height_gini = ineq(Height_year_3, type = "Gini"), .groups = "drop") %>%
  filter(!is.nan(height_gini))
site_bar(gini_ht, "height_gini", "Gini coefficient of height", "supp_gini_height_by_site.jpeg")

# Canopy Gini
gini_cc <- CC3 %>%
  group_by(Site, Treatment, Plot_number) %>%
  summarize(canopy_gini = ineq(Percent_canopy_cover, type = "Gini"), .groups = "drop") %>%
  filter(!is.nan(canopy_gini))
site_bar(gini_cc, "canopy_gini", "Gini coefficient of canopy cover", "supp_gini_canopy_by_site.jpeg")

# Survival
surv <- TreeData %>%
  group_by(Site, Treatment, Plot_number) %>%
  summarise(survival = sum(!is.na(Height_year_3)) / n() * 100, .groups = "drop")
site_bar(surv, "survival", "Survival (%)", "supp_survival_by_site.jpeg")

# Inga individual BA by site
inga_ba <- TreeData %>%
  filter(Species == "Inga edulis", !is.na(BA_year_3)) %>%
  group_by(Site, Treatment, Plot_number) %>%
  summarise(inga_BA = mean(BA_year_3, na.rm = TRUE), .groups = "drop")
site_bar(inga_ba, "inga_BA", expression("Inga edulis mean BA (" * cm^2 * ")"), "supp_inga_BA_by_site.jpeg")

# Inga individual height by site
inga_ht <- TreeData %>%
  filter(Species == "Inga edulis", !is.na(Height_year_3)) %>%
  group_by(Site, Treatment, Plot_number) %>%
  summarise(inga_height = mean(Height_year_3, na.rm = TRUE), .groups = "drop")
site_bar(inga_ht, "inga_height", "Inga edulis mean height (cm)", "supp_inga_height_by_site.jpeg")

cat("\nAll supplement figures done.\n")
