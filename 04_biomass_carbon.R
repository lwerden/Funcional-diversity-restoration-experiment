library(dplyr)
library(ggplot2)
library(lme4)
library(lmerTest)
library(emmeans)
library(multcomp)
library(BIOMASS)

source("00_preprocess.R")

trt_colors <- c("M1" = "#D55E00", "M2" = "#E69F00",
                "2SP" = "darkseagreen2", "6SP" = "darkseagreen", "12SP" = "darkolivegreen")

# Estimate AGB using Chave et al. 2014 via BIOMASS package
# WD from 5 databases (see 00_preprocess.R), DBH back-calculated from BA
bio <- TreeData %>%
  filter(!is.na(DBH_year_3), !is.na(Height_year_3), Height_year_3 > 0, !is.na(WD))

bio$AGB_Mg <- computeAGB(D = bio$DBH_year_3, WD = bio$WD, H = bio$Height_year_3 / 100)

# Plot-level carbon (Mg C/ha, 0.47 carbon fraction)
plot_area_ha <- (15 * 15) / 10000

c_plot <- bio %>%
  group_by(Site, Treatment, Plot_number) %>%
  summarise(
    n_trees = n(),
    Total_AGB_Mg = sum(AGB_Mg, na.rm = TRUE),
    AGB_Mg_ha = Total_AGB_Mg / plot_area_ha,
    C_Mg_ha = AGB_Mg_ha * 0.47,
    .groups = "drop"
  )

cat("=== Carbon storage by treatment (Mg C/ha) ===\n")
c_plot %>%
  group_by(Treatment) %>%
  summarise(
    mean_C = round(mean(C_Mg_ha), 2),
    sd_C = round(sd(C_Mg_ha), 2),
    n_plots = n(),
    .groups = "drop"
  ) %>% print()

# Model
m_c <- lmer(C_Mg_ha ~ Treatment * Site + (1 | Site:Plot_number), data = c_plot)
cat("\n=== ANOVA ===\n")
print(anova(m_c))

cld_c <- as.data.frame(cld(emmeans(m_c, ~ Treatment),
  Letters = letters, adjust = "tukey", decreasing = TRUE))

# Pooled bar plot
sm <- c_plot %>%
  group_by(Treatment) %>%
  summarise(m = mean(C_Mg_ha), s = sd(C_Mg_ha), .groups = "drop")

carbon_plot <- ggplot(sm, aes(x = Treatment, y = m, fill = Treatment)) +
  geom_bar(stat = "identity") +
  geom_jitter(data = c_plot, aes(y = C_Mg_ha),
              width = 0.2, size = 1.5, alpha = 0.5, color = "black") +
  geom_errorbar(aes(ymin = pmax(m - s, 0), ymax = m + s),
                width = 0.2, linetype = "dashed", linewidth = 0.5) +
  geom_text(data = cld_c, aes(y = Inf, label = .group), vjust = 1, size = 3) +
  scale_fill_manual(values = trt_colors, guide = "none") +
  labs(x = "Treatment", y = expression("Carbon (Mg C " * ha^-1 * ")")) +
  theme_classic()

ggsave("carbon_pooled.jpeg", plot = carbon_plot, width = 6, height = 5, dpi = 300)

# By site
carbon_site <- ggplot(c_plot %>%
    group_by(Site, Treatment) %>%
    summarise(m = mean(C_Mg_ha), s = sd(C_Mg_ha), .groups = "drop"),
  aes(x = Treatment, y = m, fill = Treatment)) +
  geom_bar(stat = "identity") +
  geom_jitter(data = c_plot, aes(y = C_Mg_ha),
              width = 0.2, size = 1.5, alpha = 0.5, color = "black") +
  geom_errorbar(aes(ymin = pmax(m - s, 0), ymax = m + s),
                width = 0.2, linetype = "dashed", linewidth = 0.5) +
  scale_fill_manual(values = trt_colors, guide = "none") +
  facet_wrap(~ Site) +
  labs(x = "Treatment", y = expression("Carbon (Mg C " * ha^-1 * ")")) +
  theme_classic()

ggsave("carbon_by_site.jpeg", plot = carbon_site, width = 10, height = 6, dpi = 300)

cat("\nSaved: carbon_pooled.jpeg, carbon_by_site.jpeg\n")
