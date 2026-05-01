library(dplyr)
library(lme4)
library(lmerTest)
library(glmmTMB)
library(emmeans)
library(ineq)
library(BIOMASS)
library(DHARMa)

source("00_preprocess.R")

CC3 <- CanopyCover %>% filter(Year == "3")

# ============================================================
# Prepare all plot-level metrics
# ============================================================
ht <- TreeData %>% filter(!is.na(Height_year_3)) %>%
  group_by(Site, Treatment, Plot_number) %>%
  summarise(mean_height = mean(Height_year_3 / 100), .groups = "drop")

cc <- CC3 %>% group_by(Site, Treatment, Plot_number) %>%
  summarise(canopy_cover = mean(Percent_canopy_cover, na.rm = TRUE), .groups = "drop")

bio <- TreeData %>%
  filter(!is.na(DBH_year_3), !is.na(Height_year_3), Height_year_3 > 0, !is.na(WD))
bio$AGB_Mg <- computeAGB(D = bio$DBH_year_3, WD = bio$WD, H = bio$Height_year_3 / 100)
bio$C_kg <- bio$AGB_Mg * 1000 * 0.47

c_dat <- bio %>% group_by(Site, Treatment, Plot_number) %>%
  summarise(carbon = sum(AGB_Mg) / 0.0225 * 0.47, .groups = "drop")

ba <- TreeData %>% filter(!is.na(BA_year_3)) %>%
  group_by(Site, Treatment, Plot_number) %>%
  summarise(total_BA = sum(BA_year_3, na.rm = TRUE), .groups = "drop")

gini_ht <- TreeData %>% filter(!is.na(Height_year_3)) %>%
  group_by(Site, Treatment, Plot_number) %>%
  summarize(height_gini = ineq(Height_year_3, type = "Gini"), .groups = "drop") %>%
  filter(!is.nan(height_gini))

gini_cc <- CC3 %>% group_by(Site, Treatment, Plot_number) %>%
  summarize(canopy_gini = ineq(Percent_canopy_cover, type = "Gini"), .groups = "drop") %>%
  filter(!is.nan(canopy_gini))

carbon_gini <- bio %>% group_by(Site, Treatment, Plot_number) %>%
  summarize(carbon_gini = ineq(C_kg, type = "Gini"), .groups = "drop") %>%
  filter(!is.nan(carbon_gini))

surv_counts <- TreeData %>% group_by(Site, Treatment, Plot_number) %>%
  summarise(n_alive = sum(!is.na(Height_year_3)), n_total = n(), .groups = "drop") %>%
  mutate(n_dead = n_total - n_alive, survival = n_alive / n_total * 100)

inga_ba <- TreeData %>% filter(Species == "Inga edulis", !is.na(BA_year_3)) %>%
  group_by(Site, Treatment, Plot_number) %>%
  summarise(inga_BA = mean(BA_year_3, na.rm = TRUE), .groups = "drop")

inga_ht <- TreeData %>% filter(Species == "Inga edulis", !is.na(Height_year_3)) %>%
  group_by(Site, Treatment, Plot_number) %>%
  summarise(inga_height = mean(Height_year_3 / 100, na.rm = TRUE), .groups = "drop")

# ============================================================
# TABLE 2: Main text ANOVA summary
# ============================================================
fit_and_extract <- function(data, val_col, metric_name, log_transform = FALSE) {
  data <- data %>% rename(val = !!val_col)
  if (log_transform) {
    m <- suppressWarnings(lmer(log(val) ~ Treatment * Site + (1 | Site:Plot_number), data = data))
  } else {
    m <- suppressWarnings(lmer(val ~ Treatment * Site + (1 | Site:Plot_number), data = data))
  }
  aov <- anova(m)
  terms <- rownames(aov)
  data.frame(
    Metric = c(metric_name, rep("", length(terms) - 1)),
    Term = terms,
    F_value = round(aov$`F value`, 2),
    NumDF = aov$NumDF,
    DenDF = round(aov$DenDF, 1),
    p_value = aov$`Pr(>F)`,
    eta2p = round(aov$`F value` * aov$NumDF / (aov$`F value` * aov$NumDF + aov$DenDF), 3),
    stringsAsFactors = FALSE
  )
}

table2 <- bind_rows(
  fit_and_extract(ht, "mean_height", "Mean height (m)"),
  fit_and_extract(cc, "canopy_cover", "Canopy cover (%)"),
  fit_and_extract(c_dat, "carbon", "Carbon (Mg C/ha)", log_transform = TRUE),
  fit_and_extract(gini_ht, "height_gini", "Height Gini"),
  fit_and_extract(gini_cc, "canopy_gini", "Canopy Gini"),
  fit_and_extract(carbon_gini, "carbon_gini", "Carbon Gini")
)

# Add survival (glmmTMB — Wald chi-square)
surv_m <- glmmTMB(cbind(n_alive, n_dead) ~ Treatment * Site + (1 | Site:Plot_number),
                  data = surv_counts, family = binomial)
surv_aov <- car::Anova(surv_m, type = 3)
surv_terms <- rownames(surv_aov)[rownames(surv_aov) != "(Intercept)"]
surv_rows <- data.frame(
  Metric = c("Survival (%)", rep("", length(surv_terms) - 1)),
  Term = surv_terms,
  F_value = round(surv_aov[surv_terms, "Chisq"] / surv_aov[surv_terms, "Df"], 2),
  NumDF = surv_aov[surv_terms, "Df"],
  DenDF = NA,
  p_value = surv_aov[surv_terms, "Pr(>Chisq)"],
  eta2p = NA,
  stringsAsFactors = FALSE
)
surv_rows$Term <- gsub("Treatment:Site", "Treatment × Site", surv_rows$Term)

table2$Term <- gsub("Treatment:Site", "Treatment × Site", table2$Term)
table2 <- bind_rows(table2, surv_rows)

table2$p_formatted <- ifelse(table2$p_value < 0.001, "< 0.001",
                              sprintf("%.3f", table2$p_value))
table2$Significance <- ifelse(table2$p_value < 0.001, "***",
                       ifelse(table2$p_value < 0.01, "**",
                       ifelse(table2$p_value < 0.05, "*", "")))

write.csv(table2, "table2_anova.csv", row.names = FALSE)
cat("Saved: table2_anova.csv\n")

# ============================================================
# TABLE S: Wood density sources
# ============================================================
wd_table <- TreeData %>%
  select(Species, WD, WD_source) %>%
  distinct() %>%
  filter(!is.na(WD)) %>%
  arrange(Species)

write.csv(wd_table, "supp_wood_density.csv", row.names = FALSE)
cat("Saved: supp_wood_density.csv\n")

# ============================================================
# TABLE S: SEM coefficients (both models)
# ============================================================
library(piecewiseSEM)

base_data <- ht %>%
  rename(mean_height_m = mean_height) %>%
  inner_join(c_dat %>% rename(C_Mg_ha = carbon), by = c("Site", "Treatment", "Plot_number")) %>%
  inner_join(gini_ht, by = c("Site", "Treatment", "Plot_number")) %>%
  left_join(cc %>% rename(canopy_cover = canopy_cover), by = c("Site", "Treatment", "Plot_number")) %>%
  left_join(gini_cc, by = c("Site", "Treatment", "Plot_number")) %>%
  filter(!is.na(canopy_gini), !is.na(C_Mg_ha)) %>%
  mutate(
    sp_richness = case_when(
      Treatment == "M1" ~ 1, Treatment == "M2" ~ 1,
      Treatment == "2SP" ~ 2, Treatment == "6SP" ~ 6, Treatment == "12SP" ~ 12),
    log_richness = log(sp_richness),
    inga_prop = case_when(
      Treatment == "M1" ~ 1, Treatment == "2SP" ~ 0.5,
      Treatment == "6SP" ~ 1/6, Treatment == "12SP" ~ 1/12, Treatment == "M2" ~ 0)
  )

# SEM B: log richness, no M2
sem_b_data <- base_data %>% filter(Treatment != "M2") %>%
  select(Site, log_richness, C_Mg_ha, canopy_cover, height_gini, canopy_gini) %>%
  mutate(Site = as.factor(Site)) %>%
  mutate(across(c(log_richness:canopy_gini), ~ as.numeric(scale(.)))) %>%
  as.data.frame()

m1b <- lmer(height_gini ~ log_richness + (1 | Site), data = sem_b_data)
m2b <- lmer(C_Mg_ha ~ log_richness + (1 | Site), data = sem_b_data)
m3b <- lmer(canopy_cover ~ log_richness + height_gini + C_Mg_ha + (1 | Site), data = sem_b_data)
m4b <- lmer(canopy_gini ~ log_richness + canopy_cover + (1 | Site), data = sem_b_data)
sem_b_summary <- summary(psem(m1b, m2b, m3b, m4b), .progressBar = FALSE)

sem_b_coefs <- sem_b_summary$coefficients
sem_b_coefs$Model <- "Main (log richness, no M2)"

# SEM A: Inga proportion, all treatments
sem_a_data <- base_data %>%
  select(Site, inga_prop, C_Mg_ha, canopy_cover, height_gini, canopy_gini) %>%
  mutate(Site = as.factor(Site)) %>%
  mutate(across(c(inga_prop:canopy_gini), ~ as.numeric(scale(.)))) %>%
  as.data.frame()

m1a <- lmer(height_gini ~ inga_prop + (1 | Site), data = sem_a_data)
m2a <- lmer(C_Mg_ha ~ inga_prop + (1 | Site), data = sem_a_data)
m3a <- lmer(canopy_cover ~ inga_prop + height_gini + C_Mg_ha + (1 | Site), data = sem_a_data)
m4a <- lmer(canopy_gini ~ inga_prop + canopy_cover + (1 | Site), data = sem_a_data)
sem_a_summary <- summary(psem(m1a, m2a, m3a, m4a), .progressBar = FALSE)

sem_a_coefs <- sem_a_summary$coefficients
sem_a_coefs$Model <- "Supplement (Inga proportion, all treatments)"

sem_table <- bind_rows(sem_b_coefs, sem_a_coefs)
write.csv(sem_table, "supp_sem_coefficients.csv", row.names = FALSE)
cat("Saved: supp_sem_coefficients.csv\n")

# ============================================================
# TABLE S: DHARMa diagnostics
# ============================================================
models_list <- list(
  list(d = ht, v = "mean_height", n = "Mean height (m)", log = FALSE),
  list(d = cc, v = "canopy_cover", n = "Canopy cover (%)", log = FALSE),
  list(d = c_dat, v = "carbon", n = "Carbon (Mg C/ha)", log = TRUE),
  list(d = ba, v = "total_BA", n = "Total BA (cm²)", log = TRUE),
  list(d = gini_ht, v = "height_gini", n = "Height Gini", log = FALSE),
  list(d = gini_cc, v = "canopy_gini", n = "Canopy Gini", log = FALSE),
  list(d = carbon_gini, v = "carbon_gini", n = "Carbon Gini", log = FALSE),
  list(d = inga_ba, v = "inga_BA", n = "Inga BA (cm²)", log = TRUE),
  list(d = inga_ht, v = "inga_height", n = "Inga height (m)", log = FALSE)
)

dharma_results <- data.frame()
for (ml in models_list) {
  d <- ml$d %>% rename(val = !!ml$v)
  if (ml$log) {
    m <- suppressWarnings(lmer(log(val) ~ Treatment * Site + (1 | Site:Plot_number), data = d))
  } else {
    m <- suppressWarnings(lmer(val ~ Treatment * Site + (1 | Site:Plot_number), data = d))
  }
  sim <- simulateResiduals(fittedModel = m, n = 1000, plot = FALSE)
  ks <- testUniformity(sim, plot = FALSE)
  disp <- testDispersion(sim, plot = FALSE)
  out <- testOutliers(sim, plot = FALSE)
  dharma_results <- rbind(dharma_results, data.frame(
    Metric = ml$n,
    Transform = ifelse(ml$log, "log", "none"),
    Family = "Gaussian",
    KS_p = round(ks$p.value, 4),
    Dispersion_p = round(disp$p.value, 4),
    Outlier_p = round(out$p.value, 4)
  ))
}

# Survival (glmmTMB)
sim_surv <- simulateResiduals(fittedModel = surv_m, n = 1000, plot = FALSE)
ks_s <- testUniformity(sim_surv, plot = FALSE)
disp_s <- testDispersion(sim_surv, plot = FALSE)
out_s <- testOutliers(sim_surv, plot = FALSE)
dharma_results <- rbind(dharma_results, data.frame(
  Metric = "Survival (%)", Transform = "logit link", Family = "Binomial",
  KS_p = round(ks_s$p.value, 4), Dispersion_p = round(disp_s$p.value, 4),
  Outlier_p = round(out_s$p.value, 4)
))

write.csv(dharma_results, "supp_dharma_diagnostics.csv", row.names = FALSE)
cat("Saved: supp_dharma_diagnostics.csv\n")

# ============================================================
# TABLE S: Per-site treatment means ± SD
# ============================================================
all_metrics <- ht %>% rename(val = mean_height) %>% mutate(Metric = "Mean height (m)") %>%
  bind_rows(cc %>% rename(val = canopy_cover) %>% mutate(Metric = "Canopy cover (%)")) %>%
  bind_rows(c_dat %>% rename(val = carbon) %>% mutate(Metric = "Carbon (Mg C/ha)")) %>%
  bind_rows(ba %>% rename(val = total_BA) %>% mutate(Metric = "Total BA (cm²)")) %>%
  bind_rows(gini_ht %>% rename(val = height_gini) %>% mutate(Metric = "Height Gini")) %>%
  bind_rows(gini_cc %>% rename(val = canopy_gini) %>% mutate(Metric = "Canopy Gini")) %>%
  bind_rows(carbon_gini %>% rename(val = carbon_gini) %>% mutate(Metric = "Carbon Gini")) %>%
  bind_rows(surv_counts %>% rename(val = survival) %>% mutate(Metric = "Survival (%)")) %>%
  bind_rows(inga_ba %>% rename(val = inga_BA) %>% mutate(Metric = "Inga BA (cm²)")) %>%
  bind_rows(inga_ht %>% rename(val = inga_height) %>% mutate(Metric = "Inga height (m)"))

site_means <- all_metrics %>%
  group_by(Metric, Site, Treatment) %>%
  summarise(Mean = round(mean(val, na.rm = TRUE), 3),
            SD = round(sd(val, na.rm = TRUE), 3),
            n = n(), .groups = "drop") %>%
  mutate(Mean_SD = paste0(Mean, " ± ", SD)) %>%
  arrange(Metric, Site, Treatment)

write.csv(site_means, "supp_site_treatment_means.csv", row.names = FALSE)
cat("Saved: supp_site_treatment_means.csv\n")

cat("\nAll tables done.\n")
