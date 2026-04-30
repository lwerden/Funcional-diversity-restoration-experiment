library(readxl)
library(dplyr)

TreeData <- read_excel("FDiv_PlantedTreeData_final_excel.xlsx", sheet = "SurvivalGrowth")

TreeData <- TreeData %>%
  mutate(
    BA_year_2 = as.numeric(BA_year_2),
    BA_year_3 = as.numeric(BA_year_3),
    Height_year_1 = as.numeric(Height_year_1),
    Height_year_2 = as.numeric(Height_year_2),
    Height_year_3 = as.numeric(Height_year_3)
  )

# Fix species spelling errors
TreeData <- TreeData %>%
  mutate(Species = case_when(
    Species == "Astronium graveolons" ~ "Astronium graveolens",
    Species == "Aspoidosperam myristicofolium" ~ "Aspidosperma myristicifolium",
    Species == "Erythrina poepiggiana" ~ "Erythrina poeppigiana",
    Species == "Platymiscium cuerense" ~ "Platymiscium curuense",
    Species == "Licania platypus" ~ "Moquilea platypus",
    TRUE ~ Species
  ))

# Fix known data entry error: FCEA, 2SP, Plot 5, Tree 27, BA=604.28 -> 60.43
TreeData <- TreeData %>%
  mutate(BA_year_3 = ifelse(
    Site == "FCEA" & Treatment == "2SP" & Plot_number == 5 & Tree_number == 27,
    BA_year_3 / 10,
    BA_year_3
  ))

# BA = 0 in living trees means diameter was too small to measure — treat as NA
TreeData <- TreeData %>%
  mutate(
    BA_year_2 = ifelse(BA_year_2 == 0, NA, BA_year_2),
    BA_year_3 = ifelse(BA_year_3 == 0, NA, BA_year_3)
  )

# Back-calculate diameter (cm) from BA (cm2): BA = pi * (d/2)^2
TreeData <- TreeData %>%
  mutate(
    DBH_year_2 = 2 * sqrt(BA_year_2 / pi),
    DBH_year_3 = 2 * sqrt(BA_year_3 / pi)
  )

# Wood density table — best estimate from up to 5 sources:
#   BIOMASS R package (Rejou-Mechain et al.), BIEN database,
#   Zanne et al. 2009 Global Wood Density Database,
#   Tervuren Xylarium (TWDD), CIRAD wood density database.
# Species-level values averaged across sources; genus-level used where no
# species-level data available from any source.
wood_density <- data.frame(
  Species = c(
    "Astronium graveolens", "Calophyllum brasiliense", "Cedrela tonduzii",
    "Citharexylum cooperi", "Cojoba arborea", "Cordia alliodora",
    "Croton draco", "Croton schiedeanus", "Dendropanax ravenii",
    "Erythrina poeppigiana", "Ficus insipida", "Ficus maxima",
    "Garcinia madruno", "Handroanthus impetiginosus", "Handroanthus ochraceus",
    "Heliocarpus appendiculatus", "Inga edulis", "Lacistema aggregatum",
    "Moquilea platypus", "Ocotea puberula", "Persea caerulea",
    "Platymiscium curuense", "Quercus insignis", "Saurauia montana",
    "Simarouba amara", "Spondias mombin", "Terminalia amazonia",
    "Theobroma simiarum", "Viburnum costaricanum", "Zygia longifolia",
    "Aspidosperma myristicifolium",
    "Chrysophyllum cainito", "Cupania rufescens",
    "Sterculia recordiana", "Vitex cooperi"),
  WD = c(
    0.861, 0.571, 0.360,
    0.667, 0.683, 0.469,
    0.510, 0.510, 0.537,
    0.302, 0.370, 0.360,
    0.735, 0.831, 0.507,
    0.185, 0.580, 0.508,
    0.621, 0.433, 0.425,
    0.630, 0.701, 0.411,
    0.374, 0.386, 0.688,
    0.532, 0.631, 0.710,
    0.750,
    0.661, 0.608,
    0.490, 0.565),
  WD_source = c(
    "BIOMASS+BIEN+Zanne", "BIOMASS+BIEN+Zanne+TWDD+CIRAD", "BIOMASS+BIEN+Zanne",
    "BIOMASS (genus)", "BIEN (species)", "BIOMASS+BIEN+Zanne+TWDD+CIRAD",
    "BIOMASS (genus)", "BIOMASS (genus)", "BIEN (species)",
    "BIOMASS+BIEN+Zanne", "BIOMASS+BIEN+Zanne+CIRAD", "BIOMASS+BIEN+Zanne",
    "BIOMASS+BIEN+Zanne", "BIEN+TWDD", "TWDD (species)",
    "BIOMASS+BIEN+Zanne", "BIOMASS+BIEN+Zanne", "BIOMASS+BIEN+Zanne",
    "BIOMASS+BIEN+Zanne", "BIOMASS+BIEN+Zanne", "BIOMASS+BIEN+Zanne",
    "Roque et al. 2026 (CR plantations)", "BIOMASS (genus)", "BIOMASS (genus)",
    "BIOMASS+BIEN+Zanne+TWDD+CIRAD", "BIOMASS+BIEN+Zanne+TWDD+CIRAD", "BIOMASS+BIEN+Zanne+TWDD+CIRAD",
    "BIOMASS (genus)", "BIOMASS (genus)", "BIOMASS+BIEN+Zanne",
    "BIOMASS (genus)",
    "BIOMASS+BIEN+Zanne", "BIOMASS (genus)",
    "BIOMASS+BIEN+Zanne", "BIOMASS (genus)"),
  stringsAsFactors = FALSE
)

TreeData <- TreeData %>%
  left_join(wood_density, by = "Species")

# Summary
cat("=== Preprocessing summary ===\n")
cat("Total rows:", nrow(TreeData), "\n")
cat("BA_year_3:", sum(!is.na(TreeData$BA_year_3)), "measured,",
    sum(is.na(TreeData$BA_year_3)), "NA\n")
cat("DBH_year_3 range:", round(min(TreeData$DBH_year_3, na.rm = TRUE), 2), "—",
    round(max(TreeData$DBH_year_3, na.rm = TRUE), 2), "cm\n")
cat("Wood density matched:", sum(!is.na(TreeData$WD)), "of", nrow(TreeData), "rows\n")
cat("Species with WD:", length(unique(TreeData$Species[!is.na(TreeData$WD)])),
    "of", length(unique(TreeData$Species)), "\n")
cat("\nFixes applied:\n")
cat("  - Species spelling: Astronium graveolons -> graveolens, Aspoidosperam -> Aspidosperma\n")
cat("  - Outlier: FCEA 2SP Plot5 Tree27 BA 604.28 -> 60.43\n")
cat("  - BA = 0 converted to NA (unmeasurable diameter)\n")

saveRDS(TreeData, "data_cleaned.rds")
cat("\nSaved: data_cleaned.rds\n")
