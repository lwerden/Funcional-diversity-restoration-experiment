library(readxl)
library(dplyr)
library(writexl)

TreeData <- read_excel("FDiv_PlantedTreeData_final_excel.xlsx", sheet = "SurvivalGrowth")

TreeData <- TreeData %>%
  mutate(
    BA_year_2 = as.numeric(BA_year_2),
    BA_year_3 = as.numeric(BA_year_3),
    Height_year_1 = as.numeric(Height_year_1),
    Height_year_2 = as.numeric(Height_year_2),
    Height_year_3 = as.numeric(Height_year_3)
  )

# Fix known data entry error: FCEA, 2SP, Plot 5, Tree 27, BA=604.28 → 60.43
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

# Back-calculate diameter (cm) from BA (cm²): BA = pi * (d/2)^2
TreeData <- TreeData %>%
  mutate(
    DBH_year_2 = 2 * sqrt(BA_year_2 / pi),
    DBH_year_3 = 2 * sqrt(BA_year_3 / pi)
  )

# Summary of changes
cat("=== Preprocessing summary ===\n")
cat("Total rows:", nrow(TreeData), "\n")
cat("BA_year_3: ", sum(!is.na(TreeData$BA_year_3)), "measured,",
    sum(is.na(TreeData$BA_year_3)), "NA\n")
cat("DBH_year_3 range:", round(min(TreeData$DBH_year_3, na.rm = TRUE), 2), "—",
    round(max(TreeData$DBH_year_3, na.rm = TRUE), 2), "cm\n")
cat("Outlier fix applied: FCEA 2SP Plot5 Tree27 BA 604.28 → 60.43\n")
cat("BA = 0 converted to NA for year 2 and year 3\n")

saveRDS(TreeData, "data_cleaned.rds")
cat("\nSaved: data_cleaned.rds\n")
