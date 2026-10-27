# =========
# CLEAN ENV
# =========
rm(list = ls())

# =========
# PACKAGES
# =========
libs <- c("tidyverse", "randomForest", "caret", "pROC")

installed_libs <- libs %in% rownames(installed.packages())
if (any(!installed_libs)) {
  install.packages(libs[!installed_libs])
}

invisible(lapply(libs, library, character.only = TRUE))

# =========
# LOAD DATA
# =========
df <- read_csv("./Data/Processed_data/combined_df.csv", show_col_types = FALSE) |>
  mutate(
    Category = recode(Category,
                      "Marsh,swamp, reedbed" = "Wetland"
    ),
    Condition = case_when(
      Condition %in% c("Good", "Fairly good") ~ "Good",
      Condition == "Moderate" ~ "Moderate",
      Condition %in% c("Poor", "Fairly poor") ~ "Poor",
      TRUE ~ NA_character_
    )
  )

# =========
# PREP DATA
# =========
analysis_df <- df |>
  select(
    Condition,
    Site,
    Category,
    ndvi_mean,
    evi_mean,
    ndvi_sd,
    evi_sd
  ) |>
  drop_na() |>
  mutate(
    Condition = factor(Condition, levels = c("Poor", "Moderate", "Good"), ordered = TRUE),
    Site = factor(Site),
    Category = factor(Category)
  )



combined_df <- analysis_df

combined_df$Site <- factor(combined_df$Site)  # Knepp / Attenborough

site_model <- train(
  Site ~ ndvi_mean + evi_mean + ndvi_sd + evi_sd,
  data = combined_df,
  method = "rf",
  trControl = trainControl(method = "cv", number = 10)
)

site_model
