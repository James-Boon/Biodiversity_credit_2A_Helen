
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

df <- read_csv("./Data/Processed_data/combined_df.csv")
head(df)

rf_df <- df |>
  select(
    Condition,
    ndvi_mean,
    evi_mean,
    ndvi_sd,
    evi_sd,
    Site,
    Category
  ) |>
  drop_na()|>
  mutate(
    Condition = case_when(
      Condition %in% c("Good", "Fairly good") ~ "Good",
      Condition == "Moderate" ~ "Moderate",
      Condition %in% c("Poor", "Fairly poor") ~ "Poor"
    )
  )


rf_df$Condition <- factor(rf_df$Condition)
rf_df$Category <- factor(rf_df$Category)
levels(rf_df$Condition)
rf_df$Condition <- factor(rf_df$Condition,
                          levels = c("Poor", "Moderate", "Good"))

set.seed(42)

# CROSS VALIDATION

control <- trainControl(
  method = "repeatedcv",
  number = 10,
  repeats = 5,
  savePredictions = "final",
  classProbs = TRUE
)

# Without habitat type
rf_no_context <- train(
  Condition ~ ndvi_mean + evi_mean + ndvi_sd + evi_sd,
  data = rf_df,
  method = "rf",
  trControl = control,
  importance = TRUE
)

# With habitat type
set.seed(42)

rf_with_habitat <- train(
  Condition ~ ndvi_mean + evi_mean + ndvi_sd + evi_sd + Category,
  data = rf_df,
  method = "rf",
  trControl = control,
  importance = TRUE
)

# Compare

rf_no_context
rf_with_habitat
results <- resamples(
  list(
    No_Context = rf_no_context,
    With_Habitat = rf_with_habitat
  )
)

summary(results)
bwplot(results, metric = "Accuracy")

# Confusion matrix heatmap
cm_df <- as.data.frame(table(
  Predicted = rf_with_habitat$pred$pred,
  Actual = rf_with_habitat$pred$obs
))

cm_df <- cm_df |>
  group_by(Actual) |>
  mutate(prop = Freq / sum(Freq))

con_plot <- ggplot(cm_df, aes(x = Actual, y = Predicted, fill = prop)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(prop, 2)), size = 5) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  theme_minimal(base_size = 14) +
  labs(
    x = "Actual Condition",
    y = "Predicted Condition",
    fill = "Proportion"
  )
con_plot

dir.create("./Plots/Confusion_matrix", showWarnings = FALSE)

ggsave(
  filename = "./Plots/Confusion_matrix/Knepp_Attenborough_confusion_matrix.png",
  plot = con_plot,
  width = 8,
  height = 6,
  dpi = 300
)


# Per habitat
lookup <- rf_df |>
  mutate(rowIndex = dplyr::row_number()) |>
  select(rowIndex, Category)

best_preds_hab <- rf_with_habitat$pred |>
  dplyr::filter(mtry == rf_with_habitat$bestTune$mtry) |>
  left_join(lookup, by = "rowIndex")

best_preds_hab$Category <- as.character(best_preds_hab$Category)

invisible(
  best_preds_hab |>
    dplyr::group_by(Category) |>
    dplyr::group_split() |>
    lapply(function(df) {
      cat("\nHabitat:", unique(df$Category), "\n")
      print(table(Predicted = df$pred, Actual = df$obs))
    })
)

# heat map

cm_df <- best_preds_hab |>
  count(Category, pred, obs) |>
  rename(Predicted = pred, Actual = obs)

cm_df <- best_preds_hab |>
  count(Category, pred, obs) |>
  rename(Predicted = pred, Actual = obs) |>
  complete(
    Category,
    Predicted = c("Poor", "Moderate", "Good"),
    Actual = c("Poor", "Moderate", "Good"),
    fill = list(n = 0)
  ) |>
  group_by(Category, Actual) |>
  mutate(prop = n / sum(n))

cm_df$Actual <- factor(cm_df$Actual, levels = c("Poor", "Moderate", "Good"))
cm_df$Predicted <- factor(cm_df$Predicted, levels = c("Poor", "Moderate", "Good"))

con_habitat_plot <- ggplot(cm_df, aes(x = Actual, y = Predicted, fill = prop)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(prop, 2)), size = 4) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  facet_wrap(~ Category) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Model Performance by Habitat Type",
    fill = "Proportion"
  ) +
  theme(
    axis.title = element_blank(),      # remove "Predicted" / "Actual"
    axis.text = element_text(size = 10),
    strip.text = element_text(size = 12, face = "bold")
  )

ggsave(
  filename = "./Plots/Confusion_matrix/Knepp_Attenborough_habitat_confusion_matrix.png",
  plot = con_habitat_plot,
  width = 8,
  height = 6,
  dpi = 300
)

# ROCs
options(scipen = 999)


best_preds <- rf_with_habitat$pred |>
  dplyr::filter(mtry == rf_with_habitat$bestTune$mtry)

roc_poor <- roc(
  response = best_preds$obs == "Poor",
  predictor = best_preds$Poor,
  quiet = TRUE
)

roc_moderate <- roc(
  response = best_preds$obs == "Moderate",
  predictor = best_preds$Moderate,
  quiet = TRUE
)

roc_good <- roc(
  response = best_preds$obs == "Good",
  predictor = best_preds$Good,
  quiet = TRUE
)

roc_df <- bind_rows(
  tibble(
    Specificity = roc_poor$specificities,
    Sensitivity = roc_poor$sensitivities,
    Class = paste0("Poor (AUC = ", round(auc(roc_poor), 2), ")")
  ),
  tibble(
    Specificity = roc_moderate$specificities,
    Sensitivity = roc_moderate$sensitivities,
    Class = paste0("Moderate (AUC = ", round(auc(roc_moderate), 2), ")")
  ),
  tibble(
    Specificity = roc_good$specificities,
    Sensitivity = roc_good$sensitivities,
    Class = paste0("Good (AUC = ", round(auc(roc_good), 2), ")")
  )
)

random_line <- data.frame(
  x = 1, y = 0,
  xend = 0, yend = 1
)

roc_plot <- ggplot(
  roc_df,
  aes(x = Specificity, y = Sensitivity, color = Class)
) +
  geom_line(linewidth = 1.2) +
  geom_segment(
    data = random_line,
    aes(x = x, y = y, xend = xend, yend = yend),
    inherit.aes = FALSE,
    color = "grey60",
    linewidth = 0.8
  ) +
  scale_x_reverse(
    limits = c(1, 0),
    breaks = seq(1, 0, by = -0.2)
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2)
  ) +
  scale_color_manual(
    values = c(
      "Poor (AUC = 0.74)" = "#D55E00",
      "Moderate (AUC = 0.64)" = "#0072B2",
      "Good (AUC = 0.78)" = "#009E73"
    )
  ) +
  labs(
    x = "Specificity",
    y = "Sensitivity",
    color = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = c(0.78, 0.18),
    legend.background = element_blank(),
    panel.grid.minor = element_blank()
  )

roc_plot

dir.create("./Plots/ROC", showWarnings = FALSE)

ggsave(
  filename = "./Plots/ROC/Knepp_att_RS_analysis_ROC.png",
  plot = roc_plot,
  width = 8,
  height = 6,
  dpi = 300
)

