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

Knepp_df <- read_csv("./Data/Processed_data/Knepp/Full_df_Knepp_LiDAR.csv")
head(Knepp_df)


rf_df <- Knepp_df |>
  select(
    Condition,
    ndvi_mean,
    evi_mean,
    lai_mean,
    Canopy_height_max_mean,
    Canopy_height_var_mean,
    Category
  ) |>
  drop_na()

rf_df$Category <- factor(rf_df$Category)
rf_df$Condition <- factor(rf_df$Condition)
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

rf_cv <- train(
  Condition ~ .,
  data = rf_df,
  method = "rf",
  trControl = control,
)

rf_cv


table(
  Predicted = rf_cv$pred$pred,
  Actual = rf_cv$pred$obs
)

# Confusion matrix
caret::confusionMatrix(
  data = rf_cv$pred$pred,
  reference = rf_cv$pred$obs
)

rf_df$Condition_binary <- ifelse(
  rf_df$Condition == "Poor",
  "Poor",
  "Not_Poor"
)

# Confusion matrix plot
cm_df <- as.data.frame(table(
  Predicted = rf_cv$pred$pred,
  Actual = rf_cv$pred$obs
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
  filename = "./Plots/Confusion_matrix/Knepp_Li_confusion_matrix.png",
  plot = con_plot,
  width = 8,
  height = 6,
  dpi = 300
)

# Confusion matrix by habitat
lookup <- rf_df |>
  mutate(rowIndex = dplyr::row_number()) |>
  select(rowIndex, Category)

best_preds_hab <- rf_cv$pred |>
  dplyr::filter(mtry == rf_cv$bestTune$mtry) |>
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

con_hab_plot <- ggplot(cm_df, aes(x = Actual, y = Predicted, fill = prop)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(prop, 2)), size = 4) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  facet_wrap(~ Category) +
  theme_minimal(base_size = 14) +
  labs(
    fill = "Proportion"
  ) +
  theme(
    axis.title = element_blank(),      # remove "Predicted" / "Actual"
    axis.text = element_text(size = 10),
    strip.text = element_text(size = 12, face = "bold")
  )
con_hab_plot
ggsave(
  filename = "./Plots/Confusion_matrix/Knepp_Li_habitat_confusion_matrix.png",
  plot = con_hab_plot,
  width = 8,
  height = 6,
  dpi = 300
)

# ROCs
options(scipen = 999)


best_preds <- rf_cv$pred |>
  dplyr::filter(mtry == rf_cv$bestTune$mtry)

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
      "Poor (AUC = 0.77)" = "#D55E00",
      "Moderate (AUC = 0.65)" = "#0072B2",
      "Good (AUC = 0.83)" = "#009E73"
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
  filename = "./Plots/ROC/Knepp_Li_ROC_curve.png",
  plot = roc_plot,
  width = 8,
  height = 6,
  dpi = 300
)

