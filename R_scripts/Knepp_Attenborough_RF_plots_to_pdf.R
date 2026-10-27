# # =========
# # CLEAN ENV
# # =========
# rm(list = ls())
# 
# # =========
# # PACKAGES
# # =========
# libs <- c("tidyverse", "randomForest", "caret", "pROC")
# 
# installed_libs <- libs %in% rownames(installed.packages())
# if (any(!installed_libs)) {
#   install.packages(libs[!installed_libs])
# # }
# 
# invisible(lapply(libs, library, character.only = TRUE))
# 
# # =========
# # SETTINGS
# # =========
# set.seed(42)
# 
# dir.create("./Plots", showWarnings = FALSE)
# dir.create("./Plots/Confusion_matrix", showWarnings = FALSE)
# dir.create("./Plots/ROC", showWarnings = FALSE)
# 
# # =========
# # HELPERS
# # =========
# clean_condition <- function(x) {
#   case_when(
#     x %in% c("Good", "Fairly good") ~ "Good",
#     x == "Moderate" ~ "Moderate",
#     x %in% c("Poor", "Fairly poor") ~ "Poor",
#     TRUE ~ NA_character_
#   )
# }
# 
# make_overall_cm_df <- function(pred_df) {
#   cm <- table(
#     Predicted = factor(pred_df$pred, levels = c("Poor", "Moderate", "Good")),
#     Actual = factor(pred_df$obs, levels = c("Poor", "Moderate", "Good"))
#   )
# 
#   cm_prop <- prop.table(cm, margin = 2)
# 
#   as.data.frame(cm_prop) |>
#     rename(prop = Freq) |>
#     mutate(
#       prop = ifelse(is.nan(prop), NA_real_, prop),
#       Actual = factor(Actual, levels = c("Poor", "Moderate", "Good")),
#       Predicted = factor(Predicted, levels = c("Poor", "Moderate", "Good"))
#     )
# }
# 
# make_habitat_cm_df <- function(pred_df) {
#   split(pred_df, pred_df$Category) |>
#     purrr::imap_dfr(function(d, cat_name) {
#       cm <- table(
#         Predicted = factor(d$pred, levels = c("Poor", "Moderate", "Good")),
#         Actual = factor(d$obs, levels = c("Poor", "Moderate", "Good"))
#       )
# 
#       cm_prop <- prop.table(cm, margin = 2)
# 
#       out <- as.data.frame(cm_prop)
#       names(out) <- c("Predicted", "Actual", "prop")
#       out$Category <- cat_name
#       out
#     }) |>
#     mutate(
#       prop = ifelse(is.nan(prop), NA_real_, prop),
#       Actual = factor(Actual, levels = c("Poor", "Moderate", "Good")),
#       Predicted = factor(Predicted, levels = c("Poor", "Moderate", "Good"))
#     )
# }
# 
# plot_overall_cm <- function(cm_df) {
#   ggplot(cm_df, aes(x = Actual, y = Predicted, fill = prop)) +
#     geom_tile(color = "white") +
#     geom_text(aes(label = ifelse(is.na(prop), "", sprintf("%.2f", prop))), size = 5) +
#     scale_fill_gradient(low = "white", high = "steelblue", limits = c(0, 1), na.value = "white") +
#     labs(
#       x = "Actual Condition",
#       y = "Predicted Condition",
#       fill = "Proportion"
#     ) +
#     theme_minimal(base_size = 14) +
#     theme(panel.grid = element_blank())
# }
# 
# plot_habitat_cm <- function(cm_df) {
#   ggplot(cm_df, aes(x = Actual, y = Predicted, fill = prop)) +
#     geom_tile(color = "white") +
#     geom_text(aes(label = ifelse(is.na(prop), "", sprintf("%.2f", prop))), size = 4) +
#     scale_fill_gradient(low = "white", high = "steelblue", limits = c(0, 1), na.value = "white") +
#     facet_wrap(~ Category) +
#     labs(fill = "Proportion") +
#     theme_minimal(base_size = 14) +
#     theme(
#       axis.title = element_blank(),
#       axis.text = element_text(size = 10),
#       strip.text = element_text(size = 12, face = "bold"),
#       panel.grid = element_blank()
#     )
# }
# 
# make_roc_plot <- function(best_preds, filename = NULL) {
#   roc_poor <- roc(best_preds$obs == "Poor", best_preds$Poor, quiet = TRUE)
#   roc_moderate <- roc(best_preds$obs == "Moderate", best_preds$Moderate, quiet = TRUE)
#   roc_good <- roc(best_preds$obs == "Good", best_preds$Good, quiet = TRUE)
# 
#   roc_df <- bind_rows(
#     tibble(
#       Specificity = roc_poor$specificities,
#       Sensitivity = roc_poor$sensitivities,
#       Class = "Poor"
#     ),
#     tibble(
#       Specificity = roc_moderate$specificities,
#       Sensitivity = roc_moderate$sensitivities,
#       Class = "Moderate"
#     ),
#     tibble(
#       Specificity = roc_good$specificities,
#       Sensitivity = roc_good$sensitivities,
#       Class = "Good"
#     )
#   )
# 
#   auc_labels <- c(
#     Poor = paste0("Poor (AUC = ", round(auc(roc_poor), 2), ")"),
#     Moderate = paste0("Moderate (AUC = ", round(auc(roc_moderate), 2), ")"),
#     Good = paste0("Good (AUC = ", round(auc(roc_good), 2), ")")
#   )
# 
#   random_line <- data.frame(x = 1, y = 0, xend = 0, yend = 1)
# 
#   p <- ggplot(roc_df, aes(x = Specificity, y = Sensitivity, color = Class)) +
#     geom_line(linewidth = 1.2) +
#     geom_segment(
#       data = random_line,
#       aes(x = x, y = y, xend = xend, yend = yend),
#       inherit.aes = FALSE,
#       color = "grey60",
#       linewidth = 0.8
#     ) +
#     scale_x_reverse(limits = c(1, 0), breaks = seq(1, 0, by = -0.2)) +
#     scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
#     scale_color_manual(
#       values = c(Poor = "#D55E00", Moderate = "#0072B2", Good = "#009E73"),
#       breaks = c("Poor", "Moderate", "Good"),
#       labels = auc_labels
#     ) +
#     labs(
#       x = "Specificity",
#       y = "Sensitivity",
#       color = NULL
#     ) +
#     theme_minimal(base_size = 14) +
#     theme(
#       legend.position = c(0.78, 0.18),
#       legend.background = element_blank(),
#       panel.grid.minor = element_blank()
#     )
# 
#   if (!is.null(filename)) {
#     ggsave(filename = filename, plot = p, width = 8, height = 6, dpi = 300)
#   }
# 
#   p
# }
# 
# # =========
# # LOAD DATA
# # =========
# df <- read_csv("./Data/Processed_data/combined_df.csv", show_col_types = FALSE)
# 
# rf_df <- df |>
#   select(
#     Condition,
#     ndvi_mean,
#     evi_mean,
#     ndvi_sd,
#     evi_sd,
#     Category
#   ) |>
#   mutate(Condition = clean_condition(as.character(Condition))) |>
#   drop_na() |>
#   mutate(
#     Condition = factor(Condition, levels = c("Poor", "Moderate", "Good")),
#     Category = factor(Category)
#   )
# 
# # =========
# # CROSS VALIDATION
# # =========
# cv_index <- createMultiFolds(rf_df$Condition, k = 10, times = 5)
# 
# control <- trainControl(
#   method = "repeatedcv",
#   number = 10,
#   repeats = 5,
#   index = cv_index,
#   savePredictions = "final",
#   classProbs = TRUE
# )
# 
# # =========
# # MODEL
# # =========
# rf_cv <- train(
#   Condition ~ ndvi_mean + evi_mean + ndvi_sd + evi_sd,
#   data = rf_df,
#   method = "rf",
#   trControl = control,
#   importance = TRUE
# )
# 
# print(rf_cv)
# 
# # =========
# # OVERALL CONFUSION MATRIX
# # =========
# best_preds_overall <- rf_cv$pred |>
#   filter(mtry == rf_cv$bestTune$mtry)
# 
# print(table(Predicted = best_preds_overall$pred, Actual = best_preds_overall$obs))
# print(confusionMatrix(data = best_preds_overall$pred, reference = best_preds_overall$obs))
# 
# cm_df_overall <- make_overall_cm_df(best_preds_overall)
# con_plot <- plot_overall_cm(cm_df_overall)
# 
# ggsave(
#   filename = "./Plots/Confusion_matrix/Knepp_Attenborough_confusion_matrix.pdf",
#   plot = con_plot,
#   width = 8,
#   height = 6,
#   dpi = 300
# )
# 
# # =========
# # CONFUSION MATRIX BY HABITAT
# # =========
# lookup <- rf_df |>
#   mutate(rowIndex = row_number()) |>
#   select(rowIndex, Category)
# 
# best_preds_hab <- rf_cv$pred |>
#   filter(mtry == rf_cv$bestTune$mtry) |>
#   left_join(lookup, by = "rowIndex") |>
#   mutate(Category = as.character(Category))
# 
# invisible(
#   best_preds_hab |>
#     group_by(Category) |>
#     group_split() |>
#     lapply(function(x) {
#       cat("\nHabitat:", unique(x$Category), "\n")
#       print(table(Predicted = x$pred, Actual = x$obs))
#     })
# )
# 
# cm_df_habitat <- make_habitat_cm_df(best_preds_hab)
# con_habitat_plot <- plot_habitat_cm(cm_df_habitat)
# 
# ggsave(
#   filename = "./Plots/Confusion_matrix/Knepp_Attenborough_habitat_confusion_matrix.pdf",
#   plot = con_habitat_plot,
#   width = 10,
#   height = 7,
#   dpi = 300
# )
# 
# # =========
# # ROC CURVES
# # =========
# roc_plot <- make_roc_plot(
#   best_preds = best_preds_overall,
#   filename = "./Plots/ROC/Knepp_Attenborough_RS_ROC.pdf"
# )

# print(roc_plot)
  
  
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
# SETTINGS
# =========
set.seed(42)

dir.create("./Plots", showWarnings = FALSE)
dir.create("./Plots/Confusion_matrix", showWarnings = FALSE)
dir.create("./Plots/Variable_importance", showWarnings = FALSE)
dir.create("./Plots/ROC", showWarnings = FALSE)

# =========
# HELPERS
# =========
clean_condition <- function(x) {
  case_when(
    x %in% c("Good", "Fairly good") ~ "Good",
    x == "Moderate" ~ "Moderate",
    x %in% c("Poor", "Fairly poor") ~ "Poor",
    TRUE ~ NA_character_
  )
}

make_overall_cm_df <- function(pred_df) {
  cm <- table(
    Predicted = factor(pred_df$pred, levels = c("Poor", "Moderate", "Good")),
    Actual = factor(pred_df$obs, levels = c("Poor", "Moderate", "Good"))
  )

  cm_prop <- prop.table(cm, margin = 2)

  as.data.frame(cm_prop) |>
    rename(prop = Freq) |>
    mutate(
      prop = ifelse(is.nan(prop), NA_real_, prop),
      Actual = factor(Actual, levels = c("Poor", "Moderate", "Good")),
      Predicted = factor(Predicted, levels = c("Poor", "Moderate", "Good"))
    )
}

make_habitat_cm_df <- function(pred_df) {
  split(pred_df, pred_df$Category) |>
    purrr::imap_dfr(function(d, cat_name) {
      cm <- table(
        Predicted = factor(d$pred, levels = c("Poor", "Moderate", "Good")),
        Actual = factor(d$obs, levels = c("Poor", "Moderate", "Good"))
      )

      cm_prop <- prop.table(cm, margin = 2)

      out <- as.data.frame(cm_prop)
      names(out) <- c("Predicted", "Actual", "prop")
      out$Category <- cat_name
      out
    }) |>
    mutate(
      prop = ifelse(is.nan(prop), NA_real_, prop),
      Actual = factor(Actual, levels = c("Poor", "Moderate", "Good")),
      Predicted = factor(Predicted, levels = c("Poor", "Moderate", "Good"))
    )
}

plot_overall_cm <- function(cm_df) {
  ggplot(cm_df, aes(x = Actual, y = Predicted, fill = prop)) +
    geom_tile(color = "white") +
    geom_text(aes(label = ifelse(is.na(prop), "", sprintf("%.2f", prop))), size = 7) +
    scale_fill_gradient(low = "white", high = "steelblue", limits = c(0, 1), na.value = "white") +
    labs(
      x = "Actual Condition",
      y = "Predicted Condition",
      fill = "Proportion"
    ) +
    theme_minimal(base_size = 18) +
    theme(panel.grid = element_blank())
}

plot_habitat_cm <- function(cm_df) {
  ggplot(cm_df, aes(x = Actual, y = Predicted, fill = prop)) +
    geom_tile(color = "white") +
    geom_text(aes(label = ifelse(is.na(prop), "", sprintf("%.2f", prop))), size = 5) +
    scale_fill_gradient(low = "white", high = "steelblue", limits = c(0, 1), na.value = "white") +
    facet_wrap(~ Category) +
    labs(fill = "Proportion") +
    theme_minimal(base_size = 14) +
    theme(
      axis.title = element_blank(),
      axis.text = element_text(size = 14),
      strip.text = element_text(size = 14, face = "bold"),
      panel.grid = element_blank()
    )
}

make_roc_plot <- function(best_preds, filename = NULL) {
  roc_poor <- roc(best_preds$obs == "Poor", best_preds$Poor, quiet = TRUE)
  roc_moderate <- roc(best_preds$obs == "Moderate", best_preds$Moderate, quiet = TRUE)
  roc_good <- roc(best_preds$obs == "Good", best_preds$Good, quiet = TRUE)
  
  roc_df <- bind_rows(
    tibble(
      Specificity = roc_poor$specificities,
      Sensitivity = roc_poor$sensitivities,
      Class = "Poor"
    ),
    tibble(
      Specificity = roc_moderate$specificities,
      Sensitivity = roc_moderate$sensitivities,
      Class = "Moderate"
    ),
    tibble(
      Specificity = roc_good$specificities,
      Sensitivity = roc_good$sensitivities,
      Class = "Good"
    )
  )
  
  auc_labels <- c(
    Poor = paste0("Poor (AUC = ", round(auc(roc_poor), 2), ")"),
    Moderate = paste0("Moderate (AUC = ", round(auc(roc_moderate), 2), ")"),
    Good = paste0("Good (AUC = ", round(auc(roc_good), 2), ")")
  )
  
  random_line <- data.frame(x = 1, y = 0, xend = 0, yend = 1)
  
  p <- ggplot(roc_df, aes(x = Specificity, y = Sensitivity, color = Class)) +
    geom_line(linewidth = 1.2) +
    geom_segment(
      data = random_line,
      aes(x = x, y = y, xend = xend, yend = yend),
      inherit.aes = FALSE,
      color = "grey60",
      linewidth = 0.8
    ) +
    scale_x_reverse(limits = c(1, 0), breaks = seq(1, 0, by = -0.2)) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
    scale_color_manual(
      values = c(Poor = "#D55E00", Moderate = "#0072B2", Good = "#009E73"),
      breaks = c("Poor", "Moderate", "Good"),
      labels = auc_labels
    ) +
    labs(
      x = "Specificity",
      y = "Sensitivity",
      color = NULL
    ) +
    theme_minimal(base_size = 18) +
    theme(
      legend.position = c(0.78, 0.18),
      legend.background = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  if (!is.null(filename)) {
    ggsave(filename = filename, plot = p, width = 8, height = 6, dpi = 1000)
  }
  
  p
}

# =========
# LOAD DATA
# =========
df <- read_csv("./Data/Processed_data/combined_df.csv", show_col_types = FALSE) %>%
  mutate(Category = recode(Category,
                           "Marsh,swamp, reedbed" = "Wetland"))

rf_df <- df |>
  select(
    Condition,
    Site,
    ndvi_mean,
    evi_mean,
    ndvi_sd,
    evi_sd,
    Category
  ) |>
  mutate(Condition = clean_condition(as.character(Condition))) |>
  drop_na() |>
  mutate(
    Condition = factor(Condition, levels = c("Poor", "Moderate", "Good")),
    Site = factor(Site),
    Category = factor(Category)
  ) |>
  group_by(Site) |>
  mutate(
    ndvi_mean_sc = as.numeric(scale(ndvi_mean)),
    evi_mean_sc  = as.numeric(scale(evi_mean)),
    ndvi_sd_sc   = as.numeric(scale(ndvi_sd)),
    evi_sd_sc    = as.numeric(scale(evi_sd))
  ) |>
  ungroup()

# =========
# CROSS VALIDATION
# =========
cv_index <- createMultiFolds(rf_df$Condition, k = 10, times = 5)

control <- trainControl(
  method = "repeatedcv",
  number = 10,
  repeats = 5,
  index = cv_index,
  savePredictions = "final",
  classProbs = TRUE
)

# =========
# MODEL
# =========
rf_cv <- train(
  Condition ~ ndvi_mean_sc + evi_mean_sc + ndvi_sd_sc + evi_sd_sc,
  data = rf_df,
  method = "rf",
  trControl = control,
  importance = TRUE
)

print(rf_cv)

# =========
# OVERALL CONFUSION MATRIX
# =========
best_preds_overall <- rf_cv$pred |>
  filter(mtry == rf_cv$bestTune$mtry)

print(table(Predicted = best_preds_overall$pred, Actual = best_preds_overall$obs))
print(confusionMatrix(data = best_preds_overall$pred, reference = best_preds_overall$obs))

cm_df_overall <- make_overall_cm_df(best_preds_overall)
con_plot <- plot_overall_cm(cm_df_overall)

ggsave(
  filename = "./Plots/Confusion_matrix/Knepp_Attenborough_confusion_matrix1.pdf",
  plot = con_plot,
  width = 8,
  height = 6
)

# =========
# CONFUSION MATRIX BY HABITAT
# =========
lookup <- rf_df |>
  mutate(rowIndex = row_number()) |>
  select(rowIndex, Category)

best_preds_hab <- rf_cv$pred |>
  filter(mtry == rf_cv$bestTune$mtry) |>
  left_join(lookup, by = "rowIndex") |>
  mutate(Category = as.character(Category))

invisible(
  best_preds_hab |>
    group_by(Category) |>
    group_split() |>
    lapply(function(x) {
      cat("\nHabitat:", unique(x$Category), "\n")
      print(table(Predicted = x$pred, Actual = x$obs))
    })
)

cm_df_habitat <- make_habitat_cm_df(best_preds_hab)
con_habitat_plot <- plot_habitat_cm(cm_df_habitat)

ggsave(
  filename = "./Plots/Confusion_matrix/Knepp_Attenborough_habitat_confusion_matrix1.pdf",
  plot = con_habitat_plot,
  width = 10,
  height = 7
)

# =========
# ROC CURVES
# =========
roc_plot <- make_roc_plot(
  best_preds = best_preds_overall,
  filename = "./Plots/ROC/Knepp_Attenborough_RS_ROC1.pdf"
)

print(roc_plot)


# Accuracy per habitat
accuracy_by_habitat_cv <- best_preds_hab |>
  mutate(correct = pred == obs) |>
  group_by(rowIndex, Category) |>
  summarise(
    mean_correct = mean(correct),
    .groups = "drop"
  ) |>
  group_by(Category) |>
  summarise(
    Accuracy = mean(mean_correct),
    n_observations = n(),
    .groups = "drop"
  ) |>
  arrange(desc(Accuracy))

print(accuracy_by_habitat_cv)

acc_plot_cv <- ggplot(accuracy_by_habitat_cv, aes(x = reorder(Category, Accuracy), y = Accuracy)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    x = NULL,
    y = "Cross-validated accuracy"
  ) +
  theme_minimal(base_size = 14) +
  theme(panel.grid.minor = element_blank())

print(acc_plot_cv)

# ggsave(
#   filename = "./Plots/Confusion_matrix/Knepp_Attenborough_RS_accuracy_by_habitat_cv1.pdf",
#   plot = acc_plot_cv,
#   width = 7,
#   height = 5
# )

# ==================
# VARIABLE IMPORTANCE
# ==================

vip <- varImp(rf_cv, scale = TRUE)

vip_df <- vip$importance |>
  tibble::rownames_to_column("Variable")

vip_df$Importance <- rowMeans(
  vip_df[, setdiff(names(vip_df), "Variable")],
  na.rm = TRUE
)

label_map <- c(
  ndvi_mean = "NDVI (mean)",
  evi_mean  = "EVI (mean)",
  ndvi_sd   = "NDVI (SD)",
  evi_sd    = "EVI (SD)"
)

vip_df <- vip_df |>
  mutate(
    Variable_label = dplyr::recode(Variable, !!!label_map)
  ) |>
  arrange(desc(Importance))

print(vip_df)

vip_plot <- ggplot(vip_df, aes(x = reorder(Variable_label, Importance), y = Importance)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(
    x = NULL,
    y = "Variable importance"
  ) +
  theme_minimal(base_size = 14) +
  theme(panel.grid.minor = element_blank())

print(vip_plot)

ggsave(
  filename = "./Plots/Variable_importance/Knepp_Attenborough_variable_importance1.pdf",
  plot = vip_plot,
  width = 7,
  height = 5
)
