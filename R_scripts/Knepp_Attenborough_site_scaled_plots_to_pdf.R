# =========
# CLEAN ENV
# =========
rm(list = ls())

# =========
# PACKAGES
# =========
libs <- c("tidyverse", "ordinal", "ggh4x")

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
# PREP DATA + WITHIN-SITE SCALING
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
# OUTPUT DIRS
# =========
dir.create("./Plots", showWarnings = FALSE)
dir.create("./Plots/BoxPlots/", showWarnings = FALSE)
dir.create("./Plots/Ordinal_Regression/", showWarnings = FALSE)
dir.create("./Plots/PCA", showWarnings = FALSE)
dir.create("./Plots/Outputs", showWarnings = FALSE)

# =======================================
# 1. GROUP DIFFERENCE ANALYSIS: MANOVA
# =======================================

manova_model <- manova(
  cbind(ndvi_mean_sc, evi_mean_sc, ndvi_sd_sc, evi_sd_sc) ~ Condition + Category,
  data = analysis_df
)

manova_pillai <- summary(manova_model, test = "Pillai")
print(manova_pillai)

manova_aov <- summary.aov(manova_model)
print(manova_aov)


# Outputs
# Overall MANOVA table
pillai_table <- as.data.frame(manova_pillai$stats, check.names = FALSE)

pillai_table$Test <- "Overall MANOVA"
pillai_table$Response <- "Multivariate"
pillai_table$Term <- rownames(pillai_table)
rownames(pillai_table) <- NULL

names(pillai_table) <- gsub("approx F", "F_value", names(pillai_table))
names(pillai_table) <- gsub("Pr\\(>F\\)", "p_value", names(pillai_table))

pillai_table <- pillai_table[, c(
  "Test", "Response", "Term", "Df", "Pillai",
  "F_value", "num Df", "den Df", "p_value"
)]

# variable-specific
manova_aov_list <- summary.aov(manova_model)

manova_table <- do.call(
  rbind,
  lapply(seq_along(manova_aov_list), function(i) {
    
    df <- as.data.frame(
      manova_aov_list[[i]],
      check.names = FALSE
    )
    
    df$Response <- names(manova_aov_list)[i]
    df$Term <- rownames(df)
    rownames(df) <- NULL
    
    names(df)[1:5] <- c("Df", "Sum_Sq", "Mean_Sq", "F_value", "p_value")
    
    df[, c("Response", "Term", "Df", "Sum_Sq", "Mean_Sq", "F_value", "p_value")]
  })
)

manova_table


# Add matching columns to univariate ANOVA table
manova_table$Test <- "Univariate ANOVA"
manova_table$Pillai <- NA
manova_table$`num Df` <- NA
manova_table$`den Df` <- NA

# Rename if needed
names(manova_table)[names(manova_table) == "F_value"] <- "F_value"
names(manova_table)[names(manova_table) == "p_value"] <- "p_value"

combined_table <- rbind(
  pillai_table,
  manova_table[, c(
    "Test", "Response", "Term", "Df", "Pillai",
    "F_value", "num Df", "den Df", "p_value"
  )]
)

combined_table


# write.csv(
#   combined_table,
#   "./Plots/Outputs/Knepp_Attenborough_MANOVA_table.csv",
#   row.names = FALSE
# )
# 
# capture.output(
#   manova_pillai,
#   file = "./Plots/Outputs/Knepp_Attenborough_MANOVA_Pillai_category_control_site_scaled.txt"
# )
# 
# capture.output(
#   manova_aov,
#   file = "./Plots/Outputs/Knepp_Attenborough_MANOVA_univariate_ANOVAs_category_control_site_scaled.txt"
# )

# =========
# SINGLE REPORT BOXPLOT
# =========

plot_df <- analysis_df |>
  pivot_longer(
    cols = c(
      ndvi_mean_sc,
      evi_mean_sc,
      ndvi_sd_sc,
      evi_sd_sc
    ),
    names_to = "Metric",
    values_to = "Value"
  ) |>
  mutate(
    Metric = recode(
      Metric,
      ndvi_mean_sc = "NDVI (mean)",
      evi_mean_sc = "EVI (mean)",
      ndvi_sd_sc = "NDVI (SD)",
      evi_sd_sc = "EVI (SD)"
    )
  )

# Remove extreme outlier for visualisation
plot_df <- plot_df |>
  filter(
    !(Category == "Lakes" &
        Condition == "Poor" &
        abs(Value) > 6)
  )

combined_boxplot <- ggplot(
  plot_df,
  aes(x = Condition, y = Value, fill = Condition)
) +
  geom_boxplot(alpha = 0.8, outlier.size = 0.5) +
  facet_grid(
    Metric ~ Category,
    scales = "free_y"
  ) +
  theme_minimal(base_size = 16) +
  labs(
    x = "Habitat condition",
    y = "Within-site scaled value"
  ) +
  theme(
    legend.position = "none",
    strip.text.y = element_text(angle = 0, face = "bold"),
    strip.background = element_rect(fill = "grey90", colour = NA),
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.spacing.y = unit(1, "lines"),
    panel.grid.major.y = element_line(size = 0.3),
    panel.grid.minor.y = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6)
  )

ggsave(
  "./Plots/BoxPlots/Knepp_att_combined_vertical_boxplot_site_scaled.pdf",
  combined_boxplot,
  width = 10,
  height = 14

)

# =======================================
# 2. ORDINAL REGRESSION
# =======================================

ord_model <- clm(
  Condition ~ ndvi_mean_sc + evi_mean_sc + ndvi_sd_sc + evi_sd_sc + Category,
  data = analysis_df,
  link = "logit"
)

ord_summary <- summary(ord_model)
print(ord_summary)

ord_coef <- as.data.frame(ord_summary$coefficients)
ord_coef$Term <- rownames(ord_coef)
rownames(ord_coef) <- NULL

ord_coef <- ord_coef |>
  rename(
    Estimate = Estimate,
    Std_Error = `Std. Error`,
    z_value = `z value`,
    p_value = `Pr(>|z|)`
  ) |>
  mutate(
    CI_low = Estimate - 1.96 * Std_Error,
    CI_high = Estimate + 1.96 * Std_Error
  )

print(ord_coef)

# write_csv(
#   ord_coef,
#   "./Plots/Outputs/Knepp_Attenborough_ordinal_regression_coefficients_category_control_site_scaled.csv"
# )

ord_plot_df <- ord_coef |>
  filter(!grepl("\\|", Term)) |>
  mutate(
    Term = recode(
      Term,
      ndvi_mean_sc = "NDVI (mean, site-scaled)",
      evi_mean_sc = "EVI (mean, site-scaled)",
      ndvi_sd_sc = "NDVI (SD, site-scaled)",
      evi_sd_sc = "EVI (SD, site-scaled)",
      `CategoryLakes` = "Habitat: Lakes",
      `CategoryScrub` = "Habitat: Scrub",
      `CategoryWetland` = "Habitat: Wetland",
      `CategoryWoodland` = "Habitat: Woodland"
    )
  )

coef_plot <- ggplot(ord_plot_df, aes(x = reorder(Term, Estimate), y = Estimate)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high), width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  coord_flip() +
  theme_minimal(base_size = 14) +
  labs(
    x = NULL,
    y = "Coefficient estimate (log-odds scale)",
    title = "Ordinal regression coefficients using within-site scaled predictors"
  )

print(coef_plot)

ggsave(
  "./Plots/Ordinal_Regression/Knepp_Attenborough_ordinal_regression_coefficients_category_control_site_scaled.pdf",
  coef_plot,
  width = 8,
  height = 6

)

# =======================================
# 3. PCA VISUALISATION
# =======================================

pca_vars <- analysis_df |>
  select(ndvi_mean_sc, evi_mean_sc, ndvi_sd_sc, evi_sd_sc)

pca_model <- prcomp(pca_vars, center = TRUE, scale. = FALSE)

pca_scores <- as.data.frame(pca_model$x) |>
  bind_cols(
    analysis_df |>
      select(Condition, Site, Category)
  )

pca_var <- summary(pca_model)$importance[2, 1:2] * 100

pca_loadings <- as.data.frame(pca_model$rotation) |>
  rownames_to_column("Variable") |>
  mutate(
    Variable = recode(
      Variable,
      ndvi_mean_sc = "NDVI mean",
      evi_mean_sc = "EVI mean",
      ndvi_sd_sc = "NDVI SD",
      evi_sd_sc = "EVI SD"
    )
  )

# write_csv(
#   pca_loadings,
#   "./Plots/Outputs/Knepp_Attenborough_PCA_loadings_category_control_site_scaled.csv"
# )

loading_multiplier <- 3

pca_biplot <- ggplot(
  pca_scores,
  aes(x = PC1, y = PC2, color = Condition, shape = Category)
) +
  geom_point(size = 2.5, alpha = 0.75) +
  stat_ellipse(aes(group = Condition), level = 0.68, linewidth = 0.8) +
  geom_segment(
    data = pca_loadings,
    aes(
      x = 0, y = 0,
      xend = PC1 * loading_multiplier,
      yend = PC2 * loading_multiplier
    ),
    inherit.aes = FALSE,
    arrow = arrow(length = grid::unit(0.2, "cm")),
    linewidth = 0.7,
    color = "grey30"
  ) +
  geom_text(
    data = pca_loadings,
    aes(
      x = PC1 * loading_multiplier,
      y = PC2 * loading_multiplier,
      label = Variable
    ),
    inherit.aes = FALSE,
    size = 4
  ) +
  theme_minimal(base_size = 14) +
  labs(
    title = "PCA biplot of within-site scaled spectral variables",
    x = paste0("PC1 (", round(pca_var[1], 1), "%)"),
    y = paste0("PC2 (", round(pca_var[2], 1), "%)")
  )

ggsave(
  "./Plots/PCA/Knepp_Attenborough_PCA_biplot_category_shape_site_scaled.pdf",
  pca_biplot,
  width = 9,
  height = 7

)



# # =========
# # CLEAN ENV
# # =========
# rm(list = ls())
# 
# # =========
# # PACKAGES
# # =========
# libs <- c("tidyverse", "ordinal")
# 
# installed_libs <- libs %in% rownames(installed.packages())
# if (any(!installed_libs)) {
#   install.packages(libs[!installed_libs])
# }
# 
# invisible(lapply(libs, library, character.only = TRUE))
# 
# # =========
# # LOAD DATA
# # =========
# df <- read_csv("./Data/Processed_data/combined_df.csv", show_col_types = FALSE) |>
#   mutate(
#     Category = recode(Category,
#                       "Marsh,swamp, reedbed" = "Wetland"
#     ),
#     Condition = case_when(
#       Condition %in% c("Good", "Fairly good") ~ "Good",
#       Condition == "Moderate" ~ "Moderate",
#       Condition %in% c("Poor", "Fairly poor") ~ "Poor",
#       TRUE ~ NA_character_
#     )
#   )
# 
# # =========
# # PREP DATA
# # =========
# analysis_df <- df |>
#   select(
#     Condition,
#     Site,
#     Category,
#     ndvi_mean,
#     evi_mean,
#     ndvi_sd,
#     evi_sd
#   ) |>
#   drop_na() |>
#   mutate(
#     Condition = factor(Condition, levels = c("Poor", "Moderate", "Good"), ordered = TRUE),
#     Site = factor(Site),
#     Category = factor(Category)
#   )
# 
# # =========
# # OUTPUT DIRS
# # =========
# dir.create("./Plots", showWarnings = FALSE)
# dir.create("./Plots/BoxPlots/", showWarnings = FALSE)
# dir.create("./Plots/Ordinal_Regression/", showWarnings = FALSE)
# dir.create("./Plots/PCA", showWarnings = FALSE)
# dir.create("./Plots/Outputs", showWarnings = FALSE)
# 
# # =======================================
# # 1. GROUP DIFFERENCE ANALYSIS: MANOVA
# # =======================================
# 
# manova_model <- manova(
#   cbind(ndvi_mean, evi_mean, ndvi_sd, evi_sd) ~ Condition + Category,
#   data = analysis_df
# )
# 
# # Multivariate test
# manova_pillai <- summary(manova_model, test = "Pillai")
# print(manova_pillai)
# 
# # Follow-up univariate ANOVAs
# manova_aov <- summary.aov(manova_model)
# print(manova_aov)
# 
# # Save MANOVA results
# capture.output(
#   manova_pillai,
#   file = "./Plots/Outputs/Knepp_Attenborough_MANOVA_Pillai_category_control.txt"
# )
# 
# capture.output(
#   manova_aov,
#   file = "./Plots/Outputs/Knepp_Attenborough_MANOVA_univariate_ANOVAs_category_control.txt"
# )
# 
# 
# # =========
# # SINGLE REPORT BOXPLOT
# # =========
# 
# plot_df <- analysis_df |>
#   pivot_longer(
#     cols = c(
#       ndvi_mean,
#       evi_mean,
#       ndvi_sd,
#       evi_sd
#     ),
#     names_to = "Metric",
#     values_to = "Value"
#   ) |>
#   mutate(
#     Metric = recode(
#       Metric,
#       ndvi_mean = "NDVI (mean)",
#       evi_mean = "EVI (mean)",
#       ndvi_sd = "NDVI (SD)",
#       evi_sd = "EVI (SD)",
#       lai_mean = "Leaf Area Index (mean)",
#       Canopy_height_max_mean = "Maximum canopy height (mean)",
#       Canopy_height_var_mean = "Canopy height variability (mean)"
#     )
#   )
# 
# combined_boxplot <- ggplot(
#   plot_df,
#   aes(x = Condition, y = Value, fill = Condition)
# ) +
#   geom_boxplot(alpha = 0.8, outlier.size = 0.5) +
#   facet_grid(
#     Metric ~ Category,
#     scales = "free_y"
#   ) +
#   facetted_pos_scales(
#     y = list(
#       Metric == "EVI (SD)" ~ scale_y_continuous(
#         breaks = c(0, 0.05, 0.10, 0.15),
#         limits = c(0, 0.16)
#       )
#     )
#   ) +
#   theme_minimal(base_size = 10) +
#   labs(
#     x = "Habitat condition",
#     y = NULL,
#   ) +
#   theme(
#     legend.position = "none",
#     strip.text.y = element_text(angle = 0, face = "bold"),
#     strip.background = element_rect(fill = "grey90", colour = NA),
#     plot.title = element_text(face = "bold"),
#     axis.text.x = element_text(angle = 45, hjust = 1),
#     
#     # 👇 this is the key part
#     panel.spacing.y = unit(1, "lines"),
#     panel.grid.major.y = element_line(size = 0.3),
#     panel.grid.minor.y = element_blank()
#   ) +
#   theme(
#     panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6),
#     panel.spacing.y = unit(0.8, "lines")
#   )
# 
# 
# ggsave(
#   "./Plots/BoxPlots/Knepp_att_combined_vertical_boxplot.pdf",
#   combined_boxplot,
#   width = 10,
#   height = 14

# )
# 
# # =======================================
# # 2. ORDINAL REGRESSION
# # =======================================
# 
# ord_model <- clm(
#   Condition ~ ndvi_mean + evi_mean + ndvi_sd + evi_sd + Category,
#   data = analysis_df,
#   link = "logit"
# )
# 
# ord_summary <- summary(ord_model)
# print(ord_summary)
# 
# # Tidy coefficients
# ord_coef <- as.data.frame(ord_summary$coefficients)
# ord_coef$Term <- rownames(ord_coef)
# rownames(ord_coef) <- NULL
# 
# ord_coef <- ord_coef |>
#   rename(
#     Estimate = Estimate,
#     Std_Error = `Std. Error`,
#     z_value = `z value`,
#     p_value = `Pr(>|z|)`
#   ) |>
#   mutate(
#     CI_low = Estimate - 1.96 * Std_Error,
#     CI_high = Estimate + 1.96 * Std_Error
#   )
# 
# print(ord_coef)
# 
# write_csv(
#   ord_coef,
#   "./Plots/Outputs/Knepp_Attenborough_ordinal_regression_coefficients_category_control.csv"
# )
# 
# # Coefficient plot on log-odds scale (excluding threshold terms)
# ord_plot_df <- ord_coef |>
#   filter(!grepl("\\|", Term)) |>
#   mutate(
#     Term = recode(Term,
#                   ndvi_mean = "NDVI (mean)",
#                   evi_mean = "EVI (mean)",
#                   ndvi_sd = "NDVI (SD)",
#                   evi_sd = "EVI (SD)",
#                   `CategoryLakes` = "Habitat: Lakes",
#                   `CategoryScrub` = "Habitat: Scrub",
#                   `CategoryWetland` = "Habitat: Wetland",
#                   `CategoryWoodland` = "Habitat: Woodland"
#     )
#   )
# 
# coef_plot <- ggplot(ord_plot_df, aes(x = reorder(Term, Estimate), y = Estimate)) +
#   geom_point(size = 3) +
#   geom_errorbar(aes(ymin = CI_low, ymax = CI_high), width = 0.2) +
#   geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
#   coord_flip() +
#   theme_minimal(base_size = 14) +
#   labs(
#     x = NULL,
#     y = "Coefficient estimate (log-odds scale)",
#     title = "Ordinal regression coefficients"
#   )
# 
# print(coef_plot)
# 
# ggsave(
#   "./Plots/Ordinal_Regression/Knepp_Attenborough_ordinal_regression_coefficients_category_control.pdf",
#   coef_plot,
#   width = 8,
#   height = 6

# )
# 
# # =======================================
# # 3. PCA VISUALISATION
# # =======================================
# 
# pca_vars <- analysis_df |>
#   select(ndvi_mean, evi_mean, ndvi_sd, evi_sd)
# 
# pca_model <- prcomp(pca_vars, center = TRUE, scale. = TRUE)
# 
# # Scores
# pca_scores <- as.data.frame(pca_model$x) |>
#   bind_cols(
#     analysis_df |>
#       select(Condition, Site, Category)
#   )
# 
# # Variance explained
# pca_var <- summary(pca_model)$importance[2, 1:2] * 100
# 
# # Loadings
# pca_loadings <- as.data.frame(pca_model$rotation) |>
#   rownames_to_column("Variable")
# 
# write_csv(
#   pca_loadings,
#   "./Plots/Outputs/Knepp_Attenborough_PCA_loadings_category_control.csv"
# )
# 
# # PCA biplot
# loading_multiplier <- 3
# 
# pca_biplot <- ggplot(
#   pca_scores,
#   aes(x = PC1, y = PC2, color = Condition, shape = Category)
# ) +
#   geom_point(size = 2.5, alpha = 0.75) +
#   stat_ellipse(aes(group = Condition), level = 0.68, linewidth = 0.8) +
#   geom_segment(
#     data = pca_loadings,
#     aes(
#       x = 0, y = 0,
#       xend = PC1 * loading_multiplier,
#       yend = PC2 * loading_multiplier
#     ),
#     inherit.aes = FALSE,
#     arrow = arrow(length = grid::unit(0.2, "cm")),
#     linewidth = 0.7,
#     color = "grey30"
#   ) +
#   geom_text(
#     data = pca_loadings,
#     aes(
#       x = PC1 * loading_multiplier,
#       y = PC2 * loading_multiplier,
#       label = Variable
#     ),
#     inherit.aes = FALSE,
#     size = 4
#   ) +
#   theme_minimal(base_size = 14) +
#   labs(
#     title = "PCA biplot of spectral variables",
#     x = paste0("PC1 (", round(pca_var[1], 1), "%)"),
#     y = paste0("PC2 (", round(pca_var[2], 1), "%)")
#   )
# 
# ggsave(
#   "./Plots/PCA/Knepp_Attenborough_PCA_biplot_category_shape.pdf",
#   pca_biplot,
#   width = 9,
#   height = 7

# )