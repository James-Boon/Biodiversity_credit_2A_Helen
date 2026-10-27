libs <- c("dplyr", "sf", "readr", "ggplot2", "randomForest", "caret")

installed_libs <- libs %in% rownames(
  installed.packages())

if (any(installed_libs == F)) {install.packages(
    libs[!installed_libs]
  )}

invisible(lapply(libs, library, character.only = T))
rm(list = ls())

folder <- "./OneDrive_1_11-03-2026/Knepp_edited_shapefiles/"

# List all .shp files
shp_files <- list.files(folder, pattern = "\\.shp$", full.names = TRUE)
RS_df <- read_csv("./OneDrive_1_11-03-2026/Pixel_level_joined.csv") %>%
  rename("ID" = "polygon_nu")

# Read and combine
habitats <- shp_files %>%
  lapply(st_read) %>%
  bind_rows()
# 
# head(habitats)
# names(habitats)

habitats_df <- st_drop_geometry(habitats)

# IDs in RS data but not in habitats
setdiff(unique(RS_df$ID), unique(habitats_df$ID))

# IDs in habitats but not in RS data
setdiff(unique(habitats_df$ID), unique(RS_df$ID))

merged_df <- merge(RS_df, habitats_df, by = "ID")

polygon_df <- merged_df %>%
  st_drop_geometry() %>%   # remove spatial column
  group_by(ID, Condition, Habitat) %>%
  summarise(
    ndvi_mean = mean(ndvi_value, na.rm = TRUE),
    ndvi_sd   = sd(ndvi_value, na.rm = TRUE),
    evi_mean  = mean(evi_value, na.rm = TRUE),
    evi_sd    = sd(evi_value, na.rm = TRUE),
    ndvi_comp_mean = mean(ndvi_comp, na.rm = TRUE),
    evi_comp_mean  = mean(evi_comp, na.rm = TRUE),
    .groups = "drop"
  )

# ================
# DATA EXPLORATION
# ================

ggplot(polygon_df, aes(x = Condition, y = ndvi_sd, fill = Condition)) +
  geom_boxplot() +
  theme_minimal()


vars <- polygon_df %>%
  select(ndvi_mean, evi_mean, ndvi_sd, evi_sd)

model <- lm(ndvi_mean ~ Condition, data = polygon_df)
hist(residuals(model), main = "ndvi_mean")
model <- lm(evi_mean ~ Condition, data = polygon_df)
hist(residuals(model), main = "evi_mean")

# Non parametric test
kruskal.test(ndvi_mean ~ Condition, data = polygon_df)

# Diagnose: habitat type major confounder
polygon_df %>%
  arrange(ndvi_mean) %>%
  head()

# Models must include habitat as a confounding variable
model <- lm(ndvi_mean ~ Condition + Habitat, data = polygon_df)
anova(model)

# Explore without ponds (likely reason for skew)
no_ponds <- polygon_df %>%
  filter(!grepl("Pond", Habitat))

ggplot(no_ponds, aes(x = Condition, y = ndvi_mean)) +
  geom_boxplot()

model <- lm(ndvi_mean ~ Condition, data = no_ponds)
hist(residuals(model), main = "ndvi_mean")
model <- lm(evi_mean ~ Condition, data = no_ponds)
hist(residuals(model), main = "evi_mean")

# Habitat specific loop

by(polygon_df, polygon_df$Habitat, function(df) {
  if(length(unique(df$Condition)) > 1) {
    summary(aov(ndvi_mean ~ Condition, data = df))
  }
})

# Issue: sample size per habitat is to small in some cases, run global model
table(polygon_df$Habitat)

# ==========
# Run models
# ==========

#### NDVI ####

# ANOVA
anova(lm(ndvi_mean ~ Condition + Habitat, data = polygon_df))
# Linear regression
summary(lm(ndvi_mean ~ Condition + Habitat, data = polygon_df))


# exclude small sample sizes (SENSITIVITY)
large_df <- polygon_df %>%
  filter(Habitat %in% c(
    "Other neutral grassland",
    "Lowland mixed deciduous woodland"
  ))

lm(ndvi_mean ~ Condition + Habitat, data = large_df)


#### EVI ####
# ANOVA
anova(lm(evi_mean ~ Condition + Habitat, data = polygon_df))
# Linear regression
summary(lm(evi_mean ~ Condition + Habitat, data = polygon_df))

# Sensitivity
lm(evi_mean ~ Condition + Habitat, data = large_df)

# ===========
# MULTI ANOVA
# ===========

vars <- polygon_df %>%
  select(ndvi_mean, evi_mean, ndvi_sd, evi_sd)

manova_model <- manova(
  cbind(ndvi_mean, evi_mean, ndvi_sd, evi_sd) ~ Condition + Habitat,
  data = polygon_df
)

summary(manova_model, test = "Pillai")
# Condition has 30% varaince in multivariate space

# =============
# RANDOM FOREST
# =============

# Keep relevant variables
rf_df <- polygon_df %>%
  select(Condition, ndvi_mean, evi_mean, ndvi_sd, evi_sd) %>%
  na.omit()

rf_df$Condition <- factor(rf_df$Condition,
                          levels = c("1. Good", "3. Moderate", "5. Poor"),
                          labels = c("Good", "Moderate", "Poor"))


# Make sure Condition is a factor
rf_df$Condition <- as.factor(rf_df$Condition)
str(rf_df)

set.seed(123)

train_index <- createDataPartition(rf_df$Condition, p = 0.8, list = FALSE)

train_data <- rf_df[train_index, ]
test_data  <- rf_df[-train_index, ]

rf_model <- randomForest(
  Condition ~ ., 
  data = train_data,
  ntree = 500,
  importance = TRUE
)

predictions <- predict(rf_model, test_data)
confusionMatrix(predictions, test_data$Condition)

# NIR 
prop.table(table(rf_df$Condition))
max(prop.table(table(rf_df$Condition)))

# Cross validation
control <- trainControl(
  method = "repeatedcv",
  number = 10,
  repeats = 20,
  classProbs = TRUE,
  savePredictions = "final"
)

set.seed(123)

rf_cv <- train(
  Condition ~ ., 
  data = rf_df,
  method = "rf",
  trControl = control,
  ntree = 500
)
head(rf_cv$pred)
rf_cv
rf_cv$results

cm_cv <- confusionMatrix(rf_cv$pred$pred, rf_cv$pred$obs)
cm_cv

results_table <- as.data.frame(cm_cv$byClass)

results_table$Class <- rownames(results_table)
rownames(results_table) <- NULL

results_table[, c("Class", "Sensitivity", "Specificity", "Pos Pred Value", "Balanced Accuracy")]

# No information rate
prop.table(table(rf_df$Condition))
max(prop.table(table(rf_df$Condition)))

summary_table <- data.frame(
  Metric = c("Accuracy", "Kappa", "No Information Rate"),
  Value = c(
    rf_cv$results$Accuracy[rf_cv$results$mtry == rf_cv$bestTune$mtry],
    rf_cv$results$Kappa[rf_cv$results$mtry == rf_cv$bestTune$mtry],
    nir
  )
)

summary_table
# =========
# PCA PLOT
# =========

vars <- polygon_df %>%
  select(ndvi_mean, evi_mean, ndvi_sd, evi_sd)

pca <- prcomp(vars, scale. = TRUE)

pca_df <- as.data.frame(pca$x)

pca_df$Condition <- polygon_df$Condition
pca_df$Habitat <- polygon_df$Habitat

p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Condition)) +
  geom_point(size = 2.5, alpha = 0.7) +
  stat_ellipse(level = 0.68, linewidth = 1) +
  theme_minimal(base_size = 14) +
  labs(
    title = "PCA of Remote Sensing Metrics",
    subtitle = "Coloured by Habitat Condition",
    x = paste0("PC1 (", round(summary(pca)$importance[2,1]*100, 1), "%)"),
    y = paste0("PC2 (", round(summary(pca)$importance[2,2]*100, 1), "%)")
  )

dir.create("./plots/PCA",showWarnings = F)
ggsave(
  filename = "./plots/PCA/PCA_plot.png",
  plot = p,
  width = 8,
  height = 6,
  dpi = 300
)
