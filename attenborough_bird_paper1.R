###############################################################################
# DATA SETS ON MY COMPUTER (JAMES)
###############################################################################

# Attenborough Audiomoth data
attenborough_am <- read.csv("ANRbird-June2024-combinedlocationresults_paper1.csv", header = TRUE)

# Attenborough Point count data June 2024
attenborough_pc <- read.csv("Attenborough_pointcounts_birds_June_24_paper1.csv", header = TRUE)

###############################################################################
# Inspect data structure
###############################################################################
str(attenborough_pc)
str(attenborough_am)

###############################################################################
# Convert columns to factors
###############################################################################
attenborough_pc$Survey.ID        <- as.factor(attenborough_pc$Survey.ID) 
attenborough_am$Location_ID      <- as.factor(attenborough_am$Location_ID) 
attenborough_pc$Scientific_Name  <- as.factor(attenborough_pc$Scientific_Name)
attenborough_am$Scientific_Name  <- as.factor(attenborough_am$Scientific_Name)

###############################################################################
# Add presence/absence columns
###############################################################################
attenborough_pc <- attenborough_pc %>%
  mutate(Presence_Absence = ifelse(Count > 0, 1, 0))    # 1 = presence

attenborough_am <- attenborough_am %>%
  mutate(Presence_Absence = ifelse(Month_ID > 0, 1, 0)) # 1 = presence

###############################################################################
# Audiomoth: presence/absence community matrix
###############################################################################
attenborough_am_wide_pa <- attenborough_am %>%
  pivot_wider(
    names_from  = Scientific_Name,          # columns = species
    values_from = Presence_Absence,         # values = presence/absence
    values_fill = 0                         # fill missing with 0
  )

attenborough_am_species_matrix <- attenborough_am_wide_pa %>%
  select(Location_ID, Tyto_alba:last_col()) %>%
  group_by(Location_ID) %>%
  summarise(across(everything(), ~ sum(. > 0)), .groups = "drop")

# Convert to numeric matrix (drop Location_ID)
attenborough_am_species_matrix <- as.matrix(attenborough_am_species_matrix %>% select(-Location_ID))

# Run species accumulation
attenborough_am_accum <- specaccum(attenborough_am_species_matrix, method = "random")

###############################################################################
# Point counts: presence/absence community matrix
###############################################################################
attenborough_pc_wide_pa <- attenborough_pc %>%
  pivot_wider(
    names_from  = Scientific_Name,          # columns = species
    values_from = Presence_Absence,         # values = presence/absence
    values_fill = 0,                        # fill missing with 0
    values_fn   = ~max(.x, na.rm = TRUE)    # handle duplicates by taking max
  )

# Step 1: Create site-by-species matrix
attenborough_pc_species_matrix <- attenborough_pc_wide_pa %>%
  select(Survey.ID, Chroicocephalus_ridibundus:last_col()) %>%
  group_by(Survey.ID) %>%
  summarise(across(everything(), ~ sum(. > 0)), .groups = "drop")

# Step 2: Remove Survey.ID but keep as data frame
attenborough_pc_species_df <- attenborough_pc_species_matrix %>% select(-Survey.ID)

# Step 3: Ensure presence/absence
attenborough_pc_species_df <- attenborough_pc_species_df %>%
  mutate(across(everything(), ~ as.integer(. > 0))) %>%
  select(where(~ all(!is.na(.))))

# Step 4: Run species accumulation
attenborough_pc_accum <- specaccum(attenborough_pc_species_df, method = "random")

###############################################################################
# Convert specaccum objects to tidy data for ggplot
###############################################################################
attenborough_accum_to_df <- function(accum_obj, label) {
  tibble(
    Sites    = seq_along(accum_obj$sites),
    Richness = accum_obj$richness,
    SD       = accum_obj$sd
  ) %>%
    mutate(
      Lower  = pmax(Richness - SD, 0),
      Upper  = Richness + SD,
      Source = label
    )
}

attenborough_am_df    <- attenborough_accum_to_df(attenborough_am_accum, "AudioMoth")
attenborough_pc_df    <- attenborough_accum_to_df(attenborough_pc_accum, "Point Count")
attenborough_accum_df <- bind_rows(attenborough_am_df, attenborough_pc_df)

###############################################################################
# ggplot: Species accumulation curves with CI
###############################################################################
library(ggplot2)

attenborough_accum_plot <- ggplot(attenborough_accum_df, 
                                  aes(x = Sites, y = Richness,
                                      colour = Source, fill = Source)) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.25, colour = NA) +
  geom_line(linewidth = 1.2) +
  scale_colour_manual(values = c("AudioMoth" = "deepskyblue3",
                                 "Point Count" = "firebrick3"),
                      labels = c("AudioMoth" = "Passive acoustic",
                                 "Point Count" = "Point count")) +
  scale_fill_manual(values   = c("AudioMoth" = "deepskyblue3",
                                 "Point Count" = "firebrick3"),
                    labels = c("AudioMoth" = "Passive acoustic",
                               "Point Count" = "Point count")) +
  scale_x_continuous(limits = c(0, 45),
                     breaks = c(10, 20, 30, 40),
                     expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 80),
                     breaks = seq(0, 75, by = 25),
                     expand = c(0, 0)) +
  labs(
    title  = NULL,
    x      = "No. of samples",
    y      = "No. of species",
    colour = NULL, fill = NULL
  ) +
  theme_bw(base_size = 20) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position  = c(0.95, 0.05),   # bottom-right inside plot
    legend.justification = c("right", "bottom"),
    legend.background = element_rect(fill = "white", colour = "black"),
    axis.title   = element_text(size = 22, face = "bold", colour = "black"),
    axis.text    = element_text(size = 20, colour = "black"),
    legend.text  = element_text(size = 12, colour = "black")
  )

# Print
attenborough_accum_plot





################################################################################
################################################################################
################################################################################

# Recreate an AM site-by-species PA data frame with the Location_ID retained
attenborough_am_species_df <- attenborough_am_wide_pa %>%
  select(Location_ID, Tyto_alba:last_col()) %>%
  group_by(Location_ID) %>%
  summarise(across(everything(), ~ sum(. > 0)), .groups = "drop")

# Convert to strict 0/1 (presence/absence)
attenborough_am_pa_df <- attenborough_am_species_df %>%
  select(-Location_ID) %>%
  mutate(across(everything(), ~ as.integer(. > 0)))

# Point-count PA data frame already exists from your pipeline
attenborough_pc_pa_df <- attenborough_pc_species_df  # already 0/1

# ---- Align species columns between the two datasets ----
add_missing_cols <- function(df, all_species) {
  missing <- setdiff(all_species, names(df))
  if (length(missing) > 0) df[missing] <- 0L
  df %>% select(all_of(all_species))
}

attenborough_common_species <- union(names(attenborough_pc_pa_df),
                                     names(attenborough_am_pa_df))

attenborough_pc_pa_df2 <- add_missing_cols(attenborough_pc_pa_df, attenborough_common_species)
attenborough_am_pa_df2 <- add_missing_cols(attenborough_am_pa_df, attenborough_common_species)

# ---- Combine with Source + IDs ----
attenborough_pc_comb <- attenborough_pc_pa_df2 %>%
  mutate(Source = "Point count",
         ID = as.character(attenborough_pc_species_matrix$Survey.ID))

attenborough_am_comb <- attenborough_am_pa_df2 %>%
  mutate(Source = "Passive acoustic",
         ID = as.character(attenborough_am_species_df$Location_ID))

attenborough_combined <- bind_rows(attenborough_pc_comb, attenborough_am_comb)

# Remove empty sites (rows with all zeros)
attenborough_combined_filt <- attenborough_combined %>%
  filter(rowSums(select(., -Source, -ID)) > 0)

# ---- Run NMDS (Jaccard on binary PA) ----
set.seed(123)
attenborough_nmds <- metaMDS(
  attenborough_combined_filt %>% select(-Source, -ID),
  distance = "jaccard",
  binary   = TRUE,
  k        = 2,
  trymax   = 400,
  autotransform = FALSE
)

# ---- Extract site scores ----
attenborough_scores <- as.data.frame(scores(attenborough_nmds, display = "sites"))
attenborough_scores$Source <- attenborough_combined_filt$Source
attenborough_scores$ID     <- attenborough_combined_filt$ID

# ---- Plot (legend inside top-right) ----
attenborough_nmds_plot <- ggplot(attenborough_scores, aes(x = NMDS1, y = NMDS2, colour = Source)) +
  geom_point(size = 4) +
  geom_polygon(
    data = attenborough_scores %>%
      group_by(Source) %>%
      slice(chull(NMDS1, NMDS2)),
    aes(fill = Source, group = Source),
    alpha = 0.15, colour = "black", linewidth = 0.8, show.legend = FALSE
  ) +
  scale_colour_manual(values = c("Passive acoustic" = "deepskyblue3",
                                 "Point count"      = "firebrick3")) +
  scale_fill_manual(values   = c("Passive acoustic" = "deepskyblue3",
                                 "Point count"      = "firebrick3")) +
  labs(x = "NMDS1", y = "NMDS2", colour = NULL, fill = NULL) +
  theme_bw(base_size = 20) +
  theme(
    panel.grid         = element_blank(),
    legend.position    = c(0.95, 0.95),      # inside top-right
    legend.justification = c("right", "top"),
    legend.background  = element_rect(fill = "white", colour = "black"),
    axis.title         = element_text(size = 22, face = "bold", colour = "black"),
    axis.text          = element_text(size = 20, colour = "black"),
    legend.text  = element_text(size = 12, colour = "black")
  )

# Print
attenborough_nmds_plot



######

# Get unique species
am_species <- unique(attenborough_am$Scientific_Name[attenborough_pc$Scientific_Name != ""])
pc_species <- unique(attenborough_pc$Scientific_Name[attenborough_pc$Scientific_Name != ""])

only_am <- setdiff(am_species, pc_species)
only_pc <- setdiff(pc_species, am_species)
both <- intersect(am_species, pc_species)

length(only_am); length(only_pc); length(both)

# Create a data frame with species and detection method
species_summary <- data.frame(
  Species = c(only_am, only_pc, both),
  Detected_By = c(
    rep("AM only", length(only_am)),
    rep("PC only", length(only_pc)),
    rep("Both", length(both))
  )
)

# View the first few rows
head(species_summary)

length(only_am)   # number of species detected by AM only
length(only_pc)   # number of species detected by PC only
length(both)      # number of species detected by both methods