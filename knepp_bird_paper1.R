## Knepp bird data first
###############################################################################
# DATA SETS ON MY COMPUTER (JAMES)
###############################################################################

# Knepp Audiomoth data
knepp_am <- read.csv("Kneppbird-June2024-combinedlocationresults.csv", header = TRUE)

# Knepp Point count data June 2024
knepp_pc <- read.csv("Knepp_pointcounts_0624.csv", header = TRUE)

###############################################################################
# HELEN'S CODE FOR HER COMPUTER 
###############################################################################
# knepp_am <- read.csv("C:/Users/helford/Documents/Postdoc paper/WP3_data/Birds/Knepp/Carbon_rewild/Kneppbird-June2024-combinedlocationresults.csv")
# knepp_pc <- read.csv("C:/Users/helford/Documents/Postdoc paper/WP3_data/Birds/Knepp/Knepp_pointcounts_0624.csv")

###############################################################################
# Inspect data structure
###############################################################################
str(knepp_pc)
str(knepp_am)

###############################################################################
# Convert columns to factors
###############################################################################
knepp_pc$Survey.ID        <- as.factor(knepp_pc$Survey.ID) 
knepp_am$Location_ID      <- as.factor(knepp_am$Location_ID) 
knepp_pc$Scientific_Name  <- as.factor(knepp_pc$Scientific_Name)
knepp_am$Scientific_Name  <- as.factor(knepp_am$Scientific_Name)

###############################################################################
# Add presence/absence columns
###############################################################################
knepp_pc <- knepp_pc %>%
  mutate(Presence_Absence = ifelse(Count > 0, 1, 0))    # 1 = presence

knepp_am <- knepp_am %>%
  mutate(Presence_Absence = ifelse(Month_ID > 0, 1, 0)) # 1 = presence

###############################################################################
# Audiomoth: presence/absence community matrix
###############################################################################
knepp_am_wide_pa <- knepp_am %>%
  pivot_wider(
    names_from  = Scientific_Name,          # columns = species
    values_from = Presence_Absence,         # values = presence/absence
    values_fill = 0                         # fill missing with 0
  )

knepp_am_species_matrix <- knepp_am_wide_pa %>%
  select(Location_ID, Tyto_alba:last_col()) %>%
  group_by(Location_ID) %>%
  summarise(across(everything(), ~ sum(. > 0)), .groups = "drop")

# Convert to numeric matrix (drop Location_ID)
knepp_am_species_matrix <- as.matrix(knepp_am_species_matrix %>% select(-Location_ID))

# Run species accumulation
knepp_am_accum <- specaccum(knepp_am_species_matrix, method = "random")

###############################################################################
# Point counts: presence/absence community matrix
###############################################################################
knepp_pc_wide_pa <- knepp_pc %>%
  pivot_wider(
    names_from  = Scientific_Name,          # columns = species
    values_from = Presence_Absence,         # values = presence/absence
    values_fill = 0,                        # fill missing with 0
    values_fn   = ~max(.x, na.rm = TRUE)    # handle duplicates by taking max
  )

# Step 1: Create site-by-species matrix
knepp_pc_species_matrix <- knepp_pc_wide_pa %>%
  select(Survey.ID, Turdus_merula:last_col()) %>%
  group_by(Survey.ID) %>%
  summarise(across(everything(), ~ sum(. > 0)), .groups = "drop")

# Step 2: Remove Survey.ID but keep as data frame
knepp_pc_species_df <- knepp_pc_species_matrix %>% select(-Survey.ID)

# Step 3: Ensure presence/absence
knepp_pc_species_df <- knepp_pc_species_df %>%
  mutate(across(everything(), ~ as.integer(. > 0))) %>%
  select(where(~ all(!is.na(.))))

# Step 4: Run species accumulation
knepp_pc_accum <- specaccum(knepp_pc_species_df, method = "random")

###############################################################################
# Convert specaccum objects to tidy data for ggplot
###############################################################################
knepp_accum_to_df <- function(accum_obj, label) {
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

knepp_am_df    <- knepp_accum_to_df(knepp_am_accum, "AudioMoth")
knepp_pc_df    <- knepp_accum_to_df(knepp_pc_accum, "Point Count")
knepp_accum_df <- bind_rows(knepp_am_df, knepp_pc_df)

###############################################################################
# ggplot: Species accumulation curves with CI
###############################################################################
library(ggplot2)

knepp_accum_plot <- ggplot(knepp_accum_df, 
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
    panel.grid.major   = element_blank(),
    panel.grid.minor   = element_blank(),
    legend.position    = c(0.95, 0.05),   # inside bottom-right
    legend.justification = c("right", "bottom"),
    legend.background  = element_rect(fill = "white", colour = "black"),
    axis.title         = element_text(size = 22, face = "bold", colour = "black"),
    axis.text          = element_text(size = 20, colour = "black"),
    legend.text        = element_text(size = 12, colour = "black")
  )

# Print
knepp_accum_plot





################################################################################
################################################################################
################################################################################


###############################################################################
# NMDS: Passive acoustic vs Point count (Knepp)
###############################################################################

# Recreate an AM site-by-species PA data frame with Location_ID retained
knepp_am_species_df_full <- knepp_am_wide_pa %>%
  select(Location_ID, Tyto_alba:last_col()) %>%
  group_by(Location_ID) %>%
  summarise(across(everything(), ~ sum(. > 0)), .groups = "drop")

# Strict 0/1 presence/absence for AM
knepp_am_pa_df <- knepp_am_species_df_full %>%
  select(-Location_ID) %>%
  mutate(across(everything(), ~ as.integer(. > 0)))

# Point-count PA data frame already exists from your pipeline
knepp_pc_pa_df <- knepp_pc_species_df  # already 0/1

# ---- Align species columns between datasets ----
add_missing_cols <- function(df, all_species) {
  missing <- setdiff(all_species, names(df))
  if (length(missing) > 0) df[missing] <- 0L
  df %>% select(all_of(all_species))
}

knepp_common_species <- union(names(knepp_pc_pa_df),
                              names(knepp_am_pa_df))

knepp_pc_pa_df2 <- add_missing_cols(knepp_pc_pa_df, knepp_common_species)
knepp_am_pa_df2 <- add_missing_cols(knepp_am_pa_df, knepp_common_species)

# ---- Combine with Source + IDs ----
knepp_pc_comb <- knepp_pc_pa_df2 %>%
  mutate(Source = "Point count",
         ID = as.character(knepp_pc_species_matrix$Survey.ID))

knepp_am_comb <- knepp_am_pa_df2 %>%
  mutate(Source = "Passive acoustic",
         ID = as.character(knepp_am_species_df_full$Location_ID))

knepp_combined <- bind_rows(knepp_pc_comb, knepp_am_comb)

# Remove empty sites (rows with all zeros)
knepp_combined_filt <- knepp_combined %>%
  filter(rowSums(select(., -Source, -ID)) > 0)

# ---- Run NMDS (Jaccard on binary PA) ----
set.seed(123)
knepp_nmds <- metaMDS(
  knepp_combined_filt %>% select(-Source, -ID),
  distance = "jaccard",
  binary   = TRUE,
  k        = 2,
  trymax   = 400,
  autotransform = FALSE
)

# ---- Extract site scores ----
knepp_scores <- as.data.frame(scores(knepp_nmds, display = "sites"))
knepp_scores$Source <- knepp_combined_filt$Source
knepp_scores$ID     <- knepp_combined_filt$ID

# Build convex hulls safely (only if group has >= 3 points)
knepp_hulls <- knepp_scores %>%
  group_by(Source) %>%
  filter(n() >= 3) %>%
  slice(chull(NMDS1, NMDS2))

# ---- Plot (match your accumulation plot styling) ----
knepp_nmds_plot <- ggplot(knepp_scores, aes(x = NMDS1, y = NMDS2, colour = Source)) +
  geom_point(size = 4) +
  geom_polygon(
    data = knepp_hulls,
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
    legend.position    = c(0.95, 0.05),      # inside bottom-right
    legend.justification = c("right", "bottom"),
    legend.background  = element_rect(fill = "white", colour = "black"),
    axis.title         = element_text(size = 22, face = "bold", colour = "black"),
    axis.text          = element_text(size = 20, colour = "black"),
    legend.text  = element_text(size = 12, colour = "black")
  )

# Print
knepp_nmds_plot







# species overlap list

# Get unique species
am_species <- unique(knepp_am$Scientific_Name[knepp_pc$Scientific_Name != ""])
pc_species <- unique(knepp_pc$Scientific_Name[knepp_pc$Scientific_Name != ""])

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



#####

################## VENN DIAGRAM PLOTS 
library(grid)  # Needed for turn clipping off
library(ggVennDiagram)

# Prepare list with custom labels
species_list <- list(
  "Passive acoustic" = am_species,
  "Point count" = pc_species
)

# Plot
p <- ggVennDiagram(species_list, label_alpha = 0, lwd = 0.7, category.names = c("Point count",
                                                                              "Passive acoustic")) +
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
  theme(
    legend.position = "none",                      # remove legend
    text = element_text(size = 14, face = "bold"), # all text bold & size 14
    axis.text = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 14, face = "bold"),
    strip.text = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 14, face = "bold")
  )

# Disable clipping so labels are not cropped
gt <- ggplotGrob(p)
gt$layout$clip[gt$layout$name == "panel"] <- "off"
grid.newpage()
grid.draw(gt)
