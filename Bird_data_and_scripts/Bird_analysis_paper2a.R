# Bird analysis
# PAM VS Point count @ Attenborough and Knepp


###############################################################################
#Load libraries
###############################################################################
library(ggplot2)
library(reshape2)
library(colorspace)
library(vegan)
library(dplyr)
library(tidyr)
library(iNEXT)
library(VennDiagram) 
library(eulerr)
library(grid)
library(stringr)
library(patchwork)
library(ggplotify)
library(cowplot)
library(tibble)
library(lubridate)
library(patchwork)
library(performance)


###############################################################################
# Import data
###############################################################################

# Attenborough passive acoustic data (i use "am" in this script which refers to Audiomoth)
attenborough_am <- read.csv("ANRbird-June2024-combinedlocationresults_paper1.csv", header = TRUE)

# Attenborough point count data
attenborough_pc <- read.csv("Attenborough_pointcounts_birds_June_24_paper1.csv", header = TRUE)

# Knepp passive acoustic data
knepp_am <- read.csv("Kneppbird-June2024-combinedlocationresults.csv", header = TRUE)

# Knepp Point count data
knepp_pc <- read.csv("Knepp_pointcounts_0624.csv", header = TRUE)

###############################################################################
# Standardise species names before any analysis
###############################################################################

attenborough_pc$Scientific_Name[attenborough_pc$Scientific_Name == "Sylvia_communis"] <- "Curruca_communis"
attenborough_am$Scientific_Name[attenborough_am$Scientific_Name == "Sylvia_communis"] <- "Curruca_communis"
knepp_pc$Scientific_Name[knepp_pc$Scientific_Name == "Sylvia_communis"]               <- "Curruca_communis"
knepp_am$Scientific_Name[knepp_am$Scientific_Name == "Sylvia_communis"]               <- "Curruca_communis"
###############################################################################
# Inspect
###############################################################################

str(attenborough_pc)
str(attenborough_am)

str(knepp_pc)
str(knepp_am)

unique(knepp_am$Location_ID)
unique(knepp_pc$Survey.ID)

###############################################################################
# Convert columns to factors and change sylvia communis to currica communis (same species)
###############################################################################

attenborough_pc$Scientific_Name <- as.character(attenborough_pc$Scientific_Name)
attenborough_am$Scientific_Name <- as.character(attenborough_am$Scientific_Name)
knepp_pc$Scientific_Name        <- as.character(knepp_pc$Scientific_Name)
knepp_am$Scientific_Name        <- as.character(knepp_am$Scientific_Name)

attenborough_pc$Scientific_Name[attenborough_pc$Scientific_Name == "Sylvia_communis"] <- "Curruca_communis"
attenborough_am$Scientific_Name[attenborough_am$Scientific_Name == "Sylvia_communis"] <- "Curruca_communis"
knepp_pc$Scientific_Name[knepp_pc$Scientific_Name == "Sylvia_communis"]               <- "Curruca_communis"
knepp_am$Scientific_Name[knepp_am$Scientific_Name == "Sylvia_communis"]               <- "Curruca_communis"

attenborough_pc$Scientific_Name <- as.factor(attenborough_pc$Scientific_Name)
attenborough_am$Scientific_Name <- as.factor(attenborough_am$Scientific_Name)
knepp_pc$Scientific_Name        <- as.factor(knepp_pc$Scientific_Name)
knepp_am$Scientific_Name        <- as.factor(knepp_am$Scientific_Name)

###############################################################################
# Look at number of species detected
###############################################################################
length(unique(attenborough_am$Scientific_Name))
length(unique(attenborough_pc$Scientific_Name))

length(unique(knepp_am$Scientific_Name))
length(unique(knepp_pc$Scientific_Name))


###############################################################################
# Add presence absence column and fill for all species in each data set
###############################################################################

attenborough_pc <- attenborough_pc %>% mutate(Presence_Absence = ifelse(Count > 0, 1, 0))    

attenborough_am <- attenborough_am %>% mutate(Presence_Absence = ifelse(Month_ID > 0, 1, 0)) 

knepp_pc <- knepp_pc %>% mutate(Presence_Absence = ifelse(Count > 0, 1, 0))    

knepp_am <- knepp_am %>% mutate(Presence_Absence = ifelse(Month_ID > 0, 1, 0)) 


###############################################################################
# Make presence absence matrix
###############################################################################

#Attenborough passive acoustic 
attenborough_am_wide_pa <- attenborough_am %>%
  pivot_wider(
    names_from  = Scientific_Name,          
    values_from = Presence_Absence,         
    values_fill = 0                         
  )

attenborough_am_species_matrix <- attenborough_am_wide_pa %>%
  select(Location_ID, Tyto_alba:last_col()) %>%
  group_by(Location_ID) %>%
  summarise(across(everything(), ~ sum(. > 0)), .groups = "drop")

attenborough_am_species_matrix <- as.matrix(attenborough_am_species_matrix %>% select(-Location_ID))

str(attenborough_am_species_matrix)# 75 columns match with 75 species 

# Attenborough point count
attenborough_pc_wide_pa <- attenborough_pc %>%
  pivot_wider(
    names_from  = Scientific_Name,          
    values_from = Presence_Absence,         
    values_fill = 0,                        
    values_fn   = ~max(.x, na.rm = TRUE)  
  )

attenborough_pc_species_matrix <- attenborough_pc_wide_pa %>%
  select(Survey.ID, Chroicocephalus_ridibundus:last_col()) %>%
  group_by(Survey.ID) %>%
  summarise(across(everything(), ~ sum(. > 0)), .groups = "drop")


attenborough_pc_species_df <- attenborough_pc_species_matrix %>% select(-Survey.ID)


attenborough_pc_species_df <- attenborough_pc_species_df %>%
  mutate(across(everything(), ~ as.integer(. > 0))) %>%
  select(where(~ all(!is.na(.))))


# knepp PAM 

knepp_am_wide_pa <- knepp_am %>%
  pivot_wider(
    names_from  = Scientific_Name,          
    values_from = Presence_Absence,         
    values_fill = 0                         
  )

knepp_am_species_matrix <- knepp_am_wide_pa %>%
  select(Location_ID, Tyto_alba:last_col()) %>%
  group_by(Location_ID) %>%
  summarise(across(everything(), ~ sum(. > 0)), .groups = "drop")

knepp_am_species_matrix <- as.matrix(knepp_am_species_matrix %>% select(-Location_ID))

# Knepp point counts 

knepp_pc_wide_pa <- knepp_pc %>%
  pivot_wider(
    names_from  = Scientific_Name,          
    values_from = Presence_Absence,         
    values_fill = 0,                       
    values_fn   = ~max(.x, na.rm = TRUE)    
  )

knepp_pc_species_matrix <- knepp_pc_wide_pa %>%
  select(Survey.ID, Turdus_merula:last_col()) %>%
  group_by(Survey.ID) %>%
  summarise(across(everything(), ~ sum(. > 0)), .groups = "drop")

knepp_pc_species_df <- knepp_pc_species_matrix %>% select(-Survey.ID)

knepp_pc_species_df <- knepp_pc_species_df %>%
  mutate(across(everything(), ~ as.integer(. > 0))) %>%
  select(where(~ all(!is.na(.))))


###############################################################################
# Species accumilation curves with iNEXT 
###############################################################################

# Attenborough first
# iNEXT needs a species by sites incidence matrix
mat_am_att <- t(as.matrix(attenborough_am_species_matrix))
mat_pc_att <- t(as.matrix(attenborough_pc_species_df))

storage.mode(mat_am_att) <- "integer"
storage.mode(mat_pc_att) <- "integer"
mat_am_att[mat_am_att > 1] <- 1L
mat_pc_att[mat_pc_att > 1] <- 1L

# Shared extrapolation endpoint of twice max effort
endpoint_att <- 2 * max(ncol(mat_am_att), ncol(mat_pc_att))

attenborough_bird_inext <- iNEXT(
  list(`Passive acoustic` = mat_am_att,
       `Point count`      = mat_pc_att),
  q        = 0,
  datatype = "incidence_raw",
  endpoint = endpoint_att,
  knots    = 80,
  se       = TRUE,
  conf     = 0.95
)


# Knepp 

knepp_am_mat_sitespp <- as.matrix(knepp_am_species_matrix)
knepp_pc_mat_sitespp <- as.matrix(knepp_pc_species_df)

storage.mode(knepp_am_mat_sitespp) <- "integer"
storage.mode(knepp_pc_mat_sitespp) <- "integer"
knepp_am_mat_sitespp[knepp_am_mat_sitespp > 1] <- 1L
knepp_pc_mat_sitespp[knepp_pc_mat_sitespp > 1] <- 1L

mat_am_knepp <- t(knepp_am_mat_sitespp)
mat_pc_knepp <- t(knepp_pc_mat_sitespp)

endpoint_knepp <- 2 * max(ncol(mat_am_knepp), ncol(mat_pc_knepp))

knepp_bird_inext <- iNEXT(
  list(`Passive acoustic` = mat_am_knepp,
       `Point count`      = mat_pc_knepp),
  q        = 0,
  datatype = "incidence_raw",
  endpoint = endpoint_knepp,
  knots    = 80,
  se       = TRUE,
  conf     = 0.95
)

###############################################################################
# Plot species accumilation curves (Figures 2a and 2b in manuscript - at the moment anyway)
###############################################################################

# Format and order stuff for both plots 

bird_cols <- c(
  `Passive acoustic` = "sienna",
  `Point count`      = "cadetblue4"
)

leg_breaks <- c("Passive acoustic", "Point count")
leg_labels <- c("PAM", "Point count")

# Attenborough plot first

attenborough_accum_plot <- ggiNEXT(attenborough_bird_inext, type = 1, color.var = "Assemblage") +
  scale_colour_manual(
    values = bird_cols,
    breaks = leg_breaks,
    labels = leg_labels
  ) +
  scale_fill_manual(values = bird_cols) +
  guides(
    fill     = "none",   
    linetype = "none",
    shape    = "none",
    colour   = guide_legend(
      override.aes = list(
        linetype = 1,
        shape    = NA,
        fill     = unname(bird_cols[leg_breaks]),
        alpha    = 0.25
      )
    )
  ) +
  labs(
    title  = NULL,
    x      = "Number of sampling events",
    y      = "Species richness",
    colour = NULL,
    fill   = NULL
  ) +
  scale_x_continuous(
    limits = c(0, 90),
    breaks = seq(0, 80, by = 20),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    limits = c(0, 90),
    breaks = seq(0, 90, by = 20),
    expand = c(0, 0)
  )+
  theme_bw(base_size = 12) +
  theme(
    panel.grid.major     = element_blank(),
    panel.grid.minor     = element_blank(),
    legend.position      = c(0.95, 0.05),
    legend.justification = c("right", "bottom"),
    legend.background    = element_rect(fill = NA, colour = NA),
    axis.title           = element_text(size = 12, colour = "black"),
    axis.text            = element_text(size = 12, colour = "black"),
    legend.text          = element_text(size = 12, colour = "black"),
    legend.title         = element_blank(),
    plot.margin          = margin(8, 8, 20, 8)
  )


# ATTENBOROUGH ACCUM PLOT (figures 2a)
attenborough_accum_plot


# Knepp plot

knepp_accum_plot <- ggiNEXT(knepp_bird_inext, type = 1, color.var = "Assemblage") +
  scale_colour_manual(
    values = bird_cols,
    breaks = leg_breaks,
    labels = leg_labels
  ) +
  scale_fill_manual(values = bird_cols) +
  guides(
    fill     = "none",   
    linetype = "none",
    shape    = "none",
    colour   = guide_legend(
      override.aes = list(
        linetype = 1,
        shape    = NA,
        # show ribbon + line in the SAME legend key
        fill     = unname(bird_cols[leg_breaks]),
        alpha    = 0.25
      )
    )
  ) +
  labs(
    title  = NULL,
    x      = "Number of sampling events",
    y      = "Species richness",
    colour = NULL,
    fill   = NULL
  ) +
  scale_x_continuous(
    limits = c(0, 82),
    breaks = seq(0, 80, by = 20),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    limits = c(0, 90),
    breaks = seq(0, 90, by = 20),
    expand = c(0, 0)
  )+
  theme_bw(base_size = 12) +
  theme(
    panel.grid.major     = element_blank(),
    panel.grid.minor     = element_blank(),
    legend.position      = "none",
    axis.title           = element_text(size = 12, colour = "black"),
    axis.text            = element_text(size = 12, colour = "black"),
    plot.margin          = margin(8, 8, 20, 8)
  )

# Knepp birds accum plot figure 2b
knepp_accum_plot



################################################################################
# NMDS plot for Attenborough
################################################################################

# Attenborough: make species-by-site PA matrix with Location_ID retained
attenborough_am_pa_df <- attenborough_am_wide_pa %>%
  select(Location_ID, Tyto_alba:last_col()) %>%
  group_by(Location_ID) %>%
  summarise(
    across(everything(), ~ as.integer(any(. > 0))),
    .groups = "drop"
  ) %>%
  select(-Location_ID)

# Point-count PA dataframe already exists
attenborough_pc_pa_df <- attenborough_pc_species_df

# Align species columns for passive acoustic and point count data sets
attenborough_common_species <- union(
  names(attenborough_pc_pa_df),
  names(attenborough_am_pa_df)
)

missing_pc <- setdiff(
  attenborough_common_species,
  names(attenborough_pc_pa_df)
)

missing_am <- setdiff(
  attenborough_common_species,
  names(attenborough_am_pa_df)
)

if (length(missing_pc) > 0) {
  attenborough_pc_pa_df[missing_pc] <- 0L
}

if (length(missing_am) > 0) {
  attenborough_am_pa_df[missing_am] <- 0L
}

attenborough_pc_pa_df <- attenborough_pc_pa_df %>%
  select(all_of(attenborough_common_species))

attenborough_am_pa_df <- attenborough_am_pa_df %>%
  select(all_of(attenborough_common_species))

# Combine datasets with source labels and IDs
attenborough_combined_filt <- bind_rows(
  
  attenborough_pc_pa_df %>%
    mutate(
      Source = "Point count",
      ID = as.character(attenborough_pc_species_matrix$Survey.ID)
    ),
  
  attenborough_am_pa_df %>%
    mutate(
      Source = "PAM",
      ID = as.character(unique(attenborough_am_wide_pa$Location_ID))
    )
  
) %>%
  filter(
    rowSums(select(., -Source, -ID)) > 0
  )

# Run  NMDS 
set.seed(123)
species_mat <- attenborough_combined_filt %>% select(-Source, -ID)  # numeric only

attenborough_nmds <- metaMDS(
  species_mat,
  distance = "jaccard",
  binary   = TRUE,
  k        = 2,
  trymax   = 400,
  autotransform = FALSE
)

# Extract site scores 
attenborough_scores <- as.data.frame(scores(attenborough_nmds, display = "sites"))
attenborough_scores$Source <- attenborough_combined_filt$Source
attenborough_scores$ID     <- attenborough_combined_filt$ID

attenborough_stress_bird <- round(attenborough_nmds$stress, 3)

#  Plot NMDS for Attenborough
attenborough_nmds_plot <- ggplot(
  attenborough_scores,
  aes(x = NMDS1, y = NMDS2, colour = Source)
) +
  geom_point(size = 0.8) +
  stat_ellipse(
    aes(fill = Source, group = Source),
    geom = "polygon",
    alpha = 0.25,
    linewidth = 0.8,
    type = "t",
    level = 0.95,
    show.legend = FALSE
  ) +
  scale_colour_manual(
    values = c(
      "PAM" = "sienna",
      "Point count" = "cadetblue4"
    )
  ) +
  scale_fill_manual(
    values = c(
      "PAM" = "sienna",
      "Point count" = "cadetblue4"
    )
  ) +
  labs(x = "NMDS1", y = "NMDS2", colour = NULL, fill = NULL) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    legend.position = c(0.999, 0.999),
    legend.justification = c("right", "top"),
    legend.background = element_rect(fill = NA, colour = NA),
    legend.key = element_rect(fill = NA, colour = NA),
    legend.key.width = unit(0.1, "cm"),
    legend.spacing.y = unit(0.01, "cm"),
    axis.title = element_text(size = 12, colour = "black"),
    axis.text = element_text(size = 12, colour = "black"),
    legend.text = element_text(size = 12, colour = "black"),
    plot.background  = element_rect(fill = NA, colour = NA),
    panel.background = element_rect(fill = NA, colour = NA)# outside panel
  ) +
  guides(colour = guide_legend(override.aes = list(size = 4, shape = 15))) +
  annotate(
    "text",
    x = -Inf, y = -Inf,
    label = paste0("Stress = ", attenborough_stress_bird),
    hjust = -0.1,
    vjust = -1,
    size = 3.3,
    colour = "black"
  )

# Attenborough NMDS plot (Figure 2c)
attenborough_nmds_plot



################################################################################
# NMDS plot for Knepp
################################################################################

knepp_am_pa_df <- knepp_am_wide_pa %>%
  select(Location_ID, Tyto_alba:last_col()) %>%
  group_by(Location_ID) %>%
  summarise(
    across(everything(), ~ as.integer(any(. > 0))),
    .groups = "drop"
  ) %>%
  select(-Location_ID)


knepp_pc_pa_df <- knepp_pc_species_df

# Align species columns between datasets
knepp_common_species <- union(
  names(knepp_pc_pa_df),
  names(knepp_am_pa_df)
)

missing_pc <- setdiff(
  knepp_common_species,
  names(knepp_pc_pa_df)
)

missing_am <- setdiff(
  knepp_common_species,
  names(knepp_am_pa_df)
)

if (length(missing_pc) > 0) {
  knepp_pc_pa_df[missing_pc] <- 0L
}

if (length(missing_am) > 0) {
  knepp_am_pa_df[missing_am] <- 0L
}

knepp_pc_pa_df <- knepp_pc_pa_df %>%
  select(all_of(knepp_common_species))

knepp_am_pa_df <- knepp_am_pa_df %>%
  select(all_of(knepp_common_species))

# Combine datasets with source labels and IDs
knepp_combined_filt <- bind_rows(
  
  knepp_pc_pa_df %>%
    mutate(
      Source = "Point count",
      ID = as.character(knepp_pc_species_matrix$Survey.ID)
    ),
  
  knepp_am_pa_df %>%
    mutate(
      Source = "Passive acoustic",
      ID = as.character(unique(knepp_am_wide_pa$Location_ID))
    )
  
) %>%
  filter(
    rowSums(select(., -Source, -ID)) > 0
  )


# Run NMDS
set.seed(123)
species_mat_knepp <- knepp_combined_filt %>% select(-Source, -ID)

knepp_nmds <- metaMDS(
  species_mat_knepp,
  distance = "jaccard",
  binary   = TRUE,
  k        = 2,
  trymax   = 400,
  autotransform = FALSE
)


# Extract  scores
knepp_scores <- as.data.frame(scores(knepp_nmds, display = "sites"))
knepp_scores$Source <- knepp_combined_filt$Source
knepp_scores$Source <- ifelse(knepp_scores$Source == "Passive acoustic", "PAM", knepp_scores$Source)
knepp_scores$Source <- factor(knepp_scores$Source, levels = c("PAM", "Point count"))

knepp_scores$ID  <- knepp_combined_filt$ID

knepp_stress_bird <- round(knepp_nmds$stress, 3)

# make plot

knepp_nmds_plot <- ggplot(
  knepp_scores,
  aes(x = NMDS1, y = NMDS2, colour = Source)
) +
  geom_point(size = 0.8) +
  stat_ellipse(
    aes(fill = Source, group = Source),
    geom = "polygon",
    alpha = 0.25,
    linewidth = 0.8,
    type = "t",
    level = 0.95,
    show.legend = FALSE
  ) +
  scale_colour_manual(
    values = c(
      "PAM" = "sienna",
      "Point count" = "cadetblue4"
    )
  ) +
  scale_fill_manual(
    values = c(
      "PAM" = "sienna",
      "Point count" = "cadetblue4"
    )
  ) +
  labs(x = "NMDS1", y = "NMDS2", colour = NULL, fill = NULL) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    legend.position = "none",
    axis.title = element_text(size = 12, colour = "black"),
    axis.text = element_text(size = 12, colour = "black"),
    plot.background  = element_rect(fill = NA, colour = NA),
    panel.background = element_rect(fill = NA, colour = NA)
  ) +
  guides(colour = guide_legend(override.aes = list(size = 4, shape = 15))) +
  annotate(
    "text",
    x = -Inf, y = -Inf,
    label = paste0("Stress = ", knepp_stress_bird),
    hjust = -0.1,
    vjust = -1,
    size = 3.3,
    colour = "black"
  )

#Knepp NMDS plot (figure 2d)
knepp_nmds_plot



################################################################################
# PERMANOVA and Multivariate dispersion test (PERMDISP) for Attenborough
################################################################################

dist_mat <- vegdist(species_mat, method = "jaccard", binary = TRUE)

env <- attenborough_combined_filt %>%
  select(Source, ID) %>%
  mutate(Source = factor(Source, levels = c("Point count", "PAM")))

set.seed(123)
adon <- adonis2(dist_mat ~ Source, data = env, permutations = 999)
print(adon)

bd <- betadisper(dist_mat, env$Source)
bd_anova <- anova(bd)
print(bd_anova)

# Extract raw p-values
attenborough_permanova_p <- adon$`Pr(>F)`[1]
attenborough_permdisp_p  <- bd_anova$`Pr(>F)`[1]

dist_df_att <- data.frame(
  Distance = bd$distances,
  Method = env$Source
)

dist_df_att_summary <- dist_df_att %>%
  dplyr::group_by(Method) %>%
  dplyr::summarise(
    mean_distance = mean(Distance, na.rm = TRUE),
    sd_distance   = sd(Distance, na.rm = TRUE),
    se_distance   = sd_distance / sqrt(n()),
    ci_lower      = mean_distance - 1.96 * se_distance,
    ci_upper      = mean_distance + 1.96 * se_distance,
    n             = n(),
    .groups = "drop"
  )

dist_df_att_summary

distance_plot_att <- ggplot(dist_df_att, aes(x = Method, y = Distance, fill = Method)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.1, alpha = 0.6) +
  theme_bw() +
  scale_colour_manual(
    values = c(
      "PAM" = "sienna",
      "Point count" = "cadetblue4"
    )
  ) +
  scale_fill_manual(
    values = c(
      "PAM" = "sienna",
      "Point count" = "cadetblue4"
    )
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    title = "Attenborough",
    y = "Distance to centroid",
    x = "Method"
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    axis.title.x = element_text(size = 14, colour = "black"),
    axis.title.y = element_text(size = 12, colour = "black"),
    axis.text.x = element_text(size = 12, colour = "black"),
    axis.text.y = element_text(size = 12, colour = "black"),
    text = element_text(size = 12)
  ) +
  annotate("text", x = 1, y = 0.85, label = "n = 40", size = 3) +
  annotate("text", x = 1, y = 0.9, label = "Mean = 0.53 (95% CI [0.50, 0.57])", size = 3) +
  annotate("text", x = 2, y = 0.6, label = "n = 43", size = 3) +
  annotate("text", x = 2, y = 0.65, label = "Mean = 0.32 (95% CI [0.29, 0.34])", size = 3)

bd_pw <- permutest(bd, pairwise = TRUE, permutations = 999)
print(bd_pw)




################################################################################
# PERMANOVA and Multivariate dispersion test (PERMDISP) for Knepp
################################################################################

dist_mat_knepp <- vegdist(species_mat_knepp, method = "jaccard", binary = TRUE)

env_knepp <- knepp_combined_filt %>%
  select(Source, ID) %>%
  mutate(
    Source = if_else(Source == "Passive acoustic", "PAM", Source),
    Source = factor(Source, levels = c("Point count", "PAM"))
  )

set.seed(123)
adon_knepp <- adonis2(dist_mat_knepp ~ Source, data = env_knepp, permutations = 999)
print(adon_knepp)

knepp_bd <- betadisper(dist_mat_knepp, env_knepp$Source)
knepp_bd_anova <- anova(knepp_bd)
print(knepp_bd_anova)

# Extract raw p-values
knepp_permanova_p <- adon_knepp$`Pr(>F)`[1]
knepp_permdisp_p  <- knepp_bd_anova$`Pr(>F)`[1]

knepp_bd_pw <- permutest(knepp_bd, pairwise = TRUE, permutations = 999)
print(knepp_bd_pw)

dist_df_knepp <- data.frame(
  Distance = knepp_bd$distances,
  Method = env_knepp$Source
)

dist_df_knepp_summary <- dist_df_knepp %>%
  dplyr::group_by(Method) %>%
  dplyr::summarise(
    mean_distance = mean(Distance, na.rm = TRUE),
    sd_distance   = sd(Distance, na.rm = TRUE),
    se_distance   = sd_distance / sqrt(n()),
    ci_lower      = mean_distance - 1.96 * se_distance,
    ci_upper      = mean_distance + 1.96 * se_distance,
    n             = n(),
    .groups = "drop"
  )

dist_df_knepp_summary

distance_plot_knepp <- ggplot(dist_df_knepp, aes(x = Method, y = Distance, fill = Method)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.1, alpha = 0.6) +
  theme_bw() +
  scale_colour_manual(
    values = c(
      "PAM" = "sienna",
      "Point count" = "cadetblue4"
    )
  ) +
  scale_fill_manual(
    values = c(
      "PAM" = "sienna",
      "Point count" = "cadetblue4"
    )
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    title = "Knepp",
    y = "Distance to centroid",
    x = "Method"
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    axis.title.x = element_text(size = 14, colour = "black"),
    axis.title.y = element_text(size = 12, colour = "black"),
    axis.text.x = element_text(size = 12, colour = "black"),
    axis.text.y = element_text(size = 12, colour = "black"),
    text = element_text(size = 12)
  ) +
  annotate("text", x = 1, y = 0.85, label = "n = 36", size = 3) +
  annotate("text", x = 1, y = 0.9, label = "Mean = 0.58 (95% CI [0.56, 0.61])", size = 3) +
  annotate("text", x = 2, y = 0.6, label = "n = 37", size = 3) +
  annotate("text", x = 2, y = 0.65, label = "Mean = 0.30 (95% CI [0.28, 0.32])", size = 3)


################################################################################
# Holm correction for PERMANOVA and PERMDISP p-values
################################################################################

raw_p_values <- c(
  Attenborough_PERMANOVA = attenborough_permanova_p,
  Knepp_PERMANOVA        = knepp_permanova_p,
  Attenborough_PERMDISP  = attenborough_permdisp_p,
  Knepp_PERMDISP         = knepp_permdisp_p
)

holm_p_values <- p.adjust(raw_p_values, method = "holm")

p_value_summary <- data.frame(
  Test = names(raw_p_values),
  Raw_p = raw_p_values,
  Holm_adjusted_p = holm_p_values,
  row.names = NULL
)

p_value_summary


################################################################################
# Combine distance plots
################################################################################

combined_plot_distance <- (distance_plot_att / distance_plot_knepp) +
  plot_annotation(
    tag_levels = "a",
    tag_prefix = "(",
    tag_suffix = ")"
  )

combined_plot_distance
library(eulerr)
library(grid)
library(patchwork)

################################################################################
# Make euler plot for Attenborough
################################################################################

am_species <- unique(attenborough_am$Scientific_Name[
  !is.na(attenborough_am$Scientific_Name) &
    attenborough_am$Scientific_Name != ""
])

pc_species <- unique(attenborough_pc$Scientific_Name[
  !is.na(attenborough_pc$Scientific_Name) &
    attenborough_pc$Scientific_Name != ""
])

only_am <- setdiff(am_species, pc_species)
only_pc <- setdiff(pc_species, am_species)
both <- intersect(am_species, pc_species)

species_summary_attenborough <- data.frame(
  Species = c(only_am, only_pc, both),
  Detected_By = c(
    rep("AM only", length(only_am)),
    rep("PC only", length(only_pc)),
    rep("Both", length(both))
  )
) %>%
  filter(!is.na(Species), Species != "")

VennDiag_Att <- euler(c(
  "PAM" = length(only_am),
  "Point count" = length(only_pc),
  "PAM&Point count" = length(both)
))

attenborough_venn_base <- plot(
  VennDiag_Att,
  quantities = list(
    type = "counts",
    fontsize = 11,
    col = "black",
    font = 2
  ),
  labels = FALSE,
  fills = list(
    fill = c("sienna", "cadetblue4"),
    alpha = 0.6
  ),
  edges = list(
    col = "black",
    lwd = 0.7
  ),
  main = NULL,
  newpage = FALSE
)

attenborough_venn_grob <- grobTree(
  attenborough_venn_base,
  
  textGrob(
    "PAM",
    x = 0.13,
    y = 0.52,
    gp = gpar(
      col = "black",
      fontsize = 11,
      fontface = "bold"
    )
  ),
  
  textGrob(
    "Point count",
    x = 0.83,
    y = 0.63,
    gp = gpar(
      col = "black",
      fontsize = 11,
      fontface = "bold"
    )
  ),
  
  textGrob(
    "Both",
    x = 0.58,
    y = 0.56,
    gp = gpar(
      col = "black",
      fontsize = 11,
      fontface = "bold"
    )
  )
)

attenborough_venn_plot <- wrap_elements(full = attenborough_venn_grob)


################################################################################
# Make euler plot for Knepp
################################################################################

am_species <- unique(knepp_am$Scientific_Name[
  !is.na(knepp_am$Scientific_Name) &
    knepp_am$Scientific_Name != ""
])

pc_species <- unique(knepp_pc$Scientific_Name[
  !is.na(knepp_pc$Scientific_Name) &
    knepp_pc$Scientific_Name != ""
])

only_am <- setdiff(am_species, pc_species)
only_pc <- setdiff(pc_species, am_species)
both <- intersect(am_species, pc_species)

species_summary_knepp <- data.frame(
  Species = c(only_am, only_pc, both),
  Detected_By = c(
    rep("AM only", length(only_am)),
    rep("PC only", length(only_pc)),
    rep("Both", length(both))
  )
) %>%
  filter(!is.na(Species), Species != "")

VennDiag_knepp <- euler(c(
  "PAM" = length(only_am),
  "Point count" = length(only_pc),
  "PAM&Point count" = length(both)
))

knepp_venn_base <- plot(
  VennDiag_knepp,
  quantities = list(
    type = "counts",
    fontsize = 11,
    col = "black",
    font = 2
  ),
  labels = FALSE,
  fills = list(
    fill = c("sienna", "#967B6D"),
    alpha = 0.6
  ),
  edges = list(
    col = "black",
    lwd = 0.7
  ),
  main = NULL,
  newpage = FALSE
)

knepp_venn_grob <- grobTree(
  knepp_venn_base,
  
  textGrob(
    "PAM",
    x = 0.09,
    y = 0.53,
    gp = gpar(
      col = "black",
      fontsize = 11,
      fontface = "bold"
    )
  ),
  
  textGrob(
    "Both",
    x = 0.57,
    y = 0.55,
    gp = gpar(
      col = "black",
      fontsize = 11,
      fontface = "bold"
    )
  )
)

knepp_venn_plot <- wrap_elements(full = knepp_venn_grob)


###########################################################################
# Multipanel (Figure 2a-f) for all bird plots in main manuscript
###########################################################################

panel_main <- (
  (attenborough_accum_plot | attenborough_venn_plot | attenborough_nmds_plot) /
    (knepp_accum_plot      | knepp_venn_plot        | knepp_nmds_plot)
) &
  theme(
    plot.margin = margin(10, 10, 10, 10)
  )

# draw everything and manual labels - att first 3 and knepp second 3
bird_panel_labeled <- ggdraw(panel_main) +
  draw_label("(a)", x = 0.03, y = 0.98, hjust = 0, vjust = 1, fontface = "bold", size = 14) +
  draw_label("(c)", x = 0.36, y = 0.98, hjust = 0, vjust = 1, fontface = "bold", size = 14) +
  draw_label("(e)", x = 0.66, y = 0.98, hjust = 0, vjust = 1, fontface = "bold", size = 14) +
  
  draw_label("(b)", x = 0.03, y = 0.51, hjust = 0, vjust = 1, fontface = "bold", size = 14) +
  draw_label("(d)", x = 0.36, y = 0.51, hjust = 0, vjust = 1, fontface = "bold", size = 14) +
  draw_label("(f)", x = 0.66, y = 0.51, hjust = 0, vjust = 1, fontface = "bold", size = 14) +
  
  draw_label("Attenborough", x = 0.015, y = 0.77, angle = 90,
             fontface = "bold", size = 18, hjust = 0.5, vjust = 0.1) +
  draw_label("Knepp", x = 0.015, y = 0.26, angle = 90,
             fontface = "bold", size = 18, hjust = 0.3, vjust = 0.1)

# final multipanel plot for birds
bird_panel_labeled

w <- dev.size("in")[1]
h <- dev.size("in")[2] 

ggsave(
  "Bird_paper1_figure.tiff",
  plot = bird_panel_labeled,
  width = w,
  height = h,
  units = "in",
  dpi = 600,
  compression = "lzw"
)




###############################################################################
# SUPPLEMENTARY ANALYSIS - SHARED 10 MINS SPECIES RICHNESS
###############################################################################

###############################################################################
# iMPORT Data
###############################################################################

# Attenborough Audiomoth
attenborough_am_time <- read.csv(
  "20240605_att_bird_PA_combined_speciestime.csv",
  header = TRUE
)

# Attenborough point counts
attenborough_pc <- read.csv(
  "Attenborough_pointcounts_birds_June_24_paper1.csv",
  header = TRUE
)

# Knepp Audiomoth
knepp_am_time <- read.csv(
  "20240626_knepp_bird_PA_combined_speciestime.csv",
  header = TRUE
)

# Knepp point counts
knepp_pc <- read.csv(
  "Knepp_pointcounts_0624.csv",
  header = TRUE
)


str(attenborough_am_time)
str(attenborough_pc)
str(knepp_am_time)
str(knepp_pc)

###############################################################################
# Standardise species names before any analysis
###############################################################################

# In the raw data Sylvia communis and Curruca communis are used, but they are the same species.
# Use Curruca communis consistently in Audiomoth data,
# and Curruca_communis consistently in point count data.

# Audiomoth data: scientific names use spaces
attenborough_am_time$Scientific.Name[
  attenborough_am_time$Scientific.Name == "Sylvia communis"
] <- "Curruca communis"

knepp_am_time$Scientific.Name[
  knepp_am_time$Scientific.Name == "Sylvia communis"
] <- "Curruca communis"


# Point count data: scientific names use underscores
attenborough_pc$Scientific_Name[
  attenborough_pc$Scientific_Name == "Sylvia_communis"
] <- "Curruca_communis"

knepp_pc$Scientific_Name[
  knepp_pc$Scientific_Name == "Sylvia_communis"
] <- "Curruca_communis"

###############################################################################
# ATTENBOROUGH ACOUSTIC: RICHNESS PER 10 MIN 
###############################################################################

att_am_10min <- attenborough_am_time %>%
  mutate(
    Date         = ymd(Date),
    Time         = hms(Time),
    datetime_utc = Date + Time,
    datetime_bst = datetime_utc + hours(1), # british  summer time means add 1 hour to this
    slot_start   = floor_date(datetime_bst, unit = "10 minutes"),
    time_of_day  = format(slot_start, "%H:%M")
  ) %>%
  group_by(Location, Date, time_of_day) %>%
  summarise(
    n_species = n_distinct(Scientific.Name),
    .groups = "drop"
  )

###############################################################################
# ATTENBOROUGH POINT COUNTS: MATCH TO 10 MIN SLOT
###############################################################################

att_pc_10min <- attenborough_pc %>%
  mutate(
    Date       = dmy(Date),
    Time_str   = sprintf("%04d", Time),
    hour_num   = as.integer(substr(Time_str, 1, 2)),
    minute_num = as.integer(substr(Time_str, 3, 4)),
    Time_posix = make_datetime(
      year(Date), month(Date), day(Date),
      hour = hour_num, min = minute_num
    ),
    slot_start = floor_date(Time_posix, unit = "10 minutes"),
    time_of_day = format(slot_start, "%H:%M")
  ) %>%
  group_by(Survey.ID, Date, time_of_day) %>%
  summarise(
    n_species = n_distinct(Scientific_Name),   
    .groups = "drop"
  )

###############################################################################
# ATTENBOROUGH: SHARED TIME SLOTS, MEAN AND SD
###############################################################################

shared_slots_att <- intersect(
  unique(att_am_10min$time_of_day),
  unique(att_pc_10min$time_of_day)
)

att_am_shared <- att_am_10min %>%
  filter(time_of_day %in% shared_slots_att)

att_pc_shared <- att_pc_10min %>%
  filter(time_of_day %in% shared_slots_att)

combined_att_10min_data <- bind_rows( att_am_shared %>%
                                        rename(ID = Location) %>%
                                        mutate(method = "Passive acoustic"),
                                      att_pc_shared %>%
                                        rename(ID = Survey.ID) %>%
                                        mutate(method = "Point count"))


combined_att_10min_data %>%
  dplyr::group_by(method) %>%
  dplyr::summarise(
    mean_species = mean(n_species, na.rm = TRUE),
    sd_species   = sd(n_species, na.rm = TRUE),
    se_species   = sd_species / sqrt(n()),
    ci_lower     = mean_species - 1.96 * se_species,
    ci_upper     = mean_species + 1.96 * se_species,
    n            = n(),
    .groups = "drop"
  )

combined_att_10min_data <- combined_att_10min_data %>%
  dplyr::mutate(
    method = factor(
      method,
      levels = c("Point count", "Passive acoustic"),
      labels = c("Point count", "PAM")
    )
  )

att_bird_10min_plot <-ggplot(combined_att_10min_data, aes(x = factor(method), y = n_species, fill = method)) +
  geom_boxplot(
    width = 0.3,
    alpha = 1,
    outlier.shape = NA ) +
  geom_jitter(
    aes(x = as.numeric(factor(method)) - 0.3),
    width = 0.1,
    height = 0,
    pch = 21,
    alpha = 0.4,
    size = 2) +
  coord_cartesian(ylim = c(0, 18)) +
  labs(
    x = NULL,
    y = "Species richness per 10 min") +
  theme_bw() +
  scale_colour_manual(
    values = c(
      "Point count" = "cadetblue4",
      "PAM" = "sienna"
      
    )
  ) +
  scale_fill_manual(
    values = c(
      "Point count" = "cadetblue4",
      "PAM" = "sienna"
    )
  ) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    axis.title.x = element_text(size = 14, colour = "black"),
    axis.title.y = element_text(size = 12, colour = "black"),
    axis.text.x = element_text(size = 12, colour = "black"),
    axis.text.y = element_text(size = 12, colour = "black"),
    text = element_text(size = 12)
  ) +
  annotate("text", x = 2, y = 16.5, label = "n = 8672", size = 3) +
  annotate("text", x = 2, y = 17.5, label = "Mean = 3.83", size = 3) +
  annotate("text", x = 1, y = 16.5, label = "n = 40", size = 3)+
  annotate("text", x = 1, y = 17.5, label = "Mean = 8.95", size = 3) 

###############################################################################
# KNEPP ACOUSTIC: RICHNESS PER 10 MIN SLOT, UTC CONVERTED TO BST
###############################################################################

knepp_am_10min <- knepp_am_time %>%
  mutate(
    Date         = ymd(Date),
    Time         = hms(Time),
    datetime_utc = Date + Time,
    datetime_bst = datetime_utc + hours(1),
    slot_start   = floor_date(datetime_bst, unit = "10 minutes"),
    time_of_day  = format(slot_start, "%H:%M")
  ) %>%
  group_by(Location, Date, time_of_day) %>%
  summarise(
    n_species = n_distinct(Scientific.Name),
    .groups = "drop"
  )

###############################################################################
# KNEPP POINT COUNTS: MATCH TO 10 MIN SLOT
###############################################################################

knepp_pc_10min <- knepp_pc %>%
  mutate(
    Date       = dmy(Date),
    Time_str   = sprintf("%04d", Time),
    hour_num   = as.integer(substr(Time_str, 1, 2)),
    minute_num = as.integer(substr(Time_str, 3, 4)),
    Time_posix = make_datetime(
      year(Date), month(Date), day(Date),
      hour = hour_num, min = minute_num
    ),
    slot_start = floor_date(Time_posix, unit = "10 minutes"),
    time_of_day = format(slot_start, "%H:%M")
  ) %>%
  group_by(Survey.ID, Date, time_of_day) %>%
  summarise(
    n_species = n_distinct(Scientific_Name),   # changed here
    .groups = "drop"
  )

###############################################################################
# KNEPP: SHARED TIME SLOTS
###############################################################################

shared_slots_knepp <- intersect(
  unique(knepp_am_10min$time_of_day),
  unique(knepp_pc_10min$time_of_day)
)

knepp_am_shared <- knepp_am_10min %>%
  filter(time_of_day %in% shared_slots_knepp)

knepp_pc_shared <- knepp_pc_10min %>%
  filter(time_of_day %in% shared_slots_knepp)

knepp_combined_10min_data <- bind_rows(
  knepp_am_shared %>%
    rename(ID = Location) %>%
    mutate(method = "Passive acoustic"),
  knepp_pc_shared %>%
    rename(ID = Survey.ID) %>%
    mutate(method = "Point count")
)

knepp_combined_10min_data %>%
  dplyr::group_by(method) %>%
  dplyr::summarise(
    mean_species = mean(n_species, na.rm = TRUE),
    sd_species   = sd(n_species, na.rm = TRUE),
    se_species   = sd_species / sqrt(n()),
    ci_lower     = mean_species - 1.96 * se_species,
    ci_upper     = mean_species + 1.96 * se_species,
    n            = n(),
    .groups = "drop"
  )

knepp_combined_10min_data <- knepp_combined_10min_data %>%
  dplyr::mutate(
    method = factor(
      method,
      levels = c("Point count", "Passive acoustic"),
      labels = c("Point count", "PAM")
    )
  )

knepp_bird_10min_plot <- ggplot(knepp_combined_10min_data, aes(x = factor(method), y = n_species, fill = method)) +
  geom_boxplot(
    width = 0.3,
    alpha = 1,
    outlier.shape = NA ) +
  geom_jitter(
    aes(x = as.numeric(factor(method)) - 0.3),
    width = 0.1,
    height = 0,
    pch = 21,
    alpha = 0.4,
    size = 2) +
  coord_cartesian(ylim = c(0, 18)) +
  labs(
    x = "Method",
    y = "Species richness per 10 min") +
  theme_bw() +
  scale_colour_manual(
    values = c(
      "PAM" = "sienna",
      "Point count" = "cadetblue4"
    )
  ) +
  scale_fill_manual(
    values = c(
      "PAM" = "sienna",
      "Point count" = "cadetblue4"
    )
  ) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    axis.title.x = element_text(size = 14, colour = "black"),
    axis.title.y = element_text(size = 12, colour = "black"),
    axis.text.x = element_text(size = 12, colour = "black"),
    axis.text.y = element_text(size = 12, colour = "black"),
    text = element_text(size = 12)
  ) +
  annotate("text", x = 2, y = 14.5, label = "n = 4833", size = 3) +
  annotate("text", x = 2, y = 15.5, label = "Mean = 3.78", size = 3) +
  annotate("text", x = 1, y = 14.5, label = "n = 38", size = 3)+
  annotate("text", x = 1, y = 15.5, label = "Mean = 5.32", size = 3) 



###############################################################################
# combine plots

att_bird_10min_plot <- att_bird_10min_plot +
  ggplot2::ggtitle("Attenborough")

knepp_bird_10min_plot <- knepp_bird_10min_plot +
  ggplot2::ggtitle("Knepp")

bird_10min_plot_supp<-(att_bird_10min_plot / knepp_bird_10min_plot) +
  plot_annotation(
    tag_levels = "a",
    tag_prefix = "(",
    tag_suffix = ")"
  )



###############################################################################
# ATTENBOROUGH: POISSON GLM
###############################################################################

att_pois <- glm(
  n_species ~ method,
  data = combined_att_10min_data,
  family = poisson
)

summary(att_pois)
check_overdispersion(att_pois)


# PAM relative to point count
exp(coef(att_pois))

# 95% confidence intervals for richness ratios
exp(confint(att_pois))

###############################################################################
# KNEPP: POISSON GLM
###############################################################################

knepp_pois <- glm(
  n_species ~ method,
  data = knepp_combined_10min_data,
  family = poisson
)

summary(knepp_pois)
check_overdispersion(knepp_pois)


# PAM relative to point count
exp(coef(knepp_pois))

# 95% confidence intervals
exp(confint(knepp_pois))