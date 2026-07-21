# Aquatic invert analysis
# Bulk, eDNA & Morphological id @ Attenborough and Knepp
# J Boon

# lines 4 to 1379 Attenborough
# lines 1390 to 2501 Knepp

# in my comments, bulk metabarcoding = meta, edna = edna and morph = morphological identification

################################################################################
# Libraries 
################################################################################
library(ggplot2)
library(reshape2)
library(colorspace)
library(vegan)
library(dplyr)
library(tidyr)
library(stringr)
library(iNEXT)
library(forcats)
library(eulerr)
library(grid)
library(patchwork)
library(ggplotify)
library(cowplot)
library(permute)

################################################################################
# Attenborough data first
################################################################################
# From metabarcoding datasets from NatureMetrics there were some non-invert taxa
# REMOVE Anseriformes FROM META DATA - AS THEYRE BIRDS!
# REMOVE Caudata FROM META DATA - SALAMANDERS!
# REMOVE Cypriniformes and Gasterosteiformes FROM META DATA - FISH!


################################################################################
#Load Attenborough morphological ID data first
################################################################################
attenborough_aq_morph <- read.csv("Morph_ID_Att_Proc.csv", header = TRUE)

str(attenborough_aq_morph)

attenborough_aq_morph$Sampling.Point <- factor(attenborough_aq_morph$Sampling.Point)
attenborough_aq_morph$Species        <- factor(attenborough_aq_morph$Species)
attenborough_aq_morph$family         <- factor(attenborough_aq_morph$family)
attenborough_aq_morph$order        <- factor(attenborough_aq_morph$order)

length(unique(attenborough_aq_morph$Sampling.Point))  # 48 sampling points

attenborough_aq_morph %>%
  count(Sampling.Point, name = "num_records")


#################################################################################
# To standardise analysis keep those only detected to species level (morph)
################################################################################

attenborough_aq_morph_species_only <- attenborough_aq_morph %>%
  filter(!is.na(Species)) %>%
  mutate(Species_trim = str_squish(as.character(Species))) %>%
  filter(str_count(Species_trim, " ") == 1)

################################################################################
# Make species-level presence/absence matrix for morph (Attenborough) (morph)
################################################################################

attenborough_aq_pa <- attenborough_aq_morph_species_only %>%
  mutate(Presence_Absence = as.integer(value > 0)) %>%   # 1 = presence, 0 = absence
  select(Sampling.Point, Species_trim, Presence_Absence)

# Community matrix: Sampling.Point x Species
attenborough_aq_wide_pa <- attenborough_aq_pa %>%
  pivot_wider(
    names_from  = Species_trim,          
    values_from = Presence_Absence,       
    values_fill = 0,                      
    values_fn   = max                     
  ) %>%
  group_by(Sampling.Point) %>%            
  summarise(across(where(is.numeric), ~ as.integer(any(. > 0))), .groups = "drop")

# Convert to numeric matrix
attenborough_aq_species_matrix <- attenborough_aq_wide_pa %>%
  arrange(Sampling.Point)

rownames(attenborough_aq_species_matrix) <- attenborough_aq_species_matrix$Sampling.Point
attenborough_aq_species_matrix <- as.matrix(
  attenborough_aq_species_matrix %>% select(-Sampling.Point)
)


################################################################################
# Make family-level presence/absence matrix for morph (Attenborough) (morph)
################################################################################

# Drop NA and blank families
attenborough_aq_morph_family_only <- attenborough_aq_morph %>%
  filter(!is.na(family) & family != "") %>%
  mutate(family_trim = str_squish(as.character(family)))

num_family <- length(unique(attenborough_aq_morph_family_only$family_trim))
family_list <- unique(attenborough_aq_morph_family_only$family_trim)

cat("\nNumber of unique families:", num_family, "\n")
print(family_list)

# Presence absence by family
attenborough_aq_pa_family <- attenborough_aq_morph_family_only %>%
  mutate(Presence_Absence = as.integer(value > 0)) %>%
  select(Sampling.Point, family_trim, Presence_Absence)

# Matrix: Sampling.Point x Family
attenborough_aq_wide_pa_family <- attenborough_aq_pa_family %>%
  pivot_wider(
    names_from  = family_trim,            
    values_from = Presence_Absence,       
    values_fill = 0,
    values_fn   = max
  ) %>%
  group_by(Sampling.Point) %>%
  summarise(across(where(is.numeric), ~ as.integer(any(. > 0))), .groups = "drop")

# Numeric
attenborough_aq_family_matrix <- attenborough_aq_wide_pa_family %>%
  arrange(Sampling.Point)

rownames(attenborough_aq_family_matrix) <- attenborough_aq_family_matrix$Sampling.Point
attenborough_aq_family_matrix <- as.matrix(
  attenborough_aq_family_matrix %>% select(-Sampling.Point)
)


################################################################################
################################################################################
################################################################################


################################################################################
# Load and format eDNA data (Attenborough)                                                                
################################################################################

attenborough_aq_edna <- read.csv("eDNA_invert_Att_Proc.csv", header = TRUE)
str(attenborough_aq_edna)

attenborough_aq_edna$Sampling.Point <- factor(attenborough_aq_edna$Sampling.Point)
attenborough_aq_edna$Species        <- factor(attenborough_aq_edna$Species)
attenborough_aq_edna$family         <- factor(attenborough_aq_edna$family)
attenborough_aq_edna$order         <- factor(attenborough_aq_edna$order)

length(unique(attenborough_aq_edna$Sampling.Point)) # 45 sampling points

attenborough_aq_edna %>%
  count(Sampling.Point, name = "num_records")

################################################################################
# No. of species, genus and orders eDNA (Attenborough)
################################################################################
num_species <- length(unique(attenborough_aq_edna$Species))
species_list <- unique(attenborough_aq_edna$Species)

num_genus <- length(unique(attenborough_aq_edna$genus))
genus_list <- unique(attenborough_aq_edna$genus)

num_order <- length(unique(attenborough_aq_edna$order))
order_list <- unique(attenborough_aq_edna$order)

################################################################################
# Remove chordata - dataset from NatureMetrics contained chordates!
################################################################################

attenborough_aq_edna <- attenborough_aq_edna %>%
  filter(phylum != "Chordata")

attenborough_aq_edna_species_only <- attenborough_aq_edna %>%
  filter(!is.na(Species)) %>%
  mutate(Species_trim = str_squish(as.character(Species))) %>%
  filter(str_count(Species_trim, " ") == 1)

num_species <- length(unique(attenborough_aq_edna_species_only$Species))

################################################################################
# Presence/absence (eDNA) at species level     (Attenborough)                               
################################################################################
attenborough_aq_edna_pa <- attenborough_aq_edna_species_only %>%
  mutate(Presence_Absence = as.integer(value > 0)) %>%   # 1 = presence, 0 = absence
  select(Sampling.Point, Species_trim, Presence_Absence)

# Community matrix Sampling.Point x Species
attenborough_aq_edna_wide_pa <- attenborough_aq_edna_pa %>%
  pivot_wider(
    names_from  = Species_trim,         
    values_from = Presence_Absence,       
    values_fill = 0,                      
    values_fn   = max                     
  ) %>%
  group_by(Sampling.Point) %>%            
  summarise(across(where(is.numeric), ~ as.integer(any(. > 0))), .groups = "drop")

# To numeric 
attenborough_edna_aq_species_matrix <- attenborough_aq_edna_wide_pa %>%
  arrange(Sampling.Point)

rownames(attenborough_edna_aq_species_matrix) <- attenborough_edna_aq_species_matrix$Sampling.Point
attenborough_edna_aq_species_matrix <- as.matrix(
  attenborough_edna_aq_species_matrix %>% select(-Sampling.Point)
)


################################################################################
# Presence/absence (eDNA) at family level    (Attenborough)                                 #
################################################################################

# Drop NA and blank families
attenborough_aq_edna_family_only <- attenborough_aq_edna %>%
  filter(!is.na(family) & family != "") %>%
  mutate(family_trim = str_squish(as.character(family)))

num_family <- length(unique(attenborough_aq_edna_family_only$family_trim))
family_list <- unique(attenborough_aq_edna_family_only$family_trim)

cat("\nNumber of unique families (eDNA):", num_family, "\n")
print(family_list)

# Presence absence by family
attenborough_aq_edna_pa_family <- attenborough_aq_edna_family_only %>%
  mutate(Presence_Absence = as.integer(value > 0)) %>%
  select(Sampling.Point, family_trim, Presence_Absence)

# Community matrix Sampling.Point x Family
attenborough_aq_edna_wide_pa_family <- attenborough_aq_edna_pa_family %>%
  pivot_wider(
    names_from  = family_trim,            
    values_from = Presence_Absence,       
    values_fill = 0,
    values_fn   = max
  ) %>%
  group_by(Sampling.Point) %>%
  summarise(across(where(is.numeric), ~ as.integer(any(. > 0))), .groups = "drop")

# To numeric
attenborough_edna_aq_family_matrix <- attenborough_aq_edna_wide_pa_family %>%
  arrange(Sampling.Point)

rownames(attenborough_edna_aq_family_matrix) <- attenborough_edna_aq_family_matrix$Sampling.Point
attenborough_edna_aq_family_matrix <- as.matrix(
  attenborough_edna_aq_family_matrix %>% select(-Sampling.Point)
)

################################################################################
# Upload and format bulk metabarcoding data (Attenborough) - called this meta in the code                                              
################################################################################

attenborough_aq_meta <- read.csv("Meta_Att_Proc.csv", header = TRUE)
str(attenborough_aq_meta)

attenborough_aq_meta$Sampling.Point <- factor(attenborough_aq_meta$Sampling.Point)
attenborough_aq_meta$Species        <- factor(attenborough_aq_meta$Species)
attenborough_aq_meta$family         <- factor(attenborough_aq_meta$family)
attenborough_aq_meta$order          <- factor(attenborough_aq_meta$order)

################################################################################
# Remove non-invert orders (birds, amphibians, fish)
################################################################################

bad_orders <- c("Anseriformes", "Caudata", "Cypriniformes", "Gasterosteiformes")

attenborough_aq_meta <- attenborough_aq_meta %>%
  mutate(order_trim = str_squish(as.character(order))) %>%
  filter(!order_trim %in% bad_orders) %>%
  select(-order_trim)

str(attenborough_aq_meta)

# check
attenborough_aq_meta %>%
  filter(str_squish(as.character(order)) %in% bad_orders) %>%
  nrow()

num_order_meta <- length(unique(attenborough_aq_meta$order))
order_list_meta <- unique(attenborough_aq_meta$order)

################################################################################
#Trim to species
################################################################################

attenborough_aq_meta_species_only <- attenborough_aq_meta %>%
  filter(!is.na(Species)) %>%
  mutate(Species_trim = str_squish(as.character(Species))) %>%
  filter(str_count(Species_trim, " ") == 1)

################################################################################
# Species presence/absence matrix for bulk metabarcoding (Attenborough)                     
################################################################################
attenborough_aq_meta_pa <- attenborough_aq_meta_species_only %>%
  mutate(Presence_Absence = as.integer(value > 0)) %>% 
  select(Sampling.Point, Species_trim, Presence_Absence)

# Community matrix: Sampling.Point x Species
attenborough_aq_meta_wide_pa <- attenborough_aq_meta_pa %>%
  pivot_wider(
    names_from  = Species_trim,           
    values_from = Presence_Absence,       
    values_fill = 0,                      
    values_fn   = max                     
  ) %>%
  group_by(Sampling.Point) %>%
  summarise(across(where(is.numeric), ~ as.integer(any(. > 0))), .groups = "drop")

# To numeric matrix
attenborough_meta_aq_species_matrix <- attenborough_aq_meta_wide_pa %>%
  arrange(Sampling.Point)

rownames(attenborough_meta_aq_species_matrix) <- attenborough_meta_aq_species_matrix$Sampling.Point
attenborough_meta_aq_species_matrix <- as.matrix(
  attenborough_meta_aq_species_matrix %>% select(-Sampling.Point)
)


################################################################################
# Family-level presence/absence matrix for bulk metabarcoding (Attenborough)                     
################################################################################

# Drop NA and not id'ed to families, trim names
attenborough_aq_meta_family_only <- attenborough_aq_meta %>%
  filter(!is.na(family) & family != "") %>%
  mutate(family_trim = str_squish(as.character(family)))

num_family_meta  <- length(unique(attenborough_aq_meta_family_only$family_trim))
family_list_meta <- unique(attenborough_aq_meta_family_only$family_trim)

cat("\nNumber of unique families (metabarcoding):", num_family_meta, "\n")
print(family_list_meta)

# Presence absence by family
attenborough_aq_meta_pa_family <- attenborough_aq_meta_family_only %>%
  mutate(Presence_Absence = as.integer(value > 0)) %>%
  select(Sampling.Point, family_trim, Presence_Absence)

# Community matrix: Sampling.Point x Family
attenborough_aq_meta_wide_pa_family <- attenborough_aq_meta_pa_family %>%
  pivot_wider(
    names_from  = family_trim,            
    values_from = Presence_Absence,       
    values_fill = 0,
    values_fn   = max
  ) %>%
  group_by(Sampling.Point) %>%
  summarise(across(where(is.numeric), ~ as.integer(any(. > 0))), .groups = "drop")

# Make numeric matrix
attenborough_meta_aq_family_matrix <- attenborough_aq_meta_wide_pa_family %>%
  arrange(Sampling.Point)

rownames(attenborough_meta_aq_family_matrix) <- attenborough_meta_aq_family_matrix$Sampling.Point
attenborough_meta_aq_family_matrix <- as.matrix(
  attenborough_meta_aq_family_matrix %>% select(-Sampling.Point)
)


################################################################################
# Species accumulation curves - inext  (Attenborough)
################################################################################

# Need to convert to species x point matrix to where rows=species and  cols=samples
mat_morph <- t(attenborough_aq_species_matrix)
mat_edna  <- t(attenborough_edna_aq_species_matrix)
mat_bulk  <- t(attenborough_meta_aq_species_matrix)

# 0/1 integers
for (m in list(mat_morph, mat_edna, mat_bulk)) {
  storage.mode(m) <- "integer"
}

mat_morph[mat_morph > 1] <- 1
mat_edna[mat_edna > 1]   <- 1
mat_bulk[mat_bulk > 1]   <- 1

# find shared endpoint (extrapolate to 2× the max sampling effort for all methods)
n_morph <- ncol(mat_morph)
n_edna  <- ncol(mat_edna)
n_bulk  <- ncol(mat_bulk)

endpoint_max <- 2 * max(n_morph, n_edna, n_bulk)

# iNEXT 

inext_three <- iNEXT(
  list(Morph = mat_morph, eDNA = mat_edna, Bulk = mat_bulk),
  q        = 0,
  datatype = "incidence_raw",
  endpoint = endpoint_max,
  knots    = 80,
  se       = TRUE,
  conf     = 0.95
)

################################################################################
# PLOT Species accumulation curves - inext  (Attenborough)
################################################################################

my_cols <- c(
  "Morph" = "#7A8F3A",
  "eDNA"  = "#4C8C8C",
  "Bulk"  = "#8A6F8F"
)

# legend order
leg_breaks <- c("Morph", "eDNA", "Bulk")

att_aquatic_accum_plot <- ggiNEXT(inext_three, type = 1, color.var = "Assemblage") +
  scale_colour_manual(
    values = my_cols,
    breaks = leg_breaks
  ) +
  scale_fill_manual(values = my_cols) +
  guides(
    fill     = "none",   
    linetype = "none",  
    shape    = "none",
    colour   = guide_legend(
      override.aes = list(
        linetype = 1,
        shape    = NA,
        # ribbon + line in same legend key
        fill     = unname(my_cols[leg_breaks]),
        alpha    = 0.25
      )
    )
  ) +
  labs(
    x = "Number of sampling events",
    y = "Species richness",
    colour = NULL,
    fill   = NULL
  ) +
  scale_y_continuous(
    breaks = c(50, 100, 150, 200, 250),
    limits = c(0, 250)
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.major     = element_blank(),
    panel.grid.minor     = element_blank(),
    axis.title           = element_text(size = 12),
    axis.text            = element_text(size = 12),
    legend.title         = element_blank(),
    legend.position      = c(0.02, 0.97),   # top left inside
    legend.justification = c("left", "top"),
    legend.background    = element_rect(fill = NA, colour = NA)
  )


################################################################################
# The actual Species accumilation curves - inext  (Attenborough) - for the paper
################################################################################

att_aquatic_accum_plot



################################################################################
# Family accumulation curves - inext  (Attenborough)
################################################################################

# family x sites
mat_morph_fam <- t(attenborough_aq_family_matrix)
mat_edna_fam  <- t(attenborough_edna_aq_family_matrix)
mat_bulk_fam  <- t(attenborough_meta_aq_family_matrix)

# integers
mat_morph_fam[mat_morph_fam > 1] <- 1
mat_edna_fam[mat_edna_fam > 1]   <- 1
mat_bulk_fam[mat_bulk_fam > 1]   <- 1

storage.mode(mat_morph_fam) <- "integer"
storage.mode(mat_edna_fam)  <- "integer"
storage.mode(mat_bulk_fam)  <- "integer"

# shared endpoint
n_morph_fam <- ncol(mat_morph_fam)
n_edna_fam  <- ncol(mat_edna_fam)
n_bulk_fam  <- ncol(mat_bulk_fam)

endpoint_max_fam <- 2 * max(n_morph_fam, n_edna_fam, n_bulk_fam)

#iNEXT  families
inext_three_fam <- iNEXT(
  list(Morph = mat_morph_fam, eDNA = mat_edna_fam, Bulk = mat_bulk_fam),
  q        = 0,
  datatype = "incidence_raw",
  endpoint = endpoint_max_fam,
  knots    = 80,
  se       = TRUE,
  conf     = 0.95
)

################################################################################
# PLOT family accumulation curves - inext  (Attenborough)
################################################################################

# Family level plot
my_cols <- c(
  "Morph" = "#7A8F3A",
  "eDNA"  = "#4C8C8C",
  "Bulk"  = "#8A6F8F"
)


leg_breaks <- c("Morph", "eDNA", "Bulk")

att_aquatic_accum_plot_family <- ggiNEXT(inext_three_fam, type = 1, color.var = "Assemblage") +
  scale_colour_manual(
    values = my_cols,
    breaks = leg_breaks
  ) +
  scale_fill_manual(values = my_cols) +
  guides(
    fill     = "none",   
    linetype = "none",   
    shape    = "none",
    colour   = guide_legend(
      override.aes = list(
        linetype = 1,
        shape    = NA,
        fill     = unname(my_cols[leg_breaks]),
        alpha    = 0.25
      )
    )
  ) +
  labs(
    x = "Number of sampling events",
    y = "Family richness",
    colour = NULL,
    fill   = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.major     = element_blank(),
    panel.grid.minor     = element_blank(),
    axis.title           = element_text(size = 12),
    axis.text            = element_text(size = 12),
    legend.title         = element_blank(),
    legend.position      = c(0.02, 0.97),  
    legend.justification = c("left", "top"),
    legend.background    = element_rect(fill = NA, colour = NA)
  )

################################################################################
# The actual family accumilation curves - inext  (Attenborough) - for the paper supp
################################################################################

att_aquatic_accum_plot_family




################################################################################
################################################################################
################################################################################



################################################################################
# NMDS species dissim and plot (Jaccard) (Attenborough)
################################################################################

# Remove ID columns and keep species matrices
mic_pa  <- attenborough_aq_wide_pa
edna_pa <- attenborough_aq_edna_wide_pa
bulk_pa <- attenborough_aq_meta_wide_pa

mic_species  <- mic_pa  %>% select(-Sampling.Point)
edna_species <- edna_pa %>% select(-Sampling.Point)
bulk_species <- bulk_pa %>% select(-Sampling.Point)

# Match species columns
all_species12 <- union(names(mic_species), names(edna_species))
all_species   <- union(all_species12, names(bulk_species))

mic_species2  <- add_missing_cols(mic_species,  all_species)
edna_species2 <- add_missing_cols(edna_species, all_species)
bulk_species2 <- add_missing_cols(bulk_species, all_species)

# Combine with Source + IDs
mic_comb <- mic_species2 %>%
  mutate(
    Source = "Morph",
    ID     = as.character(mic_pa$Sampling.Point)
  )

edna_comb <- edna_species2 %>%
  mutate(
    Source = "eDNA",
    ID     = as.character(edna_pa$Sampling.Point)
  )

bulk_comb <- bulk_species2 %>%
  mutate(
    Source = "Bulk",
    ID     = as.character(bulk_pa$Sampling.Point)
  )

attenborough_combined <- bind_rows(mic_comb, edna_comb, bulk_comb)

# Remove empty samples and add Site as the nesting/blocking variable
attenborough_combined_filt <- attenborough_combined %>%
  filter(rowSums(select(., -Source, -ID)) > 0) %>%
  mutate(
    Site   = stringr::str_extract(ID, "att\\d+"),
    Source = factor(Source, levels = c("Morph", "eDNA", "Bulk"))
  )

# Check nesting structure
table(attenborough_combined_filt$Site, attenborough_combined_filt$Source)
sum(is.na(attenborough_combined_filt$Site))

################################################################################
# Run NMDS - SPECIES
################################################################################

set.seed(123)

species_mat <- attenborough_combined_filt %>%
  select(-Source, -ID, -Site)

attenborough_nmds <- metaMDS(
  species_mat,
  distance      = "jaccard",
  binary        = TRUE,
  k             = 2,
  trymax        = 400,
  autotransform = FALSE
)

attenborough_scores <- as.data.frame(scores(attenborough_nmds, display = "sites"))
attenborough_scores$Source <- attenborough_combined_filt$Source
attenborough_scores$ID     <- attenborough_combined_filt$ID
attenborough_scores$Site   <- attenborough_combined_filt$Site

################################################################################
# Plot NMDS - SPECIES
################################################################################

stress_val <- round(attenborough_nmds$stress, 3)

attenborough_aquatic_nmds_plot <- ggplot(
  attenborough_scores,
  aes(x = NMDS1, y = NMDS2, colour = Source)
) +
  geom_point(size = 1.5) +
  stat_ellipse(
    aes(fill = Source, colour = Source, group = Source),
    geom = "polygon",
    alpha = 0.15,
    linewidth = 0.8,
    type = "t",
    level = 0.95,
    show.legend = FALSE
  ) +
  scale_colour_manual(values = c(
    "Morph" = "#7A8F3A",
    "eDNA"  = "#4C8C8C",
    "Bulk"  = "#8A6F8F"
  )) +
  scale_fill_manual(values = c(
    "Morph" = "#7A8F3A",
    "eDNA"  = "#4C8C8C",
    "Bulk"  = "#8A6F8F"
  )) +
  labs(x = "NMDS1", y = "NMDS2", colour = NULL, fill = NULL) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    legend.position = c(0.999, 0.001),
    legend.justification = c("right", "bottom"),
    legend.background = element_rect(fill = NA, colour = NA),
    legend.key = element_rect(fill = NA, colour = NA),
    legend.key.width = unit(0.1, "cm"),
    legend.spacing.y = unit(0.01, "cm"),
    axis.title = element_text(size = 12, colour = "black"),
    axis.text = element_text(size = 12, colour = "black"),
    legend.text = element_text(size = 11, colour = "black"),
    plot.background = element_rect(fill = NA, colour = NA),
    panel.background = element_rect(fill = NA, colour = NA)
  ) +
  guides(colour = guide_legend(override.aes = list(size = 3, shape = 15))) +
  annotate(
    "text",
    x = -Inf, y = Inf,
    label = paste0("Stress = ", stress_val),
    hjust = -0.1,
    vjust = 1.2,
    size = 3.2,
    colour = "black"
  )

attenborough_aquatic_nmds_plot

################################################################################
# PERMANOVA - SPECIES, blocked by Site (Attenborough)
################################################################################

dist_mat <- vegdist(species_mat, method = "jaccard", binary = TRUE)

env <- attenborough_combined_filt %>%
  select(Source, ID, Site) %>%
  mutate(
    Source = factor(Source, levels = c("Morph", "eDNA", "Bulk"))
  )

perm_species <- how(nperm = 999)
setBlocks(perm_species) <- env$Site

set.seed(123)

adon <- adonis2(
  dist_mat ~ Source,
  data = env,
  permutations = perm_species
)

print(adon)

################################################################################
# PERMDISP-style test - SPECIES
# Site-blocked permutation test on distance-to-centroid values
################################################################################

bd <- betadisper(dist_mat, env$Source)

disp_df_att_aq_sp <- data.frame(
  Distance = bd$distances,
  Source   = env$Source,
  Site     = env$Site
)

# Response must be supplied as a matrix for adonis2
disp_mat_sp <- as.matrix(disp_df_att_aq_sp$Distance)
rownames(disp_mat_sp) <- rownames(disp_df_att_aq_sp)

perm_disp_species <- how(nperm = 999)
setBlocks(perm_disp_species) <- disp_df_att_aq_sp$Site

set.seed(123)

disp_test_species <- adonis2(
  disp_mat_sp ~ Source,
  data = disp_df_att_aq_sp,
  method = "euclidean",
  permutations = perm_disp_species
)

print(disp_test_species)

################################################################################
# Pairwise PERMDISP-style tests - SPECIES, blocked by Site
################################################################################

# Morph vs eDNA
disp_sp_me <- disp_df_att_aq_sp %>%
  filter(Source %in% c("Morph", "eDNA")) %>%
  droplevels()

disp_mat_sp_me <- as.matrix(disp_sp_me$Distance)
rownames(disp_mat_sp_me) <- rownames(disp_sp_me)

perm_disp_sp_me <- how(nperm = 999)
setBlocks(perm_disp_sp_me) <- disp_sp_me$Site

set.seed(123)

disp_test_sp_me <- adonis2(
  disp_mat_sp_me ~ Source,
  data = disp_sp_me,
  method = "euclidean",
  permutations = perm_disp_sp_me
)

# Morph vs Bulk
disp_sp_mb <- disp_df_att_aq_sp %>%
  filter(Source %in% c("Morph", "Bulk")) %>%
  droplevels()

disp_mat_sp_mb <- as.matrix(disp_sp_mb$Distance)
rownames(disp_mat_sp_mb) <- rownames(disp_sp_mb)

perm_disp_sp_mb <- how(nperm = 999)
setBlocks(perm_disp_sp_mb) <- disp_sp_mb$Site

set.seed(123)

disp_test_sp_mb <- adonis2(
  disp_mat_sp_mb ~ Source,
  data = disp_sp_mb,
  method = "euclidean",
  permutations = perm_disp_sp_mb
)

# eDNA vs Bulk
disp_sp_eb <- disp_df_att_aq_sp %>%
  filter(Source %in% c("eDNA", "Bulk")) %>%
  droplevels()

disp_mat_sp_eb <- as.matrix(disp_sp_eb$Distance)
rownames(disp_mat_sp_eb) <- rownames(disp_sp_eb)

perm_disp_sp_eb <- how(nperm = 999)
setBlocks(perm_disp_sp_eb) <- disp_sp_eb$Site

set.seed(123)

disp_test_sp_eb <- adonis2(
  disp_mat_sp_eb ~ Source,
  data = disp_sp_eb,
  method = "euclidean",
  permutations = perm_disp_sp_eb
)

disp_pw_sp <- tibble(
  comparison = c("Morph vs eDNA", "Morph vs Bulk", "eDNA vs Bulk"),
  F = c(
    disp_test_sp_me$F[1],
    disp_test_sp_mb$F[1],
    disp_test_sp_eb$F[1]
  ),
  p = c(
    disp_test_sp_me$`Pr(>F)`[1],
    disp_test_sp_mb$`Pr(>F)`[1],
    disp_test_sp_eb$`Pr(>F)`[1]
  )
) %>%
  mutate(p_adj_holm = p.adjust(p, method = "holm")) %>%
  arrange(p_adj_holm)

print(disp_pw_sp)

################################################################################
# Summary and plot of dispersion distances - SPECIES
################################################################################

dist_summary_att_aq_sp <- disp_df_att_aq_sp %>%
  dplyr::group_by(Source) %>%
  dplyr::summarise(
    mean_distance = mean(Distance, na.rm = TRUE),
    sd_distance   = sd(Distance, na.rm = TRUE),
    se_distance   = sd_distance / sqrt(n()),
    ci_lower      = mean_distance - 1.96 * se_distance,
    ci_upper      = mean_distance + 1.96 * se_distance,
    n             = n(),
    .groups = "drop"
  )

print(dist_summary_att_aq_sp)

distance_plot_att_aq_sp <- ggplot(
  disp_df_att_aq_sp,
  aes(x = Source, y = Distance, fill = Source)
) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.1, alpha = 0.6) +
  scale_colour_manual(values = c(
    "Morph" = "#7A8F3A",
    "eDNA"  = "#4C8C8C",
    "Bulk"  = "#8A6F8F"
  )) +
  scale_fill_manual(values = c(
    "Morph" = "#7A8F3A",
    "eDNA"  = "#4C8C8C",
    "Bulk"  = "#8A6F8F"
  )) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    title = "Attenborough (Species)",
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
  )

distance_plot_att_aq_sp

################################################################################
# Pairwise PERMANOVA - SPECIES, blocked by Site
################################################################################

dat <- attenborough_combined_filt %>%
  select(Source, ID, Site, everything()) %>%
  mutate(
    Source = factor(Source, levels = c("Morph", "eDNA", "Bulk"))
  )

sp_cols <- setdiff(names(dat), c("Source", "ID", "Site"))

# Morph vs eDNA
dat_me <- dat %>%
  filter(Source %in% c("Morph", "eDNA"))

dist_me <- vegdist(
  dat_me %>% select(all_of(sp_cols)),
  method = "jaccard",
  binary = TRUE
)

perm_me <- how(nperm = 999)
setBlocks(perm_me) <- dat_me$Site

set.seed(123)

adon_me <- adonis2(
  dist_me ~ Source,
  data = dat_me,
  permutations = perm_me
)

# Morph vs Bulk
dat_mb <- dat %>%
  filter(Source %in% c("Morph", "Bulk"))

dist_mb <- vegdist(
  dat_mb %>% select(all_of(sp_cols)),
  method = "jaccard",
  binary = TRUE
)

perm_mb <- how(nperm = 999)
setBlocks(perm_mb) <- dat_mb$Site

set.seed(123)

adon_mb <- adonis2(
  dist_mb ~ Source,
  data = dat_mb,
  permutations = perm_mb
)

# eDNA vs Bulk
dat_eb <- dat %>%
  filter(Source %in% c("eDNA", "Bulk"))

dist_eb <- vegdist(
  dat_eb %>% select(all_of(sp_cols)),
  method = "jaccard",
  binary = TRUE
)

perm_eb <- how(nperm = 999)
setBlocks(perm_eb) <- dat_eb$Site

set.seed(123)

adon_eb <- adonis2(
  dist_eb ~ Source,
  data = dat_eb,
  permutations = perm_eb
)

################################################################################
# Combine pairwise PERMANOVA results - SPECIES
################################################################################

pw <- tibble(
  comparison = c("Morph vs eDNA", "Morph vs Bulk", "eDNA vs Bulk"),
  R2 = c(adon_me$R2[1], adon_mb$R2[1], adon_eb$R2[1]),
  F  = c(adon_me$F[1],  adon_mb$F[1],  adon_eb$F[1]),
  p  = c(
    adon_me$`Pr(>F)`[1],
    adon_mb$`Pr(>F)`[1],
    adon_eb$`Pr(>F)`[1]
  )
) %>%
  mutate(p_adj_holm = p.adjust(p, method = "holm")) %>%
  arrange(p_adj_holm)

print(pw)

################################################################################
# FAMILY-level dissim and NMDS (presence/absence Jaccard) (Attenborough)
################################################################################

mic_pa_fam  <- attenborough_aq_wide_pa_family
edna_pa_fam <- attenborough_aq_edna_wide_pa_family
bulk_pa_fam <- attenborough_aq_meta_wide_pa_family

mic_family  <- mic_pa_fam  %>% select(-Sampling.Point)
edna_family <- edna_pa_fam %>% select(-Sampling.Point)
bulk_family <- bulk_pa_fam %>% select(-Sampling.Point)

# Match family columns
all_families12 <- union(names(mic_family), names(edna_family))
all_families   <- union(all_families12, names(bulk_family))

mic_family2  <- add_missing_cols(mic_family,  all_families)
edna_family2 <- add_missing_cols(edna_family, all_families)
bulk_family2 <- add_missing_cols(bulk_family, all_families)

# Combine with Source + IDs
mic_comb_fam <- mic_family2 %>%
  mutate(
    Source = "Morph",
    ID     = as.character(mic_pa_fam$Sampling.Point)
  )

edna_comb_fam <- edna_family2 %>%
  mutate(
    Source = "eDNA",
    ID     = as.character(edna_pa_fam$Sampling.Point)
  )

bulk_comb_fam <- bulk_family2 %>%
  mutate(
    Source = "Bulk",
    ID     = as.character(bulk_pa_fam$Sampling.Point)
  )

attenborough_combined_fam <- bind_rows(
  mic_comb_fam,
  edna_comb_fam,
  bulk_comb_fam
)

# Remove empty samples and add Site as the blocking variable
attenborough_combined_fam_filt <- attenborough_combined_fam %>%
  filter(rowSums(select(., -Source, -ID)) > 0) %>%
  mutate(
    Site   = stringr::str_extract(ID, "att\\d+"),
    Source = factor(Source, levels = c("Morph", "eDNA", "Bulk"))
  )

# Check
table(attenborough_combined_fam_filt$Site, attenborough_combined_fam_filt$Source)
sum(is.na(attenborough_combined_fam_filt$Site))





################################################################################
# Run NMDS - FAMILY (Attenborough)
################################################################################

set.seed(123)

family_mat <- attenborough_combined_fam_filt %>%
  select(-Source, -ID, -Site)

attenborough_nmds_fam <- metaMDS(
  family_mat,
  distance      = "jaccard",
  binary        = TRUE,
  k             = 2,
  trymax        = 400,
  autotransform = FALSE
)

attenborough_scores_fam <- as.data.frame(scores(attenborough_nmds_fam, display = "sites"))
attenborough_scores_fam$Source <- attenborough_combined_fam_filt$Source
attenborough_scores_fam$ID     <- attenborough_combined_fam_filt$ID
attenborough_scores_fam$Site   <- attenborough_combined_fam_filt$Site

################################################################################
# Plot NMDS - FAMILY
################################################################################

stress_val_fam <- round(attenborough_nmds_fam$stress, 3)

attenborough_aquatic_nmds_plot_family <- ggplot(
  attenborough_scores_fam,
  aes(x = NMDS1, y = NMDS2, colour = Source)
) +
  geom_point(size = 1.5) +
  stat_ellipse(
    aes(fill = Source, colour = Source, group = Source),
    geom = "polygon",
    alpha = 0.15,
    linewidth = 0.8,
    type = "t",
    level = 0.95,
    show.legend = FALSE
  ) +
  scale_colour_manual(values = c(
    "Morph" = "#7A8F3A",
    "eDNA"  = "#4C8C8C",
    "Bulk"  = "#8A6F8F"
  )) +
  scale_fill_manual(values = c(
    "Morph" = "#7A8F3A",
    "eDNA"  = "#4C8C8C",
    "Bulk"  = "#8A6F8F"
  )) +
  labs(
    x = "NMDS1",
    y = "NMDS2",
    colour = NULL,
    fill = NULL
  ) +
  coord_cartesian(xlim = c(-2, NA), ylim = c(NA, 3)) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    legend.position = c(0.01, 0.99),
    legend.justification = c("left", "top"),
    legend.background = element_rect(fill = NA, colour = NA),
    legend.key = element_rect(fill = NA, colour = NA),
    legend.key.width = unit(0.1, "cm"),
    legend.spacing.y = unit(0.01, "cm"),
    axis.title = element_text(size = 12, colour = "black"),
    axis.text  = element_text(size = 12, colour = "black"),
    legend.text = element_text(size = 10, colour = "black"),
    plot.background = element_rect(fill = NA, colour = NA),
    panel.background = element_rect(fill = NA, colour = NA)
  ) +
  guides(colour = guide_legend(override.aes = list(size = 3, shape = 15))) +
  annotate(
    "text",
    x = -Inf, y = -Inf,
    label = paste0("Stress = ", stress_val_fam),
    hjust = -0.1,
    vjust = -1,
    size = 3.2,
    colour = "black"
  )

attenborough_aquatic_nmds_plot_family

################################################################################
# PERMANOVA - FAMILY, blocked by Site
################################################################################

dist_mat_fam <- vegdist(family_mat, method = "jaccard", binary = TRUE)

env_fam <- attenborough_combined_fam_filt %>%
  select(Source, ID, Site) %>%
  mutate(
    Source = factor(Source, levels = c("Morph", "eDNA", "Bulk"))
  )

perm_fam <- how(nperm = 999)
setBlocks(perm_fam) <- env_fam$Site

set.seed(123)

adon_fam <- adonis2(
  dist_mat_fam ~ Source,
  data = env_fam,
  permutations = perm_fam
)

print(adon_fam)

################################################################################
# PERMDISP-style test - FAMILY
# Site-blocked permutation test on distance-to-centroid values
################################################################################

bd_fam <- betadisper(dist_mat_fam, env_fam$Source)

disp_df_att_aq_fam <- data.frame(
  Distance = bd_fam$distances,
  Source   = env_fam$Source,
  Site     = env_fam$Site
)

# Response must be supplied as a matrix for adonis2
disp_mat_fam <- as.matrix(disp_df_att_aq_fam$Distance)
rownames(disp_mat_fam) <- rownames(disp_df_att_aq_fam)

perm_disp_fam <- how(nperm = 999)
setBlocks(perm_disp_fam) <- disp_df_att_aq_fam$Site

set.seed(123)

disp_test_fam <- adonis2(
  disp_mat_fam ~ Source,
  data = disp_df_att_aq_fam,
  method = "euclidean",
  permutations = perm_disp_fam
)

print(disp_test_fam)

################################################################################
# Pairwise PERMDISP-style tests - FAMILY, blocked by Site
################################################################################

# Morph vs eDNA
disp_fam_me <- disp_df_att_aq_fam %>%
  filter(Source %in% c("Morph", "eDNA")) %>%
  droplevels()

disp_mat_fam_me <- as.matrix(disp_fam_me$Distance)
rownames(disp_mat_fam_me) <- rownames(disp_fam_me)

perm_disp_fam_me <- how(nperm = 999)
setBlocks(perm_disp_fam_me) <- disp_fam_me$Site

set.seed(123)

disp_test_fam_me <- adonis2(
  disp_mat_fam_me ~ Source,
  data = disp_fam_me,
  method = "euclidean",
  permutations = perm_disp_fam_me
)

# Morph vs Bulk
disp_fam_mb <- disp_df_att_aq_fam %>%
  filter(Source %in% c("Morph", "Bulk")) %>%
  droplevels()

disp_mat_fam_mb <- as.matrix(disp_fam_mb$Distance)
rownames(disp_mat_fam_mb) <- rownames(disp_fam_mb)

perm_disp_fam_mb <- how(nperm = 999)
setBlocks(perm_disp_fam_mb) <- disp_fam_mb$Site

set.seed(123)

disp_test_fam_mb <- adonis2(
  disp_mat_fam_mb ~ Source,
  data = disp_fam_mb,
  method = "euclidean",
  permutations = perm_disp_fam_mb
)

# eDNA vs Bulk
disp_fam_eb <- disp_df_att_aq_fam %>%
  filter(Source %in% c("eDNA", "Bulk")) %>%
  droplevels()

disp_mat_fam_eb <- as.matrix(disp_fam_eb$Distance)
rownames(disp_mat_fam_eb) <- rownames(disp_fam_eb)

perm_disp_fam_eb <- how(nperm = 999)
setBlocks(perm_disp_fam_eb) <- disp_fam_eb$Site

set.seed(123)

disp_test_fam_eb <- adonis2(
  disp_mat_fam_eb ~ Source,
  data = disp_fam_eb,
  method = "euclidean",
  permutations = perm_disp_fam_eb
)

disp_pw_fam <- tibble(
  comparison = c("Morph vs eDNA", "Morph vs Bulk", "eDNA vs Bulk"),
  F = c(
    disp_test_fam_me$F[1],
    disp_test_fam_mb$F[1],
    disp_test_fam_eb$F[1]
  ),
  p = c(
    disp_test_fam_me$`Pr(>F)`[1],
    disp_test_fam_mb$`Pr(>F)`[1],
    disp_test_fam_eb$`Pr(>F)`[1]
  )
) %>%
  mutate(p_adj_holm = p.adjust(p, method = "holm")) %>%
  arrange(p_adj_holm)

print(disp_pw_fam)

################################################################################
# Summary and plot of dispersion distances - FAMILY
################################################################################

dist_summary_att_aq_fam <- disp_df_att_aq_fam %>%
  dplyr::group_by(Source) %>%
  dplyr::summarise(
    mean_distance = mean(Distance, na.rm = TRUE),
    sd_distance   = sd(Distance, na.rm = TRUE),
    se_distance   = sd_distance / sqrt(n()),
    ci_lower      = mean_distance - 1.96 * se_distance,
    ci_upper      = mean_distance + 1.96 * se_distance,
    n             = n(),
    .groups = "drop"
  )

print(dist_summary_att_aq_fam)

distance_plot_att_aq_fam <- ggplot(
  disp_df_att_aq_fam,
  aes(x = Source, y = Distance, fill = Source)
) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.1, alpha = 0.6) +
  scale_colour_manual(values = c(
    "Morph" = "#7A8F3A",
    "eDNA"  = "#4C8C8C",
    "Bulk"  = "#8A6F8F"
  )) +
  scale_fill_manual(values = c(
    "Morph" = "#7A8F3A",
    "eDNA"  = "#4C8C8C",
    "Bulk"  = "#8A6F8F"
  )) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    title = "Attenborough (Family)",
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
  )

distance_plot_att_aq_fam

################################################################################
# Pairwise PERMANOVA - FAMILY, blocked by Site
################################################################################

dat_fam <- attenborough_combined_fam_filt %>%
  select(Source, ID, Site, everything()) %>%
  mutate(
    Source = factor(Source, levels = c("Morph", "eDNA", "Bulk"))
  )

fam_cols <- setdiff(names(dat_fam), c("Source", "ID", "Site"))

# Morph vs eDNA
dat_fam_me <- dat_fam %>%
  filter(Source %in% c("Morph", "eDNA"))

dist_fam_me <- vegdist(
  dat_fam_me %>% select(all_of(fam_cols)),
  method = "jaccard",
  binary = TRUE
)

perm_fam_me <- how(nperm = 999)
setBlocks(perm_fam_me) <- dat_fam_me$Site

set.seed(123)

adon_fam_me <- adonis2(
  dist_fam_me ~ Source,
  data = dat_fam_me,
  permutations = perm_fam_me
)

# Morph vs Bulk
dat_fam_mb <- dat_fam %>%
  filter(Source %in% c("Morph", "Bulk"))

dist_fam_mb <- vegdist(
  dat_fam_mb %>% select(all_of(fam_cols)),
  method = "jaccard",
  binary = TRUE
)

perm_fam_mb <- how(nperm = 999)
setBlocks(perm_fam_mb) <- dat_fam_mb$Site

set.seed(123)

adon_fam_mb <- adonis2(
  dist_fam_mb ~ Source,
  data = dat_fam_mb,
  permutations = perm_fam_mb
)

# eDNA vs Bulk
dat_fam_eb <- dat_fam %>%
  filter(Source %in% c("eDNA", "Bulk"))

dist_fam_eb <- vegdist(
  dat_fam_eb %>% select(all_of(fam_cols)),
  method = "jaccard",
  binary = TRUE
)

perm_fam_eb <- how(nperm = 999)
setBlocks(perm_fam_eb) <- dat_fam_eb$Site

set.seed(123)

adon_fam_eb <- adonis2(
  dist_fam_eb ~ Source,
  data = dat_fam_eb,
  permutations = perm_fam_eb
)

################################################################################
# Combine pairwise PERMANOVA results - FAMILY
################################################################################

pw_fam <- tibble(
  comparison = c("Morph vs eDNA", "Morph vs Bulk", "eDNA vs Bulk"),
  R2 = c(adon_fam_me$R2[1], adon_fam_mb$R2[1], adon_fam_eb$R2[1]),
  F  = c(adon_fam_me$F[1],  adon_fam_mb$F[1],  adon_fam_eb$F[1]),
  p  = c(
    adon_fam_me$`Pr(>F)`[1],
    adon_fam_mb$`Pr(>F)`[1],
    adon_fam_eb$`Pr(>F)`[1]
  )
) %>%
  mutate(p_adj_holm = p.adjust(p, method = "holm")) %>%
  arrange(p_adj_holm)

print(pw_fam)

################################################################################
# Result tables for reporting
################################################################################

# Overall site-blocked tests
overall_tests_att_aq <- tibble(
  Level = c("Species", "Species", "Family", "Family"),
  Test = c("PERMANOVA", "Dispersion", "PERMANOVA", "Dispersion"),
  Term = c("Source", "Source", "Source", "Source"),
  Df = c(adon$Df[1], disp_test_species$Df[1], adon_fam$Df[1], disp_test_fam$Df[1]),
  F = c(adon$F[1], disp_test_species$F[1], adon_fam$F[1], disp_test_fam$F[1]),
  R2 = c(adon$R2[1], disp_test_species$R2[1], adon_fam$R2[1], disp_test_fam$R2[1]),
  p = c(
    adon$`Pr(>F)`[1],
    disp_test_species$`Pr(>F)`[1],
    adon_fam$`Pr(>F)`[1],
    disp_test_fam$`Pr(>F)`[1]
  )
)

print(overall_tests_att_aq)

# Pairwise site-blocked tests
pairwise_tests_att_aq <- bind_rows(
  pw %>%
    mutate(Level = "Species", Test = "PERMANOVA") %>%
    rename(p_raw = p),
  disp_pw_sp %>%
    mutate(Level = "Species", Test = "Dispersion") %>%
    rename(p_raw = p),
  pw_fam %>%
    mutate(Level = "Family", Test = "PERMANOVA") %>%
    rename(p_raw = p),
  disp_pw_fam %>%
    mutate(Level = "Family", Test = "Dispersion") %>%
    rename(p_raw = p)
) %>%
  select(Level, Test, comparison, everything())

print(pairwise_tests_att_aq)
################################################################################
################################################################################
################################################################################


################################################################################
# Euler/Venn plots for Attenborough - SPECIES 
################################################################################

# Get unique species names
mic_species  <- unique(names(attenborough_aq_wide_pa)[-1])          # exclude Sampling Point COLUMN
edna_species <- unique(names(attenborough_aq_edna_wide_pa)[-1])
bulk_species <- unique(names(attenborough_aq_meta_wide_pa)[-1])

# Species only found in one and shared  monitoring approaches
only_mic  <- setdiff(mic_species, union(edna_species, bulk_species))
only_edna <- setdiff(edna_species, union(mic_species, bulk_species))
only_bulk <- setdiff(bulk_species, union(mic_species, edna_species))

mic_edna_shared  <- intersect(mic_species, edna_species)
mic_bulk_shared  <- intersect(mic_species, bulk_species)
edna_bulk_shared <- intersect(edna_species, bulk_species)

all_three_shared <- Reduce(intersect, list(mic_species, edna_species, bulk_species))

# make df that combines them
species_summary_aq <- data.frame(
  Species = c(
    only_mic,
    only_edna,
    only_bulk,
    setdiff(mic_edna_shared, all_three_shared),
    setdiff(mic_bulk_shared, all_three_shared),
    setdiff(edna_bulk_shared, all_three_shared),
    all_three_shared
  ),
  Detected_By = c(
    rep("Microscope only", length(only_mic)),
    rep("eDNA only", length(only_edna)),
    rep("Bulk only", length(only_bulk)),
    rep("Microscope + eDNA", length(setdiff(mic_edna_shared, all_three_shared))),
    rep("Microscope + Bulk", length(setdiff(mic_bulk_shared, all_three_shared))),
    rep("eDNA + Bulk", length(setdiff(edna_bulk_shared, all_three_shared))),
    rep("All three", length(all_three_shared))
  ),
  stringsAsFactors = FALSE
)

# get numbers for each group - need to add number manually to euler plot
# (find a way to do this automatically potentially)

cat("Microscope only:", length(only_mic), "\n")
cat("eDNA only:", length(only_edna), "\n")
cat("Bulk only:", length(only_bulk), "\n")
cat("Microscope + eDNA:", length(setdiff(mic_edna_shared, all_three_shared)), "\n")
cat("Microscope + Bulk:", length(setdiff(mic_bulk_shared, all_three_shared)), "\n")
cat("eDNA + Bulk:", length(setdiff(edna_bulk_shared, all_three_shared)), "\n")
cat("All three:", length(all_three_shared), "\n")

# species from each method
attenborough_aq_morph_species_only %>% distinct(Species_trim) %>% nrow()
attenborough_aq_edna_species_only  %>% distinct(Species_trim) %>% nrow()
attenborough_aq_meta_species_only  %>% distinct(Species_trim) %>% nrow()

################################################################################
# Make Euler plots for Attenborough - SPECIES - (Attenborough)
################################################################################

# Region sizes (need to enter manually from above)

VennTri <- euler(c(
  "Morphological"              = 17,  
  "eDNA"                    = 108, 
  "Bulk"                    = 55,  
  "Morphological&eDNA"         = 4,   
  "Morphological&Bulk"         = 13,  
  "eDNA&Bulk"               = 29,  
  "Morphological&eDNA&Bulk"    = 3    
))


attenborough_euler_aquatic_plot <- plot(
  VennTri,
  quantities = list(
    type = "counts",   
    fontsize = 12,
    col = "black",
    font = 2
  ),
  labels = list(
    labels = c("Morph", "eDNA", "Bulk"),
    col = "black",
    fontsize = 12,
    font = 2,
    cex = 0.9
  ),
  fills = list(
    fill = c("#7A8F3A", "#4C8C8C", "#8A6F8F"),
    alpha = 0.6
  ),
  edges = list(
    col = "black",
    lwd = 0.7
  )
)

################################################################################
#  Euler plots for Attenborough - SPECIES - for manuscript
################################################################################

attenborough_euler_aquatic_plot




################################################################################
# Euler/Venn plots for Attenborough - Family 
################################################################################

# Get individual family names
mic_families  <- unique(names(attenborough_aq_wide_pa_family)[-1])          
edna_families <- unique(names(attenborough_aq_edna_wide_pa_family)[-1])
bulk_families <- unique(names(attenborough_aq_meta_wide_pa_family)[-1])

# Pairwise and shared comparisons
only_mic_fam  <- setdiff(mic_families, union(edna_families, bulk_families))
only_edna_fam <- setdiff(edna_families, union(mic_families, bulk_families))
only_bulk_fam <- setdiff(bulk_families, union(mic_families, edna_families))

mic_edna_shared_fam  <- intersect(mic_families, edna_families)
mic_bulk_shared_fam  <- intersect(mic_families, bulk_families)
edna_bulk_shared_fam <- intersect(edna_families, bulk_families)

all_three_shared_fam <- Reduce(intersect, list(mic_families, edna_families, bulk_families))

# make the combined summary data frame
family_summary_aq <- data.frame(
  Family = c(
    only_mic_fam,
    only_edna_fam,
    only_bulk_fam,
    setdiff(mic_edna_shared_fam, all_three_shared_fam),
    setdiff(mic_bulk_shared_fam, all_three_shared_fam),
    setdiff(edna_bulk_shared_fam, all_three_shared_fam),
    all_three_shared_fam
  ),
  Detected_By = c(
    rep("Microscope only", length(only_mic_fam)),
    rep("eDNA only", length(only_edna_fam)),
    rep("Bulk only", length(only_bulk_fam)),
    rep("Microscope + eDNA", length(setdiff(mic_edna_shared_fam, all_three_shared_fam))),
    rep("Microscope + Bulk", length(setdiff(mic_bulk_shared_fam, all_three_shared_fam))),
    rep("eDNA + Bulk", length(setdiff(edna_bulk_shared_fam, all_three_shared_fam))),
    rep("All three", length(all_three_shared_fam))
  ),
  stringsAsFactors = FALSE
)


# get numbers for each group - need to add number manually to euler plot

cat("Microscope only (families):", length(only_mic_fam), "\n")
cat("eDNA only (families):", length(only_edna_fam), "\n")
cat("Bulk only (families):", length(only_bulk_fam), "\n")
cat("Microscope + eDNA (families):", length(setdiff(mic_edna_shared_fam, all_three_shared_fam)), "\n")
cat("Microscope + Bulk (families):", length(setdiff(mic_bulk_shared_fam, all_three_shared_fam)), "\n")
cat("eDNA + Bulk (families):", length(setdiff(edna_bulk_shared_fam, all_three_shared_fam)), "\n")
cat("All three (families):", length(all_three_shared_fam), "\n")


# Region sizes for family (need to enter manually from above)
VennTri <- euler(c(
  "Morphological"              = 15,  # A only
  "eDNA"                    = 48, # B only
  "Bulk"                    = 25,  # C only
  "Morphological&eDNA"         = 12,   # A∩B only
  "Morphological&Bulk"         = 14,  # A∩C only
  "eDNA&Bulk"               = 7,  # B∩C only
  "Morphological&eDNA&Bulk"    = 18    # A∩B∩C
))


attenborough_euler_aquatic_plot_family <- plot(
  VennTri,
  quantities = list(
    type = "counts",   
    fontsize = 12,
    col = "black",
    font = 2
  ),
  labels = list(
    labels = c("Morph", "eDNA", "Bulk"),
    col = "black",
    fontsize = 12,
    font = 2,
    cex = 0.9
  ),
  fills = list(
    fill = c("#7A8F3A", "#4C8C8C", "#8A6F8F"),
    alpha = 0.6
  ),
  edges = list(
    col = "black",
    lwd = 1.5
  )
)

################################################################################
#  Euler plot for Attenborough - family - for supplementary
################################################################################

attenborough_euler_aquatic_plot_family


################################################################################
# Attenborough aquatic invert orders list
################################################################################

# morphological id
orders_morph <- attenborough_aq_morph %>%
  filter(!is.na(order), order != "") %>%
  mutate(order_trim = str_squish(as.character(order))) %>%
  pull(order_trim) %>%
  unique() %>%
  sort()

#edna
orders_edna <- attenborough_aq_edna %>%
  filter(!is.na(order), order != "") %>%
  mutate(order_trim = str_squish(as.character(order))) %>%
  pull(order_trim) %>%
  unique() %>%
  sort()

# bulk metabarocding
orders_meta <- attenborough_aq_meta %>%
  filter(!is.na(order), order != "") %>%
  mutate(order_trim = str_squish(as.character(order))) %>%
  pull(order_trim) %>%
  unique() %>%
  sort()





################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################






# KNEPP DATA NOW
################################################################################
#Load all of the Knepp data (Morpological id, eDNA and bulk (called meta in script))
################################################################################
knepp_aq_morph <- read.csv("Morph_ID_Knepp_Proc.csv", header = TRUE)
knepp_aq_edna  <- read.csv("eDNA_invert_Knepp_Proc.csv", header = TRUE)
knepp_aq_meta  <- read.csv("Meta_Knepp_Proc.csv", header = TRUE)

# have a look
str(knepp_aq_morph)
str(knepp_aq_edna)
str(knepp_aq_meta)

knepp_aq_morph$Sampling.Point <- factor(knepp_aq_morph$Sampling.Point)
knepp_aq_morph$Species        <- factor(knepp_aq_morph$Species)
knepp_aq_morph$family        <- factor(knepp_aq_morph$family)

knepp_aq_edna$Sampling.Point  <- factor(knepp_aq_edna$Sampling.Point)
knepp_aq_edna$Species         <- factor(knepp_aq_edna$Species)
knepp_aq_edna$family        <- factor(knepp_aq_edna$family)

knepp_aq_meta$Sampling.Point  <- factor(knepp_aq_meta$Sampling.Point)
knepp_aq_meta$Species         <- factor(knepp_aq_meta$Species)
knepp_aq_meta$family         <- factor(knepp_aq_meta$family)


#################################################################################
# To standardise analysis keep those only detected to species level (Knepp) (morph)
################################################################################

knepp_aq_morph_species_only <- knepp_aq_morph %>%
  filter(!is.na(Species)) %>%
  mutate(Species_trim = str_squish(as.character(Species))) %>%
  filter(str_count(Species_trim, " ") == 1)

################################################################################
# Make species-level presence/absence matrix for morph (Knepp) (morph)
################################################################################

knepp_aq_pa <- knepp_aq_morph_species_only %>%
  mutate(Presence_Absence = as.integer(value > 0)) %>%
  select(Sampling.Point, Species_trim, Presence_Absence)

knepp_aq_wide_pa <- knepp_aq_pa %>%
  pivot_wider(
    names_from  = Species_trim,
    values_from = Presence_Absence,
    values_fill = 0,
    values_fn   = max
  ) %>%
  group_by(Sampling.Point) %>%
  summarise(across(where(is.numeric), ~ as.integer(any(. > 0))), .groups = "drop")

knepp_aq_species_matrix <- knepp_aq_wide_pa %>%
  arrange(Sampling.Point)

rownames(knepp_aq_species_matrix) <- knepp_aq_species_matrix$Sampling.Point
knepp_aq_species_matrix <- as.matrix(knepp_aq_species_matrix %>% select(-Sampling.Point))


################################################################################
# Make family-level presence/absence matrix for morph (Knepp) (morph)
################################################################################

knepp_aq_morph_family_only <- knepp_aq_morph %>%
  filter(!is.na(family) & family != "") %>%
  mutate(family_trim = str_squish(as.character(family)))

knepp_aq_pa_family <- knepp_aq_morph_family_only %>%
  mutate(Presence_Absence = as.integer(value > 0)) %>%
  select(Sampling.Point, family_trim, Presence_Absence)

knepp_aq_wide_pa_family <- knepp_aq_pa_family %>%
  pivot_wider(
    names_from  = family_trim,
    values_from = Presence_Absence,
    values_fill = 0,
    values_fn   = max
  ) %>%
  group_by(Sampling.Point) %>%
  summarise(across(where(is.numeric), ~ as.integer(any(. > 0))), .groups = "drop")

knepp_aq_family_matrix <- knepp_aq_wide_pa_family %>%
  arrange(Sampling.Point)

rownames(knepp_aq_family_matrix) <- knepp_aq_family_matrix$Sampling.Point
knepp_aq_family_matrix <- as.matrix(knepp_aq_family_matrix %>% select(-Sampling.Point))



################################################################################
# Make species-level presence/absence matrix for eDNA (Knepp) (eDNA)
################################################################################

# Remove Chordata eDNA data had chordates in them (which arent aqutic inverts obviously)
knepp_aq_edna2 <- knepp_aq_edna %>%
  filter(phylum != "Chordata")

knepp_aq_edna_species_only <- knepp_aq_edna2 %>%
  filter(!is.na(Species)) %>%
  mutate(Species_trim = str_squish(as.character(Species))) %>%
  filter(str_count(Species_trim, " ") == 1)

knepp_aq_edna_pa <- knepp_aq_edna_species_only %>%
  mutate(Presence_Absence = as.integer(value > 0)) %>%
  select(Sampling.Point, Species_trim, Presence_Absence)

knepp_aq_edna_wide_pa <- knepp_aq_edna_pa %>%
  pivot_wider(
    names_from  = Species_trim,
    values_from = Presence_Absence,
    values_fill = 0,
    values_fn   = max
  ) %>%
  group_by(Sampling.Point) %>%
  summarise(across(where(is.numeric), ~ as.integer(any(. > 0))), .groups = "drop")

knepp_edna_aq_species_matrix <- knepp_aq_edna_wide_pa %>%
  arrange(Sampling.Point)

rownames(knepp_edna_aq_species_matrix) <- knepp_edna_aq_species_matrix$Sampling.Point
knepp_edna_aq_species_matrix <- as.matrix(knepp_edna_aq_species_matrix %>% select(-Sampling.Point))


################################################################################
# Make family-level presence/absence matrix for eDNA (Knepp) (eDNA)
################################################################################

knepp_aq_edna_family_only <- knepp_aq_edna2 %>%
  filter(!is.na(family) & family != "") %>%
  mutate(family_trim = str_squish(as.character(family)))

knepp_aq_edna_pa_family <- knepp_aq_edna_family_only %>%
  mutate(Presence_Absence = as.integer(value > 0)) %>%
  select(Sampling.Point, family_trim, Presence_Absence)

knepp_aq_edna_wide_pa_family <- knepp_aq_edna_pa_family %>%
  pivot_wider(
    names_from  = family_trim,
    values_from = Presence_Absence,
    values_fill = 0,
    values_fn   = max
  ) %>%
  group_by(Sampling.Point) %>%
  summarise(across(where(is.numeric), ~ as.integer(any(. > 0))), .groups = "drop")

knepp_edna_aq_family_matrix <- knepp_aq_edna_wide_pa_family %>%
  arrange(Sampling.Point)

rownames(knepp_edna_aq_family_matrix) <- knepp_edna_aq_family_matrix$Sampling.Point
knepp_edna_aq_family_matrix <- as.matrix(knepp_edna_aq_family_matrix %>% select(-Sampling.Point))


################################################################################
# Make species-level presence/absence matrix for bulk metabarcoding (Knepp) (meta)
################################################################################

# found non-aquatic inverts which i need to remove
orders_drop_meta <- c(
  "Artiodactyla",
  "Caudata",
  "Cypriniformes",
  "Gasterosteiformes",
  "Gruiformes"
)

knepp_aq_meta <- knepp_aq_meta %>%
  mutate(order_trim = str_squish(as.character(order))) %>%
  filter(!(order_trim %in% orders_drop_meta)) %>%
  select(-order_trim)


knepp_aq_meta_species_only <- knepp_aq_meta %>%
  filter(!is.na(Species)) %>%
  mutate(Species_trim = str_squish(as.character(Species))) %>%
  filter(str_count(Species_trim, " ") == 1)

knepp_aq_meta_pa <- knepp_aq_meta_species_only %>%
  mutate(Presence_Absence = as.integer(value > 0)) %>%
  select(Sampling.Point, Species_trim, Presence_Absence)

knepp_aq_meta_wide_pa <- knepp_aq_meta_pa %>%
  pivot_wider(
    names_from  = Species_trim,
    values_from = Presence_Absence,
    values_fill = 0,
    values_fn   = max
  ) %>%
  group_by(Sampling.Point) %>%
  summarise(across(where(is.numeric), ~ as.integer(any(. > 0))), .groups = "drop")

knepp_meta_aq_species_matrix <- knepp_aq_meta_wide_pa %>%
  arrange(Sampling.Point)

rownames(knepp_meta_aq_species_matrix) <- knepp_meta_aq_species_matrix$Sampling.Point
knepp_meta_aq_species_matrix <- as.matrix(knepp_meta_aq_species_matrix %>% select(-Sampling.Point))


################################################################################
# Make family-level presence/absence matrix for bulk metabarcoding (Knepp) (meta)
################################################################################


knepp_aq_meta_family_only <- knepp_aq_meta %>%
  filter(!is.na(family) & family != "") %>%
  mutate(family_trim = str_squish(as.character(family)))

knepp_aq_meta_pa_family <- knepp_aq_meta_family_only %>%
  mutate(Presence_Absence = as.integer(value > 0)) %>%
  select(Sampling.Point, family_trim, Presence_Absence)

knepp_aq_meta_wide_pa_family <- knepp_aq_meta_pa_family %>%
  pivot_wider(
    names_from  = family_trim,
    values_from = Presence_Absence,
    values_fill = 0,
    values_fn   = max
  ) %>%
  group_by(Sampling.Point) %>%
  summarise(across(where(is.numeric), ~ as.integer(any(. > 0))), .groups = "drop")

knepp_meta_aq_family_matrix <- knepp_aq_meta_wide_pa_family %>%
  arrange(Sampling.Point)

rownames(knepp_meta_aq_family_matrix) <- knepp_meta_aq_family_matrix$Sampling.Point
knepp_meta_aq_family_matrix <- as.matrix(knepp_meta_aq_family_matrix %>% select(-Sampling.Point))



################################################################################
# Species accumulation curves - inext  (Knepp)
################################################################################

# Species x sites matrices for iNEXT
mat_morph <- t(knepp_aq_species_matrix)
mat_edna  <- t(knepp_edna_aq_species_matrix)
mat_bulk  <- t(knepp_meta_aq_species_matrix)

for (m in list(mat_morph, mat_edna, mat_bulk)) storage.mode(m) <- "integer"
mat_morph[mat_morph > 1] <- 1
mat_edna[mat_edna > 1]   <- 1
mat_bulk[mat_bulk > 1]   <- 1

# 2x max effort
n_morph <- ncol(mat_morph)
n_edna  <- ncol(mat_edna)
n_bulk  <- ncol(mat_bulk)
endpoint_max <- 2 * max(n_morph, n_edna, n_bulk)

#  iNEXT
knepp_inext_three <- iNEXT(
  list(Morphological = mat_morph,
       eDNA          = mat_edna,
       Bulk          = mat_bulk),
  q        = 0,
  datatype = "incidence_raw",
  endpoint = endpoint_max,
  knots    = 80,
  se       = TRUE,
  conf     = 0.95
)

knepp_cols <- c(
  "Morphological" = "#7A8F3A",
  "eDNA"  = "#4C8C8C",
  "Bulk"  = "#8A6F8F"
)

# legend
leg_breaks <- c("Morphological", "eDNA", "Bulk")
leg_labels <- c("Morph", "eDNA", "Bulk")

knepp_aquatic_accum_plot <- ggiNEXT(knepp_inext_three, type = 1, color.var = "Assemblage") +
  scale_colour_manual(
    values = knepp_cols,
    breaks = leg_breaks,
    labels = leg_labels
  ) +
  scale_fill_manual(values = knepp_cols) +
  guides(
    fill     = "none",   
    linetype = "none",   
    shape    = "none",
    colour   = guide_legend(
      override.aes = list(
        linetype = 1,
        shape    = NA,
        fill     = unname(knepp_cols[leg_breaks]),
        alpha    = 0.25
      )
    )
  ) +
  labs(
    x = "Number of sampling events",
    y = "Species richness",
    colour = NULL,
    fill   = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.major     = element_blank(),
    panel.grid.minor     = element_blank(),
    axis.title           = element_text(size = 12),
    axis.text            = element_text(size = 12),
    legend.title         = element_blank(),
    legend.position      = "none",
    legend.justification = c("left", "top"),
    legend.background    = element_rect(fill = NA, colour = NA)
  )


################################################################################
# The actual Species accumulation curves - inext  (Knepp) - for the paper
################################################################################

knepp_aquatic_accum_plot


################################################################################
# Family accumulation curves - inext  (Knepp)
################################################################################


mat_morph_fam <- t(knepp_aq_family_matrix)
mat_edna_fam  <- t(knepp_edna_aq_family_matrix)
mat_bulk_fam  <- t(knepp_meta_aq_family_matrix)

mat_morph_fam[mat_morph_fam > 1] <- 1
mat_edna_fam[mat_edna_fam > 1]   <- 1
mat_bulk_fam[mat_bulk_fam > 1]   <- 1

storage.mode(mat_morph_fam) <- "integer"
storage.mode(mat_edna_fam)  <- "integer"
storage.mode(mat_bulk_fam)  <- "integer"

n_morph_fam <- ncol(mat_morph_fam)
n_edna_fam  <- ncol(mat_edna_fam)
n_bulk_fam  <- ncol(mat_bulk_fam)

endpoint_max_fam <- 2 * max(n_morph_fam, n_edna_fam, n_bulk_fam)

knepp_inext_three_fam <- iNEXT(
  list(Morphological = mat_morph_fam,
       eDNA          = mat_edna_fam,
       Bulk          = mat_bulk_fam),
  q        = 0,
  datatype = "incidence_raw",
  endpoint = endpoint_max_fam,
  knots    = 80,
  se       = TRUE,
  conf     = 0.95
)

knepp_cols <- c(
  "Morphological" = "#7A8F3A",
  "eDNA"  = "#4C8C8C",
  "Bulk"  = "#8A6F8F"
)

# legend 
leg_breaks <- c("Morphological", "eDNA", "Bulk")
leg_labels <- c("Morph", "eDNA", "Bulk")

knepp_aquatic_accum_plot_family <- ggiNEXT(knepp_inext_three_fam, type = 1, color.var = "Assemblage") +
  scale_colour_manual(
    values = knepp_cols,
    breaks = leg_breaks,
    labels = leg_labels
  ) +
  scale_fill_manual(values = knepp_cols) +
  guides(
    fill     = "none",   
    linetype = "none",   
    shape    = "none",
    colour   = guide_legend(
      override.aes = list(
        linetype = 1,
        shape    = NA,
        fill     = unname(knepp_cols[leg_breaks]),
        alpha    = 0.25
      )
    )
  ) +
  labs(
    x = "Number of sampling events",
    y = "Family richness",
    colour = NULL,
    fill   = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.major     = element_blank(),
    panel.grid.minor     = element_blank(),
    axis.title           = element_text(size = 12),
    axis.text            = element_text(size = 12),
    legend.title         = element_blank(),
    legend.position      = "none",
    legend.justification = c("left", "top"),
    legend.background    = element_rect(fill = NA, colour = NA)
  )

################################################################################
# The actual family accumulation curves - inext  (Knepp) - for the supplementary
################################################################################

knepp_aquatic_accum_plot_family


################################################################################
################################################################################
################################################################################




################################################################################
# NMDS species dissim and plot (Jaccard) (Knepp)
################################################################################

knepp_mic_pa  <- knepp_aq_wide_pa
knepp_edna_pa <- knepp_aq_edna_wide_pa
knepp_bulk_pa <- knepp_aq_meta_wide_pa

knepp_mic_species  <- knepp_mic_pa  %>% select(-Sampling.Point)
knepp_edna_species <- knepp_edna_pa %>% select(-Sampling.Point)
knepp_bulk_species <- knepp_bulk_pa %>% select(-Sampling.Point)

# Align species columns
knepp_all_species12 <- union(names(knepp_mic_species), names(knepp_edna_species))
knepp_all_species   <- union(knepp_all_species12, names(knepp_bulk_species))

knepp_miss_mic  <- setdiff(knepp_all_species, names(knepp_mic_species))
knepp_miss_edna <- setdiff(knepp_all_species, names(knepp_edna_species))
knepp_miss_bulk <- setdiff(knepp_all_species, names(knepp_bulk_species))

if (length(knepp_miss_mic)  > 0) knepp_mic_species[knepp_miss_mic]   <- 0L
if (length(knepp_miss_edna) > 0) knepp_edna_species[knepp_miss_edna] <- 0L
if (length(knepp_miss_bulk) > 0) knepp_bulk_species[knepp_miss_bulk] <- 0L

knepp_mic_species2  <- knepp_mic_species  %>% select(all_of(knepp_all_species))
knepp_edna_species2 <- knepp_edna_species %>% select(all_of(knepp_all_species))
knepp_bulk_species2 <- knepp_bulk_species %>% select(all_of(knepp_all_species))

# Combine with Source + IDs
knepp_mic_comb <- knepp_mic_species2 %>%
  mutate(
    Source = "Morph",
    ID     = as.character(knepp_mic_pa$Sampling.Point)
  )

knepp_edna_comb <- knepp_edna_species2 %>%
  mutate(
    Source = "eDNA",
    ID     = as.character(knepp_edna_pa$Sampling.Point)
  )

knepp_bulk_comb <- knepp_bulk_species2 %>%
  mutate(
    Source = "Bulk",
    ID     = as.character(knepp_bulk_pa$Sampling.Point)
  )

knepp_combined <- bind_rows(
  knepp_mic_comb,
  knepp_edna_comb,
  knepp_bulk_comb
)

# Remove empty samples and add Site as the nesting/blocking variable
# This extracts the middle part of IDs like x010724_site01_mzb1.
# If your Knepp IDs are site01_mzb1, it extracts site01.
knepp_combined_filt <- knepp_combined %>%
  filter(rowSums(select(., -Source, -ID)) > 0) %>%
  mutate(
    Site = if_else(
      stringr::str_detect(ID, "^[^_]+_[^_]+_"),
      stringr::str_replace(ID, "^[^_]+_([^_]+)_.*$", "\\1"),
      stringr::str_replace(ID, "_[^_]+$", "")
    ),
    Source = factor(Source, levels = c("Morph", "eDNA", "Bulk"))
  )

# Check nesting structure
table(knepp_combined_filt$Site, knepp_combined_filt$Source)
sum(is.na(knepp_combined_filt$Site))

################################################################################
# Run NMDS - SPECIES
################################################################################

set.seed(123)

knepp_species_mat <- knepp_combined_filt %>%
  select(-Source, -ID, -Site)

knepp_nmds <- metaMDS(
  knepp_species_mat,
  distance      = "jaccard",
  binary        = TRUE,
  k             = 2,
  trymax        = 400,
  autotransform = FALSE
)

knepp_scores <- as.data.frame(scores(knepp_nmds, display = "sites"))
knepp_scores$Source <- knepp_combined_filt$Source
knepp_scores$ID     <- knepp_combined_filt$ID
knepp_scores$Site   <- knepp_combined_filt$Site

################################################################################
# Plot NMDS - SPECIES
################################################################################

knepp_cols_nmds <- c(
  "Morph" = "#7A8F3A",
  "eDNA"  = "#4C8C8C",
  "Bulk"  = "#8A6F8F"
)

knepp_stress <- round(knepp_nmds$stress, 3)

knepp_aquatic_nmds_plot <- ggplot(
  knepp_scores,
  aes(x = NMDS1, y = NMDS2, colour = Source)
) +
  geom_point(size = 1.5) +
  stat_ellipse(
    aes(fill = Source, colour = Source, group = Source),
    geom = "polygon",
    alpha = 0.15,
    linewidth = 0.8,
    type = "t",
    level = 0.95,
    show.legend = FALSE
  ) +
  scale_colour_manual(values = knepp_cols_nmds) +
  scale_fill_manual(values = knepp_cols_nmds) +
  labs(x = "NMDS1", y = "NMDS2", colour = NULL, fill = NULL) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    legend.position = "none",
    axis.title = element_text(size = 12, colour = "black"),
    axis.text = element_text(size = 12, colour = "black"),
    plot.background = element_rect(fill = NA, colour = NA),
    panel.background = element_rect(fill = NA, colour = NA)
  ) +
  guides(colour = guide_legend(override.aes = list(size = 3, shape = 15))) +
  annotate(
    "text",
    x = -Inf, y = -Inf,
    label = paste0("Stress = ", knepp_stress),
    hjust = -0.1,
    vjust = -1,
    size = 3.2,
    colour = "black"
  )

knepp_aquatic_nmds_plot

################################################################################
# PERMANOVA - SPECIES, blocked by Site
################################################################################

knepp_dist_mat <- vegdist(
  knepp_species_mat,
  method = "jaccard",
  binary = TRUE
)

knepp_env <- knepp_combined_filt %>%
  select(Source, ID, Site) %>%
  mutate(
    Source = factor(Source, levels = c("Morph", "eDNA", "Bulk"))
  )

knepp_perm_species <- how(nperm = 999)
setBlocks(knepp_perm_species) <- knepp_env$Site

set.seed(123)

knepp_adon <- adonis2(
  knepp_dist_mat ~ Source,
  data = knepp_env,
  permutations = knepp_perm_species
)

print(knepp_adon)

################################################################################
# PERMDISP-style test - SPECIES
# Site-blocked permutation test on distance-to-centroid values
################################################################################

knepp_bd <- betadisper(knepp_dist_mat, knepp_env$Source)

dist_df_knepp_aq_sp <- data.frame(
  Distance = knepp_bd$distances,
  Source   = knepp_env$Source,
  Site     = knepp_env$Site
)

# Response must be supplied as a matrix for adonis2
knepp_disp_mat_sp <- as.matrix(dist_df_knepp_aq_sp$Distance)
rownames(knepp_disp_mat_sp) <- rownames(dist_df_knepp_aq_sp)

knepp_perm_disp_sp <- how(nperm = 999)
setBlocks(knepp_perm_disp_sp) <- dist_df_knepp_aq_sp$Site

set.seed(123)

knepp_disp_test_sp <- adonis2(
  knepp_disp_mat_sp ~ Source,
  data = dist_df_knepp_aq_sp,
  method = "euclidean",
  permutations = knepp_perm_disp_sp
)

print(knepp_disp_test_sp)

################################################################################
# Pairwise PERMDISP-style tests - SPECIES, blocked by Site
################################################################################

# Morph vs eDNA
knepp_disp_sp_me <- dist_df_knepp_aq_sp %>%
  filter(Source %in% c("Morph", "eDNA")) %>%
  droplevels()

knepp_disp_mat_sp_me <- as.matrix(knepp_disp_sp_me$Distance)
rownames(knepp_disp_mat_sp_me) <- rownames(knepp_disp_sp_me)

knepp_perm_disp_sp_me <- how(nperm = 999)
setBlocks(knepp_perm_disp_sp_me) <- knepp_disp_sp_me$Site

set.seed(123)

knepp_disp_test_sp_me <- adonis2(
  knepp_disp_mat_sp_me ~ Source,
  data = knepp_disp_sp_me,
  method = "euclidean",
  permutations = knepp_perm_disp_sp_me
)

# Morph vs Bulk
knepp_disp_sp_mb <- dist_df_knepp_aq_sp %>%
  filter(Source %in% c("Morph", "Bulk")) %>%
  droplevels()

knepp_disp_mat_sp_mb <- as.matrix(knepp_disp_sp_mb$Distance)
rownames(knepp_disp_mat_sp_mb) <- rownames(knepp_disp_sp_mb)

knepp_perm_disp_sp_mb <- how(nperm = 999)
setBlocks(knepp_perm_disp_sp_mb) <- knepp_disp_sp_mb$Site

set.seed(123)

knepp_disp_test_sp_mb <- adonis2(
  knepp_disp_mat_sp_mb ~ Source,
  data = knepp_disp_sp_mb,
  method = "euclidean",
  permutations = knepp_perm_disp_sp_mb
)

# eDNA vs Bulk
knepp_disp_sp_eb <- dist_df_knepp_aq_sp %>%
  filter(Source %in% c("eDNA", "Bulk")) %>%
  droplevels()

knepp_disp_mat_sp_eb <- as.matrix(knepp_disp_sp_eb$Distance)
rownames(knepp_disp_mat_sp_eb) <- rownames(knepp_disp_sp_eb)

knepp_perm_disp_sp_eb <- how(nperm = 999)
setBlocks(knepp_perm_disp_sp_eb) <- knepp_disp_sp_eb$Site

set.seed(123)

knepp_disp_test_sp_eb <- adonis2(
  knepp_disp_mat_sp_eb ~ Source,
  data = knepp_disp_sp_eb,
  method = "euclidean",
  permutations = knepp_perm_disp_sp_eb
)

knepp_disp_pw_sp <- tibble(
  comparison = c("Morph vs eDNA", "Morph vs Bulk", "eDNA vs Bulk"),
  F = c(
    knepp_disp_test_sp_me$F[1],
    knepp_disp_test_sp_mb$F[1],
    knepp_disp_test_sp_eb$F[1]
  ),
  p = c(
    knepp_disp_test_sp_me$`Pr(>F)`[1],
    knepp_disp_test_sp_mb$`Pr(>F)`[1],
    knepp_disp_test_sp_eb$`Pr(>F)`[1]
  )
) %>%
  mutate(p_adj_holm = p.adjust(p, method = "holm")) %>%
  arrange(p_adj_holm)

print(knepp_disp_pw_sp)

################################################################################
# Summary and plot of dispersion distances - SPECIES
################################################################################

dist_summary_knepp_aq_sp <- dist_df_knepp_aq_sp %>%
  dplyr::group_by(Source) %>%
  dplyr::summarise(
    mean_distance = mean(Distance, na.rm = TRUE),
    sd_distance   = sd(Distance, na.rm = TRUE),
    se_distance   = sd_distance / sqrt(n()),
    ci_lower      = mean_distance - 1.96 * se_distance,
    ci_upper      = mean_distance + 1.96 * se_distance,
    n             = n(),
    .groups = "drop"
  )

print(dist_summary_knepp_aq_sp)

distance_plot_knepp_aq_sp <- ggplot(
  dist_df_knepp_aq_sp,
  aes(x = Source, y = Distance, fill = Source)
) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.1, alpha = 0.6) +
  scale_colour_manual(values = knepp_cols_nmds) +
  scale_fill_manual(values = knepp_cols_nmds) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    title = "Knepp (Species)",
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
  )

distance_plot_knepp_aq_sp

################################################################################
# Pairwise PERMANOVA - SPECIES, blocked by Site
################################################################################

knepp_dat <- knepp_combined_filt %>%
  select(Source, ID, Site, everything()) %>%
  mutate(
    Source = factor(Source, levels = c("Morph", "eDNA", "Bulk"))
  )

knepp_sp_cols <- setdiff(names(knepp_dat), c("Source", "ID", "Site"))

# Morph vs eDNA
knepp_dat_me <- knepp_dat %>%
  filter(Source %in% c("Morph", "eDNA"))

knepp_dist_me <- vegdist(
  knepp_dat_me %>% select(all_of(knepp_sp_cols)),
  method = "jaccard",
  binary = TRUE
)

knepp_perm_me <- how(nperm = 999)
setBlocks(knepp_perm_me) <- knepp_dat_me$Site

set.seed(123)

knepp_adon_me <- adonis2(
  knepp_dist_me ~ Source,
  data = knepp_dat_me,
  permutations = knepp_perm_me
)

# Morph vs Bulk
knepp_dat_mb <- knepp_dat %>%
  filter(Source %in% c("Morph", "Bulk"))

knepp_dist_mb <- vegdist(
  knepp_dat_mb %>% select(all_of(knepp_sp_cols)),
  method = "jaccard",
  binary = TRUE
)

knepp_perm_mb <- how(nperm = 999)
setBlocks(knepp_perm_mb) <- knepp_dat_mb$Site

set.seed(123)

knepp_adon_mb <- adonis2(
  knepp_dist_mb ~ Source,
  data = knepp_dat_mb,
  permutations = knepp_perm_mb
)

# eDNA vs Bulk
knepp_dat_eb <- knepp_dat %>%
  filter(Source %in% c("eDNA", "Bulk"))

knepp_dist_eb <- vegdist(
  knepp_dat_eb %>% select(all_of(knepp_sp_cols)),
  method = "jaccard",
  binary = TRUE
)

knepp_perm_eb <- how(nperm = 999)
setBlocks(knepp_perm_eb) <- knepp_dat_eb$Site

set.seed(123)

knepp_adon_eb <- adonis2(
  knepp_dist_eb ~ Source,
  data = knepp_dat_eb,
  permutations = knepp_perm_eb
)

################################################################################
# Combine pairwise PERMANOVA results - SPECIES
################################################################################

knepp_pw <- tibble(
  comparison = c("Morph vs eDNA", "Morph vs Bulk", "eDNA vs Bulk"),
  R2 = c(knepp_adon_me$R2[1], knepp_adon_mb$R2[1], knepp_adon_eb$R2[1]),
  F  = c(knepp_adon_me$F[1],  knepp_adon_mb$F[1],  knepp_adon_eb$F[1]),
  p  = c(
    knepp_adon_me$`Pr(>F)`[1],
    knepp_adon_mb$`Pr(>F)`[1],
    knepp_adon_eb$`Pr(>F)`[1]
  )
) %>%
  mutate(p_adj_holm = p.adjust(p, method = "holm")) %>%
  arrange(p_adj_holm)

print(knepp_pw)

################################################################################
# FAMILY-level dissim and NMDS (presence/absence Jaccard) (Knepp)
################################################################################

knepp_mic_pa_fam  <- knepp_aq_wide_pa_family
knepp_edna_pa_fam <- knepp_aq_edna_wide_pa_family
knepp_bulk_pa_fam <- knepp_aq_meta_wide_pa_family

knepp_mic_family  <- knepp_mic_pa_fam  %>% select(-Sampling.Point)
knepp_edna_family <- knepp_edna_pa_fam %>% select(-Sampling.Point)
knepp_bulk_family <- knepp_bulk_pa_fam %>% select(-Sampling.Point)

# Align family columns
knepp_all_families12 <- union(names(knepp_mic_family), names(knepp_edna_family))
knepp_all_families   <- union(knepp_all_families12, names(knepp_bulk_family))

knepp_miss_mic_fam  <- setdiff(knepp_all_families, names(knepp_mic_family))
knepp_miss_edna_fam <- setdiff(knepp_all_families, names(knepp_edna_family))
knepp_miss_bulk_fam <- setdiff(knepp_all_families, names(knepp_bulk_family))

if (length(knepp_miss_mic_fam)  > 0) knepp_mic_family[knepp_miss_mic_fam]   <- 0L
if (length(knepp_miss_edna_fam) > 0) knepp_edna_family[knepp_miss_edna_fam] <- 0L
if (length(knepp_miss_bulk_fam) > 0) knepp_bulk_family[knepp_miss_bulk_fam] <- 0L

knepp_mic_family2  <- knepp_mic_family  %>% select(all_of(knepp_all_families))
knepp_edna_family2 <- knepp_edna_family %>% select(all_of(knepp_all_families))
knepp_bulk_family2 <- knepp_bulk_family %>% select(all_of(knepp_all_families))

# Combine with Source + IDs
knepp_mic_comb_fam <- knepp_mic_family2 %>%
  mutate(
    Source = "Morph",
    ID     = as.character(knepp_mic_pa_fam$Sampling.Point)
  )

knepp_edna_comb_fam <- knepp_edna_family2 %>%
  mutate(
    Source = "eDNA",
    ID     = as.character(knepp_edna_pa_fam$Sampling.Point)
  )

knepp_bulk_comb_fam <- knepp_bulk_family2 %>%
  mutate(
    Source = "Bulk",
    ID     = as.character(knepp_bulk_pa_fam$Sampling.Point)
  )

knepp_combined_fam <- bind_rows(
  knepp_mic_comb_fam,
  knepp_edna_comb_fam,
  knepp_bulk_comb_fam
)

# Remove empty samples and add Site as the nesting/blocking variable
knepp_combined_fam_filt <- knepp_combined_fam %>%
  filter(rowSums(select(., -Source, -ID)) > 0) %>%
  mutate(
    Site = if_else(
      stringr::str_detect(ID, "^[^_]+_[^_]+_"),
      stringr::str_replace(ID, "^[^_]+_([^_]+)_.*$", "\\1"),
      stringr::str_replace(ID, "_[^_]+$", "")
    ),
    Source = factor(Source, levels = c("Morph", "eDNA", "Bulk"))
  )

# Check nesting structure
table(knepp_combined_fam_filt$Site, knepp_combined_fam_filt$Source)
sum(is.na(knepp_combined_fam_filt$Site))

################################################################################
# Run NMDS - FAMILY
################################################################################

set.seed(123)

knepp_family_mat <- knepp_combined_fam_filt %>%
  select(-Source, -ID, -Site)

knepp_nmds_fam <- metaMDS(
  knepp_family_mat,
  distance      = "jaccard",
  binary        = TRUE,
  k             = 2,
  trymax        = 1000,
  autotransform = FALSE
)

knepp_scores_fam <- as.data.frame(scores(knepp_nmds_fam, display = "sites"))
knepp_scores_fam$Source <- knepp_combined_fam_filt$Source
knepp_scores_fam$ID     <- knepp_combined_fam_filt$ID
knepp_scores_fam$Site   <- knepp_combined_fam_filt$Site

################################################################################
# Plot NMDS - FAMILY
################################################################################

knepp_stress_fam <- round(knepp_nmds_fam$stress, 3)

knepp_aquatic_nmds_plot_family <- ggplot(
  knepp_scores_fam,
  aes(x = NMDS1, y = NMDS2, colour = Source)
) +
  geom_point(size = 1.5) +
  stat_ellipse(
    aes(fill = Source, colour = Source, group = Source),
    geom = "polygon",
    alpha = 0.15,
    linewidth = 0.8,
    type = "t",
    level = 0.95,
    show.legend = FALSE
  ) +
  scale_colour_manual(values = knepp_cols_nmds) +
  scale_fill_manual(values = knepp_cols_nmds) +
  labs(
    x = "NMDS1",
    y = "NMDS2",
    colour = NULL,
    fill = NULL
  ) +
  coord_cartesian(xlim = c(-1.5, NA), ylim = c(NA, 2)) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    legend.position = "none",
    axis.title = element_text(size = 12, colour = "black"),
    axis.text  = element_text(size = 12, colour = "black"),
    plot.background = element_rect(fill = NA, colour = NA),
    panel.background = element_rect(fill = NA, colour = NA)
  ) +
  guides(colour = guide_legend(override.aes = list(size = 3, shape = 15))) +
  annotate(
    "text",
    x = -Inf, y = -Inf,
    label = paste0("Stress = ", knepp_stress_fam),
    hjust = -0.1,
    vjust = -1,
    size = 3.2,
    colour = "black"
  )

knepp_aquatic_nmds_plot_family

################################################################################
# PERMANOVA - FAMILY, blocked by Site
################################################################################

knepp_dist_mat_fam <- vegdist(
  knepp_family_mat,
  method = "jaccard",
  binary = TRUE
)

knepp_env_fam <- knepp_combined_fam_filt %>%
  select(Source, ID, Site) %>%
  mutate(
    Source = factor(Source, levels = c("Morph", "eDNA", "Bulk"))
  )

knepp_perm_fam <- how(nperm = 999)
setBlocks(knepp_perm_fam) <- knepp_env_fam$Site

set.seed(123)

knepp_adon_fam <- adonis2(
  knepp_dist_mat_fam ~ Source,
  data = knepp_env_fam,
  permutations = knepp_perm_fam
)

print(knepp_adon_fam)

################################################################################
# PERMDISP-style test - FAMILY
# Site-blocked permutation test on distance-to-centroid values
################################################################################

knepp_bd_fam <- betadisper(knepp_dist_mat_fam, knepp_env_fam$Source)

dist_df_knepp_aq_fam <- data.frame(
  Distance = knepp_bd_fam$distances,
  Source   = knepp_env_fam$Source,
  Site     = knepp_env_fam$Site
)

# Response must be supplied as a matrix for adonis2
knepp_disp_mat_fam <- as.matrix(dist_df_knepp_aq_fam$Distance)
rownames(knepp_disp_mat_fam) <- rownames(dist_df_knepp_aq_fam)

knepp_perm_disp_fam <- how(nperm = 999)
setBlocks(knepp_perm_disp_fam) <- dist_df_knepp_aq_fam$Site

set.seed(123)

knepp_disp_test_fam <- adonis2(
  knepp_disp_mat_fam ~ Source,
  data = dist_df_knepp_aq_fam,
  method = "euclidean",
  permutations = knepp_perm_disp_fam
)

print(knepp_disp_test_fam)

################################################################################
# Pairwise PERMDISP-style tests - FAMILY, blocked by Site
################################################################################

# Morph vs eDNA
knepp_disp_fam_me <- dist_df_knepp_aq_fam %>%
  filter(Source %in% c("Morph", "eDNA")) %>%
  droplevels()

knepp_disp_mat_fam_me <- as.matrix(knepp_disp_fam_me$Distance)
rownames(knepp_disp_mat_fam_me) <- rownames(knepp_disp_fam_me)

knepp_perm_disp_fam_me <- how(nperm = 999)
setBlocks(knepp_perm_disp_fam_me) <- knepp_disp_fam_me$Site

set.seed(123)

knepp_disp_test_fam_me <- adonis2(
  knepp_disp_mat_fam_me ~ Source,
  data = knepp_disp_fam_me,
  method = "euclidean",
  permutations = knepp_perm_disp_fam_me
)

# Morph vs Bulk
knepp_disp_fam_mb <- dist_df_knepp_aq_fam %>%
  filter(Source %in% c("Morph", "Bulk")) %>%
  droplevels()

knepp_disp_mat_fam_mb <- as.matrix(knepp_disp_fam_mb$Distance)
rownames(knepp_disp_mat_fam_mb) <- rownames(knepp_disp_fam_mb)

knepp_perm_disp_fam_mb <- how(nperm = 999)
setBlocks(knepp_perm_disp_fam_mb) <- knepp_disp_fam_mb$Site

set.seed(123)

knepp_disp_test_fam_mb <- adonis2(
  knepp_disp_mat_fam_mb ~ Source,
  data = knepp_disp_fam_mb,
  method = "euclidean",
  permutations = knepp_perm_disp_fam_mb
)

# eDNA vs Bulk
knepp_disp_fam_eb <- dist_df_knepp_aq_fam %>%
  filter(Source %in% c("eDNA", "Bulk")) %>%
  droplevels()

knepp_disp_mat_fam_eb <- as.matrix(knepp_disp_fam_eb$Distance)
rownames(knepp_disp_mat_fam_eb) <- rownames(knepp_disp_fam_eb)

knepp_perm_disp_fam_eb <- how(nperm = 999)
setBlocks(knepp_perm_disp_fam_eb) <- knepp_disp_fam_eb$Site

set.seed(123)

knepp_disp_test_fam_eb <- adonis2(
  knepp_disp_mat_fam_eb ~ Source,
  data = knepp_disp_fam_eb,
  method = "euclidean",
  permutations = knepp_perm_disp_fam_eb
)

knepp_disp_pw_fam <- tibble(
  comparison = c("Morph vs eDNA", "Morph vs Bulk", "eDNA vs Bulk"),
  F = c(
    knepp_disp_test_fam_me$F[1],
    knepp_disp_test_fam_mb$F[1],
    knepp_disp_test_fam_eb$F[1]
  ),
  p = c(
    knepp_disp_test_fam_me$`Pr(>F)`[1],
    knepp_disp_test_fam_mb$`Pr(>F)`[1],
    knepp_disp_test_fam_eb$`Pr(>F)`[1]
  )
) %>%
  mutate(p_adj_holm = p.adjust(p, method = "holm")) %>%
  arrange(p_adj_holm)

print(knepp_disp_pw_fam)

################################################################################
# Summary and plot of dispersion distances - FAMILY
################################################################################

dist_summary_knepp_aq_fam <- dist_df_knepp_aq_fam %>%
  dplyr::group_by(Source) %>%
  dplyr::summarise(
    mean_distance = mean(Distance, na.rm = TRUE),
    sd_distance   = sd(Distance, na.rm = TRUE),
    se_distance   = sd_distance / sqrt(n()),
    ci_lower      = mean_distance - 1.96 * se_distance,
    ci_upper      = mean_distance + 1.96 * se_distance,
    n             = n(),
    .groups = "drop"
  )

print(dist_summary_knepp_aq_fam)

distance_plot_knepp_aq_fam <- ggplot(
  dist_df_knepp_aq_fam,
  aes(x = Source, y = Distance, fill = Source)
) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.1, alpha = 0.6) +
  scale_colour_manual(values = knepp_cols_nmds) +
  scale_fill_manual(values = knepp_cols_nmds) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    title = "Knepp (Family)",
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
  )

distance_plot_knepp_aq_fam

################################################################################
# Pairwise PERMANOVA - FAMILY, blocked by Site
################################################################################

knepp_dat_fam <- knepp_combined_fam_filt %>%
  select(Source, ID, Site, everything()) %>%
  mutate(
    Source = factor(Source, levels = c("Morph", "eDNA", "Bulk"))
  )

knepp_fam_cols <- setdiff(names(knepp_dat_fam), c("Source", "ID", "Site"))

# Morph vs eDNA
knepp_dat_fam_me <- knepp_dat_fam %>%
  filter(Source %in% c("Morph", "eDNA"))

knepp_dist_fam_me <- vegdist(
  knepp_dat_fam_me %>% select(all_of(knepp_fam_cols)),
  method = "jaccard",
  binary = TRUE
)

knepp_perm_fam_me <- how(nperm = 999)
setBlocks(knepp_perm_fam_me) <- knepp_dat_fam_me$Site

set.seed(123)

knepp_adon_fam_me <- adonis2(
  knepp_dist_fam_me ~ Source,
  data = knepp_dat_fam_me,
  permutations = knepp_perm_fam_me
)

# Morph vs Bulk
knepp_dat_fam_mb <- knepp_dat_fam %>%
  filter(Source %in% c("Morph", "Bulk"))

knepp_dist_fam_mb <- vegdist(
  knepp_dat_fam_mb %>% select(all_of(knepp_fam_cols)),
  method = "jaccard",
  binary = TRUE
)

knepp_perm_fam_mb <- how(nperm = 999)
setBlocks(knepp_perm_fam_mb) <- knepp_dat_fam_mb$Site

set.seed(123)

knepp_adon_fam_mb <- adonis2(
  knepp_dist_fam_mb ~ Source,
  data = knepp_dat_fam_mb,
  permutations = knepp_perm_fam_mb
)

# eDNA vs Bulk
knepp_dat_fam_eb <- knepp_dat_fam %>%
  filter(Source %in% c("eDNA", "Bulk"))

knepp_dist_fam_eb <- vegdist(
  knepp_dat_fam_eb %>% select(all_of(knepp_fam_cols)),
  method = "jaccard",
  binary = TRUE
)

knepp_perm_fam_eb <- how(nperm = 999)
setBlocks(knepp_perm_fam_eb) <- knepp_dat_fam_eb$Site

set.seed(123)

knepp_adon_fam_eb <- adonis2(
  knepp_dist_fam_eb ~ Source,
  data = knepp_dat_fam_eb,
  permutations = knepp_perm_fam_eb
)

################################################################################
# Combine pairwise PERMANOVA results - FAMILY
################################################################################

knepp_pw_fam <- tibble(
  comparison = c("Morph vs eDNA", "Morph vs Bulk", "eDNA vs Bulk"),
  R2 = c(knepp_adon_fam_me$R2[1], knepp_adon_fam_mb$R2[1], knepp_adon_fam_eb$R2[1]),
  F  = c(knepp_adon_fam_me$F[1],  knepp_adon_fam_mb$F[1],  knepp_adon_fam_eb$F[1]),
  p  = c(
    knepp_adon_fam_me$`Pr(>F)`[1],
    knepp_adon_fam_mb$`Pr(>F)`[1],
    knepp_adon_fam_eb$`Pr(>F)`[1]
  )
) %>%
  mutate(p_adj_holm = p.adjust(p, method = "holm")) %>%
  arrange(p_adj_holm)

print(knepp_pw_fam)

################################################################################
# Result tables for reporting - KNEPP
################################################################################

# Overall site-blocked tests
knepp_overall_tests <- tibble(
  Site = "Knepp",
  Level = c("Species", "Species", "Family", "Family"),
  Test = c("PERMANOVA", "Dispersion", "PERMANOVA", "Dispersion"),
  Term = c("Source", "Source", "Source", "Source"),
  Df = c(
    knepp_adon$Df[1],
    knepp_disp_test_sp$Df[1],
    knepp_adon_fam$Df[1],
    knepp_disp_test_fam$Df[1]
  ),
  F = c(
    knepp_adon$F[1],
    knepp_disp_test_sp$F[1],
    knepp_adon_fam$F[1],
    knepp_disp_test_fam$F[1]
  ),
  R2 = c(
    knepp_adon$R2[1],
    knepp_disp_test_sp$R2[1],
    knepp_adon_fam$R2[1],
    knepp_disp_test_fam$R2[1]
  ),
  p = c(
    knepp_adon$`Pr(>F)`[1],
    knepp_disp_test_sp$`Pr(>F)`[1],
    knepp_adon_fam$`Pr(>F)`[1],
    knepp_disp_test_fam$`Pr(>F)`[1]
  )
)

print(knepp_overall_tests)

# Pairwise site-blocked tests
knepp_pairwise_tests <- bind_rows(
  knepp_pw %>%
    mutate(
      Site = "Knepp",
      Level = "Species",
      Test = "PERMANOVA"
    ),
  knepp_disp_pw_sp %>%
    mutate(
      Site = "Knepp",
      Level = "Species",
      Test = "Dispersion",
      R2 = NA_real_
    ),
  knepp_pw_fam %>%
    mutate(
      Site = "Knepp",
      Level = "Family",
      Test = "PERMANOVA"
    ),
  knepp_disp_pw_fam %>%
    mutate(
      Site = "Knepp",
      Level = "Family",
      Test = "Dispersion",
      R2 = NA_real_
    )
) %>%
  rename(p_raw = p) %>%
  group_by(Level, Test) %>%
  mutate(p_adj_holm_within_test = p.adjust(p_raw, method = "holm")) %>%
  ungroup() %>%
  select(Site, Level, Test, comparison, R2, F, p_raw, p_adj_holm, p_adj_holm_within_test)

print(knepp_pairwise_tests)


################################################################################
################################################################################
################################################################################


################################################################################
# Euler/Venn plots for Knepp - SPECIES 
################################################################################

#  unique species names
knepp_mic_spp  <- unique(names(knepp_aq_wide_pa)[-1])
knepp_edna_spp <- unique(names(knepp_aq_edna_wide_pa)[-1])
knepp_bulk_spp <- unique(names(knepp_aq_meta_wide_pa)[-1])

# different and shared sets
knepp_only_mic    <- setdiff(knepp_mic_spp,  union(knepp_edna_spp, knepp_bulk_spp))
knepp_only_edna   <- setdiff(knepp_edna_spp, union(knepp_mic_spp,  knepp_bulk_spp))
knepp_only_bulk   <- setdiff(knepp_bulk_spp, union(knepp_mic_spp,  knepp_edna_spp))

knepp_mic_edna_sh  <- intersect(knepp_mic_spp,  knepp_edna_spp)
knepp_mic_bulk_sh  <- intersect(knepp_mic_spp,  knepp_bulk_spp)
knepp_edna_bulk_sh <- intersect(knepp_edna_spp, knepp_bulk_spp)

knepp_all_three_sh <- Reduce(intersect, list(knepp_mic_spp, knepp_edna_spp, knepp_bulk_spp))

knepp_species_summary_aq <- data.frame(
  Species = c(
    knepp_only_mic,
    knepp_only_edna,
    knepp_only_bulk,
    setdiff(knepp_mic_edna_sh,  knepp_all_three_sh),
    setdiff(knepp_mic_bulk_sh,  knepp_all_three_sh),
    setdiff(knepp_edna_bulk_sh, knepp_all_three_sh),
    knepp_all_three_sh
  ),
  Detected_By = c(
    rep("Microscope only", length(knepp_only_mic)),
    rep("eDNA only",       length(knepp_only_edna)),
    rep("Bulk only",       length(knepp_only_bulk)),
    rep("Microscope + eDNA", length(setdiff(knepp_mic_edna_sh,  knepp_all_three_sh))),
    rep("Microscope + Bulk", length(setdiff(knepp_mic_bulk_sh,  knepp_all_three_sh))),
    rep("eDNA + Bulk",       length(setdiff(knepp_edna_bulk_sh, knepp_all_three_sh))),
    rep("All three",         length(knepp_all_three_sh))
  ),
  stringsAsFactors = FALSE
)

# count the above
cat("Knepp - Microscope only:", length(knepp_only_mic), "\n")
cat("Knepp - eDNA only:",       length(knepp_only_edna), "\n")
cat("Knepp - Bulk only:",       length(knepp_only_bulk), "\n")
cat("Knepp - Microscope + eDNA:", length(setdiff(knepp_mic_edna_sh,  knepp_all_three_sh)), "\n")
cat("Knepp - Microscope + Bulk:", length(setdiff(knepp_mic_bulk_sh,  knepp_all_three_sh)), "\n")
cat("Knepp - eDNA + Bulk:",       length(setdiff(knepp_edna_bulk_sh, knepp_all_three_sh)), "\n")
cat("Knepp - All three:",         length(knepp_all_three_sh), "\n")

################################################################################
# Make euler plot Knepp species level
################################################################################

# take lengths from above (saves doing it manually like for Att)
knepp_counts <- list(
  morph_only  = length(knepp_only_mic),
  edna_only   = length(knepp_only_edna),
  bulk_only   = length(knepp_only_bulk),
  morph_edna  = length(setdiff(knepp_mic_edna_sh,  knepp_all_three_sh)),
  morph_bulk  = length(setdiff(knepp_mic_bulk_sh,  knepp_all_three_sh)),
  edna_bulk   = length(setdiff(knepp_edna_bulk_sh, knepp_all_three_sh)),
  all_three   = length(knepp_all_three_sh)
)

knepp_VennTri <- euler(c(
  "Morphological"           = knepp_counts$morph_only,
  "eDNA"                    = knepp_counts$edna_only,
  "Bulk"                    = knepp_counts$bulk_only,
  "Morphological&eDNA"      = knepp_counts$morph_edna,
  "Morphological&Bulk"      = knepp_counts$morph_bulk,
  "eDNA&Bulk"               = knepp_counts$edna_bulk,
  "Morphological&eDNA&Bulk" = knepp_counts$all_three
))

knepp_euler_aquatic_plot <- plot(
  knepp_VennTri,
  quantities = list(type = "counts", fontsize = 12, col = "black", font = 2),
  labels     = list(labels = c("Morph", "eDNA", "Bulk"),
                    col = "black",fontsize = 12, font = 2, cex = 0.9),
  fills      = list(fill = c("#7A8F3A", "#4C8C8C", "#8A6F8F"), alpha = 0.6),
  edges      = list(col = "black", lwd = 0.7)
)


################################################################################
# Knepp euler plot for SPECIES - Manuscript 
################################################################################

knepp_euler_aquatic_plot




################################################################################
# Euler/Venn plots for Knepp - FAMILY 
################################################################################

knepp_mic_fam  <- unique(names(knepp_aq_wide_pa_family)[-1])
knepp_edna_fam <- unique(names(knepp_aq_edna_wide_pa_family)[-1])
knepp_bulk_fam <- unique(names(knepp_aq_meta_wide_pa_family)[-1])

knepp_only_mic_fam  <- setdiff(knepp_mic_fam,  union(knepp_edna_fam, knepp_bulk_fam))
knepp_only_edna_fam <- setdiff(knepp_edna_fam, union(knepp_mic_fam,  knepp_bulk_fam))
knepp_only_bulk_fam <- setdiff(knepp_bulk_fam, union(knepp_mic_fam,  knepp_edna_fam))

knepp_mic_edna_sh_fam  <- intersect(knepp_mic_fam,  knepp_edna_fam)
knepp_mic_bulk_sh_fam  <- intersect(knepp_mic_fam,  knepp_bulk_fam)
knepp_edna_bulk_sh_fam <- intersect(knepp_edna_fam, knepp_bulk_fam)

knepp_all_three_sh_fam <- Reduce(intersect, list(knepp_mic_fam, knepp_edna_fam, knepp_bulk_fam))

knepp_family_summary_aq <- data.frame(
  Family = c(
    knepp_only_mic_fam,
    knepp_only_edna_fam,
    knepp_only_bulk_fam,
    setdiff(knepp_mic_edna_sh_fam,  knepp_all_three_sh_fam),
    setdiff(knepp_mic_bulk_sh_fam,  knepp_all_three_sh_fam),
    setdiff(knepp_edna_bulk_sh_fam, knepp_all_three_sh_fam),
    knepp_all_three_sh_fam
  ),
  Detected_By = c(
    rep("Microscope only", length(knepp_only_mic_fam)),
    rep("eDNA only",       length(knepp_only_edna_fam)),
    rep("Bulk only",       length(knepp_only_bulk_fam)),
    rep("Microscope + eDNA", length(setdiff(knepp_mic_edna_sh_fam,  knepp_all_three_sh_fam))),
    rep("Microscope + Bulk", length(setdiff(knepp_mic_bulk_sh_fam,  knepp_all_three_sh_fam))),
    rep("eDNA + Bulk",       length(setdiff(knepp_edna_bulk_sh_fam, knepp_all_three_sh_fam))),
    rep("All three",         length(knepp_all_three_sh_fam))
  ),
  stringsAsFactors = FALSE
)

cat("Knepp families - Microscope only:",     length(knepp_only_mic_fam),  "\n")
cat("Knepp families - eDNA only:",           length(knepp_only_edna_fam), "\n")
cat("Knepp families - Bulk only:",           length(knepp_only_bulk_fam), "\n")
cat("Knepp families - Microscope + eDNA:",   length(setdiff(knepp_mic_edna_sh_fam,  knepp_all_three_sh_fam)), "\n")
cat("Knepp families - Microscope + Bulk:",   length(setdiff(knepp_mic_bulk_sh_fam,  knepp_all_three_sh_fam)), "\n")
cat("Knepp families - eDNA + Bulk:",         length(setdiff(knepp_edna_bulk_sh_fam, knepp_all_three_sh_fam)), "\n")
cat("Knepp families - All three:",           length(knepp_all_three_sh_fam), "\n")



################################################################################
# Euler plot for families Knepp                                                   #
################################################################################

knepp_counts_fam <- list(
  morph_only  = length(knepp_only_mic_fam),
  edna_only   = length(knepp_only_edna_fam),
  bulk_only   = length(knepp_only_bulk_fam),
  morph_edna  = length(setdiff(knepp_mic_edna_sh_fam,  knepp_all_three_sh_fam)),
  morph_bulk  = length(setdiff(knepp_mic_bulk_sh_fam,  knepp_all_three_sh_fam)),
  edna_bulk   = length(setdiff(knepp_edna_bulk_sh_fam, knepp_all_three_sh_fam)),
  all_three   = length(knepp_all_three_sh_fam)
)

knepp_VennTri_fam <- euler(c(
  "Morphological"           = knepp_counts_fam$morph_only,
  "eDNA"                    = knepp_counts_fam$edna_only,
  "Bulk"                    = knepp_counts_fam$bulk_only,
  "Morphological&eDNA"      = knepp_counts_fam$morph_edna,
  "Morphological&Bulk"      = knepp_counts_fam$morph_bulk,
  "eDNA&Bulk"               = knepp_counts_fam$edna_bulk,
  "Morphological&eDNA&Bulk" = knepp_counts_fam$all_three
))

knepp_euler_aquatic_plot_family <- plot(
  knepp_VennTri_fam,
  quantities = list(type = "counts", fontsize = 12, col = "black", font = 2),
  labels     = list(labels = c("Morph", "eDNA", "Bulk"),
                    col = "black", font = 2,fontsize = 12, cex = 0.9),
  fills      = list(fill = c("#7A8F3A", "#4C8C8C", "#8A6F8F"), alpha = 0.6),
  edges      = list(col = "black", lwd = 1.5)
)


################################################################################
# Knepp euler plot for FAMILY - for supplementary 
################################################################################

knepp_euler_aquatic_plot_family



################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################



################################################################################
#  Combine all plots for the manuscript multifigure (Figure 2a-f) - Species first
################################################################################

# look at all plots 
attenborough_euler_aquatic_plot
att_aquatic_accum_plot
attenborough_aquatic_nmds_plot
knepp_aquatic_accum_plot
knepp_aquatic_nmds_plot
knepp_euler_aquatic_plot

# euler package doesnt make managable figures, so make as ggplot object
att_venn <- ggplotify::as.ggplot(attenborough_euler_aquatic_plot)
kne_venn <- ggplotify::as.ggplot(knepp_euler_aquatic_plot)


# Make aquatic invert species plot then i can add labels etc. later
panel_main_aquatic <- (
  (att_aquatic_accum_plot | att_venn | attenborough_aquatic_nmds_plot) /
    (knepp_aquatic_accum_plot         | kne_venn | knepp_aquatic_nmds_plot)
)&
  theme(
    plot.margin = margin(10, 10, 10, 10)
  )


# edit species multipanel figure 
panel_labeled_aquatic <- ggdraw(panel_main_aquatic) +
  draw_label("(a)", x = 0.03, y = 0.98, hjust = 0, vjust = 1, fontface = "bold", size = 14) +
  draw_label("(c)", x = 0.36, y = 0.98, hjust = 0, vjust = 1, fontface = "bold", size = 14) +
  draw_label("(e)", x = 0.66, y = 0.98, hjust = 0, vjust = 1, fontface = "bold", size = 14) +
  draw_label("(b)", x = 0.03, y = 0.51, hjust = 0, vjust = 1, fontface = "bold", size = 14) +
  draw_label("(d)", x = 0.36, y = 0.51, hjust = 0, vjust = 1, fontface = "bold", size = 14) +
  draw_label("(f)", x = 0.66, y = 0.51, hjust = 0, vjust = 1, fontface = "bold", size = 14) +
  draw_label("Attenborough", x = 0.015, y = 0.77, angle = 90,
             fontface = "bold", size = 18, hjust = 0.5, vjust = 0.1) +
  draw_label("Knepp",        x = 0.015, y = 0.26, angle = 90,
             fontface = "bold", size = 18, hjust = 0.3, vjust = 0.1)


################################################################################
# Species multifigure panel aquatic invert (Figure 2)
################################################################################

panel_labeled_aquatic

################################################################################
#save it
################################################################################

w <- dev.size("in")[1]
h <- dev.size("in")[2] 

ggsave(
  "aquatic_invert_species_paper1_figure.tiff",
  plot = panel_labeled_aquatic,
  width = w,
  height = h,
  units = "in",
  dpi = 600,
  compression = "lzw"
)






################################################################################
################################################################################
################################################################################


################################################################################
#  Combine all family aquatic invert plots for the supp multifigure - families 
################################################################################


# look at all plots 
attenborough_euler_aquatic_plot_family
attenborough_aquatic_nmds_plot_family
att_aquatic_accum_plot_family
knepp_aquatic_accum_plot_family
knepp_aquatic_nmds_plot_family
knepp_euler_aquatic_plot_family

# euler package doesnt make managable figures, so make as ggplot object
att_venn_family <- ggplotify::as.ggplot(attenborough_euler_aquatic_plot_family)
kne_venn_family <- ggplotify::as.ggplot(knepp_euler_aquatic_plot_family)

# Make aquatic invert family plot then  add labels after
panel_main_aquatic_family <- (
  (att_aquatic_accum_plot_family | att_venn_family | attenborough_aquatic_nmds_plot_family) /
    (knepp_aquatic_accum_plot_family | kne_venn_family | knepp_aquatic_nmds_plot_family)
) &
  theme(
    plot.margin = margin(10, 10, 10, 10)
  )


# edit species multipanel figure 
panel_labeled_aquatic_family <- ggdraw(panel_main_aquatic_family) +
  draw_label("(a)", x = 0.03, y = 0.98, hjust = 0, vjust = 1, fontface = "bold", size = 14) +
  draw_label("(c)", x = 0.36, y = 0.98, hjust = 0, vjust = 1, fontface = "bold", size = 14) +
  draw_label("(e)", x = 0.66, y = 0.98, hjust = 0, vjust = 1, fontface = "bold", size = 14) +
  draw_label("(b)", x = 0.03, y = 0.51, hjust = 0, vjust = 1, fontface = "bold", size = 14) +
  draw_label("(d)", x = 0.36, y = 0.51, hjust = 0, vjust = 1, fontface = "bold", size = 14) +
  draw_label("(f)", x = 0.66, y = 0.51, hjust = 0, vjust = 1, fontface = "bold", size = 14) +
  draw_label("Attenborough", x = 0.015, y = 0.77, angle = 90,
             fontface = "bold", size = 18, hjust = 0.5, vjust = 0.1) +
  draw_label("Knepp",        x = 0.015, y = 0.27, angle = 90,
             fontface = "bold", size = 18, hjust = 0.3, vjust = 0.1)


################################################################################
# Family multifigure panel aquatic invert for supplementary material
################################################################################

panel_labeled_aquatic_family


################################################################################
#save it
################################################################################

w <- dev.size("in")[1]
h <- dev.size("in")[2] 

ggsave(
  "aquatic_invert_family_paper1_figure.tiff",
  plot = panel_labeled_aquatic_family,
  width = w,
  height = h,
  units = "in",
  dpi = 600,
  compression = "lzw"
)


################################################################################
################################################################################
################################################################################



################################################################################
# Combined NMDS distance plots for supplementary Attenborough and Knepp
################################################################################

#Species aquatic ivnert distance combined plot
combined_plot_distance_aq_sp <- (distance_plot_att_aq_sp / distance_plot_knepp_aq_sp) +
  plot_annotation(
    tag_levels = "a",
    tag_prefix = "(",
    tag_suffix = ")"
  ) 

combined_plot_distance_aq_sp

#Families aquatic invert distance combined plot

combined_plot_distance_aq_fam <- (distance_plot_att_aq_fam / distance_plot_knepp_aq_fam) +
  plot_annotation(
    tag_levels = "a",
    tag_prefix = "(",
    tag_suffix = ")"
  ) 

combined_plot_distance_aq_fam

################################################################################
################################################################################
################################################################################

