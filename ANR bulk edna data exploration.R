#Explore ANR freshwater datasheets, runs species accumulation curves and display taxonomic data

#Read in libraries
library("ggplot2")
library("reshape2")
library(colorspace)
library(vegan)
library(dplyr)
library(tidyr)


## Read in eDNA data for insects and bulk invertebrates data from sampling at Attenborough in July 2024.

edna<-read.csv("C:\\Users\\helford\\Documents\\Postdoc paper\\WP3_data\\Freshwater\\2024\\NatureMetrics\\0724_eDNA_ANR_percentage.csv")

bulk<-read.csv("C:\\Users\\helford\\Documents\\Postdoc paper\\WP3_data\\Freshwater\\2024\\NatureMetrics\\0724_Bulk_invert_ANR_percentage.csv")

morph<-read.csv("C:\\Users\\helford\\Documents\\Postdoc paper\\WP3_data\\Freshwater\\2024\\Morphology_Att_MZB-ID-Feb24.csv")


str(edna)

# FOR METABARCODING Convert to long format
bulk_long <- bulk %>%
  pivot_longer(
    cols = starts_with("X"),  # Select all sample-specific columns
    names_to = "Sample",        # New column for sample identifiers
    values_to = "Perc"       # New column for percentage of total read of sample
  ) %>%
  filter(Perc > 0)           # Filter out rows with 0 count (if needed)
# Remove Count column if no longer needed

#Add column for site - Remove ".MZB" and any following digit(s)
bulk_long <- bulk_long %>%
mutate(Site = sub("\\.MZB[0-9]+$", "", Sample),
         Site = factor(Site) )
str(bulk_long)

# Add a new column for presence/absence based on the Count column
bulk_long <- bulk_long %>%
  mutate(Presence_Absence = as.numeric(ifelse(Perc > 0, 1, 0)))  # 1 for presence, 0 for absence

# FOR EDNA Convert to long format
edna_long <- edna %>%
  pivot_longer(
    cols = starts_with("X"),  # Select all sample-specific columns
    names_to = "Sample",      # New column for sample identifiers
    values_to = "Perc"        # New column for percentage of total read of sample
  ) %>%
  filter(Perc > 0) %>%        # Filter out rows with 0 counts
  mutate(
    Site = sub("\\.[Ss][0-9]+$", "", Sample),  # remove ".s" replicate suffix
    Site = factor(Site)                        # make Site a factor
  )

table(edna_long$Sample)

str(edna_long)
#Remove Chordates from Phylum
edna_long <- edna_long %>%
  filter(Phylum != "Chordata")

# Add a new column for presence/absence based on the Count column
edna_long <- edna_long %>%
  mutate(Presence_Absence = as.numeric(ifelse(Perc > 0, 1, 0)))  # 1 for presence, 0 for absence

# Convert to long format
morph_long <- morph %>%
  pivot_longer(
    cols = starts_with("X"),  # Select all sample-specific columns
    names_to = "Sample",        # New column for site identifiers
    values_to = "Count"       # New column for count of sample
  )

# Add a new column for presence/absence based on the Count column
morph_long <- morph_long %>%
  mutate(Presence_Absence = as.numeric(ifelse(Count > 0, 1, 0)), # 1 for presence, 0 for absence
         Site = sub("\\_MZB[0-9]+$", "", Sample),
         Site = factor(Site) ) 

morph_long [is.na(morph_long )] <- 0


# View structure of the modified data
str(bulk_long)
head(bulk_long)
summary(morph_long)
str(morph_long)
str(edna_long)

#Prep for community matrix)

# For species accumulation curves we don't need to subset to only insect groups here as we are just interested in the abilities of the methods to grab species richness?

set.seed(123)  # for reproducibility

edna_presence_matrix <- edna_long %>%
  mutate(Species = ifelse(Species == "" | is.na(Species), "Unknown", Species)) %>%
  group_by(Site, Sample, Species) %>%
  summarise(Presence_Absence = max(Presence_Absence), .groups = "drop") %>% 
  group_by(Site) %>%
  slice_sample(n = 1) %>%   # randomly select ONE replicate sample per Site
  ungroup() %>%
  select(Site, Species, Presence_Absence) %>%
  pivot_wider(
    names_from = Species,
    values_from = Presence_Absence,
    values_fill = 0
  )


library(tibble)

edna_presence_matrix_numeric <- edna_presence_matrix %>%
  column_to_rownames("Site") %>%
  as.matrix()

# View the first few rows of the matrix
names(edna_presence_matrix)
rownames(edna_presence_matrix_numeric)

# View the numeric matrix
head(edna_presence_matrix_numeric)


# Run species accumulation curve
edna_accum <- specaccum(edna_presence_matrix_numeric, method = "random")

# Plot the species accumulation curve
plot(
  edna_accum,
  ci.type = "polygon",   # Confidence interval shaded
  col = "blue",          # Curve color
  lwd = 2,               # Line width
  xlab = "Samples",
  ylab = "Species Richness",
  main = "Species Accumulation Curve for eDNA Data"
)

# Prepare bulk invertebrates data

names(bulk_long)

bulk_presence_matrix <- bulk_long %>%
  mutate(Species = ifelse(Species == "" | is.na(Species), "Unknown", Species)) %>%
  group_by(Site, Sample, Species) %>%
  summarise(Presence_Absence = max(Presence_Absence), .groups = "drop") %>% 
  group_by(Site) %>%
  slice_sample(n = 1) %>%   # randomly select ONE replicate sample per Site
  ungroup() %>%
  select(Site, Species, Presence_Absence) %>%
  pivot_wider(
    names_from = Species,
    values_from = Presence_Absence,
    values_fill = 0
  )

bulk_presence_matrix_numeric <- bulk_presence_matrix %>%
  column_to_rownames("Site") %>%
  as.matrix()

# Species accumulation - BULK
bulk_accum <- specaccum(bulk_presence_matrix_numeric, method = "random")

# Morphology data set up
morph_presence_matrix <- morph_long %>%
  mutate(Species = ifelse(Species == "" | is.na(Species), "Unknown", Species)) %>%
  group_by(Site, Sample, Species) %>%
  summarise(Presence_Absence = max(Presence_Absence), .groups = "drop") %>% 
  group_by(Site) %>%
  slice_sample(n = 1) %>%   # randomly select ONE replicate sample per Site
  ungroup() %>%
  select(Site, Species, Presence_Absence) %>%
  pivot_wider(
    names_from = Species,
    values_from = Presence_Absence,
    values_fill = 0
  )

morph_presence_matrix_numeric <- morph_presence_matrix %>%
  column_to_rownames("Site") %>%
  as.matrix()

# Now run specaccum
morph_accum <- specaccum(morph_presence_matrix_numeric, method = "random")

nrow(morph_presence_matrix_numeric)  # should now be 48 (your real sites)

# Plot both curves for comparison
plot(
  edna_accum,
  ci.type = "polygon",
  col = "blue",
  lwd = 2,
  xlab = "Samples",
  ylab = "Species Richness",
  main = "Species Accumulation Curve"
)
lines(bulk_accum, col = "red", lwd = 2)
legend("bottomright", legend = c("eDNA", "Bulk Invertebrates"), col = c("blue", "red"), lwd = 2)
lines(morph_accum, col = "yellow", lwd = 2)
legend("bottomright", legend = c("eDNA", "Bulk Invertebrates", "Morphological ID"), col = c("blue","red", "yellow"), lwd = 2)

str(edna_accum)
str(bulk_accum)
str(morph_accum)

# Extracting richness values
edna_richness <- edna_accum$richness
bulk_richness <- bulk_accum$richness
morph_richness <- morph_accum$richness

# Improved plot
plot(
  edna_richness,
  ci.type = "polygon",
  col = "deepskyblue3",  # Blue color for eDNA
  lwd = 2,
  xlab = "Samples",
  ylab = "Species Count",
  main = "Species Accumulation Curve",
  cex.lab = 1.5,  # Larger axis labels
  cex.main = 2,   # Larger title
  cex.axis = 1.2, # Larger axis numbers
  font.main = 2,  # Bold title
  las = 1,        # Horizontal axis labels
  col.main = "black",
  xlim = c(0, max(c(length(edna_richness), length(bulk_richness), length(morph_richness)))), # Adjust X-axis limit
  ylim = c(0, max(c(max(edna_richness), max(bulk_richness), max(morph_richness)))) # Adjust Y-axis limit
)
# Add bulk invertebrates curve
lines(bulk_richness, col = "firebrick3", lwd = 2, lty = 2)  # Red color, dashed line
# Add morphological ID curve
lines(morph_richness, col = "goldenrod", lwd = 2, lty = 3)  # Yellow color, dotted line
legend(
  "bottomright",
  legend = c("eDNA", "Bulk Invertebrates", "Morphological ID"),
  col = c("deepskyblue3", "firebrick3", "goldenrod"),
  lwd = 2,
  lty = c(1, 2, 3),  # Line types
  cex = 1.2,         # Increase legend size
  bty = "n"          # Remove box around legend
)
# Add grid lines for readability
grid(lwd = 0.5, col = "gray90")


# ----- Function to plot species accumulation with confidence intervals -----
plot_species_accum <- function(accum_obj, curve_col, ci_col, curve_label, lty = 1) {
  x_vals <- 1:length(accum_obj$sites)
  upper_ci <- accum_obj$richness + accum_obj$sd
  lower_ci <- accum_obj$richness - accum_obj$sd
  lower_ci[lower_ci < 0] <- 0  # Avoid negatives
  
  polygon(
    c(x_vals, rev(x_vals)),
    c(upper_ci, rev(lower_ci)),
    col = adjustcolor(ci_col, alpha.f = 0.3),
    border = NA
  )
  
  lines(
    x_vals, accum_obj$richness,
    col = curve_col,
    lwd = 2,
    lty = lty
  )
  
  text(
    x = max(x_vals),
    y = tail(accum_obj$richness, 1),
    labels = curve_label,
    col = curve_col,
    pos = 4,
    cex = 1.1,
    font = 2
  )
}

# ----- Plot all three species accumulation curves with CIs -----
plot(
  1, type = "n",
  xlim = c(1, max(length(edna_accum$sites), length(bulk_accum$sites), length(morph_accum$sites))),
  ylim = c(0, max(
    edna_accum$richness + edna_accum$sd,
    bulk_accum$richness + bulk_accum$sd,
    morph_accum$richness + morph_accum$sd
  )),
  xlab = "Number of Samples",
  ylab = "Species Richness",
  main = "Species Accumulation Curves with 95% Confidence Intervals",
  cex.lab = 1.5,
  cex.main = 1.8,
  cex.axis = 1.2,
  las = 1
)

# Plot each curve
plot_species_accum(edna_accum, curve_col = "deepskyblue3", ci_col = "deepskyblue3", curve_label = NULL)
plot_species_accum(bulk_accum, curve_col = "firebrick3", ci_col = "firebrick3", curve_label = NULL, lty = 2)
plot_species_accum(morph_accum, curve_col = "goldenrod", ci_col = "goldenrod", curve_label = NULL, lty = 3)

# Add legend
legend(
  "bottomright",
  legend = c("eDNA", "Bulk Invertebrates", "Morphological ID"),
  col = c("deepskyblue3", "firebrick3", "goldenrod"),
  lwd = 2,
  lty = c(1, 2, 3),
  cex = 1.2,
  bty = "n"
)

# Add grid
grid(lwd = 0.5, col = "gray85")




#Venn Diagrams


#Randomly sample one replicate for fairer comparison between methods

# eDNA: pick 1 random replicate per Site
edna_subset <- edna_long %>%
  group_by(Site) %>%
  slice_sample(n = 1) %>%
  ungroup()

bulk_subset <- bulk_long %>%
  group_by(Site) %>%
  slice_sample(n = 1) %>%
  ungroup()

morph_subset <- morph_long %>%
  group_by(Site) %>%
  slice_sample(n = 1) %>%
  ungroup()



# Get unique class for each method
edna_Class <- unique(edna_subset$Class)
bulk_Class <- unique(bulk_subset$Class)
morph_Class <- unique(morph_subset$Class)

# Remove any empty or unknown categories
edna_Class  <- edna_Class[edna_Class != "" & !is.na(edna_Class)]
bulk_Class  <- bulk_Class[bulk_Class != "" & !is.na(bulk_Class)]
morph_Class <- morph_Class[morph_Class != "" & !is.na(morph_Class)]


library(ggVennDiagram) #https://gaospecial.github.io/ggVennDiagram/articles/using-ggVennDiagram.html

# Create a list of species by method
class_list <- list(eDNA = edna_Class, Bulk = bulk_Class, Morph=morph_Class)

# Generate Venn diagram
venn <- ggVennDiagram(class_list,  category.names = names(class_list), label="both", edge_size = 0.7) +
  scale_fill_distiller(palette = "OrRd")+  # Color gradient
  theme_minimal(base_size = 14) +                          # Minimal theme
  labs(title = "Venn Diagram of Classes Detected by eDNA, Bulk and Morphological ID Methods") +
  theme(plot.title = element_text(hjust = 0.5, size = 16)) # Centered title

# Display the diagram
print(venn)



# Subset rows where Class is "Insecta"
edna_insecta <- edna_subset %>%
  filter(Class == "Insecta")

bulk_insecta <- bulk_subset %>%
  filter(Class == "Insecta")

morph_insecta <- morph_subset %>%
  filter(Class == "Insecta")



# Get unique species for each method
edna_insect <- unique(edna_insecta$Species)
bulk_insect <- unique(bulk_insecta$Species)
morph_insect <- unique(morph_insecta$Species)


# Remove any empty or unknown categories

# Remove any empty or unknown categories
edna_insect  <- edna_insect[edna_insect != "" & !is.na(edna_insect)]
bulk_insect  <- bulk_insect[bulk_insect != "" & !is.na(bulk_insect)]
morph_insect <- morph_insect[morph_insect != "" & !is.na(morph_insect)]



# Create a list of species by method
species_list <- list(eDNA = edna_insect, Bulk = bulk_insect, Morph=morph_insect)
species_list <- list(
  eDNA = edna_insect[edna_insect != "Unknown"],
  Bulk = bulk_insect[bulk_insect != "Unknown"],
  Morphological = morph_insect[morph_insect != "Unknown"]
)

# Generate Venn diagram
 ggVennDiagram(species_list,  category.names = names(species_list), label="both", edge_size = 0.7) +
  scale_fill_distiller(palette = "BuGn")+  # Color gradient
  theme_minimal(base_size = 14) +                          # Minimal theme
  labs(title = "Venn Diagram of Insect species Detected by eDNA, Bulk and Morphological ID Methods") +
  theme(plot.title = element_text(hjust = 0.5, size = 16)) # Centered title



# Combine both datasets into a single data frame for plotting
combined_insecta <- bind_rows(
  edna_insecta %>% mutate(Source = "eDNA"),
  bulk_insecta %>% mutate(Source = "Bulk"),
  morph_insecta %>% mutate(Source = "Morph")
)

ggplot(combined_insecta, aes(x = Source, y = Presence_Absence, fill = Order)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_minimal(base_size = 14) +
  labs(
    title = "Number of MOTUs/Morphospecies in each Insect Order",
    x = "Family",
    y = "Number of MOTUs or Morphospecies"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  ) +
  scale_fill_brewer(palette = "Paired")  # Use a color palette for families


combined_insecta_proportional <- combined_insecta %>%
  group_by(Source) %>%
  mutate(Proportional_Presence = Presence_Absence / sum(Presence_Absence)) %>%
  ungroup()

ggplot(combined_insecta_proportional, aes(x = Source, y = Proportional_Presence, fill = Order)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_minimal(base_size = 14) +
  labs(
    title = "Proportional Number of MOTUs by Order for Insecta",
    x = "Source",
    y = "Proportional MOTUs"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  ) +
  scale_fill_brewer(palette = "Paired")  # Use a color palette for orders
