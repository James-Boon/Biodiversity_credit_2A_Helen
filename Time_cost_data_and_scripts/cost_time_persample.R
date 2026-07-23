# Cost and time results: Paper 1 (Helen's analysis)

################################################################################
# PACKAGES
################################################################################

library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)
library(gt)
library(readr)
library(grid)


################################################################################
# ANALYSIS NOTE
################################################################################

# Time and cost per sampling event are identical between sites for bird and
# aquatic invertebrate surveys, so these results are not separated by site.


################################################################################
# READ DATA
################################################################################

bird_surveys_costs_time <- read.csv(
  "bird_cost_time_paper2a.csv",
  header = TRUE
)

ai_surveys_costs_time <- read.csv(
  "aquaticinvert_costs_paper1_helen.csv",
  header = TRUE
)

veg_surveys_costs_time <- read.csv(
  "veg_costs_paper1_helen.csv",
  header = TRUE
)


################################################################################
# THEMES AND COLOURS
################################################################################

time_cols <- c(
  "Field time" = "#173a38",
  "Processing time" = "#7fb7ad"
)

cost_cols <- c(
  "Field labour" = "#173a38",
  "Processing labour" = "#7fb7ad",
  "Equipment" = "#502b2d",
  "Additional processing" = "#dcbdac"
)

base_theme <- theme_classic() +
  theme(
    text = element_text(colour = "black"),
    axis.text = element_text(colour = "black"),
    axis.title = element_text(colour = "black"),
    legend.text = element_text(colour = "black"),
    legend.title = element_text(colour = "black"),
    legend.position = "top"
  )

legend_style <- theme(
  legend.key.size = unit(0.35, "cm"),
  legend.spacing.x = unit(0.15, "cm"),
  legend.spacing.y = unit(0.1, "cm"),
  legend.text = element_text(size = 9),
  legend.title = element_text(size = 8),
  legend.margin = margin(0, 0, 0, 0),
  legend.box.margin = margin(0, 0, 0, 0)
)


################################################################################
# BIRD SURVEYS
################################################################################

# Retain one row for each unique combination of method, time, and cost values.
bird_unique <- bird_surveys_costs_time %>%
  select(
    Method,
    Field.time.per.sample,
    Processing.time.per.sample,
    Consumable.cost.per.sample,
    Additional.processing.cost.per.sample
  ) %>%
  distinct()


# ------------------------------------------------------------------------------
# Time per sampling event
# ------------------------------------------------------------------------------

bird_time_long <- bird_unique %>%
  pivot_longer(
    cols = c(
      Field.time.per.sample,
      Processing.time.per.sample
    ),
    names_to = "Time_component",
    values_to = "Hours"
  ) %>%
  mutate(
    Time_component = recode(
      Time_component,
      "Field.time.per.sample" = "Field time",
      "Processing.time.per.sample" = "Processing time"
    ),
    Time_component = factor(
      Time_component,
      levels = c("Field time", "Processing time")
    )
  )

bird_time_plot <- ggplot(
  bird_time_long,
  aes(x = Method, y = Hours, fill = Time_component)
) +
  geom_col(position = position_stack(reverse = TRUE)) +
  coord_flip() +
  scale_fill_manual(values = time_cols) +
  labs(
    x = NULL,
    y = "Time per sampling event (h)",
    fill = NULL
  ) +
  base_theme +
  legend_style +
  theme(
    plot.margin = margin(0, 1, 0, 1)
  )

bird_time_plot


# ------------------------------------------------------------------------------
# Cost per sampling event
# ------------------------------------------------------------------------------

bird_cost_long <- bird_unique %>%
  transmute(
    Method,
    `Field labour` = Field.time.per.sample * 35,
    `Processing labour` = Processing.time.per.sample * 35,
    Equipment = Consumable.cost.per.sample,
    `Additional processing` = Additional.processing.cost.per.sample
  ) %>%
  pivot_longer(
    -Method,
    names_to = "Cost_component",
    values_to = "Cost"
  ) %>%
  mutate(
    Cost_component = factor(
      Cost_component,
      levels = c(
        "Field labour",
        "Processing labour",
        "Equipment",
        "Additional processing"
      )
    )
  )

bird_cost_plot <- ggplot(
  bird_cost_long,
  aes(x = Method, y = Cost, fill = Cost_component)
) +
  geom_col(position = position_stack(reverse = TRUE)) +
  coord_flip() +
  scale_fill_manual(values = cost_cols) +
  labs(
    x = NULL,
    y = "Cost per sampling event (£)",
    fill = NULL
  ) +
  base_theme +
  legend_style +
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    plot.margin = margin(0, 1, 0, 1)
  )

bird_cost_plot


################################################################################
# AQUATIC INVERTEBRATE SURVEYS
################################################################################

# Retain one row for each unique combination of method, time, and cost values.
ai_unique <- ai_surveys_costs_time %>%
  select(
    Survey.method,
    Field.time.per.sample..h.,
    Processing.time.per.sample..h.,
    Consumable.cost.per.sample
  ) %>%
  distinct()


# ------------------------------------------------------------------------------
# Time per sampling event
# ------------------------------------------------------------------------------

ai_time_long <- ai_unique %>%
  transmute(
    Survey.method,
    `Field time` = Field.time.per.sample..h.,
    `Processing time` = Processing.time.per.sample..h.
  ) %>%
  pivot_longer(
    cols = c(`Field time`, `Processing time`),
    names_to = "Time_component",
    values_to = "Hours"
  ) %>%
  mutate(
    Time_component = factor(
      Time_component,
      levels = c("Field time", "Processing time")
    )
  )

ai_time_plot <- ggplot(
  ai_time_long,
  aes(x = Survey.method, y = Hours, fill = Time_component)
) +
  geom_col(position = position_stack(reverse = TRUE)) +
  coord_flip() +
  scale_fill_manual(values = time_cols) +
  labs(
    x = NULL,
    y = "Time per sampling event (h)",
    fill = NULL
  ) +
  base_theme +
  legend_style +
  theme(
    plot.margin = margin(0, 1, 0, 1)
  )


# ------------------------------------------------------------------------------
# Cost per sampling event
# ------------------------------------------------------------------------------

ai_cost_long <- ai_unique %>%
  transmute(
    Survey.method,
    `Field labour` = Field.time.per.sample..h. * 35,
    `Processing labour` = Processing.time.per.sample..h. * 35,
    Equipment = Consumable.cost.per.sample
  ) %>%
  pivot_longer(
    -Survey.method,
    names_to = "Cost_component",
    values_to = "Cost"
  ) %>%
  mutate(
    Cost_component = factor(
      Cost_component,
      levels = c(
        "Field labour",
        "Processing labour",
        "Equipment"
      )
    )
  )

ai_cost_plot <- ggplot(
  ai_cost_long,
  aes(x = Survey.method, y = Cost, fill = Cost_component)
) +
  geom_col(position = position_stack(reverse = TRUE)) +
  coord_flip() +
  scale_fill_manual(
    values = cost_cols[c(
      "Field labour",
      "Processing labour",
      "Equipment"
    )]
  ) +
  labs(
    x = NULL,
    y = "Cost per sampling event (£)",
    fill = NULL
  ) +
  base_theme +
  legend_style +
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    plot.margin = margin(0, 1, 0, 1)
  )


################################################################################
# VEGETATION AND HABITAT-CONDITION SURVEYS
################################################################################

# ------------------------------------------------------------------------------
# Total survey time
# ------------------------------------------------------------------------------

# Convert field and processing times to long format.
veg_time_long <- veg_surveys_costs_time %>%
  transmute(
    Site,
    Survey.method,
    Field.time,
    Processing.time
  ) %>%
  pivot_longer(
    cols = c(Field.time, Processing.time),
    names_to = "Time_component",
    values_to = "Hours"
  ) %>%
  mutate(
    Time_component = recode(
      Time_component,
      "Field.time" = "Field time",
      "Processing.time" = "Processing time"
    ),
    Time_component = factor(
      Time_component,
      levels = c("Field time", "Processing time")
    )
  )

# Remove the Attenborough UAV survey because it was not used in the analysis.
veg_time_long <- veg_time_long %>%
  filter(!(Site == "Attenborough" & Survey.method == "UAV"))

# Calculate a common upper limit with 5% headroom.
veg_time_common_cost_max <- veg_time_long %>%
  group_by(Site, Survey.method) %>%
  summarise(total = sum(Hours, na.rm = TRUE), .groups = "drop") %>%
  summarise(mx = max(total, na.rm = TRUE) * 1.05) %>%
  pull(mx)

# Prepare a combined method-and-site label for the plot axis.
veg_time_combined <- veg_time_long %>%
  filter(!(Site == "Attenborough" & Survey.method == "UAV")) %>%
  mutate(
    # Remove trailing spaces from survey-method names.
    Survey.method = stringr::str_trim(Survey.method),
    
    # Combine survey method and site in a single axis label.
    Axis_Label = paste0(Survey.method, "\n(", Site, ")"),
    
    # Set the plotted order from bottom to top after coord_flip().
    Axis_Label = factor(
      Axis_Label,
      levels = c(
        "DEFRA metric\n(Knepp)",
        "Satellite\n(Knepp)",
        "UAV\n(Knepp)",
        "DEFRA metric\n(Attenborough)",
        "Satellite\n(Attenborough)"
      )
    )
  )

# Create the stacked time plot.
final_plot_time_veg <- ggplot(
  veg_time_combined,
  aes(x = Axis_Label, y = Hours, fill = Time_component)
) +
  geom_col(position = position_stack(reverse = TRUE)) +
  coord_flip() +
  scale_fill_manual(
    values = time_cols,
    breaks = c("Field time", "Processing time")
  ) +
  labs(
    x = NULL,
    y = "Total time (h)",
    fill = NULL
  ) +
  base_theme +
  legend_style


# ------------------------------------------------------------------------------
# Total survey cost
# ------------------------------------------------------------------------------

# Convert cost components to long format.
veg_cost_long <- veg_surveys_costs_time %>%
  transmute(
    Site,
    Survey.method,
    Equipment = Consumable.cost,
    `Field labour` = Total.field.labour.cost,
    `Additional processing` = Data.access.cost,
    `Processing labour` = Data.processing.cost,
  ) %>%
  pivot_longer(
    -c(Site, Survey.method,),
    names_to = "Cost_component",
    values_to = "Cost"
  ) %>%
  mutate(
    Cost_component = factor(
      Cost_component,
      levels = c(
        "Field labour",
        "Processing labour",
        "Additional processing",
        "Equipment"
      )
    )
  )

# Remove the Attenborough UAV survey because it was not used in the analysis.
veg_cost_long <- veg_cost_long %>%
  filter(!(Site == "Attenborough" & Survey.method == "UAV"))

# Calculate a common upper limit with 5% headroom.
veg_common_cost_max <- veg_cost_long %>%
  group_by(Site, Survey.method) %>%
  summarise(total = sum(Cost, na.rm = TRUE), .groups = "drop") %>%
  summarise(mx = max(total, na.rm = TRUE) * 1.05) %>%
  pull(mx)

# Prepare the cost data and combined method-and-site axis labels.
veg_cost_combined <- veg_surveys_costs_time %>%
  transmute(
    Site,
    Survey.method,
    Equipment = Consumable.cost,
    `Field labour` = Total.field.labour.cost,
    `Additional processing` = Data.access.cost,
    `Processing labour` = Data.processing.cost,
  ) %>%
  pivot_longer(
    -c(Site, Survey.method),
    names_to = "Cost_component",
    values_to = "Cost"
  ) %>%
  mutate(
    Cost_component = factor(
      Cost_component,
      levels = c(
        "Field labour",
        "Processing labour",
        "Additional processing",
        "Equipment"
      )
    )
  ) %>%
  
  # Remove the Attenborough UAV survey because it was not used in the analysis.
  filter(!(Site == "Attenborough" & Survey.method == "UAV")) %>%
  mutate(
    # Remove trailing spaces from survey-method names.
    Survey.method = stringr::str_trim(Survey.method),
    
    # Combine survey method and site in a single axis label.
    Axis_Label = paste0(Survey.method, "\n(", Site, ")"),
    
    # Match the plotted order used for the time plot.
    Axis_Label = factor(
      Axis_Label,
      levels = c(
        "DEFRA metric\n(Knepp)",
        "Satellite\n(Knepp)",
        "UAV\n(Knepp)",
        "DEFRA metric\n(Attenborough)",
        "Satellite\n(Attenborough)"
      )
    )
  )

# Create the stacked cost plot.
veg_final_cost_plot <- ggplot(
  veg_cost_combined,
  aes(x = Axis_Label, y = Cost, fill = Cost_component)
) +
  geom_col(position = position_stack(reverse = TRUE)) +
  coord_flip() +
  scale_fill_manual(values = cost_cols) +
  labs(
    x = NULL,
    y = "Total cost (£)",
    fill = NULL
  ) +
  base_theme +
  legend_style +
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    plot.margin = margin(0, 1, 0, 1)
  )



################################################################################
# MULTIPANEL FIGURE 
################################################################################
# look at all plots

bird_time_plot
bird_cost_plot

ai_time_plot
ai_cost_plot

final_plot_time_veg
veg_final_cost_plot

current_rel_widths <- c(1.05, 1) 

################################################################################
# EXTRACT LEGENDS
################################################################################

legend_time <- get_legend(
  bird_time_plot +
    theme(
      legend.position = "top",
      legend.direction = "horizontal"
    )
)

legend_cost <- get_legend(
  bird_cost_plot +
    theme(
      legend.position = "top",
      legend.direction = "horizontal"
    )
)

################################################################################
# REMOVE LEGENDS, RESET MARGINS, & FIX INNER AXIS CLIPPING
################################################################################

# Left plots (Time): Clean legends/margins, and expand the right side slightly (0.08)
bird_time_clean <- bird_time_plot +
  theme(legend.position = "none", plot.margin = margin()) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.08)))

ai_time_clean <- ai_time_plot +
  theme(legend.position = "none", plot.margin = margin()) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.08)))

veg_time_clean <- final_plot_time_veg +
  theme(legend.position = "none", plot.margin = margin()) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.08)))


# Right plots (Cost): Clean legends/margins, and expand the left side slightly (0.08)
bird_cost_clean <- bird_cost_plot +
  theme(legend.position = "none", plot.margin = margin()) +
  scale_y_continuous(expand = expansion(mult = c(0.08, 0.02)))

ai_cost_clean <- ai_cost_plot +
  theme(legend.position = "none", plot.margin = margin()) +
  scale_y_continuous(expand = expansion(mult = c(0.08, 0.02)))

veg_cost_clean <- veg_final_cost_plot +
  theme(legend.position = "none", plot.margin = margin()) +
  scale_y_continuous(expand = expansion(mult = c(0.08, 0.02)))
################################################################################
# ALIGN TIME PANELS
################################################################################

time_plots <- align_plots(
  bird_time_clean,
  ai_time_clean,
  veg_time_clean,
  align = "v",
  axis = "l"
)

bird_time_aligned <- time_plots[[1]]
ai_time_aligned   <- time_plots[[2]]
veg_time_aligned  <- time_plots[[3]]

################################################################################
# ALIGN COST PANELS
################################################################################

cost_plots <- align_plots(
  bird_cost_clean,
  ai_cost_clean,
  veg_cost_clean,
  align = "v",
  axis = "l"
)

bird_cost_aligned <- cost_plots[[1]]
ai_cost_aligned   <- cost_plots[[2]]
veg_cost_aligned  <- cost_plots[[3]]

################################################################################
# BUILD ROWS (Using your tight rel_widths)
################################################################################

row_a <- plot_grid(
  bird_time_aligned,
  bird_cost_aligned,
  nrow = 1,
  rel_widths = current_rel_widths
)

row_b <- plot_grid(
  ai_time_aligned,
  ai_cost_aligned,
  nrow = 1,
  rel_widths = current_rel_widths
)

row_c <- plot_grid(
  veg_time_aligned,
  veg_cost_aligned,
  nrow = 1,
  rel_widths = current_rel_widths
)

################################################################################
# STACK ROWS
################################################################################

panel_plots <- plot_grid(
  row_a,
  row_b,
  row_c,
  ncol = 1,
  align = "v",
  axis = "l",
  rel_heights = c(1, 1, 1)
)

################################################################################
# LEGEND ROW
################################################################################

panel_legends <- plot_grid(
  legend_time,
  legend_cost,
  nrow = 1,
  rel_widths = current_rel_widths
)

################################################################################
# COMBINE
################################################################################

panel_base <- plot_grid(
  panel_legends,
  panel_plots,
  ncol = 1,
  rel_heights = c(0.08, 1)
)

################################################################################
# PANEL LABELS
################################################################################

panel_final <- ggdraw(panel_base) +
  draw_label(
    "(a)",
    x = 0.01,
    y = 0.92,
    hjust = 0,
    vjust = 1,
    fontface = "bold",
    size = 12
  ) +
  draw_label(
    "(b)",
    x = 0.01,
    y = 0.64,
    hjust = 0,
    vjust = 1,
    fontface = "bold",
    size = 12
  ) +
  draw_label(
    "(c)",
    x = 0.01,
    y = 0.36,
    hjust = 0,
    vjust = 1,
    fontface = "bold",
    size = 12
  )

panel_final


# final multipanel plot for birds!
panel_final

w <- dev.size("in")[1]
h <- dev.size("in")[2] 

# 1. Add a solid white background theme to the final panel object
panel_final_fixed <- panel_final + 
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )

ggsave(
  "costtimefigure.tiff",
  plot = panel_final_fixed,
  width = w,
  height = h,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"        # <-- Explicitly forces transparency to render as white
)
