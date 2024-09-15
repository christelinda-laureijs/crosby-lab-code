# Theme Options


treatment_names_and_colours <- data.frame(
  treatment = c(
    "WT",
    "Cre"
  ),
  display_names = c(
    "WT",
    "Cre"
  ),
  colours = c(
    "#6600cc",
    "#e86c00"
  ),
  very_pale_colours = c(
    "#d6b8f5",
    "#ffc38f"
  )
)


# treatment_names_and_colours <- data.frame(
#   treatment = c(
#     "Control",
#     "HNMPA",
#     "PPP",
#     "PPP_and_HNMPA"
#   ),
#   display_names = c(
#     "Control",
#     "HNMPA",
#     "PPP",
#     "PPP\n&\nHNMPA"
#   ),
#   colours = c(
#     "#6600cc",
#     "#e86c00",
#     "#0093fb",
#     "#411900"
#   ),
#   very_pale_colours = c(
#     "#d6b8f5",
#     "#ffc38f",
#     "#8fd0ff",
#     "#bf9b84"
#   )
# )


# Set colours here only for consistency

gray_shading_colour <- "#dcdcdc"
line_col <- "#333333" # Sets colour of x and y-axes

# Rectangle highlights intervals from 5-10 min and 15-20 min
rectangle_shading_colour <- "#f6f6f6"

# Custom fonts may cause issues depending on what fonts you have in your system
# If it does not work, you could always delete 'family = plot_font_family' in the ggplot theme set below
plot_font_family <- "Segoe UI"
plot_light_font_family <- "Segoe UI Light"
significance_stars_font <- plot_font_family


available_sEPSC_parameters <- data.frame(available_parameters = c("amplitude", "raw_amplitude", "frequency", "raw_frequency"))

sEPSC_parameters_and_treatments <- expand.grid(
  treatment_names_and_colours$treatment,
  available_sEPSC_parameters$available_parameters
) %>%
  rename(treatment = Var1, parameter = Var2)

sEPSC_parameters_and_treatments

# Shapes
male_shape <- 16
female_shape <- 17

list_of_significance_stars <- c(
  "***" = 0.001,
  "**" = 0.01,
  "*" = 0.05,
  "ns" = 2
)

# Set default colours and point sizes for action potential summary plots
mean_point_colour <- "#000000"
baseline_group_colour <- "#727a85"
insulin_group_colour <- "#6600cc"
connecting_line_width <- 0.1
connecting_line_colour_aps <- "#bfbfbf"
connecting_line_width_PPR <- 0.2

mean_point_size <- 0.5
geom_sina_size <- 3
geom_signif_text_size <- 8

AP_trace_size <- 0.7
scale_bar_shift_y <- 5
scale_bar_shift_x <- 30

# These sizes work better for the multi-plot figure in the PDF output
if (knitr::is_latex_output()) {
  geom_sina_size <- 2
  geom_signif_text_size <- 4
  AP_trace_size <- 0.6
  scale_bar_shift_y <- 13
  scale_bar_shift_x <- 200
  connecting_line_width <- 0.04
  connecting_line_width_PPR <- 0.04
}

# A consistent y-axis enables comparison across multiple experiments, treatments, etc.
# This is not applied to the raw plots of eEPSCs vs. time for individual cells
# If you change this from 175, make sure to regenerate all summary plots so that you can compare across the same y-axis.

y_axis_limit <- 175
PPR_y_min <- 0
PPR_y_max <- 5

# Y-axis for spontaneous currents
custom_y_axis_spont_amplitude <- c(70, 120)
custom_y_axis_spont_frequency <- c(0, 25)


# Requires the extrafont() package (loaded in the load-libraries chunk) for custom font

modified_theme_classic <- theme_classic() +
  theme(
    text = element_text(family = plot_font_family),
    plot.title = element_text(
      color = "black",
      size = 20,
      family = plot_font_family,
      # face = "plain",
      margin = margin(b = 25),
      hjust = 0.5
    ),
    plot.margin = margin(25, 25, 25, 25),
    plot.caption = element_text(
      hjust = 0,
      family = plot_font_family,
      size = 12
    ),
    plot.caption.position = "plot",
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(
      margin = margin(r = 25),
      angle = 90,
      vjust = 0.5
    ),
    axis.title.x = element_text(margin = margin(b = 25, t = 20)),
    axis.ticks = element_blank(),
    strip.background = element_rect(color = NA, fill = NA),
    strip.text = element_text(size = 20)
  )


theme_set(modified_theme_classic)


modified_facet_theme <- modified_theme_classic +
  theme(
    legend.position = "none",
    axis.line.x = element_line(color = "gray"),
    axis.line.y = element_line(color = "gray"),
    panel.spacing = unit(2, "lines"),
    plot.title = element_text(family = plot_light_font_family, size = 40),
    plot.subtitle = element_text(
      size = 20,
      hjust = 0.5,
      margin = margin(b = 50)
    )
  )