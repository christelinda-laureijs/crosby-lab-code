# Theme Options

# Set colours here only for consistency

line_col <- "#333333" # Sets colour of x and y-axes
my_colours <- c("#6600cc", "#0093fb", "#55b323","#ffe70f","#e86c00","#333333", "#411900","#e11584")
my_colours_pale <- c("#b080e0","#5cb9fa", "#92d46e","#ebdf05","#f58c31","#b1b1b1","#a16b4a","#f57dbe")

# Rectangle highlights intervals from 5-10 min and 15-20 min
rectangle_shading_colour <- "#f6f6f6"

# Required for better contrast between the sexes in summary plots
# Do NOT use these for any of the raw P1 vs. Time plots because they are too pale and/or too dark
my_colours_very_dark <- c("#4d0299","#026bb5","#398511","#a69502","#994700","#000000","#785138","#910150")
my_colours_very_pale <- c("#d6b8f5", "#8fd0ff", "#aee691","#ebdf05","#ffc38f","#dcdcdc","#bf9b84","#eba2c9")


# Custom fonts may cause issues depending on what fonts you have in your system
# Troubleshooting steps: 
# Try changing the font to one that you have on your computer
# If it does not work, you could always delete 'family = plot_font_family' in the ggplot theme set below
plot_font_family <- "Segoe UI"
plot_light_font_family <- "Segoe UI Light"
significance_stars_font <- plot_font_family

treatment_shapes <- c(16, 17, 15, 18)
list_of_treatments <- c("Control", "Fasting", "HNMPA", "PPP", "PPP_and_HNMPA")
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

# ----------- Set ggplot theme ----------------------------------------------------------------------------------
# Formatting changes like increasing font size & removing gray background
# Requires the extrafont() package (loaded in the load-libraries chunk) for custom font

theme_set(
  theme_classic() %+replace%
    theme(
      text = element_text(
        family = plot_font_family
      ),
      plot.title = element_text(
        color = "black",
        size = 20,
        family = plot_font_family,
        #face = "plain",
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
      axis.text = element_text(
        size = 12,
        color = "black"
      ),
      axis.title = element_text(
        size = 16,
        face = "bold"
      ),
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
)