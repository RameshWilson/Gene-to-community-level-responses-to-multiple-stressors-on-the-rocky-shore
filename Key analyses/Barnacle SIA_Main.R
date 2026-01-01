########## 0) Configuration ##########

# install.packages("here") # Install once if needed
library(here)

# CSV
SIA_CSV <- here("Dataframes", "SIAdata.csv")

# Treatment colours
colors <- c(
  "Non-polluted + Ambient" = "grey75",
  "Non-polluted + Warm"    = "orange",
  "Polluted + Ambient"     = "darkgreen",
  "Polluted + Warm"        = "#CC79A7"
)

########## 1) Packages ##########

suppressPackageStartupMessages({
  library(tidyverse)  
  library(patchwork)
  library(cowplot)
  library(ggthemes)
  library(DT) # For HTML viewing      
})

########## 2) Load data ##########

sia_data <- readr::read_csv(SIA_CSV, show_col_types = FALSE)
head(sia_data)
# Expect columns: `15N (Sam_DeltaAir)`, `13C (Sam_DeltaPDB)`, `C (Sam_ug)`, `N (Sam_ug)`

########## 3) Cleaning and factor assignment ##########

# Assign Group from sample Name
sia_clean <- sia_data %>%
  filter(!is.na(`15N (Sam_DeltaAir)`), !is.na(`13C (Sam_DeltaPDB)`)) %>%
  mutate(Group = case_when(
    grepl("SC", Name)       ~ "Standards (Seal Collagen)",
    grepl("S-BP", Name)     ~ "Polluted + Warm",
    grepl("S-WP", Name)     ~ "Polluted + Ambient",
    grepl("S-BNP", Name)    ~ "Non-polluted + Warm",
    grepl("S-WNP", Name)    ~ "Non-polluted + Ambient",
    grepl("Alanine", Name)  ~ "Standards (Alanine)",
    grepl("ch-6", Name, ignore.case = TRUE) | grepl("CH6", Name) ~ "Standards (IAEA-CH-6)",
    TRUE ~ NA_character_
  ))

# Include only specified groups
groups_to_include <- c(
  "Standards (Seal Collagen)", "Standards (Alanine)", "Standards (IAEA-CH-6)",
  "Polluted + Warm", "Polluted + Ambient", "Non-polluted + Warm", "Non-polluted + Ambient"
)

sia_filtered <- sia_clean %>%
  filter(Group %in% groups_to_include)

# Derive Pollution and Warming factors from Group label
sia_filtered <- sia_filtered %>%
  mutate(
    Pollution = case_when(
      grepl("Polluted", Group)      ~ "Polluted",
      grepl("Non-polluted", Group)  ~ "Non-Polluted",
      TRUE                          ~ NA_character_
    ),
    Warming = case_when(
      grepl("Warm", Group)    ~ "Warmed",
      grepl("Ambient", Group) ~ "Ambient",
      TRUE                    ~ NA_character_
    ),
    Pollution = factor(Pollution, levels = c("Non-Polluted", "Polluted")),
    Warming   = factor(Warming,   levels = c("Ambient", "Warmed")),
    CN_ratio  = `C (Sam_ug)` / `N (Sam_ug)`
  )

print(sia_filtered, n=50) # check full dataset

########## 4) Descriptive statistics ##########

summary_stats <- sia_filtered %>%
  group_by(Group) %>%
  summarise(
    Mean_δ15N = mean(`15N (Sam_DeltaAir)`, na.rm = TRUE),
    SD_δ15N   = sd(`15N (Sam_DeltaAir)`,   na.rm = TRUE),
    Min_δ15N  = min(`15N (Sam_DeltaAir)`,  na.rm = TRUE),
    Max_δ15N  = max(`15N (Sam_DeltaAir)`,  na.rm = TRUE),
    
    Mean_δ13C = mean(`13C (Sam_DeltaPDB)`, na.rm = TRUE),
    SD_δ13C   = sd(`13C (Sam_DeltaPDB)`,   na.rm = TRUE),
    Min_δ13C  = min(`13C (Sam_DeltaPDB)`,  na.rm = TRUE),
    Max_δ13C  = max(`13C (Sam_DeltaPDB)`,  na.rm = TRUE),
    
    Mean_CN   = mean(CN_ratio, na.rm = TRUE),
    SD_CN     = sd(CN_ratio,   na.rm = TRUE),
    .groups   = "drop"
  )

print(summary_stats)

########## 5) Two-way ANOVAs (δ15N, δ13C, C:N) ##########

anova_δ15N <- aov(`15N (Sam_DeltaAir)` ~ Pollution * Warming, data = sia_filtered)
anova_δ13C <- aov(`13C (Sam_DeltaPDB)` ~ Pollution * Warming, data = sia_filtered)
anova_CN   <- aov(CN_ratio ~ Pollution * Warming, data = sia_filtered)

print(summary(anova_δ15N))
print(summary(anova_δ13C))
print(summary(anova_CN))

########## 6) Main plot: boxplots ##########

# Match GLM ordering
LEGEND_ORDER <- c(
  "Polluted + Ambient",
  "Polluted + Warm",
  "Non-polluted + Ambient",
  "Non-polluted + Warm"
)

# Helper to make darker outline for jitter points
darken_hex <- function(cols, factor = 0.7) {
  out <- vapply(cols, function(col) {
    rgb_mat <- grDevices::col2rgb(col)
    rgb_new <- pmax(0, pmin(255, round(rgb_mat * factor)))
    grDevices::rgb(rgb_new[1], rgb_new[2], rgb_new[3], maxColorValue = 255)
  }, character(1))
  names(out) <- names(cols)
  out
}
stroke_colors <- darken_hex(colors, 0.7)

# Filter out standards; fix factor order
sia_filtered_no_standards <- sia_filtered %>%
  filter(!is.na(Group)) %>%
  filter(!Group %in% c("Standards (Seal Collagen)", "Standards (Alanine)", "Standards (IAEA-CH-6)")) %>%
  mutate(Group = factor(Group, levels = LEGEND_ORDER))

create_boxplot <- function(data, y_variable, y_label) {
  ggplot(data, aes(x = Group, y = .data[[y_variable]], fill = Group)) +
    geom_boxplot(alpha = 0.6, outlier.shape = NA) +
    geom_jitter(
      aes(color = Group),
      shape = 21, width = 0.1, size = 1.9, stroke = 0.6,
      alpha = 0.9, show.legend = FALSE
    ) +
    labs(x = "Treatment", y = y_label, fill = "Treatment") +
    scale_x_discrete(limits = LEGEND_ORDER) +
    scale_fill_manual(values = colors, limits = LEGEND_ORDER) +
    scale_color_manual(values = stroke_colors, limits = LEGEND_ORDER, guide = "none") +
    theme_minimal() +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      axis.text.x  = element_text(size = 8, angle = 45, hjust = 1),
      axis.text.y  = element_text(size = 8),
      axis.title.x = element_text(size = 10),
      axis.title.y = element_text(size = 10),
      plot.margin  = margin(5, 5, 5, 5)
    )
}

boxplot_δ15N <- create_boxplot(sia_filtered_no_standards, "15N (Sam_DeltaAir)", "δ15N value (‰)")

# Manual ticks for δ13C
boxplot_δ13C <- create_boxplot(sia_filtered_no_standards, "13C (Sam_DeltaPDB)", "δ13C value (‰)") +
  scale_y_continuous(breaks = c(-18.5, -18, -17.5, -17, -16.5))

boxplot_CN   <- create_boxplot(sia_filtered_no_standards, "CN_ratio", "C:N mass ratio")
gap_rel <- 0.05

# Panel plot
SIA_boxplots <- (
  boxplot_δ15N | plot_spacer() | boxplot_δ13C | plot_spacer() | boxplot_CN
) +
  plot_layout(widths = c(1, gap_rel, 1, gap_rel, 1), guides = "collect") &
  theme(
    legend.position       = "right",
    legend.justification  = "center",
    legend.text           = element_text(size = 10),
    legend.title          = element_text(face = "bold", size = 11),
    legend.key.size       = grid::unit(10, "mm"),
    legend.key.width      = grid::unit(10, "mm"),
    legend.spacing.y      = grid::unit(4, "mm"),
    legend.background     = element_blank(),
    legend.box.background = element_blank()
  ) &
  guides(fill = guide_legend(ncol = 1))   # single-column legend

print(SIA_boxplots)

########## 7) Additional plots ##########

# Scatter plot with elipses
  ELLIPSE_FILL_LVL    <- 0.68  # 68% core
  ELLIPSE_OUTLINE_LVL <- 0.95  # 95% outline
  
  # Compute 2D centroids
  centroids <- sia_filtered_no_standards %>%
    dplyr::group_by(Group) %>%
    dplyr::summarise(
      cx = mean(`13C (Sam_DeltaPDB)`, na.rm = TRUE),
      cy = mean(`15N (Sam_DeltaAir)`, na.rm = TRUE),
      .groups = "drop"
    )
  
  SIA_scatterplot <- ggplot(
    sia_filtered_no_standards,
    aes(x = `13C (Sam_DeltaPDB)`, y = `15N (Sam_DeltaAir)`)
  ) +
    # Filled translucent 68% ellipses
    stat_ellipse(
      aes(fill = Group, group = Group),
      type = "norm", level = ELLIPSE_FILL_LVL,
      alpha = 0.12, geom = "polygon", color = NA
    ) +
    # Sample points
    geom_point(aes(color = Group), size = 3, alpha = 1) +
    # Dashed outline 95% ellipses 
    stat_ellipse(
      aes(color = Group, group = Group),
      type = "norm", level = ELLIPSE_OUTLINE_LVL,
      linetype = 5, geom = "path", linewidth = 0.3, show.legend = FALSE
    ) +
    # Group centroids
    geom_point(
      data = centroids,
      aes(x = cx, y = cy, fill = Group),
      shape = 21, size = 4.0, stroke = 2.2, color = "black", show.legend = FALSE
    ) +
    labs(
      title = NULL,
      x = "δ13C Value (‰)",
      y = "δ15N Value (‰)",
      color = "Treatment",
      fill  = "Treatment"
    ) +
    scale_x_reverse() +
    theme_light() +
    scale_color_manual(values = colors) +  
    scale_fill_manual(values = colors) +  
    theme(
      axis.text.x = element_text(size = 8),
      axis.text.y = element_text(size = 8),
      axis.title.x = element_text(size = 10),
      axis.title.y = element_text(size = 10),
      legend.title = element_text(size = 11, face = "bold"),
      legend.text  = element_text(size = 10),
      legend.background = element_rect(color = NA, fill = "white"),
      plot.margin = margin(10, 10, 10, 10)
    )
  
  print(SIA_scatterplot)
  
# Density plot
  # Long format
  sia_long <- sia_filtered_no_standards %>%
    tidyr::pivot_longer(
      cols = c(`15N (Sam_DeltaAir)`, `13C (Sam_DeltaPDB)`),
      names_to = "Isotope",
      values_to = "Value"
    ) %>%
    dplyr::mutate(
      Isotope = dplyr::case_when(
        Isotope == "15N (Sam_DeltaAir)" ~ "δ15N",
        Isotope == "13C (Sam_DeltaPDB)" ~ "δ13C",
        TRUE ~ Isotope
      )
    )
  
  # Split by isotope to reverse δ13C cleanly
  sia_c13 <- dplyr::filter(sia_long, Isotope == "δ13C")
  sia_n15 <- dplyr::filter(sia_long, Isotope == "δ15N")
  
  # Density
  p_den_c13 <- ggplot(sia_c13, aes(x = Value, fill = Group)) +
    geom_density(alpha = 0.5) +
    scale_x_reverse() +
    labs(x = "Content (‰)", y = "Density", fill = "Treatment") +
    scale_fill_manual(values = colors) +
    theme_minimal() +
    theme(
      plot.title        = element_text(size = 14, face = "bold"),
      legend.title      = element_text(size = 9, face = "bold"),
      legend.background = element_rect(color = NA, fill = "white"),
      legend.text       = element_text(size = 10),
      axis.title.x      = element_text(margin = margin(t = 10))
    ) +
    ggtitle("δ13C")
  
  p_den_n15 <- ggplot(sia_n15, aes(x = Value, fill = Group)) +
    geom_density(alpha = 0.5) +
    labs(x = "Content (‰)", y = "Density", fill = "Treatment") +
    scale_fill_manual(values = colors) +
    theme_minimal() +
    theme(
      plot.title        = element_text(size = 14, face = "bold"),
      legend.title      = element_text(size = 9, face = "bold"),
      legend.background = element_rect(color = NA, fill = "white"),
      legend.text       = element_text(size = 10),
      axis.title.x      = element_text(margin = margin(t = 10))
    ) +
    ggtitle("δ15N")
  
  density_plot <- (p_den_c13 | p_den_n15) + plot_layout(guides = "collect")
  density_plot <- density_plot & theme(legend.position = "right")
  density_plot <- density_plot + plot_annotation(
    title = "Density plots of δ15N and δ13C by treatment group",
    theme = theme(plot.title = element_text(size = 14, face = "bold"))
  )
  print(density_plot)

# Histograms
  p_hist_c13 <- ggplot(sia_c13, aes(x = Value, fill = Group)) +
    geom_histogram(bins = 30, alpha = 0.7, position = "identity", color = "black") +
    scale_x_reverse() +
    labs(x = "Content (‰)", y = "Count", fill = "Treatment") +
    scale_fill_manual(values = colors) +
    theme_minimal() +
    theme(
      plot.title        = element_text(size = 14, face = "bold"),
      legend.title      = element_text(size = 9, face = "bold"),
      legend.background = element_rect(color = NA, fill = "white"),
      legend.text       = element_text(size = 10),
      axis.title.x      = element_text(margin = margin(t = 10))
    ) +
    ggtitle("δ13C")
  
  p_hist_n15 <- ggplot(sia_n15, aes(x = Value, fill = Group)) +
    geom_histogram(bins = 30, alpha = 0.7, position = "identity", color = "black") +
    labs(x = "Content (‰)", y = "Count", fill = "Treatment") +
    scale_fill_manual(values = colors) +
    theme_minimal() +
    theme(
      plot.title        = element_text(size = 14, face = "bold"),
      legend.title      = element_text(size = 9, face = "bold"),
      legend.background = element_rect(color = NA, fill = "white"),
      legend.text       = element_text(size = 10),
      axis.title.x      = element_text(margin = margin(t = 10))
    ) +
    ggtitle("δ15N")
  
  histogram_plot <- (p_hist_c13 | p_hist_n15) + plot_layout(guides = "collect")
  histogram_plot <- histogram_plot & theme(legend.position = "right")
  histogram_plot <- histogram_plot + plot_annotation(
    title = "Histograms of δ15N and δ13C by treatment group",
    theme = theme(plot.title = element_text(size = 14, face = "bold"))
  )
  print(histogram_plot)