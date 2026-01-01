# SI figure script: Peacehaven/Newhaven precipitation vs sewage releases
#   This script generates a visual (non-inferential) figure for Supplementary
#   Information, overlaying:
#     (1) daily precipitation totals (from Visual Crossing; user-downloaded), and
#     (2) sewage release event durations (public-access dataset; included in repo).
#
#   The raw precipitation CSV sourced from Visual Crossing is not stored or
#   distributed in this GitHub repository. Visual Crossing’s license terms place
#   restrictions on publishing/making their service components and/or returned
#   data available for others to copy, and on storage beyond what a license allows.
#
#   To reproduce this figure, one must obtain their own VC precipitation export
#   if desired (via your own account/subscription), and place it locally at:
#     Dataframes/Precipitation data.csv
#
#   Expected minimum columns in the VC precipitation CSV:
#     - datetime  (date or datetime; will be coerced to Date)
#     - precip    (daily total; units specified below)
#
#   VC license terms:
#     https://www.visualcrossing.com/weather-services-terms/
#
#   Because the VC precipitation file is not included in our repo, the below script
#   is provided for illustrative means to demonstrate the method for SI figure generation.

# ==============================================================================

########## 0) Configuration ##########

# install.packages("here") # Install once if needed
library(here)

# File paths (relative to the RStudio project root)
PRECIP_CSV  <- here("Dataframes", "Precipitation data.csv")
RELEASE_CSV <- here("Dataframes", "Newhaven releases.csv")

# Units in precip file: "in" (inches) or "mm"
PRECIP_UNITS <- "in"

# X-axis ticks
selected_dates <- as.Date(c(
  "2023-06-01", "2023-06-20", "2023-07-06", "2023-07-20",
  "2023-08-01", "2023-08-17", "2023-09-01", "2023-09-15"
))
custom_labels <- c("01/06/2023", "June 1", "July 1", "July 2",
                   "August 1", "August 2", "September 1", "September 2")

# Packages
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(lubridate)
  library(readr)
  library(stringr)
})

########## 1) Load & prepare precipitation data ##########

precip_data <- read_csv(PRECIP_CSV, show_col_types = FALSE)

precip_data <- precip_data %>%
  mutate(
    datetime = as.Date(datetime),
    precip   = if (tolower(PRECIP_UNITS) == "in") precip * 25.4 else precip  # inches -> mm
  )
head(precip_data)

########## 2) Load & prepare sewage release data ##########

releases_raw <- read_csv(RELEASE_CSV, show_col_types = FALSE)
head(releases_raw)

# Standardise column names so spaces become dots
names(releases_raw) <- make.names(names(releases_raw))

# Contingency parser if durations are textual
parse_duration_to_minutes <- function(x) {
  x2 <- tolower(trimws(x))
  hrs  <- as.numeric(str_match(x2, "(\\d+)\\s*hour")[, 2])
  mins <- as.numeric(str_match(x2, "(\\d+)\\s*minute")[, 2])
  hrs[is.na(hrs)]   <- 0
  mins[is.na(mins)] <- 0
  60 * hrs + mins
}

releases_sum <- releases_raw %>%
  mutate(
    # Try both "mdy HMS" and "mdy HM"
    release_date = as.Date(parse_date_time(Start, orders = c("mdy HMS", "mdy HM"))),
    minutes_from_decimals = if ("Duration.in.Decimals" %in% names(.)) {
      60 * suppressWarnings(as.numeric(.data$Duration.in.Decimals))
    } else NA_real_,
    minutes_from_text = if ("Duration" %in% names(.)) {
      parse_duration_to_minutes(.data$Duration)
    } else NA_real_,
    total_minutes_row = dplyr::coalesce(minutes_from_decimals, minutes_from_text)
  ) %>%
  group_by(release_date) %>%
  summarise(total_duration = sum(total_minutes_row, na.rm = TRUE), .groups = "drop")

# Filter out zero-duration rows
releases_nonzero <- filter(releases_sum, total_duration > 0)

########## 3) Scaling ratio for sewage bars ##########

max_precip  <- max(precip_data$precip, na.rm = TRUE)
max_release <- ifelse(nrow(releases_nonzero) > 0,
                      max(releases_nonzero$total_duration, na.rm = TRUE), 0)

# tallest sewage bar = 50% of precip peak
ratio <- if (max_release > 0) 0.5 * (max_precip / max_release) else 0

lower_limit <- -ratio * max_release
upper_limit <- max_precip

########## 4) Plot ##########

rainandsewagefigure <- ggplot() +
  # Precipitation (above y=0)
  geom_line(
    data = precip_data,
    aes(x = datetime, y = precip, color = "Precipitation (daily total)"),
    linetype = "longdash",
    linewidth = 0.4
  ) +
  geom_point(
    data = precip_data,
    size = 1,
    aes(x = datetime, y = precip, color = "Precipitation (daily total)")
  ) +
  
  # Sewage bars (below y=0)
  geom_col(
    data = releases_nonzero,
    aes(x = release_date, y = -ratio * total_duration, fill = total_duration),
    width = 1
  ) +
  
  # Reference line
  geom_hline(yintercept = 0, color = "black") +
  
  # Scales
  scale_x_date(
    breaks = selected_dates,
    labels = custom_labels,
    expand = expansion(mult = 0.02)
  ) +
  scale_y_continuous(
    breaks = seq(0, 25, 5),
    labels = seq(0, 25, 5),
    expand = c(0.1, 0.1)
  ) +
  coord_cartesian(ylim = c(lower_limit, upper_limit)) +
  
  # Colours & legends
  scale_color_manual(
    name   = "Rainfall",
    values = c("Precipitation (daily total)" = "blue"),
    labels = c("Precipitation (daily total)" = "Total precipitation (mm)"),
    guide  = guide_legend(order = 1, label.theme = element_text(face = "plain"))
  ) +
  scale_fill_gradient(
    name = "Total sewage duration (minutes)",
    low  = "lightgreen",
    high = "darkgreen",
    guide = guide_colorbar(order = 2, label.theme = element_text(face = "plain"))
  ) +
  
  # Labels & theme
  labs(
    title = "Daily precipitation and sewage releases",
    x     = "Sampling date",
    y     = "Precipitation (mm)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
    legend.position  = "right",
    legend.direction = "vertical",
    legend.title     = element_text(size = 10, face = "bold"),
    
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10, face = "plain"),
    axis.text.x  = element_text(size = 8, angle = 45, hjust = 1, margin = margin(t = -10)),
    axis.text.y  = element_text(size = 8),
    
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA)
  )

print(rainandsewagefigure)