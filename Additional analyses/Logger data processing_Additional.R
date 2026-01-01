########## 0) Configuration ##########

# install.packages("here") # Install once if needed
library(here)

# Directory of logger CSVs
LOG_DIR  <- here("Dataframes", "Logger data", "Loggers")

# Direct CSV files
ID_CSV   <- here("Dataframes", "Logger data", "Logger ID list.csv")
TIDE_CSV <- here("Dataframes", "Brighton tide heights.csv")

# Timezone & analysis windows
TZ              <- "Europe/London"
INSTALL_FROM    <- as.Date("2022-12-06")
INSTALL_TO      <- as.Date("2023-09-16")
SUMMER_FROM     <- as.Date("2023-06-20")
SUMMER_TO       <- as.Date("2023-09-16")

# Tide detection
TIDE_HOUR_MIN   <- 10   # 10:00
TIDE_HOUR_MAX   <- 14   # 14:00
SMOOTH_K        <- 3    # rollmean window for tide smoothing

# Plate IDs as characters
WHITE <- c('10','7','5','6','2','12','11','3','8','1','9','4')   # Ambient
BLACK <- c('21','20','17','19','13','15','23','18')               # Warmed

# Collapse any identical plate_id+time readings
DEDUP_RULE <- "max"

WRITE_CSV <- FALSE
OUT_DIR   <- here("Outputs")

# Positive y-limit (start at 0) for °C difference plots
pos_temp_limit <- function(ymax, step = 0.5) {
  if (!is.finite(ymax)) ymax <- step
  m <- ceiling(ymax / step) * step
  if (m <= 0) m <- step
  c(0, m)
}

########## 1) Packages ##########

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(patchwork)
  library(zoo)
})
options(dplyr.summarise.inform = FALSE)

########## 2) Read & combine logger data ##########

csv_files <- list.files(path = LOG_DIR, recursive = TRUE, pattern = "\\.csv$", full.names = TRUE)
if (!length(csv_files)) stop("No logger CSVs found in: ", LOG_DIR)

# Read all CSVs, skipping 20 header lines of metadata; attach file base as logger_name
logger_data <- purrr::map_dfr(
  csv_files,
  \(f) readr::read_csv(f, skip = 20, show_col_types = FALSE) %>%
    mutate(logger_name = tools::file_path_sans_ext(basename(f))),
  .id = NULL
)

# Basic check
stopifnot(all(c("logger_name","time","temp") %in% names(logger_data)))

########## 3) Join Logger ID list & derive time/labels ##########

id_map <- readr::read_csv(ID_CSV, show_col_types = FALSE)
head(id_map) # Expected columns: logger_name, Plate_ID

logger_data <- logger_data %>%
  left_join(id_map, by = "logger_name") %>%
  mutate(
    time     = lubridate::ymd_hms(time, tz = TZ, quiet = TRUE),
    plate_id = as.character(Plate_ID),
    group    = dplyr::case_when(
      plate_id %in% WHITE ~ "Ambient",
      plate_id %in% BLACK ~ "Warmed",
      TRUE ~ NA_character_
    )
  )

########## 4) QC (Plate_ID + time; duplicates) ##########

# Diagnostic: duplicates by identical plate_id + time (should loggers have recorded twice at a given timepoint)
dups <- logger_data %>%
  group_by(plate_id, time) %>%
  filter(n() > 1) %>%
  arrange(plate_id, time)

if (nrow(dups) > 0) {
  message("Duplicate readings detected at identical plate_id + time: n=", nrow(dups))
  if (WRITE_CSV) {
    if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)
    dups_qc <- dups %>%
      summarise(
        n_dups = n(),
        temps  = paste(round(temp, 2), collapse = ","),
        sources= paste(unique(logger_name), collapse = "|"),
        .groups = "drop"
      )
    readr::write_csv(dups_qc, file.path(OUT_DIR, "qc_duplicates_by_plate_time.csv"))
  }
}

dedup_by_rule <- function(data, method = c("mean","max","first","last")) {
  method <- match.arg(method)
  data %>%
    arrange(plate_id, time) %>%
    group_by(plate_id, time) %>%
    summarise(
      temp = switch(
        method,
        mean  = mean(temp, na.rm = TRUE),
        max   = max(temp, na.rm = TRUE),
        first = dplyr::first(temp),
        last  = dplyr::last(temp)
      ),
      logger_name = dplyr::first(logger_name),
      Plate_ID    = dplyr::first(Plate_ID),
      group       = dplyr::first(group),
      .groups = "drop"
    )
}

logger_data <- dedup_by_rule(logger_data, DEDUP_RULE)

# Derive date parts after dedup for consistency
logger_data <- logger_data %>%
  mutate(
    date = as.Date(time),
    hour = lubridate::hour(time),
    min  = lubridate::minute(time),
    sec  = lubridate::second(time),
    plate_id = as.character(Plate_ID)
  )

# Temperature range check (check for any anomalous readings)
if (any(!dplyr::between(logger_data$temp, -5, 55), na.rm = TRUE)) {
  message("Warning: temps outside sanity range [-5, 55] detected.")
}

########## 6) Subsets ##########

data_entire       <- logger_data
data_installation <- logger_data %>% filter(date > INSTALL_FROM, date < INSTALL_TO)
data_ADM          <- data_installation
data_summer       <- logger_data %>% filter(date > SUMMER_FROM, date < SUMMER_TO)

########## 7) Tide-adjusted days (days in which low tide occured during peak zenith hours of 10:00–14:00) ##########

# Return Date vector of days with low tide within window, after smoothing
low_tide_days <- function(tide_csv, from, to, hour_min = 10, hour_max = 14, k = 3, tz = "Europe/London") {
  
  # Read as character to avoid readr guessing/parsing issues
  tide <- readr::read_csv(
    tide_csv,
    skip = 11,
    col_names = FALSE,
    col_types = readr::cols(.default = readr::col_character()),
    trim_ws = TRUE,
    show_col_types = FALSE,
    na = c("", "NA")
  )
  
  if (ncol(tide) < 5) stop("Tide file has fewer than 5 columns after skipping headers.")
  names(tide)[1:5] <- c("Row", "Date", "Time", "Height", "Residual")
  tide <- tide %>% dplyr::select(1:5)
  
  # Clean Date strings
  clean_date <- function(x) {
    x %>%
      stringr::str_trim() %>%
      stringr::str_replace("^X", "") %>%                         # drop leading "X" if present
      stringr::str_extract("\\d{2}[/.]\\d{2}[/.]\\d{4}") %>%     # pull a DD/MM/YYYY or DD.MM.YYYY token
      stringr::str_replace_all("\\.", "/") %>%
      dplyr::na_if("")
  }
  
  # Clean Time strings
  clean_time <- function(x) {
    x %>%
      stringr::str_trim() %>%
      stringr::str_replace("^X", "") %>%              # drop leading "X" if present
      stringr::str_replace_all("\\.", ":") %>%        
      dplyr::na_if("") %>%
      dplyr::if_else(stringr::str_detect(., "^\\d{1,2}:\\d{2}$"), paste0(., ":00"), .) %>%  # add seconds if missing
      dplyr::if_else(stringr::str_detect(., "^\\d{1,2}$"),       paste0(., ":00:00"), .)    # hour only -> add mm:ss
  }
  
  tide <- tide %>%
    mutate(
      Date_clean = clean_date(.data$Date),
      Time_clean = clean_time(.data$Time),
      
      # Parse in UTC first, then convert to tz
      DateTime_utc = lubridate::dmy_hms(paste(.data$Date_clean, .data$Time_clean), tz = "UTC", quiet = TRUE),
      DateTime     = lubridate::with_tz(.data$DateTime_utc, tzone = tz),
      
      Height = readr::parse_number(.data$Height)
    )
  
  # Inclusive date window
  from_dt <- as.POSIXct(from, tz = tz)
  to_dt   <- as.POSIXct(to + 1, tz = tz)
  
  # Show any unparsed rows
  bad_dt <- sum(is.na(tide$DateTime))
  bad_h  <- sum(is.na(tide$Height))
  if (bad_dt > 0 || bad_h > 0) {
    message("Tide parse NAs -> DateTime: ", bad_dt, " | Height: ", bad_h)
    print(tide %>% dplyr::filter(is.na(DateTime) | is.na(Height)) %>% dplyr::slice_head(n = 12))
  }
  
  tide <- tide %>%
    filter(!is.na(DateTime), DateTime >= from_dt, DateTime < to_dt) %>%
    arrange(DateTime) %>%
    mutate(
      Height_smooth = zoo::rollmean(Height, k = k, fill = NA, align = "center"),
      lag_s         = dplyr::lag(Height_smooth),
      lead_s        = dplyr::lead(Height_smooth),
      is_trough     = !is.na(lag_s) & !is.na(lead_s) & lag_s > Height_smooth & lead_s > Height_smooth
    )
  
  troughs <- tide %>%
    filter(is_trough) %>%
    mutate(hr = lubridate::hour(DateTime)) %>%
    filter(hr >= hour_min, hr <= hour_max)
  
  unique(as.Date(troughs$DateTime))
}

summer_low_tide_days <- low_tide_days(
  tide_csv = TIDE_CSV,
  from     = SUMMER_FROM,
  to       = SUMMER_TO,
  hour_min = TIDE_HOUR_MIN,
  hour_max = TIDE_HOUR_MAX,
  k        = SMOOTH_K,
  tz       = TZ
)

data_tide <- data_summer %>% filter(date %in% summer_low_tide_days)

########## 8) Summary helper functions ##########

max_or_na <- function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)

# Mean of pad-level maxima over period, with group mean + SE
summarize_max_with_se <- function(data, metric_name = "Mean_Max") {
  pad_max <- data %>%
    filter(!is.na(plate_id), !is.na(group)) %>%
    group_by(plate_id) %>%
    summarise(pad_max = max_or_na(temp), group = dplyr::first(group), .groups = "drop")
  
  pad_max %>%
    group_by(group) %>%
    summarise(
      mean_val = mean(pad_max, na.rm = TRUE),
      se_val   = sd(pad_max,   na.rm = TRUE) / sqrt(dplyr::n()),
      .groups  = "drop"
    ) %>%
    mutate(metric = metric_name)
}

# ADM (Average of Daily Maxima): mean across days per plate, then group mean + SE
summarize_ADM_with_se <- function(data, metric_name = "ADM") {
  pad_daily <- data %>%
    filter(!is.na(plate_id), !is.na(group)) %>%
    group_by(plate_id, date) %>%
    summarise(daily_max = max_or_na(temp), .groups = "drop")
  
  pad_ADM <- pad_daily %>%
    group_by(plate_id) %>%
    summarise(ADM = mean(daily_max, na.rm = TRUE), .groups = "drop") %>%
    left_join(data %>% distinct(plate_id, group), by = "plate_id")
  
  pad_ADM %>%
    group_by(group) %>%
    summarise(
      mean_val = mean(ADM, na.rm = TRUE),
      se_val   = sd(ADM,  na.rm = TRUE) / sqrt(dplyr::n()),
      .groups  = "drop"
    ) %>%
    mutate(metric = metric_name)
}

# Difference helper (Warmed - Ambient)
diff_of_means <- function(df_group_stats) {
  wide <- df_group_stats %>%
    dplyr::select(metric, group, mean_val, se_val) %>%
    tidyr::pivot_wider(
      names_from  = group,
      values_from = c(mean_val, se_val),
      names_sep   = "_"
    )
  
  wide %>%
    dplyr::mutate(
      delta_mean = .data$mean_val_Warmed - .data$mean_val_Ambient,
      se_diff    = sqrt((.data$se_val_Warmed)^2 + (.data$se_val_Ambient)^2),
      ci_lwr     = delta_mean - 1.96 * se_diff,
      ci_upr     = delta_mean + 1.96 * se_diff
    ) %>%
    dplyr::select(metric, delta_mean, se_diff, ci_lwr, ci_upr)
}

########## 9) Calculate summaries ##########

# Group summaries
summary_mean_max_stats <- summarize_max_with_se(data_installation, "Mean_Max")
summary_ADM_JanSep     <- summarize_ADM_with_se(data_ADM,         "ADM_Jan_Sep")
summary_ADM_Summer     <- summarize_ADM_with_se(data_summer,      "ADM_Jun_Sep")
summary_ADM_Tide       <- summarize_ADM_with_se(data_tide,        "ADM_TideAdj_Jun_Sep")

# Check: both groups present
check_groups <- function(x, name = "object") {
  need <- c("Ambient", "Warmed")
  have <- unique(x$group)
  miss <- setdiff(need, have)
  if (length(miss) > 0) {
    warning(name, " is missing group(s): ", paste(miss, collapse = ", "),
            ". Δ(Warmed - Ambient) may be NA.")
  }
}
check_groups(summary_mean_max_stats, "summary_mean_max_stats")
check_groups(summary_ADM_Tide,       "summary_ADM_Tide")

# Compute Δ(Warmed - Ambient) with 95% CI
diff_mean_max <- diff_of_means(summary_mean_max_stats) %>%
  dplyr::mutate(label = "Mean Max")

diff_tide_ADM <- diff_of_means(summary_ADM_Tide) %>%
  dplyr::mutate(label = "Average Daily Mean (Tide-adjusted)")

# Print
all_summaries <- dplyr::bind_rows(
  summary_mean_max_stats,
  summary_ADM_JanSep,
  summary_ADM_Summer,
  summary_ADM_Tide
)

all_deltas <- dplyr::bind_rows(
  diff_mean_max,
  diff_tide_ADM
)

print(all_summaries)
print(all_deltas)

########## 10) Plotting ##########

# Recompute CI to guarantee it is the same rule for both panels
recalc_quick_ci <- function(df) {
  df %>%
    dplyr::mutate(
      ci_lwr = delta_mean - 1.96 * se_diff,
      ci_upr = delta_mean + 1.96 * se_diff
    )
}

diff_mean_max_ci <- recalc_quick_ci(diff_mean_max)
diff_tide_ADM_ci <- recalc_quick_ci(diff_tide_ADM)

# Y-limits
upper_mean_max_ci <- pos_temp_limit(
  max(diff_mean_max_ci$ci_upr, diff_mean_max_ci$delta_mean, na.rm = TRUE),
  step = 0.5
)
upper_tide_adm_ci <- pos_temp_limit(
  max(diff_tide_ADM_ci$ci_upr, diff_tide_ADM_ci$delta_mean, na.rm = TRUE),
  step = 0.5
)

# Δ Mean of pad-level maxima
p_mean_max_ci <- ggplot(diff_mean_max_ci, aes(x = label, y = delta_mean)) +
  geom_col(width = 0.6, fill = "#FFA500") +
  geom_errorbar(
    aes(ymin = ci_lwr, ymax = ci_upr),
    width = 0.15, linewidth = 0.5, color = "black"
  ) +
  scale_y_continuous(limits = upper_mean_max_ci) +
  labs(
    title = "Mean max",
    x = NULL,
    y = "Absolute warming increase (°C)"
  ) +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_blank())

# Δ ADM (tide-adjusted)
p_tide_adm_ci <- ggplot(diff_tide_ADM_ci, aes(x = label, y = delta_mean)) +
  geom_col(width = 0.6, fill = "#FFA500") +
  geom_errorbar(
    aes(ymin = ci_lwr, ymax = ci_upr),
    width = 0.15, linewidth = 0.5, color = "black"
  ) +
  scale_y_continuous(limits = upper_tide_adm_ci) +
  labs(
    title = "Average daily mean",
    x = NULL,
    y = "Absolute warming increase (°C)"
  ) +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_blank())

final_logger_panel_ci <- (p_mean_max_ci + p_tide_adm_ci) &
  theme(plot.title = element_text(hjust = 0.5))

print(final_logger_panel_ci)