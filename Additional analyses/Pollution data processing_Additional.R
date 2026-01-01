########## 0) Configuration ##########

# install.packages("here") # Install once if needed
library(here)

# Pre-experiment pollution CSV
POLLUTION_CSV <- here("Dataframes", "Pre-experiment pollution.csv")

# During-experiment pollution
EXP_POLLUTION_CSV <- here("Dataframes", "Experiment pollution.csv")

# Plot aesthetics
BASE_SIZE <- 10
pad_fill  <- c("Polluted" = "darkgreen", "Non-polluted" = "grey75") #No non-polluted for these figures, but just kept for consistency

# Date order for pre-experiment figure
DATE_LEVELS <- c("Sep-22", "Oct-22", "Nov-22")

# Date order + labels for during-experiment figure
EXP_DATE_LEVELS <- c("Jun1-23","Jul1-23","Jul2-23","Aug1-23","Aug2-23","Sep1-23","Sep2-23")
EXP_DATE_LABELS <- c("20/06/23","06/07/23","20/07/23","01/08/23","17/08/23","01/09/23","15/09/23")
EXP_LABEL_MAP   <- setNames(EXP_DATE_LABELS, EXP_DATE_LEVELS)

# Positive y-limit (start at 0) for % plots
pos_percent_limit <- function(ymax, step = 25) {
  if (!is.finite(ymax)) ymax <- step
  m <- ceiling(ymax / step) * step
  if (m <= 0) m <- step
  c(0, m)
}

########## 1) Packages ##########

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
})
options(dplyr.summarise.inform = FALSE)

########## 2) Load & normalise (pre-experiment) ##########

df <- readr::read_csv(POLLUTION_CSV, show_col_types = FALSE)

# Normalise column names for 'Non-polluted' vs 'Nonpolluted'
need_cols <- c("Date", "Pollutant", "Polluted")
if (!all(need_cols %in% names(df))) {
  stop("CSV must contain columns: ", paste(need_cols, collapse = ", "))
}
if (!("Non-polluted" %in% names(df) | "Nonpolluted" %in% names(df))) {
  stop("CSV must contain either 'Non-polluted' or 'Nonpolluted'.")
}
if ("Nonpolluted" %in% names(df)) {
  df <- dplyr::rename(df, `Non-polluted` = `Nonpolluted`)
}

########## 3) Tidy to long + factor levels (pre-experiment) ##########

df_long <- df %>%
  tidyr::pivot_longer(
    cols = c(Polluted, `Non-polluted`),
    names_to = "Site",
    values_to = "Measurement"
  ) %>%
  dplyr::mutate(
    Site = factor(Site, levels = c("Polluted", "Non-polluted")),
    Date = factor(Date, levels = intersect(DATE_LEVELS, unique(Date)))
  )

########## 4) Summaries (pre-experiment) ##########

# Per-date means and SE
df_summary <- df_long %>%
  dplyr::group_by(Date, Pollutant, Site) %>%
  dplyr::summarise(
    n        = dplyr::n(),
    mean_val = mean(Measurement, na.rm = TRUE),
    se       = sd(Measurement,   na.rm = TRUE) / sqrt(n),
    .groups  = "drop"
  )

# Overall (across dates)
df_agg <- df_long %>%
  dplyr::group_by(Pollutant, Site) %>%
  dplyr::summarise(
    n        = dplyr::n(),
    mean_val = mean(Measurement, na.rm = TRUE),
    se       = sd(Measurement,   na.rm = TRUE) / sqrt(n),
    .groups  = "drop"
  ) %>%
  dplyr::mutate(Site = factor(Site, levels = c("Polluted","Non-polluted")))

print(df_summary)
print(df_agg)

# Proportional changes
# Uses log-ratio delta method for CIs; baseline = 0% when Polluted == Non-polluted

# Per-date % difference (Polluted vs Non-polluted)
pre_prop_by_date <- df_summary %>%
  dplyr::select(Date, Pollutant, Site, mean_val, se, n) %>%
  tidyr::pivot_wider(
    names_from = Site,
    values_from = c(mean_val, se, n),
    names_sep = "_"
  ) %>%
  dplyr::mutate(
    invalid = is.na(mean_val_Polluted) |
      is.na(`mean_val_Non-polluted`) |
      mean_val_Polluted <= 0 |
      `mean_val_Non-polluted` <= 0,
    logR  = ifelse(invalid, NA_real_,
                   log(mean_val_Polluted) - log(`mean_val_Non-polluted`)),
    seLog = ifelse(invalid, NA_real_,
                   sqrt((se_Polluted/mean_val_Polluted)^2 +
                          (`se_Non-polluted`/`mean_val_Non-polluted`)^2)),
    R_est   = exp(logR),
    R_lwr   = exp(logR - 1.96*seLog),
    R_upr   = exp(logR + 1.96*seLog),
    pct_est = 100*(R_est - 1),
    pct_lwr = 100*(R_lwr - 1),
    pct_upr = 100*(R_upr - 1)
  ) %>%
  dplyr::select(Date, Pollutant, pct_est, pct_lwr, pct_upr)

# Overall (across dates) % difference
pre_prop_overall <- df_agg %>%
  dplyr::select(Pollutant, Site, mean_val, se, n) %>%
  tidyr::pivot_wider(
    names_from = Site,
    values_from = c(mean_val, se, n),
    names_sep = "_"
  ) %>%
  dplyr::mutate(
    invalid = is.na(mean_val_Polluted) |
      is.na(`mean_val_Non-polluted`) |
      mean_val_Polluted <= 0 |
      `mean_val_Non-polluted` <= 0,
    logR  = ifelse(invalid, NA_real_,
                   log(mean_val_Polluted) - log(`mean_val_Non-polluted`)),
    seLog = ifelse(invalid, NA_real_,
                   sqrt((se_Polluted/mean_val_Polluted)^2 +
                          (`se_Non-polluted`/`mean_val_Non-polluted`)^2)),
    R_est   = exp(logR),
    R_lwr   = exp(logR - 1.96*seLog),
    R_upr   = exp(logR + 1.96*seLog),
    pct_est = 100*(R_est - 1),
    pct_lwr = 100*(R_lwr - 1),
    pct_upr = 100*(R_upr - 1)
  ) %>%
  dplyr::select(Pollutant, pct_est, pct_lwr, pct_upr)

print(pre_prop_by_date)
print(pre_prop_overall)

########## 5) Plots; per-date (pre-experiment) ##########

# Positive y-limits by pollutant
upper_pre <- pre_prop_by_date %>%
  dplyr::group_by(Pollutant) %>%
  dplyr::summarise(
    ymax = suppressWarnings(max(pct_upr, pct_est, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  dplyr::mutate(lims = purrr::map(ymax, ~pos_percent_limit(.x))) %>%
  {setNames(.$lims, .$Pollutant)}

# Nitrate
p_nitrate <- pre_prop_by_date %>%
  dplyr::filter(Pollutant == "Nitrate") %>%
  ggplot2::ggplot(ggplot2::aes(x = Date, y = pct_est)) +
  ggplot2::geom_col(width = 0.8, fill = "darkgreen") +
  ggplot2::geom_errorbar(ggplot2::aes(ymin = pct_lwr, ymax = pct_upr),
                         width = 0.2, size = 0.5, color = "black") +
  ggplot2::scale_y_continuous(limits = unlist(upper_pre[["Nitrate"]])) +
  ggplot2::labs(title = "Nitrate", x = "Date",
                y = "Proportional increase (%)") +
  ggplot2::theme_minimal(base_size = BASE_SIZE) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

# Phosphate
p_phosphate <- pre_prop_by_date %>%
  dplyr::filter(Pollutant == "Phosphate") %>%
  ggplot2::ggplot(ggplot2::aes(x = Date, y = pct_est)) +
  ggplot2::geom_col(width = 0.8, fill = "darkgreen") +
  ggplot2::geom_errorbar(ggplot2::aes(ymin = pct_lwr, ymax = pct_upr),
                         width = 0.2, size = 0.5, color = "black") +
  ggplot2::scale_y_continuous(limits = unlist(upper_pre[["Phosphate"]])) +
  ggplot2::labs(title = "Phosphate", x = "Date",
                y = "Proportional increase (%)") +
  ggplot2::theme_minimal(base_size = BASE_SIZE) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

final_pollution_panel <- (p_nitrate + p_phosphate) &
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
print(final_pollution_panel)

########## 6) Plots; overall (aggregated, pre-experiment) ##########

# Positive y-limits for aggregated
upper_pre_agg <- pre_prop_overall %>%
  dplyr::group_by(Pollutant) %>%
  dplyr::summarise(
    ymax = suppressWarnings(max(pct_upr, pct_est, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  dplyr::mutate(lims = purrr::map(ymax, ~pos_percent_limit(.x))) %>%
  {setNames(.$lims, .$Pollutant)}

p_nitrate_agg <- pre_prop_overall %>%
  dplyr::filter(Pollutant == "Nitrate") %>%
  ggplot2::ggplot(ggplot2::aes(x = Pollutant, y = pct_est)) +
  ggplot2::geom_col(width = 0.6, fill = "darkgreen") +
  ggplot2::geom_errorbar(ggplot2::aes(ymin = pct_lwr, ymax = pct_upr),
                         width = 0.15, size = 0.5, color = "black") +
  ggplot2::scale_y_continuous(limits = unlist(upper_pre_agg[["Nitrate"]])) +
  ggplot2::labs(title = "Nitrate", x = NULL, y = "Proportional increase (%)") +
  ggplot2::theme_minimal(base_size = BASE_SIZE) +
  ggplot2::theme(axis.text.x = ggplot2::element_blank())

p_phosphate_agg <- pre_prop_overall %>%
  dplyr::filter(Pollutant == "Phosphate") %>%
  ggplot2::ggplot(ggplot2::aes(x = Pollutant, y = pct_est)) +
  ggplot2::geom_col(width = 0.6, fill = "darkgreen") +
  ggplot2::geom_errorbar(ggplot2::aes(ymin = pct_lwr, ymax = pct_upr),
                         width = 0.15, size = 0.5, color = "black") +
  ggplot2::scale_y_continuous(limits = unlist(upper_pre_agg[["Phosphate"]])) +
  ggplot2::labs(title = "Phosphate", x = NULL, y = "Proportional increase (%)") +
  ggplot2::theme_minimal(base_size = BASE_SIZE) +
  ggplot2::theme(axis.text.x = ggplot2::element_blank())

final_pollution_panel_agg <- (p_nitrate_agg + p_phosphate_agg) &
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
print(final_pollution_panel_agg)

########## 7) During-experiment pollution; load & tidy ##########

df_exp <- readr::read_csv(EXP_POLLUTION_CSV, show_col_types = FALSE)

df_exp_long <- df_exp %>%
  tidyr::pivot_longer(
    cols = c(Polluted, `Non-polluted`),
    names_to = "Site",
    values_to = "Measurement"
  ) %>%
  dplyr::mutate(
    Site = factor(Site, levels = c("Polluted", "Non-polluted")),
    Date = factor(Date, levels = EXP_DATE_LEVELS)   # keep empty Jun1-23 for nitrate if present
  )

########## 8) During-experiment; per-date summaries ##########

df_exp_summary <- df_exp_long %>%
  dplyr::group_by(Date, Pollutant, Site) %>%
  dplyr::summarise(
    n_non_na = sum(!is.na(Measurement)),
    mean_val = mean(Measurement, na.rm = TRUE),
    se       = ifelse(n_non_na > 0, sd(Measurement, na.rm = TRUE) / sqrt(n_non_na), NA_real_),
    .groups  = "drop"
  )

# Proportional changes
exp_prop_by_date <- df_exp_summary %>%
  dplyr::select(Date, Pollutant, Site, mean_val, se) %>%
  tidyr::pivot_wider(
    names_from = Site,
    values_from = c(mean_val, se),
    names_sep = "_"
  ) %>%
  dplyr::mutate(
    invalid = is.na(mean_val_Polluted) |
      is.na(`mean_val_Non-polluted`) |
      mean_val_Polluted <= 0 |
      `mean_val_Non-polluted` <= 0,
    logR  = ifelse(invalid, NA_real_,
                   log(mean_val_Polluted) - log(`mean_val_Non-polluted`)),
    seLog = ifelse(invalid, NA_real_,
                   sqrt((se_Polluted/mean_val_Polluted)^2 +
                          (`se_Non-polluted`/`mean_val_Non-polluted`)^2)),
    R_est   = exp(logR),
    R_lwr   = exp(logR - 1.96*seLog),
    R_upr   = exp(logR + 1.96*seLog),
    pct_est = 100*(R_est - 1),
    pct_lwr = 100*(R_lwr - 1),
    pct_upr = 100*(R_upr - 1)
  ) %>%
  dplyr::select(Date, Pollutant, pct_est, pct_lwr, pct_upr) %>%
  dplyr::mutate(Date = factor(Date, levels = EXP_DATE_LEVELS))

########## 9) During-experiment; per-date plots ##########

# Positive y-limits by pollutant (during experiment)
upper_exp <- exp_prop_by_date %>%
  dplyr::group_by(Pollutant) %>%
  dplyr::summarise(
    ymax = suppressWarnings(max(pct_upr, pct_est, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  dplyr::mutate(lims = purrr::map(ymax, ~pos_percent_limit(.x))) %>%
  {setNames(.$lims, .$Pollutant)}

# Nitrate
dfN <- exp_prop_by_date %>% dplyr::filter(Pollutant == "Nitrate")
dfN_draw <- dfN %>% dplyr::filter(!is.na(pct_est))

p_nitrate_exp <- ggplot2::ggplot(dfN, ggplot2::aes(x = Date, y = pct_est)) +
  ggplot2::geom_col(data = dfN_draw, width = 0.8, fill = "darkgreen") +
  ggplot2::geom_errorbar(
    data = dfN_draw,
    ggplot2::aes(ymin = pct_lwr, ymax = pct_upr),
    width = 0.2, size = 0.5, color = "black"
  ) +
  ggplot2::scale_x_discrete(limits = EXP_DATE_LEVELS, labels = EXP_LABEL_MAP) +
  ggplot2::scale_y_continuous(limits = unlist(upper_exp[["Nitrate"]])) +
  ggplot2::labs(title = "Nitrate", x = "Date",
                y = "Proportional increase (%)") +
  ggplot2::theme_minimal(base_size = BASE_SIZE) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

# Phosphate
p_phosphate_exp <- exp_prop_by_date %>%
  dplyr::filter(Pollutant == "Phosphate") %>%
  ggplot2::ggplot(ggplot2::aes(x = Date, y = pct_est)) +
  ggplot2::geom_col(width = 0.8, fill = "darkgreen") +
  ggplot2::geom_errorbar(ggplot2::aes(ymin = pct_lwr, ymax = pct_upr),
                         width = 0.2, size = 0.5, color = "black") +
  ggplot2::scale_x_discrete(limits = EXP_DATE_LEVELS, labels = EXP_LABEL_MAP) +
  ggplot2::scale_y_continuous(limits = unlist(upper_exp[["Phosphate"]])) +
  ggplot2::labs(title = "Phosphate", x = "Date",
                y = "Proportional increase (%)") +
  ggplot2::theme_minimal(base_size = BASE_SIZE) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

final_experiment_panel <- (p_nitrate_exp / p_phosphate_exp) &
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
print(final_experiment_panel)

########## 10) Summary tables ##########

# Pre-experiment per-date (%)
preexp_by_date_tbl <- pre_prop_by_date %>%
  dplyr::arrange(Pollutant, Date)
print(preexp_by_date_tbl, n = 12)

# Pre-experiment aggregated (%)
preexp_agg_tbl <- pre_prop_overall %>%
  dplyr::arrange(Pollutant)
print(preexp_agg_tbl)

# During-experiment per-date (%)
exp_by_date_tbl <- exp_prop_by_date %>%
  dplyr::arrange(Pollutant, Date)
print(exp_by_date_tbl, n = 28)