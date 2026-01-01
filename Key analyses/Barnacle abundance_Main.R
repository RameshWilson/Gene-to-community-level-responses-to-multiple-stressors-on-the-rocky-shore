# Barnacle abundance_Main.R
#
#   Brighton experiment (2×2): Warming (Pad_Colour) × Pollution (Pad_Region)
#   Response: Barnacle_Count (counts)
#
#   Outputs:
#     1) Per-date Negative Binomial GLMs (MASS::glm.nb)
#     2) GLM interaction classification table
#     3) Whole-summer GLM-style visualisation (ggeffects predictions + raw overlay)
#     4) Conditional per-date panels (pollution-stratified lines/ribbons across pad colour)
#     5) Whole-season NB GAM (mgcv::gam) with separate smooth per 2×2 treatment
#     6) GAM interaction classification table (parametric fixed effects)
#     7) Whole-season GAM visualisation (smooth predictions + raw overlay)
#
#   Evidence/diagnostics (Poisson vs NB, dispersion tests, GAM checks) are kept in:
#     Barnacle abundance_Additional.R

########## 0) Configuration ##########

library(here)

DATA_CSV <- here("Dataframes", "Full Brighton dataframe.csv")

# Sampling dates (as Date)
D_JUN  <- as.Date("2023-06-20")
D_JUL1 <- as.Date("2023-07-06")
D_JUL2 <- as.Date("2023-07-20")
D_AUG1 <- as.Date("2023-08-01")
D_AUG2 <- as.Date("2023-08-17")
D_SEP1 <- as.Date("2023-09-01")
D_SEP2 <- as.Date("2023-09-15")

# Labels used across plots/tables
DATE_LABELS <- c("June", "July 1", "July 2", "August 1", "August 2", "September 1", "September 2")

# Date label lookup
DATE_LOOKUP <- setNames(
  DATE_LABELS,
  as.character(c(D_JUN, D_JUL1, D_JUL2, D_AUG1, D_AUG2, D_SEP1, D_SEP2))
)

# Interaction classification constants (used across GLM + GAM)
#   - Treat exp(beta_int) as the multiplicative interaction ratio on the response scale
#   - "Additive" if CI overlaps the ±5% equivalence band [0.95, 1.05]
EQUIV_LO <- 0.95
EQUIV_HI <- 1.05
P_ALPHA  <- 0.05

# Colours for Treatment
TREATMENT_COLORS <- c(
  "Non-Polluted + Ambient" = "grey75",
  "Non-Polluted + Warm"    = "orange",
  "Polluted + Ambient"     = "darkgreen",
  "Polluted + Warm"        = "#CC79A7"
)

# Legend order
LEGEND_ORDER <- c(
  "Polluted + Ambient",
  "Polluted + Warm",
  "Non-Polluted + Ambient",
  "Non-Polluted + Warm"
)

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(lubridate)
  library(MASS)
  library(ggeffects)
  library(colorspace)
  library(broom)
  library(cowplot)
  library(mgcv)
  library(DT)
  library(scales)
  library(stringr)
  library(tibble)
})

########## 1) Load & prep data ##########

brighton_data <- readr::read_csv(DATA_CSV, show_col_types = FALSE)

# Parse Sampling_Date: try dmy(); if too many NAs, fall back to ymd()
if (!inherits(brighton_data$Sampling_Date, "Date")) {
  tmp_dmy <- suppressWarnings(lubridate::dmy(brighton_data$Sampling_Date))
  frac_na <- mean(is.na(tmp_dmy))
  brighton_data$Sampling_Date <- if (frac_na < 0.5) tmp_dmy else suppressWarnings(lubridate::ymd(brighton_data$Sampling_Date))
}

# Recode factor baselines (controls first):
#   Pad_Colour: W = Ambient (control); B = Warm (treatment)
#   Pad_Region: N = Non-Polluted (control); P = Polluted (treatment)
brighton_data <- brighton_data %>%
  mutate(
    Pad_Colour = factor(Pad_Colour, levels = c("W", "B")),
    Pad_Region = factor(Pad_Region, levels = c("N", "P"))
  )

# Helper: ensure correct reference levels within any subset
prep_subset <- function(df) {
  df$Pad_Colour <- relevel(as.factor(df$Pad_Colour), ref = "W")
  df$Pad_Region <- relevel(as.factor(df$Pad_Region), ref = "N")
  df
}

########## 2) Individual-date NB GLMs (final models) ##########

# Subsets per date
JUN   <- prep_subset(dplyr::filter(brighton_data, Sampling_Date == D_JUN))
JUL_1 <- prep_subset(dplyr::filter(brighton_data, Sampling_Date == D_JUL1))
JUL_2 <- prep_subset(dplyr::filter(brighton_data, Sampling_Date == D_JUL2))
AUG_1 <- prep_subset(dplyr::filter(brighton_data, Sampling_Date == D_AUG1))
AUG_2 <- prep_subset(dplyr::filter(brighton_data, Sampling_Date == D_AUG2))
SEP_1 <- prep_subset(dplyr::filter(brighton_data, Sampling_Date == D_SEP1))
SEP_2 <- prep_subset(dplyr::filter(brighton_data, Sampling_Date == D_SEP2))

# Final NB models
glm.nb_JUN   <- MASS::glm.nb(Barnacle_Count ~ Pad_Colour * Pad_Region, data = JUN)
glm.nb_JUL_1 <- MASS::glm.nb(Barnacle_Count ~ Pad_Colour * Pad_Region, data = JUL_1)
glm.nb_JUL_2 <- MASS::glm.nb(Barnacle_Count ~ Pad_Colour * Pad_Region, data = JUL_2)
glm.nb_AUG_1 <- MASS::glm.nb(Barnacle_Count ~ Pad_Colour * Pad_Region, data = AUG_1)
glm.nb_AUG_2 <- MASS::glm.nb(Barnacle_Count ~ Pad_Colour * Pad_Region, data = AUG_2)
glm.nb_SEP_1 <- MASS::glm.nb(Barnacle_Count ~ Pad_Colour * Pad_Region, data = SEP_1)
glm.nb_SEP_2 <- MASS::glm.nb(Barnacle_Count ~ Pad_Colour * Pad_Region, data = SEP_2)

########## 3) GLM interaction classification + summary table ##########

extract_barnacle_abundance_model_stats <- function(model, date_label) {
  
  # Coefficient table (log link scale)
  #   α    = baseline (Ambient + Non-Polluted)
  #   βw   = warming effect in Non-Polluted
  #   βp   = pollution effect under Ambient
  #   βint = interaction (deviation beyond additivity)
  coef_summary <- broom::tidy(model, conf.int = TRUE, conf.level = 0.95)
  
  # If conf.int was not returned, construct Wald intervals
  if (!all(c("conf.low", "conf.high") %in% names(coef_summary))) {
    coef_summary <- coef_summary %>%
      mutate(
        conf.low  = estimate - 1.96 * std.error,
        conf.high = estimate + 1.96 * std.error
      )
  }
  
  α_row   <- coef_summary %>% filter(term == "(Intercept)")
  w_row   <- coef_summary %>% filter(term == "Pad_ColourB")
  p_row   <- coef_summary %>% filter(term == "Pad_RegionP")
  int_row <- coef_summary %>% filter(term == "Pad_ColourB:Pad_RegionP")
  
  α    <- α_row$estimate
  βw   <- w_row$estimate
  βp   <- p_row$estimate
  βint <- int_row$estimate
  
  # Interaction p-value (fallback to Wald z if missing)
  se_int <- int_row$std.error
  p_int  <- int_row$p.value
  if (is.null(p_int) || is.na(p_int)) {
    z_int <- βint / se_int
    p_int <- 2 * pnorm(abs(z_int), lower.tail = FALSE)
  }
  
  # Additive (null) vs observed (with interaction), on the link scale
  Δadd_log <- βw + βp
  Δobs_log <- Δadd_log + βint
  
  # % change relative to baseline on response scale
  pct_add <- (exp(Δadd_log) - 1) * 100
  pct_obs <- (exp(Δobs_log) - 1) * 100
  
  # Expected means (response scale)
  null_pred     <- exp(α + Δadd_log)
  observed_pred <- exp(α + Δobs_log)
  
  # Interaction ratio on response scale
  # exp(βint) > 1: combined effect larger than additive expectation
  # exp(βint) < 1: combined effect smaller than additive expectation
  int_ratio    <- exp(βint)
  int_ratio_lo <- exp(int_row$conf.low)
  int_ratio_hi <- exp(int_row$conf.high)
  
  # Directional label (only used if CI fully outside equivalence band AND p < alpha)
  int_class_point <- case_when(
    abs(int_ratio - 1) < 0.05                                ~ "Additive (|Δobs–Δadd|<5%)",
    Δadd_log >  0 & Δobs_log <  0                            ~ "Reversal (Pos to Neg)",
    Δadd_log <  0 & Δobs_log >  0                            ~ "Reversal (Neg to Pos)",
    βw > 0 & βp > 0 & βint >  0                              ~ "Synergism (>1)",
    βw > 0 & βp > 0 & βint <  0                              ~ "Antagonism (<1)",
    βw < 0 & βp < 0 & βint <  0                              ~ "Synergism (<1)",
    βw < 0 & βp < 0 & βint >  0                              ~ "Antagonism (>1)",
    sign(Δadd_log) == sign(βint)                             ~ "Synergism (same–sign)",
    sign(Δadd_log) != sign(βint)                             ~ "Antagonism (opp–sign)",
    TRUE                                                     ~ "Unclassified"
  )
  
  # Equivalence-band logic using CI of interaction ratio
  ci_overlaps_band <- !is.na(int_ratio_lo) &&
    !(int_ratio_hi < EQUIV_LO || int_ratio_lo > EQUIV_HI)
  
  ci_outside_band <- !is.na(int_ratio_lo) &&
    (int_ratio_hi < EQUIV_LO || int_ratio_lo > EQUIV_HI)
  
  # Final interaction class (p-value gate first; then CI band)
  int_class_final <- case_when(
    is.na(p_int) | p_int >= P_ALPHA     ~ "Unclassified (p≥0.05)",
    ci_overlaps_band                    ~ "Additive (CI overlaps ±5%)",
    ci_outside_band & p_int < P_ALPHA   ~ int_class_point,
    TRUE                                ~ "Unclassified"
  )
  
  # Perspective metrics under each stressor stratum
  pollution_pct               <- (exp(βp)        - 1) * 100
  pollution_under_warming_pct <- (exp(βp + βint) - 1) * 100
  warming_pct                 <- (exp(βw)        - 1) * 100
  warming_under_pollution_pct <- (exp(βw + βint) - 1) * 100
  
  # Assemble table rows (one row per term, with interaction-only columns populated)
  coef_summary %>%
    filter(term %in% c("Pad_ColourB", "Pad_RegionP", "Pad_ColourB:Pad_RegionP")) %>%
    transmute(
      Date                         = date_label,
      Term                         = case_when(
        term == "Pad_ColourB"             ~ "Warming",
        term == "Pad_RegionP"             ~ "Pollution",
        term == "Pad_ColourB:Pad_RegionP" ~ "Warming × Pollution"
      ),
      Estimate                     = round(estimate, 4),
      `Std. Error`                 = round(std.error, 4),
      `p-value`                    = format.pval(p.value, digits = 3, eps = 0.001),
      Pct_Add                      = if_else(term == "Pad_ColourB:Pad_RegionP", paste0(round(pct_add, 1), "%"), NA_character_),
      Pct_Obs                      = if_else(term == "Pad_ColourB:Pad_RegionP", paste0(round(pct_obs, 1), "%"), NA_character_),
      `Pollution -> under Warming`  = if_else(
        term == "Pad_ColourB:Pad_RegionP",
        sprintf("Pollution’s %+0.1f%% to %+0.1f%%", pollution_pct, pollution_under_warming_pct),
        NA_character_
      ),
      `Warming -> under Pollution`  = if_else(
        term == "Pad_ColourB:Pad_RegionP",
        sprintf("Warming’s %+0.1f%% to %+0.1f%%", warming_pct, warming_under_pollution_pct),
        NA_character_
      ),
      `Interaction Class`          = if_else(term == "Pad_ColourB:Pad_RegionP", int_class_final, NA_character_),
      `Null Prediction`            = if_else(term == "Pad_ColourB:Pad_RegionP", round(null_pred, 2), NA_real_),
      `Observed Prediction`        = if_else(term == "Pad_ColourB:Pad_RegionP", round(observed_pred, 2), NA_real_),
      `Interaction Ratio`          = if_else(term == "Pad_ColourB:Pad_RegionP", round(int_ratio, 3), NA_real_),
      `Interaction Ratio CI`       = if_else(term == "Pad_ColourB:Pad_RegionP",
                                             sprintf("[%.3f, %.3f]", int_ratio_lo, int_ratio_hi),
                                             NA_character_)
    )
}

# Collect models in date order
glm_nb_models <- list(
  glm.nb_JUN, glm.nb_JUL_1, glm.nb_JUL_2,
  glm.nb_AUG_1, glm.nb_AUG_2,
  glm.nb_SEP_1, glm.nb_SEP_2
)

dates_vec <- DATE_LABELS

barnacle_abundance_results_table <- purrr::map2_df(
  glm_nb_models, dates_vec, extract_barnacle_abundance_model_stats
)

DT::datatable(
  barnacle_abundance_results_table,
  caption = "Brighton barnacles: NB-GLM fixed effects & interaction summary",
  options = list(pageLength = 7, scrollX = TRUE)
)

########## 4) Whole-summer GLM-style visualisation ##########

# Predictions per date
pred_list <- vector("list", length(glm_nb_models))

for (i in seq_along(glm_nb_models)) {
  
  pr <- ggeffects::ggpredict(glm_nb_models[[i]], terms = c("Pad_Colour", "Pad_Region")) %>%
    mutate(
      Sampling_Date = dates_vec[i],
      Pad_Region    = factor(group, levels = c("N", "P"), labels = c("Non-Polluted", "Polluted")),
      Pad_Colour    = factor(x,     levels = c("W", "B"), labels = c("Ambient", "Warm")),
      Treatment     = paste(Pad_Region, Pad_Colour, sep = " + ")
    )
  
  pred_list[[i]] <- pr
}

all_preds <- bind_rows(pred_list) %>%
  filter(!is.na(Treatment))

all_preds$Sampling_Date <- factor(all_preds$Sampling_Date, levels = dates_vec)

# Harmonise raw data labels for overlay
brighton_data_plot <- brighton_data %>%
  mutate(
    Pad_Region = factor(Pad_Region, levels = c("N", "P"), labels = c("Non-Polluted", "Polluted")),
    Pad_Colour = factor(Pad_Colour, levels = c("W", "B"), labels = c("Ambient", "Warm")),
    Treatment  = paste(Pad_Region, Pad_Colour, sep = " + "),
    Sampling_Date = unname(DATE_LOOKUP[as.character(Sampling_Date)])
  ) %>%
  mutate(Sampling_Date = factor(Sampling_Date, levels = levels(all_preds$Sampling_Date))) %>%
  filter(!is.na(Sampling_Date))

# Plot
line_plot_barnacles <- ggplot(all_preds, aes(x = Sampling_Date, y = predicted, group = Treatment)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = Treatment), alpha = 0.1, color = NA) +
  geom_line(
    aes(color = Treatment),
    linewidth = 0.5,
    alpha = 0.30,
    linetype = "dashed"
  ) +
  geom_point(
    aes(fill = Treatment, colour = after_scale(colorspace::darken(fill, 0.35))),
    shape = 21,
    size  = 2.7,  
    stroke = 0.9, 
    show.legend = FALSE
  ) +
  geom_jitter(
    data = brighton_data_plot,
    aes(x = Sampling_Date, y = Barnacle_Count, color = Treatment),
    shape = 20, width = 0.2, alpha = 0.7, size = 1.5, inherit.aes = FALSE
  ) +
  labs(
    title = "Barnacle abundance",
    x = "Sampling date",
    y = expression(Barnacle~abundance~"(log"[10]*" scale)")
  ) +
  theme_minimal() +
  theme(
    plot.title   = element_text(face = "bold", hjust = 0.5, size = 12),
    axis.title   = element_text(size = 10),
    axis.text    = element_text(size = 7),
    legend.title = element_text(face = "bold", size = 9),
    legend.text  = element_text(size = 8),
    legend.box.background = element_blank(),
    axis.title.x = element_text(margin = margin(t = 10))
  ) +
  scale_y_continuous(trans = "log10", labels = scales::comma) +
  scale_color_manual(values = TREATMENT_COLORS, breaks = LEGEND_ORDER) +
  scale_fill_manual(values = TREATMENT_COLORS, breaks = LEGEND_ORDER)

print(line_plot_barnacles)

########## 5) Conditional plot panels ##########

# Within each date: predicted means across Pad_Colour, stratified by Pad_Region
REGION_LEVELS <- c("Non-Polluted", "Polluted")
REGION_VALUES <- unname(c(
  TREATMENT_COLORS["Non-Polluted + Ambient"],
  TREATMENT_COLORS["Polluted + Ambient"]
))

# Extract the interaction p-value from a fitted GLM
get_interaction_pvalue <- function(model) {
  sm <- broom::tidy(model)
  r  <- dplyr::filter(sm, term == "Pad_ColourB:Pad_RegionP")
  if (nrow(r) > 0) r$p.value[1] else NA_real_
}

create_conditional_plot <- function(model, date_label, include_legend = FALSE, italicize_p = FALSE) {
  
  pr <- ggeffects::ggpredict(model, terms = c("Pad_Colour", "Pad_Region")) %>%
    dplyr::mutate(
      Pad_Region = factor(group, levels = c("N", "P"), labels = REGION_LEVELS) |> droplevels(),
      x          = factor(x,     levels = c("W", "B"), labels = c("Ambient", "Warm"))
    ) %>%
    dplyr::filter(!is.na(x), !is.na(Pad_Region))
  
  p_val <- get_interaction_pvalue(model)
  p_lab <- ifelse(is.na(p_val), "p = NA", paste0("p = ", signif(p_val, 4)))
  
  ggplot(pr, aes(x = .data$x, y = .data$predicted, group = .data$Pad_Region)) +
    geom_ribbon(
      aes(ymin = .data$conf.low, ymax = .data$conf.high, fill = .data$Pad_Region),
      alpha = 0.2, colour = NA
    ) +
    geom_line(aes(color = .data$Pad_Region), linewidth = 0.8) +
    geom_point(aes(color = .data$Pad_Region), size = 2) +
    labs(
      title    = date_label,
      subtitle = p_lab,
      x = "Pad colour",
      y = expression(Barnacle~abundance~"(log"[10]*" scale)")
    ) +
    theme_minimal() +
    theme(
      plot.title    = element_text(face = "bold", hjust = 0.5, size = 12, margin = margin(b = 1)),
      plot.subtitle = element_text(face = if (italicize_p) "italic" else "plain",
                                   hjust = 0.5, size = 9, margin = margin(b = 4)),
      plot.margin   = margin(5, 5, 5, 5),
      axis.title    = element_text(size = 10),
      axis.text     = element_text(size = 7),
      legend.title  = element_text(face = "bold", size = 9),
      legend.text   = element_text(size = 8),
      legend.position = if (include_legend) "right" else "none",
      legend.box.background = element_blank()
    ) +
    scale_y_continuous(trans = "log10", labels = scales::comma) +
    scale_color_manual(
      name   = "Pad region",
      breaks = REGION_LEVELS,
      values = REGION_VALUES
    ) +
    scale_fill_manual(
      name   = "Pad region",
      breaks = REGION_LEVELS,
      values = REGION_VALUES,
      guide  = "none"
    ) +
    guides(color = guide_legend(override.aes = list(shape = 16, size = 4, fill = NA)))
}

plot_JUN   <- create_conditional_plot(glm.nb_JUN,   "June")
plot_JUL_1 <- create_conditional_plot(glm.nb_JUL_1, "July 1")
plot_JUL_2 <- create_conditional_plot(glm.nb_JUL_2, "July 2")
plot_AUG_1 <- create_conditional_plot(glm.nb_AUG_1, "August 1")
plot_AUG_2 <- create_conditional_plot(glm.nb_AUG_2, "August 2")
plot_SEP_1 <- create_conditional_plot(glm.nb_SEP_1, "September 1")
plot_SEP_2 <- create_conditional_plot(glm.nb_SEP_2, "September 2", include_legend = TRUE)

# Keep just one legend for the panel plot
legs <- cowplot::get_plot_component(plot_SEP_2, "guide-box", return_all = TRUE)
legend_g <- legs[[1]]

panel_plot <- cowplot::ggdraw() +
  cowplot::draw_plot(plot_JUN,   x = 0.00, y = 0.55, width = 0.25, height = 0.45) +
  cowplot::draw_plot(plot_JUL_1, x = 0.25, y = 0.55, width = 0.25, height = 0.45) +
  cowplot::draw_plot(plot_JUL_2, x = 0.50, y = 0.55, width = 0.25, height = 0.45) +
  cowplot::draw_plot(plot_AUG_1, x = 0.75, y = 0.55, width = 0.25, height = 0.45) +
  cowplot::draw_plot(plot_AUG_2, x = 0.00, y = 0.05, width = 0.25, height = 0.45) +
  cowplot::draw_plot(plot_SEP_1, x = 0.25, y = 0.05, width = 0.25, height = 0.45) +
  cowplot::draw_plot(plot_SEP_2 + theme(legend.position = "none"),
                     x = 0.50, y = 0.05, width = 0.25, height = 0.45) +
  cowplot::draw_grob(legend_g, x = 0.75, y = 0.20, width = 0.25, height = 0.40)

print(panel_plot)

########## 6) Whole-season GAM ##########

# Use dplyr::select() explicitly to avoid select() masking
brighton_data_GAM <- brighton_data %>%
  dplyr::select(Sampling_Date, Pad_Colour, Pad_Region, Barnacle_Count, Pad_ID) %>%
  mutate(
    Sampling_Date    = as.Date(Sampling_Date),
    Days_Since_Start = as.numeric(Sampling_Date - min(Sampling_Date, na.rm = TRUE)),
    Pad_ID           = factor(Pad_ID)
  ) %>%
  filter(!is.na(Sampling_Date))

# Whole-season GAM (NB log link):
barnacle_gam_model <- mgcv::gam(
  Barnacle_Count ~
    s(Days_Since_Start, by = interaction(Pad_Colour, Pad_Region), k = 7) +
    Pad_Colour * Pad_Region +
    s(Pad_ID, bs = "re"),
  data   = brighton_data_GAM,
  family = nb(link = "log"),
  method = "REML"
)

summary(barnacle_gam_model)

########## 7) GAM interaction classification + summary table ##########

extract_barnacle_gam_model_stats <- function(model, date_label = "Whole Season") {
  
  # Same equivalence-band logic as per-date GLMs:
  # classify interaction on the link scale (log link), interpret exp(beta_int) as multiplicative interaction ratio.
  ptab <- as.data.frame(summary(model)$p.table) %>%
    tibble::rownames_to_column("term")
  
  # Normalise column names (robust across mgcv versions)
  names(ptab) <- tolower(gsub("[[:punct:] ]+", "_", names(ptab)))
  pcol <- grep("^pr|p_value", names(ptab), value = TRUE)[1]
  ptab <- dplyr::rename(ptab, p_value = all_of(pcol))
  
  # Term matching
  term_lc <- tolower(ptab$term)
  i_int   <- which(term_lc %in% "(intercept)")
  i_w     <- grep("pad[._ ]*colourb$", term_lc)
  i_p     <- grep("pad[._ ]*regionp$", term_lc)
  i_ix    <- grep("pad[._ ]*colourb[:*]pad[._ ]*regionp", term_lc)
  
  getv <- function(col, idx) if (length(idx)) ptab[[col]][idx[1]] else NA_real_
  
  alpha_hat <- getv("estimate",  i_int)
  beta_w    <- getv("estimate",  i_w)
  beta_p    <- getv("estimate",  i_p)
  beta_int  <- getv("estimate",  i_ix)
  
  se_w      <- getv("std_error", i_w)
  se_p      <- getv("std_error", i_p)
  se_i      <- getv("std_error", i_ix)
  p_int     <- getv("p_value",   i_ix)
  
  # Additive vs observed
  d_add_log <- beta_w + beta_p
  d_obs_log <- d_add_log + beta_int
  null_pred     <- exp(alpha_hat + d_add_log)
  observed_pred <- exp(alpha_hat + d_obs_log)
  
  pct_add <- (exp(d_add_log) - 1) * 100
  pct_obs <- (exp(d_obs_log) - 1) * 100
  
  warming_pct                 <- (exp(beta_w)           - 1) * 100
  warming_under_pollution_pct <- (exp(beta_w + beta_int) - 1) * 100
  pollution_pct               <- (exp(beta_p)           - 1) * 100
  pollution_under_warming_pct <- (exp(beta_p + beta_int) - 1) * 100
  
  int_ratio    <- exp(beta_int)
  int_ratio_lo <- exp(beta_int - 1.96 * se_i)
  int_ratio_hi <- exp(beta_int + 1.96 * se_i)
  
  int_class_point <- dplyr::case_when(
    is.na(int_ratio)                                      ~ NA_character_,
    abs(int_ratio - 1) < 0.05                             ~ "Additive (|Δobs–Δadd|<5%)",
    (d_add_log >  0 & d_obs_log <  0)                     ~ "Reversal (Pos to Neg)",
    (d_add_log <  0 & d_obs_log >  0)                     ~ "Reversal (Neg to Pos)",
    (beta_w > 0 & beta_p > 0 & beta_int >  0)             ~ "Synergism (>1)",
    (beta_w > 0 & beta_p > 0 & beta_int <  0)             ~ "Antagonism (<1)",
    (beta_w < 0 & beta_p < 0 & beta_int <  0)             ~ "Synergism (<1)",
    (beta_w < 0 & beta_p < 0 & beta_int >  0)             ~ "Antagonism (>1)",
    sign(d_add_log) == sign(beta_int)                     ~ "Synergism (same–sign)",
    sign(d_add_log) != sign(beta_int)                     ~ "Antagonism (opp–sign)",
    TRUE                                                  ~ "Unclassified"
  )
  
  ci_overlaps_band <- !(int_ratio_hi < EQUIV_LO || int_ratio_lo > EQUIV_HI)
  ci_outside_band  <-  (int_ratio_hi < EQUIV_LO || int_ratio_lo > EQUIV_HI)
  
  int_class_final <- dplyr::case_when(
    is.na(p_int) | is.na(int_ratio) | p_int >= P_ALPHA ~ "Unclassified (p≥0.05)",
    ci_overlaps_band                                   ~ "Additive (CI overlaps ±5%)",
    ci_outside_band & p_int < P_ALPHA                  ~ int_class_point,
    TRUE                                               ~ "Unclassified"
  )
  
  out <- tibble::tibble(
    Date                        = date_label,
    Term                        = c("Warming", "Pollution", "Warming × Pollution"),
    Estimate                    = c(beta_w, beta_p, beta_int),
    `Std. Error`                = c(se_w,   se_p,   se_i),
    `p-value`                   = c(
      ifelse(!is.na(se_w) & se_w > 0, 2 * pnorm(abs(beta_w / se_w), lower.tail = FALSE), NA_real_),
      ifelse(!is.na(se_p) & se_p > 0, 2 * pnorm(abs(beta_p / se_p), lower.tail = FALSE), NA_real_),
      p_int
    ),
    Pct_Add                     = c(NA, NA, ifelse(is.na(pct_add), NA, paste0(round(pct_add, 1), "%"))),
    Pct_Obs                     = c(NA, NA, ifelse(is.na(pct_obs), NA, paste0(round(pct_obs, 1), "%"))),
    `Pollution -> under Warming` = c(NA, NA,
                                     ifelse(any(is.na(c(pollution_pct, pollution_under_warming_pct))), NA,
                                            sprintf("Pollution’s %+0.1f%% to %+0.1f%%",
                                                    pollution_pct, pollution_under_warming_pct))),
    `Warming -> under Pollution` = c(NA, NA,
                                     ifelse(any(is.na(c(warming_pct, warming_under_pollution_pct))), NA,
                                            sprintf("Warming’s %+0.1f%% to %+0.1f%%",
                                                    warming_pct, warming_under_pollution_pct))),
    `Interaction Class`         = c(NA, NA, int_class_final),
    `Null Prediction`           = c(NA, NA, round(null_pred, 2)),
    `Observed Prediction`       = c(NA, NA, round(observed_pred, 2)),
    `Interaction Ratio`         = c(NA, NA, round(int_ratio, 3)),
    `Interaction Ratio CI`      = c(NA, NA,
                                    ifelse(any(is.na(c(int_ratio_lo, int_ratio_hi))), NA,
                                           sprintf("[%.3f, %.3f]", int_ratio_lo, int_ratio_hi)))
  )
  
  out$`p-value` <- ifelse(is.na(out$`p-value`), NA, format.pval(out$`p-value`, digits = 3, eps = 0.001))
  out
}

barnacle_gam_results_table <- extract_barnacle_gam_model_stats(barnacle_gam_model, "Whole Season")

DT::datatable(
  barnacle_gam_results_table,
  caption = "GAM (NB link): parametric fixed effects & interaction summary",
  options = list(pageLength = 10, scrollX = TRUE)
)

########## 8) Whole-season GAM visualisation ##########

# Map calendar dates to Days_Since_Start axis used by the GAM
start_date  <- min(brighton_data_GAM$Sampling_Date, na.rm = TRUE)
date_breaks <- as.numeric(as.Date(names(DATE_LOOKUP)) - start_date)

# Prediction grid across the full season
pred_grid <- tidyr::expand_grid(
  Days_Since_Start = seq(min(brighton_data_GAM$Days_Since_Start, na.rm = TRUE),
                         max(brighton_data_GAM$Days_Since_Start, na.rm = TRUE),
                         length.out = 300),
  Pad_Colour = factor(c("W", "B"), levels = c("W", "B")),
  Pad_Region = factor(c("N", "P"), levels = c("N", "P"))
)

# Because the model contains s(Pad_ID), predict.gam requires a Pad_ID in newdata
pad_levels <- levels(brighton_data_GAM$Pad_ID)

pred_grid <- pred_grid %>%
  mutate(Pad_ID = factor(if (length(pad_levels)) pad_levels[1] else NA_character_,
                         levels = pad_levels))

# Predict on the link scale and back-transform to counts
pr_link <- predict(
  barnacle_gam_model,
  newdata  = pred_grid,
  type     = "link",
  se.fit   = TRUE,
  exclude  = "s(Pad_ID)",
  newdata.guaranteed = TRUE
)

gam_preds <- pred_grid %>%
  mutate(
    fit       = pr_link$fit,
    se        = pr_link$se.fit,
    predicted = exp(fit),
    conf.low  = exp(fit - 1.96 * se),
    conf.high = exp(fit + 1.96 * se),
    Pad_Region_lbl = factor(Pad_Region, levels = c("N", "P"), labels = c("Non-Polluted", "Polluted")),
    Pad_Colour_lbl = factor(Pad_Colour, levels = c("W", "B"), labels = c("Ambient", "Warm")),
    Treatment      = paste(Pad_Region_lbl, Pad_Colour_lbl, sep = " + ")
  )

# Date-point predictions (at sampling dates only)
datepoint_grid <- tidyr::expand_grid(
  Days_Since_Start = date_breaks,
  Pad_Colour = factor(c("W", "B"), levels = c("W", "B")),
  Pad_Region = factor(c("N", "P"), levels = c("N", "P"))
) %>%
  mutate(Pad_ID = factor(if (length(pad_levels)) pad_levels[1] else NA_character_,
                         levels = pad_levels))

dp_link <- predict(
  barnacle_gam_model,
  newdata  = datepoint_grid,
  type     = "link",
  se.fit   = TRUE,
  exclude  = "s(Pad_ID)",
  newdata.guaranteed = TRUE
)

gam_datepoints <- datepoint_grid %>%
  mutate(
    fit       = dp_link$fit,
    se        = dp_link$se.fit,
    predicted = exp(fit),
    Pad_Region_lbl = factor(Pad_Region, levels = c("N", "P"), labels = c("Non-Polluted", "Polluted")),
    Pad_Colour_lbl = factor(Pad_Colour, levels = c("W", "B"), labels = c("Ambient", "Warm")),
    Treatment      = paste(Pad_Region_lbl, Pad_Colour_lbl, sep = " + ")
  )

brighton_data_GAM_plot <- brighton_data_GAM %>%
  mutate(
    Pad_Region_lbl = factor(Pad_Region, levels = c("N", "P"), labels = c("Non-Polluted", "Polluted")),
    Pad_Colour_lbl = factor(Pad_Colour, levels = c("W", "B"), labels = c("Ambient", "Warm")),
    Treatment      = paste(Pad_Region_lbl, Pad_Colour_lbl, sep = " + ")
  )

# Whole-season GAM plot
gam_time_plot <- ggplot(gam_preds, aes(x = Days_Since_Start, y = predicted, group = Treatment)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = Treatment), alpha = 0.15, color = NA) +
  geom_line(aes(color = Treatment), linewidth = 0.7) +
  geom_point(
    data = gam_datepoints,
    aes(x = Days_Since_Start, y = predicted, fill = Treatment,
        colour = after_scale(colorspace::darken(fill, 0.35))),
    shape = 21, size = 2.7, stroke = 0.9, show.legend = FALSE
  ) +
  geom_jitter(
    data = brighton_data_GAM_plot,
    aes(x = Days_Since_Start, y = Barnacle_Count, color = Treatment),
    shape = 20, width = 0.2, alpha = 0.7, size = 1.5, inherit.aes = FALSE
  ) +
  scale_x_continuous(
    breaks = date_breaks,
    labels = DATE_LABELS,
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  scale_y_continuous(trans = "log10", labels = scales::comma) +
  scale_color_manual(values = TREATMENT_COLORS, breaks = LEGEND_ORDER) +
  scale_fill_manual(values = TREATMENT_COLORS,  breaks = LEGEND_ORDER) +
  labs(
    title = "Whole-season GAM: Barnacle abundance",
    x = "Sampling date",
    y = expression(Barnacle~abundance~"(log"[10]*" scale)")
  ) +
  theme_minimal() +
  theme(
    plot.title   = element_text(face = "bold", hjust = 0.5, size = 12),
    axis.title   = element_text(size = 10),
    axis.text    = element_text(size = 7),
    legend.title = element_text(face = "bold", size = 9),
    legend.text  = element_text(size = 8),
    legend.box.background = element_blank(),
    axis.title.x = element_text(margin = margin(t = 10))
  )

print(gam_time_plot)