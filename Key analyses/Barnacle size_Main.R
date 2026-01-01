# Barnacle size_Main.R
#
#   Brighton experiment (2×2): Warming (Pad_Colour) × Pollution (Pad_Region)
#   Response: AVG_Barnacle_Size (mm; continuous)
#
#   Outputs:
#     1) Per-date Gaussian GLMs (identity)
#     2) GLM interaction classification table
#     3) Whole-summer GLM-style visualisation (ggeffects predictions + raw overlay)
#     4) Conditional per-date panels (pollution-stratified lines/ribbons across pad colour)
#     5) Whole-season Gaussian GAM (mgcv::gam) with time smooths varying by Region only
#     6) GAM interaction classification table (parametric fixed effects)
#     7) Whole-season GAM visualisation (smooth predictions + raw overlay)
#
#   Evidence/diagnostics (residual checks, alternative structures) are kept in:
#     Barnacle size_Additional.R

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
#   - For size (identity link), beta_int is in mm: the interaction is a deviation from additivity on the mm scale.
#   - We treat interactions as "additive" if the interaction CI overlaps a ±5% band around the additive (null) prediction.
P_ALPHA  <- 0.05
EQUIV_PC <- 0.05  # ±5%

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
  library(ggeffects)
  library(colorspace)
  library(broom)
  library(cowplot)
  library(mgcv)
  library(DT)
  library(scales)
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
    Pad_Region = factor(Pad_Region, levels = c("N", "P")),
    Pad_ID     = factor(Pad_ID)
  )

# Helper: ensure correct reference levels within any subset
prep_subset <- function(df) {
  df$Pad_Colour <- relevel(as.factor(df$Pad_Colour), ref = "W")
  df$Pad_Region <- relevel(as.factor(df$Pad_Region), ref = "N")
  df
}

########## 2) Individual-date GLMs (Gaussian identity) ##########

# Subsets per date
JUN   <- prep_subset(dplyr::filter(brighton_data, Sampling_Date == D_JUN))
JUL_1 <- prep_subset(dplyr::filter(brighton_data, Sampling_Date == D_JUL1))
JUL_2 <- prep_subset(dplyr::filter(brighton_data, Sampling_Date == D_JUL2))
AUG_1 <- prep_subset(dplyr::filter(brighton_data, Sampling_Date == D_AUG1))
AUG_2 <- prep_subset(dplyr::filter(brighton_data, Sampling_Date == D_AUG2))
SEP_1 <- prep_subset(dplyr::filter(brighton_data, Sampling_Date == D_SEP1))
SEP_2 <- prep_subset(dplyr::filter(brighton_data, Sampling_Date == D_SEP2))

# Final per-date models:
#   Identity link (mm scale): coefficients are in mm differences from the baseline cell.
glm_JUN   <- glm(AVG_Barnacle_Size ~ Pad_Colour * Pad_Region, data = JUN,   family = gaussian(link = "identity"))
glm_JUL_1 <- glm(AVG_Barnacle_Size ~ Pad_Colour * Pad_Region, data = JUL_1, family = gaussian(link = "identity"))
glm_JUL_2 <- glm(AVG_Barnacle_Size ~ Pad_Colour * Pad_Region, data = JUL_2, family = gaussian(link = "identity"))
glm_AUG_1 <- glm(AVG_Barnacle_Size ~ Pad_Colour * Pad_Region, data = AUG_1, family = gaussian(link = "identity"))
glm_AUG_2 <- glm(AVG_Barnacle_Size ~ Pad_Colour * Pad_Region, data = AUG_2, family = gaussian(link = "identity"))
glm_SEP_1 <- glm(AVG_Barnacle_Size ~ Pad_Colour * Pad_Region, data = SEP_1, family = gaussian(link = "identity"))
glm_SEP_2 <- glm(AVG_Barnacle_Size ~ Pad_Colour * Pad_Region, data = SEP_2, family = gaussian(link = "identity"))

glm_models <- list(glm_JUN, glm_JUL_1, glm_JUL_2, glm_AUG_1, glm_AUG_2, glm_SEP_1, glm_SEP_2)

########## 3) GLM interaction classification + summary table ##########

extract_barnacle_size_glm_stats <- function(model, date_label) {
  
  # Coefficient table (identity scale, mm)
  #   α    = baseline mean size (Ambient + Non-Polluted)
  #   βw   = warming effect in Non-Polluted
  #   βp   = pollution effect under Ambient
  #   βint = interaction (mm deviation beyond additivity)
  fe <- broom::tidy(model, conf.int = TRUE, conf.level = 0.95)
  
  # If conf.int was not returned, construct Wald intervals
  if (!all(c("conf.low", "conf.high") %in% names(fe))) {
    fe <- fe %>%
      mutate(
        conf.low  = estimate - 1.96 * std.error,
        conf.high = estimate + 1.96 * std.error
      )
  }
  
  ic  <- dplyr::filter(fe, term == "(Intercept)")
  w   <- dplyr::filter(fe, term == "Pad_ColourB")
  p   <- dplyr::filter(fe, term == "Pad_RegionP")
  int <- dplyr::filter(fe, term == "Pad_ColourB:Pad_RegionP")
  
  α     <- ic$estimate
  βw    <- w$estimate
  βp    <- p$estimate
  βint  <- int$estimate
  seint <- int$std.error
  pint  <- int$p.value
  
  # Interaction p-value (fallback to Wald z if missing)
  if (is.null(pint) || is.na(pint)) pint <- 2 * pnorm(abs(βint / seint), lower.tail = FALSE)
  
  # Additive (null) vs observed (identity scale, mm)
  Δadd      <- βw + βp
  Δobs      <- Δadd + βint
  null_pred <- α + Δadd
  obs_pred  <- α + Δobs
  
  # ±5% equivalence band around the additive (null) prediction.
  # We treat the interaction as "additive" if the CI for βint overlaps this band.
  tol     <- EQUIV_PC * abs(null_pred)
  band_lo <- -tol
  band_hi <-  tol
  
  # 95% CI for βint (mm) from broom
  int_lo <- int$conf.low
  int_hi <- int$conf.high
  
  # Direction label (only used if CI fully outside band AND p < alpha)
  dir_label <- dplyr::case_when(
    (Δadd < 0 & Δobs > 0)               ~ "Reversal (Neg to Pos)",
    (Δadd > 0 & Δobs < 0)               ~ "Reversal (Pos to Neg)",
    (Δadd > 0 &  βint > 0)              ~ "Synergism",
    (Δadd > 0 &  βint < 0)              ~ "Antagonism",
    (Δadd < 0 &  βint < 0)              ~ "Synergism",
    (Δadd < 0 &  βint > 0)              ~ "Antagonism",
    sign(Δadd) == sign(βint)            ~ "Synergism (same–sign)",
    TRUE                                ~ "Antagonism (opp–sign)"
  )
  
  # Equivalence-band logic based on CI for βint (mm)
  ci_overlaps_band <- !(int_hi < band_lo || int_lo > band_hi)
  ci_outside_band  <-  (int_hi < band_lo || int_lo > band_hi)
  
  # Final interaction class (p-value gate first; then CI band)
  int_class <- dplyr::case_when(
    is.na(pint) | pint >= P_ALPHA       ~ "Unclassified (p≥0.05)",
    ci_overlaps_band                    ~ "Additive (CI overlaps ±5% band)",
    ci_outside_band & pint < P_ALPHA    ~ dir_label,
    TRUE                                ~ "Unclassified"
  )
  
  # Perspective metrics under each stressor stratum
  warming_under_pollution <- βw + βint
  pollution_under_warming <- βp + βint
  
  # Assemble table rows (one row per term, with interaction-only columns populated)
  fe %>%
    dplyr::filter(term %in% c("Pad_ColourB", "Pad_RegionP", "Pad_ColourB:Pad_RegionP")) %>%
    dplyr::transmute(
      Date      = date_label,
      Term      = dplyr::case_when(
        term == "Pad_ColourB"             ~ "Warming",
        term == "Pad_RegionP"             ~ "Pollution",
        term == "Pad_ColourB:Pad_RegionP" ~ "Warming × Pollution",
        TRUE                              ~ term
      ),
      Estimate     = round(estimate, 4),
      `Std. Error` = round(std.error, 4),
      `p-value`    = format.pval(p.value, digits = 3, eps = 0.001),
      `ΔAdd (mm)`                = if_else(Term == "Warming × Pollution", round(Δadd, 4), NA_real_),
      `ΔObs (mm)`                = if_else(Term == "Warming × Pollution", round(Δobs, 4), NA_real_),
      `Null Prediction (mm)`     = if_else(Term == "Warming × Pollution", round(null_pred, 4), NA_real_),
      `Observed Prediction (mm)` = if_else(Term == "Warming × Pollution", round(obs_pred, 4), NA_real_),
      `Interaction 95% CI (mm)`  = if_else(Term == "Warming × Pollution", sprintf("[%.4f, %.4f]", int_lo, int_hi), NA_character_),
      `±5% Band (mm)`            = if_else(Term == "Warming × Pollution", sprintf("[%.4f, %.4f]", band_lo, band_hi), NA_character_),
      `Interaction Class`        = if_else(Term == "Warming × Pollution", int_class, NA_character_),
      `Warming -> under Pollution (mm)` =
        if_else(Term == "Warming × Pollution", sprintf("%+.4f to %+.4f", βw, warming_under_pollution), NA_character_),
      `Pollution -> under Warming (mm)` =
        if_else(Term == "Warming × Pollution", sprintf("%+.4f to %+.4f", βp, pollution_under_warming), NA_character_)
    )
}

barnacle_size_glm_table <- purrr::map2_df(glm_models, DATE_LABELS, extract_barnacle_size_glm_stats)

DT::datatable(
  barnacle_size_glm_table,
  caption = "Barnacle size (Gaussian GLM, identity): fixed effects & interaction summary",
  options = list(pageLength = 10, scrollX = TRUE)
)

########## 4) Whole-summer GLM-style visualisation ##########

# Predictions per date
pred_list <- vector("list", length(glm_models))

for (i in seq_along(glm_models)) {
  
  pr <- ggeffects::ggpredict(glm_models[[i]], terms = c("Pad_Colour", "Pad_Region")) %>%
    dplyr::mutate(
      Sampling_Date = DATE_LABELS[i],
      Pad_Region    = factor(group, levels = c("N", "P"), labels = c("Non-Polluted", "Polluted")),
      Pad_Colour    = factor(x,     levels = c("W", "B"), labels = c("Ambient", "Warm")),
      Treatment     = paste(Pad_Region, Pad_Colour, sep = " + ")
    )
  
  pred_list[[i]] <- pr
}

all_preds <- bind_rows(pred_list) %>%
  filter(!is.na(Treatment))

all_preds$Sampling_Date <- factor(all_preds$Sampling_Date, levels = DATE_LABELS)

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
barnacle_size_line_plot <- ggplot(all_preds, aes(x = Sampling_Date, y = predicted, group = Treatment)) +
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
    aes(x = Sampling_Date, y = AVG_Barnacle_Size, color = Treatment),
    shape = 20, width = 0.25, height = 0, alpha = 0.7, size = 1.5, inherit.aes = FALSE
  ) +
  labs(
    title = "Barnacle size",
    x = "Sampling date",
    y = "Barnacle size (mm)"
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
  scale_color_manual(values = TREATMENT_COLORS, breaks = LEGEND_ORDER) +
  scale_fill_manual(values = TREATMENT_COLORS, breaks = LEGEND_ORDER)

print(barnacle_size_line_plot)

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
      y = "Barnacle size (mm)"
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

plot_JUN   <- create_conditional_plot(glm_JUN,   "June")
plot_JUL_1 <- create_conditional_plot(glm_JUL_1, "July 1")
plot_JUL_2 <- create_conditional_plot(glm_JUL_2, "July 2")
plot_AUG_1 <- create_conditional_plot(glm_AUG_1, "August 1")
plot_AUG_2 <- create_conditional_plot(glm_AUG_2, "August 2")
plot_SEP_1 <- create_conditional_plot(glm_SEP_1, "September 1")
plot_SEP_2 <- create_conditional_plot(glm_SEP_2, "September 2", include_legend = TRUE)

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
  dplyr::select(Sampling_Date, Pad_Colour, Pad_Region, AVG_Barnacle_Size, Pad_ID) %>%
  mutate(
    Sampling_Date    = as.Date(Sampling_Date),
    Days_Since_Start = as.numeric(Sampling_Date - min(Sampling_Date, na.rm = TRUE)),
    Pad_ID           = factor(Pad_ID)
  ) %>%
  filter(!is.na(Sampling_Date))

# Whole-season GAM (Gaussian identity):
barnacle_size_gam <- mgcv::gam(
  AVG_Barnacle_Size ~
    s(Days_Since_Start, by = interaction(Pad_Colour, Pad_Region), k = 7) +
    Pad_Colour * Pad_Region +
    s(Pad_ID, bs = "re"),
  data   = brighton_data_GAM,
  family = gaussian(link = "identity"),
  method = "REML"
)

print(summary(barnacle_size_gam))

########## 7) GAM parametric effects & interaction summary ##########

extract_barnacle_size_gam_stats <- function(model, date_label = "Whole Season") {
  
  # Same equivalence-band logic as per-date GLMs:
  # For size (identity link), treat beta_int (mm) as "additive" if its CI overlaps
  # a ±5% band around the additive (null) prediction on the mm scale.
  ptab <- as.data.frame(summary(model)$p.table) %>%
    tibble::rownames_to_column("term")
  
  # Normalise column names (robust across mgcv versions)
  names(ptab) <- tolower(gsub("[[:punct:] ]+", "_", names(ptab)))
  pcol <- grep("^pr|p_value", names(ptab), value = TRUE)[1]
  ptab <- dplyr::rename(ptab, p_value = all_of(pcol))
  
  # Term matching
  term_lc <- tolower(ptab$term)
  i_int <- which(term_lc %in% "(intercept)")
  i_w   <- grep("pad[._ ]*colourb$", term_lc)
  i_p   <- grep("pad[._ ]*regionp$", term_lc)
  i_ix  <- grep("pad[._ ]*colourb[:*]pad[._ ]*regionp", term_lc)
  
  gv <- function(col, idx) if (length(idx)) ptab[[col]][idx[1]] else NA_real_
  
  α     <- gv("estimate",  i_int)
  βw    <- gv("estimate",  i_w)
  βp    <- gv("estimate",  i_p)
  βint  <- gv("estimate",  i_ix)
  
  se_w  <- gv("std_error", i_w)
  se_p  <- gv("std_error", i_p)
  se_i  <- gv("std_error", i_ix)
  
  pv_w  <- gv("p_value",   i_w)
  pv_p  <- gv("p_value",   i_p)
  pv_i  <- gv("p_value",   i_ix)
  
  # Additive vs observed (identity scale, mm)
  Δadd      <- βw + βp
  Δobs      <- Δadd + βint
  null_pred <- α + Δadd
  obs_pred  <- α + Δobs
  
  # Wald CI for beta_int (mm)
  int_lo <- βint - 1.96 * se_i
  int_hi <- βint + 1.96 * se_i
  
  # ±5% equivalence band around the additive (null) prediction
  tol     <- EQUIV_PC * abs(null_pred)
  band_lo <- -tol
  band_hi <-  tol
  
  dir_label <- dplyr::case_when(
    (Δadd < 0 & Δobs > 0)               ~ "Reversal (Neg to Pos)",
    (Δadd > 0 & Δobs < 0)               ~ "Reversal (Pos to Neg)",
    (Δadd > 0 &  βint > 0)              ~ "Synergism",
    (Δadd > 0 &  βint < 0)              ~ "Antagonism",
    (Δadd < 0 &  βint < 0)              ~ "Synergism",
    (Δadd < 0 &  βint > 0)              ~ "Antagonism",
    sign(Δadd) == sign(βint)            ~ "Synergism (same–sign)",
    TRUE                                ~ "Antagonism (opp–sign)"
  )
  
  ci_overlaps_band <- !(int_hi < band_lo || int_lo > band_hi)
  ci_outside_band  <-  (int_hi < band_lo || int_lo > band_hi)
  
  int_class <- dplyr::case_when(
    is.na(pv_i) | pv_i >= P_ALPHA       ~ "Unclassified (p≥0.05)",
    ci_overlaps_band                    ~ "Additive (CI overlaps ±5% band)",
    ci_outside_band & pv_i < P_ALPHA    ~ dir_label,
    TRUE                                ~ "Unclassified"
  )
  
  out <- tibble::tibble(
    Date                        = date_label,
    Term                        = c("Warming", "Pollution", "Warming × Pollution"),
    Estimate                    = c(βw, βp, βint),
    `Std. Error`                = c(se_w, se_p, se_i),
    `p-value`                   = c(pv_w, pv_p, pv_i),
    `ΔAdd (mm)`                 = c(NA, NA, Δadd),
    `ΔObs (mm)`                 = c(NA, NA, Δobs),
    `Null Prediction (mm)`      = c(NA, NA, null_pred),
    `Observed Prediction (mm)`  = c(NA, NA, obs_pred),
    `Interaction 95% CI (mm)`   = c(NA, NA, ifelse(any(is.na(c(int_lo, int_hi))), NA, sprintf("[%.4f, %.4f]", int_lo, int_hi))),
    `±5% Band (mm)`             = c(NA, NA, ifelse(is.na(tol), NA, sprintf("[%.4f, %.4f]", band_lo, band_hi))),
    `Interaction Class`         = c(NA, NA, int_class)
  )
  
  # Formatting
  out$Estimate     <- round(out$Estimate, 4)
  out$`Std. Error` <- round(out$`Std. Error`, 4)
  out$`p-value`    <- ifelse(is.na(out$`p-value`), NA, format.pval(out$`p-value`, digits = 3, eps = 0.001))
  
  out$`ΔAdd (mm)`                <- ifelse(is.na(out$`ΔAdd (mm)`), NA, round(out$`ΔAdd (mm)`, 4))
  out$`ΔObs (mm)`                <- ifelse(is.na(out$`ΔObs (mm)`), NA, round(out$`ΔObs (mm)`, 4))
  out$`Null Prediction (mm)`     <- ifelse(is.na(out$`Null Prediction (mm)`), NA, round(out$`Null Prediction (mm)`, 4))
  out$`Observed Prediction (mm)` <- ifelse(is.na(out$`Observed Prediction (mm)`), NA, round(out$`Observed Prediction (mm)`, 4))
  
  out
}

barnacle_size_gam_table <- extract_barnacle_size_gam_stats(barnacle_size_gam, "Whole Season")

DT::datatable(
  barnacle_size_gam_table,
  caption = "Barnacle size GAM (Gaussian identity): parametric effects & interaction summary",
  options = list(pageLength = 10, scrollX = TRUE)
)

########## 8) Whole-season GAM visualisation ##########

# Map calendar dates to Days_Since_Start axis used by the GAM
start_date  <- min(brighton_data_GAM$Sampling_Date, na.rm = TRUE)
date_breaks <- as.numeric(as.Date(names(DATE_LOOKUP)) - start_date)

# Because the model contains s(Pad_ID), predict.gam requires a Pad_ID in newdata
pad_levels <- levels(brighton_data_GAM$Pad_ID)
dummy_pad  <- if (length(pad_levels)) pad_levels[1] else NA_character_

# Prediction grid across the full season
pred_grid <- tidyr::expand_grid(
  Days_Since_Start = seq(min(brighton_data_GAM$Days_Since_Start, na.rm = TRUE),
                         max(brighton_data_GAM$Days_Since_Start, na.rm = TRUE),
                         length.out = 300),
  Pad_Colour = factor(c("W", "B"), levels = c("W", "B")),
  Pad_Region = factor(c("N", "P"), levels = c("N", "P"))
) %>%
  mutate(Pad_ID = factor(dummy_pad, levels = pad_levels))

# Predict on the link scale (identity here, so link == response)
pr_link <- predict(
  barnacle_size_gam,
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
    predicted = fit,
    conf.low  = fit - 1.96 * se,
    conf.high = fit + 1.96 * se,
    Pad_Region_lbl = factor(Pad_Region, levels = c("N", "P"), labels = c("Non-Polluted", "Polluted")),
    Pad_Colour_lbl = factor(Pad_Colour, levels = c("W", "B"), labels = c("Ambient", "Warm")),
    Treatment      = paste(Pad_Region_lbl, Pad_Colour_lbl, sep = " + ")
  )

# Date-point predictions at sampling dates
datepoint_grid <- tidyr::expand_grid(
  Days_Since_Start = date_breaks,
  Pad_Colour = factor(c("W", "B"), levels = c("W", "B")),
  Pad_Region = factor(c("N", "P"), levels = c("N", "P"))
) %>%
  mutate(Pad_ID = factor(dummy_pad, levels = pad_levels))

dp_link <- predict(
  barnacle_size_gam,
  newdata  = datepoint_grid,
  type     = "link",
  se.fit   = TRUE,
  exclude  = "s(Pad_ID)",
  newdata.guaranteed = TRUE
)

gam_datepoints <- datepoint_grid %>%
  mutate(
    predicted = dp_link$fit,
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

# Whole-season GAM plot (identity scale)
barnacle_size_gam_timeplot <- ggplot(gam_preds, aes(x = Days_Since_Start, y = predicted, group = Treatment)) +
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
    aes(x = Days_Since_Start, y = AVG_Barnacle_Size, color = Treatment),
    shape = 20, width = 0.8, height = 0, alpha = 0.7, size = 1.5, inherit.aes = FALSE
  ) +
  scale_x_continuous(breaks = date_breaks, labels = DATE_LABELS, expand = expansion(mult = c(0.02, 0.02))) +
  scale_color_manual(values = TREATMENT_COLORS, breaks = LEGEND_ORDER) +
  scale_fill_manual(values = TREATMENT_COLORS,  breaks = LEGEND_ORDER) +
  labs(
    title = "Whole-season GAM: Barnacle size",
    x = "Sampling date",
    y = "Barnacle size (mm)"
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

print(barnacle_size_gam_timeplot)