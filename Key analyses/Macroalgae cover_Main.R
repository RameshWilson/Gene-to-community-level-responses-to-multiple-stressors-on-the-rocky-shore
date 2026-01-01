# Macroalgae cover_Main.R
#
#   Brighton experiment (2×2): Warming (Pad_Colour) × Pollution (Pad_Region)
#   Response: macroalgae cover proportion (Algae_Cover / 100), modelled with Beta regression / Beta GAM (logit)
#
#   Outputs:
#     1) Per-date Beta GLMs (betareg; logit)
#     2) GLM interaction classification table (±5% equivalence band on exp(beta_int); odds scale)
#     3) Whole-summer GLM-style visualisation (ggeffects predictions + raw overlay)
#     4) Conditional per-date panels (pollution-stratified lines/ribbons across pad colour)
#     5) Whole-season Beta GAM (mgcv::gam family=betar(link="logit")) + interaction table
#     6) Whole-season GAM visualisation (smooth predictions + raw overlay)
#
#   Evidence/diagnostics (0/1 inflation handling, alternative families/structures) live in:
#     Macroalgae cover_Additional.R

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
#   - Beta-logit models operate on the logit(mean cover proportion); exp(beta_int) is an 'interaction odds ratio'
#   for the mean cover (not an OR of binary occurrence).
P_ALPHA  <- 0.05
EQUIV_LO <- 0.95
EQUIV_HI <- 1.05

# Beta regression requires y in (0,1); clamp any 0/1 values to (EPS, 1-EPS).
EPS <- 1e-4

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
  library(betareg)
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

# Macroalgae cover as proportion; clamp to (0,1) for Beta models
brighton_data <- brighton_data %>%
  mutate(
    Algae_Cover_prop = pmax(EPS, pmin(1 - EPS, Algae_Cover / 100))
  )

# Helper: ensure correct reference levels within any subset
prep_subset <- function(df) {
  df$Pad_Colour <- stats::relevel(as.factor(df$Pad_Colour), ref = "W")
  df$Pad_Region <- stats::relevel(as.factor(df$Pad_Region), ref = "N")
  df
}

########## 2) Individual-date GLMs (Beta–logit) ##########

# Subsets per date
JUN   <- prep_subset(dplyr::filter(brighton_data, Sampling_Date == D_JUN))
JUL_1 <- prep_subset(dplyr::filter(brighton_data, Sampling_Date == D_JUL1))
JUL_2 <- prep_subset(dplyr::filter(brighton_data, Sampling_Date == D_JUL2))
AUG_1 <- prep_subset(dplyr::filter(brighton_data, Sampling_Date == D_AUG1))
AUG_2 <- prep_subset(dplyr::filter(brighton_data, Sampling_Date == D_AUG2))
SEP_1 <- prep_subset(dplyr::filter(brighton_data, Sampling_Date == D_SEP1))
SEP_2 <- prep_subset(dplyr::filter(brighton_data, Sampling_Date == D_SEP2))

# Final per-date models:
#   - Beta regressions (0-1)
glm_JUN_beta   <- betareg::betareg(Algae_Cover_prop ~ Pad_Colour * Pad_Region, data = JUN,   link = "logit")
glm_JUL_1_beta <- betareg::betareg(Algae_Cover_prop ~ Pad_Colour * Pad_Region, data = JUL_1, link = "logit")
glm_JUL_2_beta <- betareg::betareg(Algae_Cover_prop ~ Pad_Colour * Pad_Region, data = JUL_2, link = "logit")
glm_AUG_1_beta <- betareg::betareg(Algae_Cover_prop ~ Pad_Colour * Pad_Region, data = AUG_1, link = "logit")
glm_AUG_2_beta <- betareg::betareg(Algae_Cover_prop ~ Pad_Colour * Pad_Region, data = AUG_2, link = "logit")
glm_SEP_1_beta <- betareg::betareg(Algae_Cover_prop ~ Pad_Colour * Pad_Region, data = SEP_1, link = "logit")
glm_SEP_2_beta <- betareg::betareg(Algae_Cover_prop ~ Pad_Colour * Pad_Region, data = SEP_2, link = "logit")

glm_models <- list(glm_JUN_beta, glm_JUL_1_beta, glm_JUL_2_beta, glm_AUG_1_beta, glm_AUG_2_beta, glm_SEP_1_beta, glm_SEP_2_beta)

########## 3) GLM interaction classification + summary table ##########

extract_macro_glm_stats <- function(model, date_label) {
  
  # betareg has mean + precision components; we classify interactions on the mean (conditional) component only
  cs <- broom::tidy(model, conf.int = TRUE, conf.level = 0.95)
  if ("component" %in% names(cs)) cs <- dplyr::filter(cs, component == "mean")
  
  # If conf.int wasn't returned, construct Wald intervals
  if (!all(c("conf.low", "conf.high") %in% names(cs))) {
    cs <- cs %>%
      mutate(
        conf.low  = estimate - 1.96 * std.error,
        conf.high = estimate + 1.96 * std.error
      )
  }
  
  ic  <- dplyr::filter(cs, term == "(Intercept)")
  w   <- dplyr::filter(cs, term == "Pad_ColourB")
  p   <- dplyr::filter(cs, term == "Pad_RegionP")
  int <- dplyr::filter(cs, term == "Pad_ColourB:Pad_RegionP")
  
  # Coefficient table (logit scale):
  #   α    = baseline OR (Ambient + Non-Polluted)
  #   βw   = warming effect (OR change) in Non-Polluted
  #   βp   = pollution effect (OR change) under Ambient
  #   βint = interaction (OR deviation beyond additivity)
  α     <- ic$estimate
  βw    <- w$estimate
  βp    <- p$estimate
  βint  <- int$estimate
  pint  <- int$p.value
  
  # Additive vs observed (logit scale)
  Δadd <- βw + βp
  Δobs <- Δadd + βint
  
  # Back-transform to mean cover proportion
  null_prop <- plogis(α + Δadd)
  obs_prop  <- plogis(α + Δobs)
  
  # Interaction odds ratio on logit(mean proportion) scale + CI
  int_or    <- exp(βint)
  int_or_lo <- exp(int$conf.low)
  int_or_hi <- exp(int$conf.high)
  
  #   # Perspective metrics under each stressor stratum (% change in odds of mean cover proportion)
  warming_pct                 <- (exp(βw)      - 1) * 100
  warming_under_pollution_pct <- (exp(βw + βint) - 1) * 100
  pollution_pct               <- (exp(βp)      - 1) * 100
  pollution_under_warming_pct <- (exp(βp + βint) - 1) * 100
  
  # Equivalence-band logic on exp(beta_int)
  ci_overlaps_band <- !(int_or_hi < EQUIV_LO || int_or_lo > EQUIV_HI)
  ci_outside_band  <-  (int_or_hi < EQUIV_LO || int_or_lo > EQUIV_HI)
  
  # Direction labels (only used if CI fully outside band AND p < alpha)
  dir_label <- dplyr::case_when(
    (Δadd > 0 & Δobs < 0) ~ "Reversal (Pos to Neg)",
    (Δadd < 0 & Δobs > 0) ~ "Reversal (Neg to Pos)",
    (βw > 0 & βp > 0 & βint > 0) ~ "Synergism (>1)",
    (βw > 0 & βp > 0 & βint < 0) ~ "Antagonism (<1)",
    (βw < 0 & βp < 0 & βint < 0) ~ "Synergism (<1)",
    (βw < 0 & βp < 0 & βint > 0) ~ "Antagonism (>1)",
    sign(Δadd) == sign(βint) ~ "Synergism (same–sign)",
    TRUE ~ "Antagonism (opp–sign)"
  )
  
  # Final interaction class (p-value gate first; then CI band)
  int_class <- dplyr::case_when(
    is.na(pint) | pint >= P_ALPHA        ~ "Unclassified (p≥0.05)",
    ci_overlaps_band                     ~ "Additive (CI overlaps ±5%)",
    ci_outside_band & pint < P_ALPHA     ~ dir_label,
    TRUE                                 ~ "Unclassified"
  )
  
  # Assemble table rows (one row per term, with interaction-only columns populated)
  cs %>%
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
      `ΔAdd (logit)`        = if_else(Term == "Warming × Pollution", round(Δadd, 4), NA_real_),
      `ΔObs (logit)`        = if_else(Term == "Warming × Pollution", round(Δobs, 4), NA_real_),
      `Null Proportion`     = if_else(Term == "Warming × Pollution", round(null_prop, 4), NA_real_),
      `Observed Proportion` = if_else(Term == "Warming × Pollution", round(obs_prop,  4), NA_real_),
      `Odds Ratio`          = if_else(Term == "Warming × Pollution", round(int_or, 4), NA_real_),
      `Interaction Ratio CI` = if_else(
        Term == "Warming × Pollution",
        sprintf("[%.3f, %.3f]", int_or_lo, int_or_hi),
        NA_character_
      ),
      `Interaction Class` = if_else(Term == "Warming × Pollution", int_class, NA_character_),
      `Warming -> under Pollution` = if_else(
        Term == "Warming × Pollution",
        sprintf("%+.1f%% to %+.1f%%", warming_pct, warming_under_pollution_pct),
        NA_character_
      ),
      `Pollution -> under Warming` = if_else(
        Term == "Warming × Pollution",
        sprintf("%+.1f%% to %+.1f%%", pollution_pct, pollution_under_warming_pct),
        NA_character_
      )
    )
}

macro_glm_table <- purrr::map2_df(glm_models, DATE_LABELS, extract_macro_glm_stats)

DT::datatable(
  macro_glm_table,
  caption = "Macroalgae cover (Beta GLM, logit link): fixed effects & interaction summary",
  options = list(pageLength = 10, scrollX = TRUE)
)

########## 4) Whole-summer GLM-style visualisation ##########

# Predictions per date
pred_list <- vector("list", length(glm_models))

for (i in seq_along(glm_models)) {
  
  pr <- ggeffects::ggpredict(glm_models[[i]], terms = c("Pad_Colour", "Pad_Region")) %>%
    mutate(
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

# Upper y-limit (data-driven; cap to keep plots readable)
y_upper_glm <- min(0.30, max(c(all_preds$conf.high, brighton_data_plot$Algae_Cover_prop), na.rm = TRUE) * 1.10)

# Plot
macro_line_plot <- ggplot(all_preds, aes(x = Sampling_Date, y = predicted, group = Treatment)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = Treatment), alpha = 0.1, color = NA) +
  geom_line(
    aes(color = Treatment),
    linewidth = 0.5,
    alpha = 0.3,
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
    aes(x = Sampling_Date, y = Algae_Cover_prop, color = Treatment),
    shape = 20, width = 0.25, height = 0.01, alpha = 0.7, size = 1.5, inherit.aes = FALSE
  ) +
  labs(
    title = "Macroalgae cover",
    x = "Sampling date",
    y = "Macroalgae cover (%)"
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
  scale_y_continuous(labels = scales::percent, limits = c(0, y_upper_glm), oob = scales::squish) +
  scale_color_manual(values = TREATMENT_COLORS, breaks = LEGEND_ORDER) +
  scale_fill_manual(values = TREATMENT_COLORS,  breaks = LEGEND_ORDER)

print(macro_line_plot)

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
  if ("component" %in% names(sm)) sm <- dplyr::filter(sm, component == "mean")
  r <- dplyr::filter(sm, term == "Pad_ColourB:Pad_RegionP")
  if (nrow(r) > 0) r$p.value[1] else NA_real_
}

create_conditional_plot <- function(model, date_label, include_legend = FALSE) {
  
  # ggpredict gives predicted mean value + CI at Ambient/Warm for each region
  pr <- ggeffects::ggpredict(model, terms = c("Pad_Colour", "Pad_Region")) %>%
    mutate(
      Pad_Region = factor(group, levels = c("N", "P"), labels = REGION_LEVELS) |> droplevels(),
      x          = factor(x,     levels = c("W", "B"), labels = c("Ambient", "Warm"))
    ) %>%
    filter(!is.na(x), !is.na(Pad_Region))
  
  p_val <- get_interaction_pvalue(model)
  p_lab <- ifelse(is.na(p_val), "p = NA", paste0("p = ", signif(p_val, 4)))
  
  ggplot(pr, aes(x = .data$x, y = .data$predicted, group = .data$Pad_Region)) +
    geom_ribbon(aes(ymin = .data$conf.low, ymax = .data$conf.high, fill = .data$Pad_Region),
                alpha = 0.2, colour = NA) +
    geom_line(aes(color = .data$Pad_Region), linewidth = 0.8) +
    geom_point(aes(color = .data$Pad_Region), size = 2) +
    labs(
      title    = date_label,
      subtitle = p_lab,
      x = "Pad colour",
      y = "Macroalgae cover (%)"
    ) +
    theme_minimal() +
    theme(
      plot.title    = element_text(face = "bold", hjust = 0.5, size = 12, margin = margin(b = 1)),
      plot.subtitle = element_text(face = "italic", hjust = 0.5, size = 9, margin = margin(b = 4)),
      plot.margin   = margin(5, 5, 5, 5),
      axis.title    = element_text(size = 10),
      axis.text     = element_text(size = 7),
      legend.title  = element_text(face = "bold", size = 9),
      legend.text   = element_text(size = 8),
      legend.position = if (include_legend) "right" else "none",
      legend.box.background = element_blank()
    ) +
    scale_y_continuous(labels = scales::percent, oob = scales::squish) +
    scale_color_manual(name = "Pad region", breaks = REGION_LEVELS, values = REGION_VALUES) +
    scale_fill_manual( name = "Pad region", breaks = REGION_LEVELS, values = REGION_VALUES, guide = "none") +
    guides(color = guide_legend(override.aes = list(shape = 16, size = 4, fill = NA)))
}

plot_JUN   <- create_conditional_plot(glm_JUN_beta,   "June")
plot_JUL_1 <- create_conditional_plot(glm_JUL_1_beta, "July 1")
plot_JUL_2 <- create_conditional_plot(glm_JUL_2_beta, "July 2")
plot_AUG_1 <- create_conditional_plot(glm_AUG_1_beta, "August 1")
plot_AUG_2 <- create_conditional_plot(glm_AUG_2_beta, "August 2")
plot_SEP_1 <- create_conditional_plot(glm_SEP_1_beta, "September 1")
plot_SEP_2 <- create_conditional_plot(glm_SEP_2_beta, "September 2", include_legend = TRUE)

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

########## 6) Whole-season GAM (Beta–logit) ##########

d <- brighton_data %>%
  mutate(
    Days_Since_Start = as.numeric(Sampling_Date - min(Sampling_Date, na.rm = TRUE))
  ) %>%
  filter(!is.na(Sampling_Date))

# Whole-season GAM:
macroalgae_gam_final <- mgcv::gam(
  Algae_Cover_prop ~
    s(Days_Since_Start, by = interaction(Pad_Colour, Pad_Region), k = 7) +
    Pad_Colour * Pad_Region +
    s(Pad_ID, bs = "re"),
  data   = d,
  family = mgcv::betar(link = "logit"),
  method = "REML"
)

print(summary(macroalgae_gam_final))

########## 7) GAM interaction classification + summary table ##########

extract_macroalgae_gam_stats <- function(model, date_label = "Whole Season") {
  
  # Same equivalence-band logic as per-date GLMs
  # Parametric coefficient table from mgcv
  ptab <- as.data.frame(summary(model)$p.table) %>%
    tibble::rownames_to_column("term")
  
  names(ptab) <- tolower(gsub("[[:punct:] ]+", "_", names(ptab)))
  p_candidates <- grep("^pr", names(ptab), value = TRUE)
  if (!length(p_candidates)) p_candidates <- grep("p_value|p$", names(ptab), value = TRUE)
  pcol <- p_candidates[1]
  ptab <- dplyr::rename(ptab, p_value = all_of(pcol))
  
  # Term matching
  tlc <- tolower(ptab$term)
  i_int <- which(tlc %in% "(intercept)")
  i_w   <- which(tlc %in% "pad_colourb")
  i_p   <- which(tlc %in% "pad_regionp")
  i_ix  <- which(tlc %in% "pad_colourb:pad_regionp")
  
  getv <- function(col, i) if (length(i)) ptab[[col]][i[1]] else NA_real_
  
  α    <- getv("estimate",  i_int)
  βw   <- getv("estimate",  i_w)
  βp   <- getv("estimate",  i_p)
  βint <- getv("estimate",  i_ix)
  
  sew  <- getv("std_error", i_w)
  sep  <- getv("std_error", i_p)
  sei  <- getv("std_error", i_ix)
  
  pv_w <- getv("p_value",   i_w)
  pv_p <- getv("p_value",   i_p)
  pv_i <- getv("p_value",   i_ix)
  
  # Additive vs observed (logit)
  Δadd <- βw + βp
  Δobs <- Δadd + βint
  
  null_prop <- plogis(α + Δadd)
  obs_prop  <- plogis(α + Δobs)
  
  # Interaction OR + CI (Wald on logit)
  int_or    <- exp(βint)
  int_or_lo <- exp(βint - 1.96 * sei)
  int_or_hi <- exp(βint + 1.96 * sei)
  
  # Perspective metrics (% change)
  w_pct  <- (exp(βw)        - 1) * 100
  wp_pct <- (exp(βw + βint) - 1) * 100
  p_pct  <- (exp(βp)        - 1) * 100
  pw_pct <- (exp(βp + βint) - 1) * 100
  
  ci_overlaps_band <- !(int_or_hi < EQUIV_LO || int_or_lo > EQUIV_HI)
  ci_outside_band  <-  (int_or_hi < EQUIV_LO || int_or_lo > EQUIV_HI)
  
  dir_lab <- dplyr::case_when(
    (Δadd > 0 & Δobs < 0) ~ "Reversal (Pos to Neg)",
    (Δadd < 0 & Δobs > 0) ~ "Reversal (Neg to Pos)",
    (βw > 0 & βp > 0 & βint > 0) ~ "Synergism (>1)",
    (βw > 0 & βp > 0 & βint < 0) ~ "Antagonism (<1)",
    (βw < 0 & βp < 0 & βint < 0) ~ "Synergism (<1)",
    (βw < 0 & βp < 0 & βint > 0) ~ "Antagonism (>1)",
    sign(Δadd) == sign(βint) ~ "Synergism (same–sign)",
    TRUE ~ "Antagonism (opp–sign)"
  )
  
  int_class <- dplyr::case_when(
    is.na(pv_i) | pv_i >= P_ALPHA         ~ "Unclassified (p≥0.05)",
    ci_overlaps_band                      ~ "Additive (CI overlaps ±5%)",
    ci_outside_band & pv_i < P_ALPHA      ~ dir_lab,
    TRUE                                  ~ "Unclassified"
  )
  
  out <- tibble::tibble(
    Date                         = date_label,
    Term                         = c("Warming", "Pollution", "Warming × Pollution"),
    Estimate                     = c(βw, βp, βint),
    `Std. Error`                 = c(sew, sep, sei),
    `p-value`                    = c(pv_w, pv_p, pv_i),
    `ΔAdd (logit)`               = c(NA, NA, Δadd),
    `ΔObs (logit)`               = c(NA, NA, Δobs),
    `Null Proportion`            = c(NA, NA, null_prop),
    `Observed Proportion`        = c(NA, NA, obs_prop),
    `Odds Ratio`                 = c(NA, NA, int_or),
    `Interaction Ratio CI`       = c(NA, NA, sprintf("[%.3f, %.3f]", int_or_lo, int_or_hi)),
    `Interaction Class`          = c(NA, NA, int_class),
    `Warming -> under Pollution`  = c(NA, NA, sprintf("%+.1f%% to %+.1f%%", w_pct,  wp_pct)),
    `Pollution -> under Warming`  = c(NA, NA, sprintf("%+.1f%% to %+.1f%%", p_pct,  pw_pct))
  )
  # Formatting
  out$Estimate     <- round(out$Estimate, 4)
  out$`Std. Error` <- round(out$`Std. Error`, 4)
  out$`p-value`    <- ifelse(is.na(out$`p-value`), NA, format.pval(out$`p-value`, digits = 3, eps = 0.001))
  out$`ΔAdd (logit)`        <- ifelse(is.na(out$`ΔAdd (logit)`), NA, round(out$`ΔAdd (logit)`, 4))
  out$`ΔObs (logit)`        <- ifelse(is.na(out$`ΔObs (logit)`), NA, round(out$`ΔObs (logit)`, 4))
  out$`Null Proportion`     <- ifelse(is.na(out$`Null Proportion`), NA, round(out$`Null Proportion`, 4))
  out$`Observed Proportion` <- ifelse(is.na(out$`Observed Proportion`), NA, round(out$`Observed Proportion`, 4))
  
  out
}

macroalgae_gam_table <- extract_macroalgae_gam_stats(macroalgae_gam_final, "Whole Season")

DT::datatable(
  macroalgae_gam_table,
  caption = "Macroalgae cover GAM (beta–logit; 4 smooths by treatment; + pad RE): parametric fixed effects & interaction summary",
  options = list(pageLength = 10, scrollX = TRUE)
)

########## 8) Whole-season GAM visualisation ##########

# Map calendar dates to Days_Since_Start axis used by the GAM
start_date  <- min(d$Sampling_Date, na.rm = TRUE)
date_breaks <- as.numeric(as.Date(names(DATE_LOOKUP)) - start_date)

# Because the model contains s(Pad_ID), predict.gam requires a Pad_ID in newdata
pad_levels <- levels(d$Pad_ID)
pad_ref    <- if (length(pad_levels)) pad_levels[1] else NA_character_

# Prediction grid across the full season
pred_grid <- tidyr::expand_grid(
  Days_Since_Start = seq(min(d$Days_Since_Start, na.rm = TRUE),
                         max(d$Days_Since_Start, na.rm = TRUE),
                         length.out = 300),
  Pad_Colour = factor(c("W", "B"), levels = c("W", "B")),
  Pad_Region = factor(c("N", "P"), levels = c("N", "P"))
) %>%
  mutate(Pad_ID = factor(pad_ref, levels = pad_levels))

# Predict on the link scale (logit), then back-transform
pr_link <- predict(
  macroalgae_gam_final,
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
    predicted = plogis(fit),
    conf.low  = plogis(fit - 1.96 * se),
    conf.high = plogis(fit + 1.96 * se),
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
  mutate(Pad_ID = factor(pad_ref, levels = pad_levels))

dp_link <- predict(
  macroalgae_gam_final,
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
    predicted = plogis(fit),
    Pad_Region_lbl = factor(Pad_Region, levels = c("N", "P"), labels = c("Non-Polluted", "Polluted")),
    Pad_Colour_lbl = factor(Pad_Colour, levels = c("W", "B"), labels = c("Ambient", "Warm")),
    Treatment      = paste(Pad_Region_lbl, Pad_Colour_lbl, sep = " + ")
  )

brighton_data_GAM_plot <- d %>%
  mutate(
    Pad_Region_lbl = factor(Pad_Region, levels = c("N", "P"), labels = c("Non-Polluted", "Polluted")),
    Pad_Colour_lbl = factor(Pad_Colour, levels = c("W", "B"), labels = c("Ambient", "Warm")),
    Treatment      = paste(Pad_Region_lbl, Pad_Colour_lbl, sep = " + ")
  )

y_upper_gam <- min(0.30, max(c(gam_preds$conf.high, brighton_data_GAM_plot$Algae_Cover_prop), na.rm = TRUE) * 1.10)

# Whole-season GAM plot (logit scale)
macroalgae_gam_time_plot <- ggplot(gam_preds, aes(x = Days_Since_Start, y = predicted, group = Treatment)) +
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
    aes(x = Days_Since_Start, y = Algae_Cover_prop, color = Treatment),
    shape = 20, width = 0.8, height = 0.01, alpha = 0.7, size = 1.5, inherit.aes = FALSE
  ) +
  scale_x_continuous(breaks = date_breaks, labels = DATE_LABELS, expand = expansion(mult = c(0.02, 0.02))) +
  scale_y_continuous(labels = scales::percent, limits = c(0, y_upper_gam), oob = scales::squish) +
  scale_color_manual(values = TREATMENT_COLORS, breaks = LEGEND_ORDER) +
  scale_fill_manual(values = TREATMENT_COLORS,  breaks = LEGEND_ORDER) +
  labs(
    title = "Whole-season GAM: macroalgae cover",
    x = "Sampling date",
    y = "Macroalgae cover (%)"
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

print(macroalgae_gam_time_plot)