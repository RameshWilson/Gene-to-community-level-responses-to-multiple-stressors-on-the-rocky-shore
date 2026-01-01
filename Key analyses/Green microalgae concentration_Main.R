# Green microalgae concentration_Main.R
#
#   Brighton experiment (2×2): Warming (Pad_Colour) × Pollution (Pad_Region)
#   Response: green microalgae concentration
#
#   Outputs:
#     1) Per-date Gaussian GLMs(identity; June only)
#     2) GLM interaction classification table
#     3) Whole-summer GLM-style visualisation (ggeffects predictions + raw overlay)
#     4) Conditional per-date panels (pollution-stratified lines/ribbons across pad colour)
#
#   Evidence/diagnostics (residual checks, alternative structures) are kept in:
#     Green microalgae concentration_Additional.R

########## 0) Configuration ##########

library(here)

DATA_CSV <- here("Dataframes", "Full Brighton dataframe.csv")

# Sampling dates (fixed)
D_JUN  <- as.Date("2023-06-20")
D_JUL1 <- as.Date("2023-07-06")
D_JUL2 <- as.Date("2023-07-20")
D_AUG1 <- as.Date("2023-08-01")
D_AUG2 <- as.Date("2023-08-17")
D_SEP1 <- as.Date("2023-09-01")
D_SEP2 <- as.Date("2023-09-15")

DATE_LOOKUP <- tibble::tibble(
  Sampling_Date = c(D_JUN, D_JUL1, D_JUL2, D_AUG1, D_AUG2, D_SEP1, D_SEP2),
  DATE_LABEL    = c("June", "July 1", "July 2", "August 1", "August 2", "September 1", "September 2")
)

DATE_LABELS <- DATE_LOOKUP$DATE_LABEL

# Treatment colors
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

# Interaction-classification gate + equivalence band (MAIN convention)
P_ALPHA   <- 0.05
EQUIV_PCT <- 0.05

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(lubridate)
  library(ggeffects)
  library(colorspace)
  library(broom)
  library(DT)
  library(scales)
  library(tibble)
})

# Helper: ensure correct reference levels within any subset
parse_sampling_date <- function(x) {
  if (inherits(x, "Date")) return(x)
  tmp_dmy <- suppressWarnings(lubridate::dmy(x))
  frac_na <- mean(is.na(tmp_dmy))
  if (is.finite(frac_na) && frac_na < 0.5) {
    return(tmp_dmy)
  } else {
    return(suppressWarnings(lubridate::ymd(x)))
  }
}

# Prep a per-date subset with guaranteed reference levels (W and N).
prep_subset <- function(df, date_value) {
  out <- df %>%
    dplyr::filter(.data$Sampling_Date == date_value)
  
  out$Pad_Colour <- relevel(out$Pad_Colour, ref = "W")
  out$Pad_Region <- relevel(out$Pad_Region, ref = "N")
  
  out
}

########## 1) Load & prep data ##########

brighton_data <- readr::read_csv(DATA_CSV, show_col_types = FALSE)

# Parse Sampling_Date
brighton_data$Sampling_Date <- parse_sampling_date(brighton_data$Sampling_Date)

# Ensure factors are coded with controls first:
#   Pad_Colour: W=Ambient, B=Warm
#   Pad_Region: N=Non-Polluted, P=Polluted
brighton_data <- brighton_data %>%
  mutate(
    Pad_Colour = factor(Pad_Colour, levels = c("W", "B")),
    Pad_Region = factor(Pad_Region, levels = c("N", "P"))
  )

# Subsets per date (June is modelled; others used descriptively)
JUN   <- prep_subset(brighton_data, D_JUN)
JUL_1 <- brighton_data %>% dplyr::filter(Sampling_Date == D_JUL1)
JUL_2 <- brighton_data %>% dplyr::filter(Sampling_Date == D_JUL2)
AUG_1 <- brighton_data %>% dplyr::filter(Sampling_Date == D_AUG1)
AUG_2 <- brighton_data %>% dplyr::filter(Sampling_Date == D_AUG2)
SEP_1 <- brighton_data %>% dplyr::filter(Sampling_Date == D_SEP1)
SEP_2 <- brighton_data %>% dplyr::filter(Sampling_Date == D_SEP2)

########## 2) June GLM (Gamma, identity) ##########

glm_green_microalgae_JUN <- glm(
  Green_Algae ~ Pad_Colour * Pad_Region,
  data   = JUN,
  family = Gamma(link = "identity")
)

print(summary(glm_green_microalgae_JUN))

########## 3) June interaction classification + summary table ##########

extract_green_glm_stats <- function(model, date_label = "June", model_label = "Gamma(identity)") {
  
  fe <- broom::tidy(model, conf.int = TRUE, conf.level = 0.95)
  
  # If conf.int was not returned, construct Wald intervals
  if (!all(c("conf.low", "conf.high") %in% names(fe))) {
    fe <- fe %>%
      mutate(
        conf.low  = estimate - 1.96 * std.error,
        conf.high = estimate + 1.96 * std.error
      )
  }
  
  a_row   <- fe %>% dplyr::filter(term == "(Intercept)")
  w_row   <- fe %>% dplyr::filter(term == "Pad_ColourB")
  p_row   <- fe %>% dplyr::filter(term == "Pad_RegionP")
  int_row <- fe %>% dplyr::filter(term == "Pad_ColourB:Pad_RegionP")
  
  # Extract key terms (identity scale)
  α    <- a_row$estimate
  βw   <- w_row$estimate
  βp   <- p_row$estimate
  βint <- int_row$estimate
  pint <- int_row$p.value
  
  # Additivity vs observed (identity scale)
  Δadd      <- βw + βp
  Δobs      <- Δadd + βint
  null_pred <- α + Δadd
  obs_pred  <- α + Δobs
  
  # ±5% equivalence band around additivity on identity scale
  tol     <- EQUIV_PCT * abs(null_pred)
  band_lo <- -tol
  band_hi <-  tol
  
  # 95% CI for interaction term
  int_lo <- int_row$conf.low
  int_hi <- int_row$conf.high
  
  dir_label <- dplyr::case_when(
    (Δadd < 0 & Δobs > 0) ~ "Reversal (Neg to Pos)",
    (Δadd > 0 & Δobs < 0) ~ "Reversal (Pos to Neg)",
    (Δadd > 0 &  βint > 0) ~ "Synergism",
    (Δadd > 0 &  βint < 0) ~ "Antagonism",
    (Δadd < 0 &  βint < 0) ~ "Synergism",
    (Δadd < 0 &  βint > 0) ~ "Antagonism",
    sign(Δadd) == sign(βint) ~ "Synergism (same–sign)",
    TRUE                     ~ "Antagonism (opp–sign)"
  )
  
  ci_overlaps_band <- !(int_hi < band_lo || int_lo > band_hi)
  
  int_class <- dplyr::case_when(
    is.na(pint) | pint >= P_ALPHA ~ "Unclassified (p≥0.05)",
    ci_overlaps_band               ~ "Additive (CI overlaps ±5% band)",
    TRUE                           ~ dir_label
  )
  
  # Perspective contrasts (identity units)
  warming_under_pollution <- βw + βint
  pollution_under_warming <- βp + βint
  
  fe %>%
    dplyr::filter(term %in% c("Pad_ColourB", "Pad_RegionP", "Pad_ColourB:Pad_RegionP")) %>%
    dplyr::transmute(
      Date   = date_label,
      Model  = model_label,
      Term   = dplyr::case_when(
        term == "Pad_ColourB"             ~ "Warming",
        term == "Pad_RegionP"             ~ "Pollution",
        term == "Pad_ColourB:Pad_RegionP" ~ "Warming × Pollution"
      ),
      Estimate     = formatC(estimate,  format = "f", digits = 4),
      `Std. Error` = formatC(std.error, format = "f", digits = 4),
      `p-value`    = ifelse(p.value < 1e-4,
                            formatC(p.value, format = "e", digits = 2),
                            formatC(p.value, format = "f", digits = 3)),
      `ΔAdd (identity)` = if_else(Term == "Warming × Pollution", round(Δadd,      4), NA_real_),
      `ΔObs (identity)` = if_else(Term == "Warming × Pollution", round(Δobs,      4), NA_real_),
      `Null Prediction` = if_else(Term == "Warming × Pollution", round(null_pred, 4), NA_real_),
      `Obs Prediction`  = if_else(Term == "Warming × Pollution", round(obs_pred,  4), NA_real_),
      `Interaction 95% CI (identity)` =
        if_else(Term == "Warming × Pollution",
                sprintf("[%.4f, %.4f]", int_lo, int_hi), NA_character_),
      `±5% Band (identity)` =
        if_else(Term == "Warming × Pollution",
                sprintf("[%.4f, %.4f]", band_lo, band_hi), NA_character_),
      `Interaction Class` =
        if_else(Term == "Warming × Pollution", int_class, NA_character_),
      `Warming -> under Pollution (identity)` =
        if_else(Term == "Warming × Pollution",
                sprintf("%+.4f to %+.4f", βw, warming_under_pollution), NA_character_),
      `Pollution -> under Warming (identity)` =
        if_else(Term == "Warming × Pollution",
                sprintf("%+.4f to %+.4f", βp, pollution_under_warming), NA_character_)
    )
}

green_microalgae_results <- extract_green_glm_stats(glm_green_microalgae_JUN, "June", "Gamma(identity)")

DT::datatable(
  green_microalgae_results,
  caption = "Green microalgae (June GLM; Gamma identity): effects, predictions & interaction",
  options = list(pageLength = 10, scrollX = TRUE)
)

########## 4) Whole summer visualisation (June = GLM; others = descriptive) ##########

# June: model-based predictions (Gamma identity)
pred_June <- ggeffects::ggpredict(glm_green_microalgae_JUN, terms = c("Pad_Colour", "Pad_Region")) %>%
  tibble::as_tibble() %>%
  mutate(
    Sampling_Date = "June",
    Pad_Region = factor(group, levels = c("N", "P"), labels = c("Non-Polluted", "Polluted")),
    Pad_Colour = factor(x,     levels = c("W", "B"), labels = c("Ambient", "Warm")),
    Treatment  = paste(Pad_Region, Pad_Colour, sep = " + ")
  ) %>%
  dplyr::select(Sampling_Date, Treatment, predicted, conf.low, conf.high)

# Other dates: descriptive mean ± 95% CI (t-based) from raw data
other_dates <- bind_rows(JUL_1, JUL_2, AUG_1, AUG_2, SEP_1, SEP_2) %>%
  dplyr::left_join(DATE_LOOKUP, by = "Sampling_Date") %>%
  mutate(
    Sampling_Date = factor(DATE_LABEL, levels = DATE_LABELS),
    Pad_Region = factor(Pad_Region, levels = c("N", "P"), labels = c("Non-Polluted", "Polluted")),
    Pad_Colour = factor(Pad_Colour, levels = c("W", "B"), labels = c("Ambient", "Warm")),
    Treatment  = paste(Pad_Region, Pad_Colour, sep = " + ")
  ) %>%
  dplyr::filter(Sampling_Date != "June")

green_summary <- other_dates %>%
  group_by(Sampling_Date, Treatment) %>%
  summarise(
    predicted = mean(Green_Algae, na.rm = TRUE),
    n  = dplyr::n(),
    sd = sd(Green_Algae, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    se       = sd / sqrt(pmax(n, 1)),
    tcrit    = qt(0.975, df = pmax(n - 1, 1)),
    conf.low = predicted - tcrit * se,
    conf.high= predicted + tcrit * se
  ) %>%
  dplyr::mutate(Sampling_Date = factor(Sampling_Date, levels = DATE_LABELS))

# Combine June model + later descriptive summaries
all_preds_green <- bind_rows(
  pred_June %>% mutate(Sampling_Date = factor(Sampling_Date, levels = DATE_LABELS)),
  green_summary
) %>%
  dplyr::filter(!is.na(Treatment))

# Raw points for overlay (all dates)
brighton_data_plot <- brighton_data %>%
  dplyr::left_join(DATE_LOOKUP, by = "Sampling_Date") %>%
  mutate(
    Sampling_Date = factor(DATE_LABEL, levels = DATE_LABELS),
    Pad_Region = factor(Pad_Region, levels = c("N", "P"), labels = c("Non-Polluted", "Polluted")),
    Pad_Colour = factor(Pad_Colour, levels = c("W", "B"), labels = c("Ambient", "Warm")),
    Treatment  = paste(Pad_Region, Pad_Colour, sep = " + ")
  ) %>%
  dplyr::filter(!is.na(Treatment), !is.na(Sampling_Date))

# Plot
line_plot_green <- ggplot(all_preds_green, aes(x = Sampling_Date, y = predicted, group = Treatment)) +
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
    data  = brighton_data_plot,
    aes(x = Sampling_Date, y = Green_Algae, color = Treatment),
    shape = 20, width = 0.25, height = 0, alpha = 0.7, size = 1.5, inherit.aes = FALSE
  ) +
  labs(
    title = "Green microalgae concentration",
    x = "Sampling date",
    y = expression("Green microalgae concentration (" ~ plain("mg") * "/" * plain("cm")^2 * ")")
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
  scale_fill_manual(values  = TREATMENT_COLORS, breaks = LEGEND_ORDER)

print(line_plot_green)

########## 5) Conditional plot (June only) ##########

REGION_LEVELS <- c("Non-Polluted", "Polluted")
REGION_VALUES <- unname(c(
  TREATMENT_COLORS["Non-Polluted + Ambient"],
  TREATMENT_COLORS["Polluted + Ambient"]
))

# Interaction p-value (June model)
sm <- broom::tidy(glm_green_microalgae_JUN)
r  <- dplyr::filter(sm, term == "Pad_ColourB:Pad_RegionP")
p_val <- if (nrow(r) > 0) r$p.value[1] else NA_real_

p_lab <- ifelse(is.na(p_val), "p = NA", paste0("p = ", signif(p_val, 4)))
face  <- if (!is.na(p_val) && p_val < P_ALPHA) "plain" else "italic"

pr_JUN <- ggeffects::ggpredict(glm_green_microalgae_JUN, terms = c("Pad_Colour", "Pad_Region")) %>%
  tibble::as_tibble() %>%
  mutate(
    Pad_Region = factor(group, levels = c("N", "P"), labels = REGION_LEVELS) |> droplevels(),
    x          = factor(x,     levels = c("W", "B"), labels = c("Ambient", "Warm"))
  ) %>%
  dplyr::filter(!is.na(x), !is.na(Pad_Region))

plot_JUN_green <- ggplot(pr_JUN, aes(x = x, y = predicted, group = Pad_Region)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = Pad_Region), alpha = 0.2, colour = NA) +
  geom_line(aes(color = Pad_Region), linewidth = 0.8) +
  geom_point(aes(color = Pad_Region), size = 2) +
  labs(
    title    = "June",
    subtitle = p_lab,
    x = "Pad colour",
    y = expression("Green microalgae concentration (" ~ plain("mg") * "/" * plain("cm")^2 * ")")
  ) +
  theme_minimal() +
  theme(
    plot.title    = element_text(face = "bold", hjust = 0.5, size = 12, margin = margin(b = 1)),
    plot.subtitle = element_text(face = face,   hjust = 0.5, size = 9,  margin = margin(b = 4)),
    axis.title    = element_text(size = 10),
    axis.text     = element_text(size = 7),
    legend.title  = element_text(face = "bold", size = 9),
    legend.text   = element_text(size = 8),
    legend.box.background = element_blank()
  ) +
  scale_color_manual(name = "Pad region", breaks = REGION_LEVELS, values = REGION_VALUES) +
  scale_fill_manual(name = "Pad region", breaks = REGION_LEVELS, values = REGION_VALUES, guide = "none") +
  guides(color = guide_legend(override.aes = list(shape = 16, size = 4, fill = NA)))

print(plot_JUN_green)