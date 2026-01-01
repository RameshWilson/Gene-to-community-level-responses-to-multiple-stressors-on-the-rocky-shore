# Cyanobacteria concentration_Additional
#
#   Evidence/diagnostics for Brighton cyanobacteria concentration.

########## 0) Configuration ##########

library(here)

DATA_CSV <- here("Dataframes", "Full Brighton dataframe.csv")

D_JUN  <- as.Date("2023-06-20")
D_JUL1 <- as.Date("2023-07-06")
D_JUL2 <- as.Date("2023-07-20")
D_AUG1 <- as.Date("2023-08-01")
D_AUG2 <- as.Date("2023-08-17")
D_SEP1 <- as.Date("2023-09-01")
D_SEP2 <- as.Date("2023-09-15")

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(lubridate)
  library(broom)
  library(mgcv)
  library(statmod)
})

# Per-date GLM family:
#   statmod::tweedie(var.power=1.5, link.power=1) gives Tweedie with identity link.
#   var.power≈1.5 is a mid-range default for "close to zero + continuous" ecological concentration data.
#   For GLMs, setting consistently to 1.5 ensures comparability across season.
TW_FAM <- statmod::tweedie(var.power = 1.5, link.power = 1)

########## 1) Load & prep data ##########

brighton_data <- readr::read_csv(DATA_CSV, show_col_types = FALSE)

# Parse dates: try dmy first; fallback to ymd
if (!inherits(brighton_data$Sampling_Date, "Date")) {
  tmp_dmy <- suppressWarnings(lubridate::dmy(brighton_data$Sampling_Date))
  frac_na <- mean(is.na(tmp_dmy))
  brighton_data$Sampling_Date <- if (frac_na < 0.5) tmp_dmy else
    suppressWarnings(lubridate::ymd(brighton_data$Sampling_Date))
}

# Recode factors with controls first
brighton_data <- brighton_data %>%
  mutate(
    Pad_Colour = factor(Pad_Colour, levels = c("W","B")), # W = Ambient, B = Warm
    Pad_Region = factor(Pad_Region, levels = c("N","P"))  # N = Non-Polluted, P = Polluted
  )

# Split by sampling date (per-date GLMs).
JUN   <- filter(brighton_data, Sampling_Date == D_JUN)
JUL_1 <- filter(brighton_data, Sampling_Date == D_JUL1)
JUL_2 <- filter(brighton_data, Sampling_Date == D_JUL2)
AUG_1 <- filter(brighton_data, Sampling_Date == D_AUG1)
AUG_2 <- filter(brighton_data, Sampling_Date == D_AUG2)
SEP_1 <- filter(brighton_data, Sampling_Date == D_SEP1)
SEP_2 <- filter(brighton_data, Sampling_Date == D_SEP2)

# Relevel within each subset (belt-and-braces so the reference is always W and N).
JUN$Pad_Colour   <- relevel(JUN$Pad_Colour,   ref = "W"); JUN$Pad_Region   <- relevel(JUN$Pad_Region,   ref = "N")
JUL_1$Pad_Colour <- relevel(JUL_1$Pad_Colour, ref = "W"); JUL_1$Pad_Region <- relevel(JUL_1$Pad_Region, ref = "N")
JUL_2$Pad_Colour <- relevel(JUL_2$Pad_Colour, ref = "W"); JUL_2$Pad_Region <- relevel(JUL_2$Pad_Region, ref = "N")
AUG_1$Pad_Colour <- relevel(AUG_1$Pad_Colour, ref = "W"); AUG_1$Pad_Region <- relevel(AUG_1$Pad_Region, ref = "N")
AUG_2$Pad_Colour <- relevel(AUG_2$Pad_Colour, ref = "W"); AUG_2$Pad_Region <- relevel(AUG_2$Pad_Region, ref = "N")
SEP_1$Pad_Colour <- relevel(SEP_1$Pad_Colour, ref = "W"); SEP_1$Pad_Region <- relevel(SEP_1$Pad_Region, ref = "N")
SEP_2$Pad_Colour <- relevel(SEP_2$Pad_Colour, ref = "W"); SEP_2$Pad_Region <- relevel(SEP_2$Pad_Region, ref = "N")

# Diagnostics helper:
# Residual diagnostics (Pearson residual histogram + QQ plot)
diag_plots <- function(m, main_prefix = "") {
  op <- par(mfrow = c(1,2)); on.exit(par(op), add = TRUE)
  rp <- residuals(m, type = "pearson")
  hist(rp, main = paste0(main_prefix, "Pearson residuals (hist)"), xlab = "Pearson residuals")
  qqnorm(rp, main = paste0(main_prefix, "Pearson residuals (QQ)")); qqline(rp)
  invisible(NULL)
}

########## 2) Per-date GLMs: Tweedie(identity) vs Gaussian(identity) + diagnostics ##########
# AIC is not used as Tweedie GLM AIC is NA here (statmod Tweedie family).
# We rely on residual diagnostics (hist + QQ) across dates

glm_JUN_tw <- glm(Cyanobacteria ~ Pad_Colour * Pad_Region, data = JUN, family = TW_FAM)
glm_JUN_ga <- glm(Cyanobacteria ~ Pad_Colour * Pad_Region, data = JUN, family = gaussian())
print(summary(glm_JUN_tw)); print(summary(glm_JUN_ga))
diag_plots(glm_JUN_tw, "June Tweedie  "); diag_plots(glm_JUN_ga, "June Gaussian ")

glm_JUL1_tw <- glm(Cyanobacteria ~ Pad_Colour * Pad_Region, data = JUL_1, family = TW_FAM)
glm_JUL1_ga <- glm(Cyanobacteria ~ Pad_Colour * Pad_Region, data = JUL_1, family = gaussian())
print(summary(glm_JUL1_tw)); print(summary(glm_JUL1_ga))
diag_plots(glm_JUL1_tw, "July1 Tweedie  "); diag_plots(glm_JUL1_ga, "July1 Gaussian ")

glm_JUL2_tw <- glm(Cyanobacteria ~ Pad_Colour * Pad_Region, data = JUL_2, family = TW_FAM)
glm_JUL2_ga <- glm(Cyanobacteria ~ Pad_Colour * Pad_Region, data = JUL_2, family = gaussian())
print(summary(glm_JUL2_tw)); print(summary(glm_JUL2_ga))
diag_plots(glm_JUL2_tw, "July2 Tweedie  "); diag_plots(glm_JUL2_ga, "July2 Gaussian ")

glm_AUG1_tw <- glm(Cyanobacteria ~ Pad_Colour * Pad_Region, data = AUG_1, family = TW_FAM)
glm_AUG1_ga <- glm(Cyanobacteria ~ Pad_Colour * Pad_Region, data = AUG_1, family = gaussian())
print(summary(glm_AUG1_tw)); print(summary(glm_AUG1_ga))
diag_plots(glm_AUG1_tw, "Aug1 Tweedie  "); diag_plots(glm_AUG1_ga, "Aug1 Gaussian ")

glm_AUG2_tw <- glm(Cyanobacteria ~ Pad_Colour * Pad_Region, data = AUG_2, family = TW_FAM)
glm_AUG2_ga <- glm(Cyanobacteria ~ Pad_Colour * Pad_Region, data = AUG_2, family = gaussian())
print(summary(glm_AUG2_tw)); print(summary(glm_AUG2_ga))
diag_plots(glm_AUG2_tw, "Aug2 Tweedie  "); diag_plots(glm_AUG2_ga, "Aug2 Gaussian ")

glm_SEP1_tw <- glm(Cyanobacteria ~ Pad_Colour * Pad_Region, data = SEP_1, family = TW_FAM)
glm_SEP1_ga <- glm(Cyanobacteria ~ Pad_Colour * Pad_Region, data = SEP_1, family = gaussian())
print(summary(glm_SEP1_tw)); print(summary(glm_SEP1_ga))
diag_plots(glm_SEP1_tw, "Sep1 Tweedie  "); diag_plots(glm_SEP1_ga, "Sep1 Gaussian ")

glm_SEP2_tw <- glm(Cyanobacteria ~ Pad_Colour * Pad_Region, data = SEP_2, family = TW_FAM)
glm_SEP2_ga <- glm(Cyanobacteria ~ Pad_Colour * Pad_Region, data = SEP_2, family = gaussian())
print(summary(glm_SEP2_tw)); print(summary(glm_SEP2_ga))
diag_plots(glm_SEP2_tw, "Sep2 Tweedie  "); diag_plots(glm_SEP2_ga, "Sep2 Gaussian ")

# Across dates, Tweedie(identity) residual histograms/QQs are comparatively more consistent and clean
# than Gaussian(identity), which tends to show heavier tails / clear deviation from normality for this response.
# This supports using Tweedie for the main per-date inference; used across all dates for consistency and comparability

########## 3) Whole-season GAM data  ##########

brighton_data_GAM <- readr::read_csv(DATA_CSV, show_col_types = FALSE)

if (!inherits(brighton_data_GAM$Sampling_Date, "Date")) {
  tmp_dmy <- suppressWarnings(lubridate::dmy(brighton_data_GAM$Sampling_Date))
  frac_na <- mean(is.na(tmp_dmy))
  brighton_data_GAM$Sampling_Date <- if (frac_na < 0.5) tmp_dmy else
    suppressWarnings(lubridate::ymd(brighton_data_GAM$Sampling_Date))
}

brighton_data_GAM <- brighton_data_GAM %>%
  mutate(
    Pad_Colour = factor(Pad_Colour, levels = c("W","B")),
    Pad_Region = factor(Pad_Region, levels = c("N","P")),
    Days_Since_Start = as.numeric(Sampling_Date - min(Sampling_Date, na.rm = TRUE))
  )

########## 4) Whole-season GAM + diagnostics ##########

cyano_gam_model <- mgcv::gam(
  Cyanobacteria ~
    s(Days_Since_Start, by = interaction(Pad_Colour, Pad_Region), k = 7) +
    Pad_Colour * Pad_Region +
    s(Pad_ID, bs = "re"),
  data   = brighton_data_GAM,
  family = mgcv::tw(link = "identity"),
  method = "REML"
)

# About the warning:
#   - "NA/Inf replaced by maximum positive value" originates inside mgcv's Tweedie p-search: tries candidate p values and sometimes the likelihood evaluation briefly will blow up.
#   - It is a warning about numerical evaluation during the internal search for p, not that the fitted model is invalid.
#   - The fitted model converged cleanly, Hessian is positive definite, fitted values are non-negative, and residuals are acceptable. 
#   - Given diagnostics, we treat the notification as an optimisation warning, but retain the model. 

print(summary(cyano_gam_model))
gam.check(cyano_gam_model)

# Try by fixing theta explicitly at mgcv estimated output (1.45)
#   - The same model as the original, just skipping the p-search decision.
#   - Reproduces the same parametric/smooth inference and fit metrics.
#   - Still shows the same warning, which reinforces that the warning is tied to evaluation/search internals rather than a failure to converge to a stable solution.
cyano_gam_model_1.45 <- mgcv::gam(
  Cyanobacteria ~
    s(Days_Since_Start, by = interaction(Pad_Colour, Pad_Region), k = 7) +
    Pad_Colour * Pad_Region +
    s(Pad_ID, bs = "re"),
  data   = brighton_data_GAM,
  family = mgcv::tw(theta = 1.45, link = "identity"),
  method = "REML"
)

print(summary(cyano_gam_model_1.45))
gam.check(cyano_gam_model_1.45)