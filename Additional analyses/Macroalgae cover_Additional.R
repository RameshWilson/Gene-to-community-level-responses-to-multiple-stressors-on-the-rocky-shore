# Macroalgae cover_Additional
#
#   Evidence/diagnostics for Brighton macroalgae cover.

########## 0) Configuration ##########

library(here)

DATA_CSV <- here("Dataframes", "Full Brighton dataframe.csv")

# Sampling dates
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
  library(betareg)
  library(mgcv)
})

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
  dplyr::mutate(
    Pad_Colour = factor(Pad_Colour, levels = c("W","B")),  # W=Ambient, B=Warm
    Pad_Region = factor(Pad_Region, levels = c("N","P")),  # N=Non-Polluted, P=Polluted
    Pad_ID     = factor(Pad_ID),
    Algae_Cover_prop = pmax(1e-4, pmin(0.9999, Algae_Cover / 100)) # Bind 0-1
  ) %>%
  dplyr::filter(!is.na(Sampling_Date))

########## 2) Per-date Beta GLMs (logit) + diagnostics ##########

# Diagnostics helper:
# Residual diagnostics (Quantile residual histogram + QQ plot)
diag_plots_beta <- function(m, main_prefix = "") {
  op <- par(mfrow = c(1,2)); on.exit(par(op), add = TRUE)
  rq <- residuals(m, type = "quantile")
  hist(rq, main = paste0(main_prefix, "Quantile residuals (hist)"), xlab = "Quantile residuals")
  qqnorm(rq, main = paste0(main_prefix, "Quantile residuals (QQ)")); qqline(rq)
  invisible(NULL)
}

# June
JUN <- dplyr::filter(brighton_data, Sampling_Date == D_JUN)
JUN$Pad_Colour <- stats::relevel(JUN$Pad_Colour, ref = "W")
JUN$Pad_Region <- stats::relevel(JUN$Pad_Region, ref = "N")
print(JUN %>% dplyr::group_by(Pad_Colour, Pad_Region) %>% dplyr::summarise(n=dplyr::n(), mean_prop=mean(Algae_Cover_prop), .groups="drop"))

glm_JUN_beta <- betareg::betareg(Algae_Cover_prop ~ Pad_Colour * Pad_Region, data = JUN, link = "logit")
print(summary(glm_JUN_beta))
diag_plots_beta(glm_JUN_beta, "June  ")

# July 1
JUL_1 <- dplyr::filter(brighton_data, Sampling_Date == D_JUL1)
JUL_1$Pad_Colour <- stats::relevel(JUL_1$Pad_Colour, ref = "W")
JUL_1$Pad_Region <- stats::relevel(JUL_1$Pad_Region, ref = "N")
print(JUL_1 %>% dplyr::group_by(Pad_Colour, Pad_Region) %>% dplyr::summarise(n=dplyr::n(), mean_prop=mean(Algae_Cover_prop), .groups="drop"))

glm_JUL_1_beta <- betareg::betareg(Algae_Cover_prop ~ Pad_Colour * Pad_Region, data = JUL_1, link = "logit")
print(summary(glm_JUL_1_beta))
diag_plots_beta(glm_JUL_1_beta, "July 1")

# July 2
JUL_2 <- dplyr::filter(brighton_data, Sampling_Date == D_JUL2)
JUL_2$Pad_Colour <- stats::relevel(JUL_2$Pad_Colour, ref = "W")
JUL_2$Pad_Region <- stats::relevel(JUL_2$Pad_Region, ref = "N")
print(JUL_2 %>% dplyr::group_by(Pad_Colour, Pad_Region) %>% dplyr::summarise(n=dplyr::n(), mean_prop=mean(Algae_Cover_prop), .groups="drop"))

glm_JUL_2_beta <- betareg::betareg(Algae_Cover_prop ~ Pad_Colour * Pad_Region, data = JUL_2, link = "logit")
print(summary(glm_JUL_2_beta))
diag_plots_beta(glm_JUL_2_beta, "July 2")

# August 1
AUG_1 <- dplyr::filter(brighton_data, Sampling_Date == D_AUG1)
AUG_1$Pad_Colour <- stats::relevel(AUG_1$Pad_Colour, ref = "W")
AUG_1$Pad_Region <- stats::relevel(AUG_1$Pad_Region, ref = "N")
print(AUG_1 %>% dplyr::group_by(Pad_Colour, Pad_Region) %>% dplyr::summarise(n=dplyr::n(), mean_prop=mean(Algae_Cover_prop), .groups="drop"))

glm_AUG_1_beta <- betareg::betareg(Algae_Cover_prop ~ Pad_Colour * Pad_Region, data = AUG_1, link = "logit")
print(summary(glm_AUG_1_beta))
diag_plots_beta(glm_AUG_1_beta, "Aug 1 ")

# August 2
AUG_2 <- dplyr::filter(brighton_data, Sampling_Date == D_AUG2)
AUG_2$Pad_Colour <- stats::relevel(AUG_2$Pad_Colour, ref = "W")
AUG_2$Pad_Region <- stats::relevel(AUG_2$Pad_Region, ref = "N")
print(AUG_2 %>% dplyr::group_by(Pad_Colour, Pad_Region) %>% dplyr::summarise(n=dplyr::n(), mean_prop=mean(Algae_Cover_prop), .groups="drop"))

glm_AUG_2_beta <- betareg::betareg(Algae_Cover_prop ~ Pad_Colour * Pad_Region, data = AUG_2, link = "logit")
print(summary(glm_AUG_2_beta))
diag_plots_beta(glm_AUG_2_beta, "Aug 2 ")

# September 1
SEP_1 <- dplyr::filter(brighton_data, Sampling_Date == D_SEP1)
SEP_1$Pad_Colour <- stats::relevel(SEP_1$Pad_Colour, ref = "W")
SEP_1$Pad_Region <- stats::relevel(SEP_1$Pad_Region, ref = "N")
print(SEP_1 %>% dplyr::group_by(Pad_Colour, Pad_Region) %>% dplyr::summarise(n=dplyr::n(), mean_prop=mean(Algae_Cover_prop), .groups="drop"))

glm_SEP_1_beta <- betareg::betareg(Algae_Cover_prop ~ Pad_Colour * Pad_Region, data = SEP_1, link = "logit")
print(summary(glm_SEP_1_beta))
diag_plots_beta(glm_SEP_1_beta, "Sep 1 ")

# September 2
SEP_2 <- dplyr::filter(brighton_data, Sampling_Date == D_SEP2)
SEP_2$Pad_Colour <- stats::relevel(SEP_2$Pad_Colour, ref = "W")
SEP_2$Pad_Region <- stats::relevel(SEP_2$Pad_Region, ref = "N")
print(SEP_2 %>% dplyr::group_by(Pad_Colour, Pad_Region) %>% dplyr::summarise(n=dplyr::n(), mean_prop=mean(Algae_Cover_prop), .groups="drop"))

glm_SEP_2_beta <- betareg::betareg(Algae_Cover_prop ~ Pad_Colour * Pad_Region, data = SEP_2, link = "logit")
print(summary(glm_SEP_2_beta))
diag_plots_beta(glm_SEP_2_beta, "Sep 2 ")

########## 3) Whole-season GAM diagnostics  ##########

brighton_data_GAM <- brighton_data %>%
  dplyr::mutate(
    Days_Since_Start = as.numeric(Sampling_Date - min(Sampling_Date, na.rm = TRUE))
  )

# Whole-season GAM:
macroalgae_gam_final <- mgcv::gam(
  Algae_Cover_prop ~
    s(Days_Since_Start, by = interaction(Pad_Colour, Pad_Region), k = 7) +
    Pad_Colour * Pad_Region +
    s(Pad_ID, bs = "re"),
  data   = brighton_data_GAM,
  family = mgcv::betar(link = "logit"),
  method = "REML"
)

print(summary(macroalgae_gam_final))
mgcv::gam.check(macroalgae_gam_final)