# Grazer abundance_Additional
#
#   Evidence/diagnostics for Brighton grazer abundance.

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

DATE_LABELS <- c("June", "July 1", "July 2", "August 1", "August 2", "September 1", "September 2")

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(lubridate)
  library(brglm2)
  library(mgcv)
})

########## 1) Load & prep data ##########

brighton_data <- readr::read_csv(DATA_CSV, show_col_types = FALSE)

# Parse dates: try dmy first; fallback to ymd if needed
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
    Pad_Region = factor(Pad_Region, levels = c("N","P"))   # N=Non-Polluted, P=Polluted
  )

# Aggregate grazer abundance across pad and convert to binary presence/absence
brighton_data <- brighton_data %>%
  dplyr::mutate(
    Gastropod_Count_TOTAL = Gastropod_Count_ON + Gastropod_Count_SIDE,
    Gastropod_Binary      = as.integer(Gastropod_Count_TOTAL > 0)
  )

########## 2) Per-date GLMs (BRGLM) + diagnostics ##########

# Residual diagnostics (deviance residual histogram + QQ plot)
diag_plots_binom <- function(m, main_prefix = "") {
  op <- par(mfrow = c(1, 2)); on.exit(par(op), add = TRUE)
  r <- residuals(m, type = "deviance")
  hist(r, main = paste0(main_prefix, "Residuals (deviance)"), xlab = "Deviance residuals")
  qqnorm(r, main = paste0(main_prefix, "QQ plot")); qqline(r)
  invisible(NULL)
}

# Relevel factors within each subset (controls as reference)
relevel_ctrl <- function(df) {
  df$Pad_Colour <- stats::relevel(as.factor(df$Pad_Colour), ref = "W")
  df$Pad_Region <- stats::relevel(as.factor(df$Pad_Region), ref = "N")
  df
}

# Subsets per date (explicit)
JUN   <- relevel_ctrl(dplyr::filter(brighton_data, Sampling_Date == D_JUN))
JUL_1 <- relevel_ctrl(dplyr::filter(brighton_data, Sampling_Date == D_JUL1))
JUL_2 <- relevel_ctrl(dplyr::filter(brighton_data, Sampling_Date == D_JUL2))
AUG_1 <- relevel_ctrl(dplyr::filter(brighton_data, Sampling_Date == D_AUG1))
AUG_2 <- relevel_ctrl(dplyr::filter(brighton_data, Sampling_Date == D_AUG2))
SEP_1 <- relevel_ctrl(dplyr::filter(brighton_data, Sampling_Date == D_SEP1))
SEP_2 <- relevel_ctrl(dplyr::filter(brighton_data, Sampling_Date == D_SEP2))

# Each time-point GLM  shows:
#   Cell counts by treatment (0/1): 0 or 1 indicates complete separation risk (exploratory)
#   Summary and diagnostic plots for GLM
#   Comparative non-biased-reduced GLM and diagnostics for comparison

# June
print(
  JUN %>%
    dplyr::group_by(Pad_Colour, Pad_Region) %>%
    dplyr::summarise(n = dplyr::n(),
                     prop_present = mean(Gastropod_Binary),
                     .groups = "drop")
)
glm_JUN_bi <- glm(Gastropod_Binary ~ Pad_Colour * Pad_Region,
                  data = JUN, family = binomial("logit"), method = "brglmFit")
print(summary(glm_JUN_bi))
diag_plots_binom(glm_JUN_bi, "June (BRGLM)  ")

glm_JUN_mle <- glm(Gastropod_Binary ~ Pad_Colour * Pad_Region,
                   data = JUN, family = binomial("logit"))
print(summary(glm_JUN_mle))
diag_plots_binom(glm_JUN_mle, "June (MLE)    ")

# July 1
print(with(JUL_1, table(Pad_Colour, Pad_Region, Gastropod_Binary)))
print(
  JUL_1 %>%
    dplyr::group_by(Pad_Colour, Pad_Region) %>%
    dplyr::summarise(n = dplyr::n(),
                     prop_present = mean(Gastropod_Binary),
                     .groups = "drop")
)
glm_JUL_1_bi <- glm(Gastropod_Binary ~ Pad_Colour * Pad_Region,
                    data = JUL_1, family = binomial("logit"), method = "brglmFit")
print(summary(glm_JUL_1_bi))
diag_plots_binom(glm_JUL_1_bi, "July 1 (BRGLM)")

glm_JUL_1_mle <- glm(Gastropod_Binary ~ Pad_Colour * Pad_Region,
                     data = JUL_1, family = binomial("logit"))
print(summary(glm_JUL_1_mle))
diag_plots_binom(glm_JUL_1_mle, "July 1 (MLE)  ")

# July 2
print(with(JUL_2, table(Pad_Colour, Pad_Region, Gastropod_Binary)))
print(
  JUL_2 %>%
    dplyr::group_by(Pad_Colour, Pad_Region) %>%
    dplyr::summarise(n = dplyr::n(),
                     prop_present = mean(Gastropod_Binary),
                     .groups = "drop")
)
glm_JUL_2_bi <- glm(Gastropod_Binary ~ Pad_Colour * Pad_Region,
                    data = JUL_2, family = binomial("logit"), method = "brglmFit")
print(summary(glm_JUL_2_bi))
diag_plots_binom(glm_JUL_2_bi, "July 2 (BRGLM)")
glm_JUL_2_mle <- glm(Gastropod_Binary ~ Pad_Colour * Pad_Region,
                     data = JUL_2, family = binomial("logit"))
print(summary(glm_JUL_2_mle))
diag_plots_binom(glm_JUL_2_mle, "July 2 (MLE)  ")

# August 1
print(with(AUG_1, table(Pad_Colour, Pad_Region, Gastropod_Binary)))
print(
  AUG_1 %>%
    dplyr::group_by(Pad_Colour, Pad_Region) %>%
    dplyr::summarise(n = dplyr::n(),
                     prop_present = mean(Gastropod_Binary),
                     .groups = "drop")
)
glm_AUG_1_bi <- glm(Gastropod_Binary ~ Pad_Colour * Pad_Region,
                    data = AUG_1, family = binomial("logit"), method = "brglmFit")
print(summary(glm_AUG_1_bi))
diag_plots_binom(glm_AUG_1_bi, "Aug 1 (BRGLM) ")
glm_AUG_1_mle <- glm(Gastropod_Binary ~ Pad_Colour * Pad_Region,
                     data = AUG_1, family = binomial("logit"))
print(summary(glm_AUG_1_mle))
diag_plots_binom(glm_AUG_1_mle, "Aug 1 (MLE)   ")

# August 2
print(with(AUG_2, table(Pad_Colour, Pad_Region, Gastropod_Binary)))
print(
  AUG_2 %>%
    dplyr::group_by(Pad_Colour, Pad_Region) %>%
    dplyr::summarise(n = dplyr::n(),
                     prop_present = mean(Gastropod_Binary),
                     .groups = "drop")
)
glm_AUG_2_bi <- glm(Gastropod_Binary ~ Pad_Colour * Pad_Region,
                    data = AUG_2, family = binomial("logit"), method = "brglmFit")
print(summary(glm_AUG_2_bi))
diag_plots_binom(glm_AUG_2_bi, "Aug 2 (BRGLM) ")
glm_AUG_2_mle <- glm(Gastropod_Binary ~ Pad_Colour * Pad_Region,
                     data = AUG_2, family = binomial("logit"))
print(summary(glm_AUG_2_mle))
diag_plots_binom(glm_AUG_2_mle, "Aug 2 (MLE)   ")

# September 1
print(with(SEP_1, table(Pad_Colour, Pad_Region, Gastropod_Binary)))
print(
  SEP_1 %>%
    dplyr::group_by(Pad_Colour, Pad_Region) %>%
    dplyr::summarise(n = dplyr::n(),
                     prop_present = mean(Gastropod_Binary),
                     .groups = "drop")
)
glm_SEP_1_bi <- glm(Gastropod_Binary ~ Pad_Colour * Pad_Region,
                    data = SEP_1, family = binomial("logit"), method = "brglmFit")
print(summary(glm_SEP_1_bi))
diag_plots_binom(glm_SEP_1_bi, "Sep 1 (BRGLM) ")
glm_SEP_1_mle <- glm(Gastropod_Binary ~ Pad_Colour * Pad_Region,
                     data = SEP_1, family = binomial("logit"))
print(summary(glm_SEP_1_mle))
diag_plots_binom(glm_SEP_1_mle, "Sep 1 (MLE)   ")

# September 2
print(with(SEP_2, table(Pad_Colour, Pad_Region, Gastropod_Binary)))
print(
  SEP_2 %>%
    dplyr::group_by(Pad_Colour, Pad_Region) %>%
    dplyr::summarise(n = dplyr::n(),
                     prop_present = mean(Gastropod_Binary),
                     .groups = "drop")
)
glm_SEP_2_bi <- glm(Gastropod_Binary ~ Pad_Colour * Pad_Region,
                    data = SEP_2, family = binomial("logit"), method = "brglmFit")
print(summary(glm_SEP_2_bi))
diag_plots_binom(glm_SEP_2_bi, "Sep 2 (BRGLM) ")
glm_SEP_2_mle <- glm(Gastropod_Binary ~ Pad_Colour * Pad_Region,
                     data = SEP_2, family = binomial("logit"))
print(summary(glm_SEP_2_mle))
diag_plots_binom(glm_SEP_2_mle, "Sep 2 (MLE)   ")

########## 3) Whole-season GAM diagnostics  ##########

d_gam <- brighton_data %>%
  dplyr::mutate(
    Days_Since_Start = as.numeric(Sampling_Date - min(Sampling_Date, na.rm = TRUE)),
    Pad_ID = factor(Pad_ID)
  ) %>%
  dplyr::filter(!is.na(Sampling_Date))

# Binomial logit + 4 by-smooths
grazer_gam_final <- mgcv::gam(
  Gastropod_Binary ~
    s(Days_Since_Start, by = interaction(Pad_Colour, Pad_Region), k = 7) +
    Pad_Colour * Pad_Region +
    s(Pad_ID, bs = "re"),
  data   = d_gam,
  family = binomial(link = "logit"),
  method = "REML",
  select = TRUE
)

summary(grazer_gam_final)
gam.check(grazer_gam_final)