# Barnacle abundance_Additional
#
#   Evidence/diagnostics for Brighton barnacle abundance.

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
  library(MASS)
  library(AER)
  library(pscl)
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

# Parse dates: try dmy first; fallback to ymd if needed
if (!inherits(brighton_data$Sampling_Date, "Date")) {
  tmp_dmy <- suppressWarnings(lubridate::dmy(brighton_data$Sampling_Date))
  frac_na <- mean(is.na(tmp_dmy))
  brighton_data$Sampling_Date <- if (frac_na < 0.5) tmp_dmy else suppressWarnings(lubridate::ymd(brighton_data$Sampling_Date))
}

# Recode factors with controls first
brighton_data <- brighton_data %>%
  mutate(
    Pad_Colour = factor(Pad_Colour, levels = c("W", "B")),
    Pad_Region = factor(Pad_Region, levels = c("N", "P"))
  )

########## 2) Per-date GLMs with diagnostics (Poisson vs NB) ##########

# Residual diagnostics (Pearson residual histogram + QQ plot)
diag_plots <- function(m, main_prefix = "") {
  op <- par(mfrow = c(1, 2)); on.exit(par(op), add = TRUE)
  r <- residuals(m, type = "pearson")
  hist(r, main = paste(main_prefix, "Histogram of Residuals"), xlab = "Pearson Residuals")
  qqnorm(r, main = paste(main_prefix, "QQ Plot of Residuals")); qqline(r)
  invisible(NULL)
}

# Subsets per date
JUN   <- dplyr::filter(brighton_data, Sampling_Date == D_JUN)
JUL_1 <- dplyr::filter(brighton_data, Sampling_Date == D_JUL1)
JUL_2 <- dplyr::filter(brighton_data, Sampling_Date == D_JUL2)
AUG_1 <- dplyr::filter(brighton_data, Sampling_Date == D_AUG1)
AUG_2 <- dplyr::filter(brighton_data, Sampling_Date == D_AUG2)
SEP_1 <- dplyr::filter(brighton_data, Sampling_Date == D_SEP1)
SEP_2 <- dplyr::filter(brighton_data, Sampling_Date == D_SEP2)

# Relevel within each subset
JUN$Pad_Colour   <- relevel(as.factor(JUN$Pad_Colour), ref = "W")
JUN$Pad_Region   <- relevel(as.factor(JUN$Pad_Region), ref = "N")

JUL_1$Pad_Colour <- relevel(as.factor(JUL_1$Pad_Colour), ref = "W")
JUL_1$Pad_Region <- relevel(as.factor(JUL_1$Pad_Region), ref = "N")

JUL_2$Pad_Colour <- relevel(as.factor(JUL_2$Pad_Colour), ref = "W")
JUL_2$Pad_Region <- relevel(as.factor(JUL_2$Pad_Region), ref = "N")

AUG_1$Pad_Colour <- relevel(as.factor(AUG_1$Pad_Colour), ref = "W")
AUG_1$Pad_Region <- relevel(as.factor(AUG_1$Pad_Region), ref = "N")

AUG_2$Pad_Colour <- relevel(as.factor(AUG_2$Pad_Colour), ref = "W")
AUG_2$Pad_Region <- relevel(as.factor(AUG_2$Pad_Region), ref = "N")

SEP_1$Pad_Colour <- relevel(as.factor(SEP_1$Pad_Colour), ref = "W")
SEP_1$Pad_Region <- relevel(as.factor(SEP_1$Pad_Region), ref = "N")

SEP_2$Pad_Colour <- relevel(as.factor(SEP_2$Pad_Colour), ref = "W")
SEP_2$Pad_Region <- relevel(as.factor(SEP_2$Pad_Region), ref = "N")

# Fit Poisson and NB for each date; run diagnostics
glm_JUN    <- glm(Barnacle_Count ~ Pad_Colour * Pad_Region, data = JUN, family = poisson())
summary(glm_JUN); AER::dispersiontest(glm_JUN); diag_plots(glm_JUN, "June (Poisson)  ")

glm.nb_JUN <- MASS::glm.nb(Barnacle_Count ~ Pad_Colour * Pad_Region, data = JUN)
summary(glm.nb_JUN); pscl::odTest(glm.nb_JUN); diag_plots(glm.nb_JUN, "June (NB)       ")

glm_JUL_1    <- glm(Barnacle_Count ~ Pad_Colour * Pad_Region, data = JUL_1, family = poisson())
summary(glm_JUL_1); AER::dispersiontest(glm_JUL_1); diag_plots(glm_JUL_1, "July 1 (Poisson)")

glm.nb_JUL_1 <- MASS::glm.nb(Barnacle_Count ~ Pad_Colour * Pad_Region, data = JUL_1)
summary(glm.nb_JUL_1); pscl::odTest(glm.nb_JUL_1); diag_plots(glm.nb_JUL_1, "July 1 (NB)     ")

glm_JUL_2    <- glm(Barnacle_Count ~ Pad_Colour * Pad_Region, data = JUL_2, family = poisson())
summary(glm_JUL_2); AER::dispersiontest(glm_JUL_2); diag_plots(glm_JUL_2, "July 2 (Poisson)")

glm.nb_JUL_2 <- MASS::glm.nb(Barnacle_Count ~ Pad_Colour * Pad_Region, data = JUL_2)
summary(glm.nb_JUL_2); pscl::odTest(glm.nb_JUL_2); diag_plots(glm.nb_JUL_2, "July 2 (NB)     ")

glm_AUG_1    <- glm(Barnacle_Count ~ Pad_Colour * Pad_Region, data = AUG_1, family = poisson())
summary(glm_AUG_1); AER::dispersiontest(glm_AUG_1); diag_plots(glm_AUG_1, "Aug 1 (Poisson) ")

glm.nb_AUG_1 <- MASS::glm.nb(Barnacle_Count ~ Pad_Colour * Pad_Region, data = AUG_1)
summary(glm.nb_AUG_1); pscl::odTest(glm.nb_AUG_1); diag_plots(glm.nb_AUG_1, "Aug 1 (NB)      ")

glm_AUG_2    <- glm(Barnacle_Count ~ Pad_Colour * Pad_Region, data = AUG_2, family = poisson())
summary(glm_AUG_2); AER::dispersiontest(glm_AUG_2); diag_plots(glm_AUG_2, "Aug 2 (Poisson) ")

glm.nb_AUG_2 <- MASS::glm.nb(Barnacle_Count ~ Pad_Colour * Pad_Region, data = AUG_2)
summary(glm.nb_AUG_2); pscl::odTest(glm.nb_AUG_2); diag_plots(glm.nb_AUG_2, "Aug 2 (NB)      ")

glm_SEP_1    <- glm(Barnacle_Count ~ Pad_Colour * Pad_Region, data = SEP_1, family = poisson())
summary(glm_SEP_1); AER::dispersiontest(glm_SEP_1); diag_plots(glm_SEP_1, "Sep 1 (Poisson) ")

glm.nb_SEP_1 <- MASS::glm.nb(Barnacle_Count ~ Pad_Colour * Pad_Region, data = SEP_1)
summary(glm.nb_SEP_1); pscl::odTest(glm.nb_SEP_1); diag_plots(glm.nb_SEP_1, "Sep 1 (NB)      ")

glm_SEP_2    <- glm(Barnacle_Count ~ Pad_Colour * Pad_Region, data = SEP_2, family = poisson())
summary(glm_SEP_2); AER::dispersiontest(glm_SEP_2); diag_plots(glm_SEP_2, "Sep 2 (Poisson) ")

glm.nb_SEP_2 <- MASS::glm.nb(Barnacle_Count ~ Pad_Colour * Pad_Region, data = SEP_2)
summary(glm.nb_SEP_2); pscl::odTest(glm.nb_SEP_2); diag_plots(glm.nb_SEP_2, "Sep 2 (NB)      ")

# AIC comparison table
aic_tbl <- tibble::tibble(
  Date = DATE_LABELS,
  AIC_Poisson = c(AIC(glm_JUN), AIC(glm_JUL_1), AIC(glm_JUL_2), AIC(glm_AUG_1), AIC(glm_AUG_2), AIC(glm_SEP_1), AIC(glm_SEP_2)),
  AIC_NB      = c(AIC(glm.nb_JUN), AIC(glm.nb_JUL_1), AIC(glm.nb_JUL_2), AIC(glm.nb_AUG_1), AIC(glm.nb_AUG_2), AIC(glm.nb_SEP_1), AIC(glm.nb_SEP_2))
) %>%
  mutate(Delta_AIC = AIC_Poisson - AIC_NB)

print(aic_tbl)

########## 3) Whole-season GAM diagnostics  ##########

# Whole-season dataset for GAM fitting
brighton_data_GAM <- brighton_data %>%
  dplyr::select(Sampling_Date, Pad_Colour, Pad_Region, Barnacle_Count, Pad_ID) %>%  # explicit to avoid masking
  mutate(
    Sampling_Date    = as.Date(Sampling_Date),
    Days_Since_Start = as.numeric(Sampling_Date - min(Sampling_Date, na.rm = TRUE)),
    Pad_ID           = factor(Pad_ID)
  ) %>%
  filter(!is.na(Sampling_Date))

# Full GAM model
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
gam.check(barnacle_gam_model)