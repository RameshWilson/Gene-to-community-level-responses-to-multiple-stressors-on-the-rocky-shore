# Diatoms concentration_Additional
#
#   Evidence/diagnostics for Brighton cyanobacteria concentration.

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

# Ensure factors are coded with controls first
brighton_data <- brighton_data %>%
  mutate(
    Pad_Colour = factor(Pad_Colour, levels = c("W","B")),  # W=Ambient, B=Warm
    Pad_Region = factor(Pad_Region, levels = c("N","P"))   # N=Non-Polluted, P=Polluted
  )

# Gamma GLMs/GAMs assume y > 0 (strictly positive). Diatoms includes zeros in this dataset.
cat("Min Diatoms:", min(brighton_data$Diatoms, na.rm = TRUE), "\n")
cat("N(Diatoms == 0):", sum(brighton_data$Diatoms == 0, na.rm = TRUE), "\n")
# We do not add a pseudo-count (e.g., y + 0.001) as that changes the scale and can bias inference at low concentrations.
# We therefore proceed with diagnostic checks for Gaussian (identity)

########## 3) Per-date GLMs: Gaussian(identity) + residual diagnostics ##########

# Subsets per date
JUN   <- filter(brighton_data, Sampling_Date == D_JUN)
JUL_1 <- filter(brighton_data, Sampling_Date == D_JUL1)
JUL_2 <- filter(brighton_data, Sampling_Date == D_JUL2)
AUG_1 <- filter(brighton_data, Sampling_Date == D_AUG1)
AUG_2 <- filter(brighton_data, Sampling_Date == D_AUG2)
SEP_1 <- filter(brighton_data, Sampling_Date == D_SEP1)
SEP_2 <- filter(brighton_data, Sampling_Date == D_SEP2)

# Relevel within each subset (control first)
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
  hist(rp, main = paste0(main_prefix, "Pearson residuals (hist)"), xlab = "Pearson residuals",
       col = "grey80", border = "white")
  qqnorm(rp, main = paste0(main_prefix, "Pearson residuals (QQ)")); qqline(rp)
  invisible(NULL)
}

glm_JUN <- glm(Diatoms ~ Pad_Colour * Pad_Region, data = JUN, family = gaussian(link = "identity"))
print(summary(glm_JUN))
diag_plots(glm_JUN, "June  ")

glm_JUL1 <- glm(Diatoms ~ Pad_Colour * Pad_Region, data = JUL_1, family = gaussian(link = "identity"))
print(summary(glm_JUL1))
diag_plots(glm_JUL1, "July1 ")

glm_JUL2 <- glm(Diatoms ~ Pad_Colour * Pad_Region, data = JUL_2, family = gaussian(link = "identity"))
print(summary(glm_JUL2))
diag_plots(glm_JUL2, "July2 ")

glm_AUG1 <- glm(Diatoms ~ Pad_Colour * Pad_Region, data = AUG_1, family = gaussian(link = "identity"))
print(summary(glm_AUG1))
diag_plots(glm_AUG1, "Aug1  ")

glm_AUG2 <- glm(Diatoms ~ Pad_Colour * Pad_Region, data = AUG_2, family = gaussian(link = "identity"))
print(summary(glm_AUG2))
diag_plots(glm_AUG2, "Aug2  ")

glm_SEP1 <- glm(Diatoms ~ Pad_Colour * Pad_Region, data = SEP_1, family = gaussian(link = "identity"))
print(summary(glm_SEP1))
diag_plots(glm_SEP1, "Sep1  ")

glm_SEP2 <- glm(Diatoms ~ Pad_Colour * Pad_Region, data = SEP_2, family = gaussian(link = "identity"))
print(summary(glm_SEP2))
diag_plots(glm_SEP2, "Sep2  ")

########## 4) Whole-season GAM data prep ##########

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

########## 5) Whole-season GAM (Gaussian identity) + diagnostics ##########

diatoms_gam_model <- mgcv::gam(
  Diatoms ~
    s(Days_Since_Start, by = interaction(Pad_Colour, Pad_Region), k = 7) +
    Pad_Colour * Pad_Region +
    s(Pad_ID, bs = "re"),
  data   = brighton_data_GAM,
  family = gaussian(link = "identity"),
  method = "REML"
)

summary(diatoms_gam_model)
gam.check(diatoms_gam_model)