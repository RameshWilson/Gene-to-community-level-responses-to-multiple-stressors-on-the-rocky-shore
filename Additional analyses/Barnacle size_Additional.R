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

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(lubridate)
  library(broom)
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
  mutate(
    Pad_Colour = factor(Pad_Colour, levels = c("W","B")),  # W=Ambient, B=Warmed
    Pad_Region = factor(Pad_Region, levels = c("N","P"))   # N=Non-polluted, P=Polluted
  )

########## 2) Per-date GLM diagnostics (Gaussian identity) ##########
# Residual diagnostics (Pearson residual histogram + QQ plot)

JUN   <- dplyr::filter(brighton_data, Sampling_Date == D_JUN)
JUL_1 <- dplyr::filter(brighton_data, Sampling_Date == D_JUL1)
JUL_2 <- dplyr::filter(brighton_data, Sampling_Date == D_JUL2)
AUG_1 <- dplyr::filter(brighton_data, Sampling_Date == D_AUG1)
AUG_2 <- dplyr::filter(brighton_data, Sampling_Date == D_AUG2)
SEP_1 <- dplyr::filter(brighton_data, Sampling_Date == D_SEP1)
SEP_2 <- dplyr::filter(brighton_data, Sampling_Date == D_SEP2)

# Ensure control references within each subset
JUN$Pad_Colour   <- relevel(JUN$Pad_Colour,   ref = "W"); JUN$Pad_Region   <- relevel(JUN$Pad_Region,   ref = "N")
JUL_1$Pad_Colour <- relevel(JUL_1$Pad_Colour, ref = "W"); JUL_1$Pad_Region <- relevel(JUL_1$Pad_Region, ref = "N")
JUL_2$Pad_Colour <- relevel(JUL_2$Pad_Colour, ref = "W"); JUL_2$Pad_Region <- relevel(JUL_2$Pad_Region, ref = "N")
AUG_1$Pad_Colour <- relevel(AUG_1$Pad_Colour, ref = "W"); AUG_1$Pad_Region <- relevel(AUG_1$Pad_Region, ref = "N")
AUG_2$Pad_Colour <- relevel(AUG_2$Pad_Colour, ref = "W"); AUG_2$Pad_Region <- relevel(AUG_2$Pad_Region, ref = "N")
SEP_1$Pad_Colour <- relevel(SEP_1$Pad_Colour, ref = "W"); SEP_1$Pad_Region <- relevel(SEP_1$Pad_Region, ref = "N")
SEP_2$Pad_Colour <- relevel(SEP_2$Pad_Colour, ref = "W"); SEP_2$Pad_Region <- relevel(SEP_2$Pad_Region, ref = "N")

glm_JUN   <- glm(AVG_Barnacle_Size ~ Pad_Colour * Pad_Region, data = JUN,   family = gaussian())
glm_JUL_1 <- glm(AVG_Barnacle_Size ~ Pad_Colour * Pad_Region, data = JUL_1, family = gaussian())
glm_JUL_2 <- glm(AVG_Barnacle_Size ~ Pad_Colour * Pad_Region, data = JUL_2, family = gaussian())
glm_AUG_1 <- glm(AVG_Barnacle_Size ~ Pad_Colour * Pad_Region, data = AUG_1, family = gaussian())
glm_AUG_2 <- glm(AVG_Barnacle_Size ~ Pad_Colour * Pad_Region, data = AUG_2, family = gaussian())
glm_SEP_1 <- glm(AVG_Barnacle_Size ~ Pad_Colour * Pad_Region, data = SEP_1, family = gaussian())
glm_SEP_2 <- glm(AVG_Barnacle_Size ~ Pad_Colour * Pad_Region, data = SEP_2, family = gaussian())

par(mfrow = c(1,2))
r <- residuals(glm_JUN, type="pearson");   hist(r, main="June: Residual histogram",       xlab="Pearson residuals"); qqnorm(r, main="June: QQ");       qqline(r)

par(mfrow = c(1,2))
r <- residuals(glm_JUL_1, type="pearson"); hist(r, main="July 1: Residual histogram",    xlab="Pearson residuals"); qqnorm(r, main="July 1: QQ");    qqline(r)

par(mfrow = c(1,2))
r <- residuals(glm_JUL_2, type="pearson"); hist(r, main="July 2: Residual histogram",    xlab="Pearson residuals"); qqnorm(r, main="July 2: QQ");    qqline(r)

par(mfrow = c(1,2))
r <- residuals(glm_AUG_1, type="pearson"); hist(r, main="August 1: Residual histogram",  xlab="Pearson residuals"); qqnorm(r, main="August 1: QQ");  qqline(r)

par(mfrow = c(1,2))
r <- residuals(glm_AUG_2, type="pearson"); hist(r, main="August 2: Residual histogram",  xlab="Pearson residuals"); qqnorm(r, main="August 2: QQ");  qqline(r)

par(mfrow = c(1,2))
r <- residuals(glm_SEP_1, type="pearson"); hist(r, main="September 1: Residual histogram", xlab="Pearson residuals"); qqnorm(r, main="September 1: QQ"); qqline(r)

par(mfrow = c(1,2))
r <- residuals(glm_SEP_2, type="pearson"); hist(r, main="September 2: Residual histogram", xlab="Pearson residuals"); qqnorm(r, main="September 2: QQ"); qqline(r)

par(mfrow = c(1,1))

########## 3) Whole-season GAM diagnostics  ##########

# Whole-season dataset for GAM fitting
brighton_data_GAM <- brighton_data %>%
  dplyr::select(Sampling_Date, Pad_Colour, Pad_Region, AVG_Barnacle_Size, Pad_ID) %>%
  mutate(
    Sampling_Date    = as.Date(Sampling_Date, tryFormats = c("%d/%m/%Y", "%Y-%m-%d")),
    Pad_Colour       = factor(Pad_Colour, levels = c("W","B")),
    Pad_Region       = factor(Pad_Region, levels = c("N","P")),
    Pad_ID           = factor(Pad_ID),
    Days_Since_Start = as.numeric(Sampling_Date - min(Sampling_Date, na.rm = TRUE))
  ) %>%
  filter(!is.na(Sampling_Date))

# Full GAM model
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
gam.check(barnacle_size_gam)