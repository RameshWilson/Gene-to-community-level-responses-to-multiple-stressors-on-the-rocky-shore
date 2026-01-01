# Green microalgae concentration_Additional
#
#   Evidence/diagnostics for Brighton green microalgae concentration.

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

brighton_data <- brighton_data %>%
  mutate(
    Pad_Colour = factor(Pad_Colour, levels = c("W","B")),
    Pad_Region = factor(Pad_Region, levels = c("N","P"))
  )

########## 2) June-only rationale ##########

# Subsequent dates with 0 values
zero_tbl <- brighton_data %>%
  filter(Sampling_Date %in% c(D_JUN, D_JUL1, D_JUL2, D_AUG1, D_AUG2, D_SEP1, D_SEP2)) %>%
  mutate(
    Date = case_when(
      Sampling_Date == D_JUN  ~ "June",
      Sampling_Date == D_JUL1 ~ "July 1",
      Sampling_Date == D_JUL2 ~ "July 2",
      Sampling_Date == D_AUG1 ~ "August 1",
      Sampling_Date == D_AUG2 ~ "August 2",
      Sampling_Date == D_SEP1 ~ "September 1",
      Sampling_Date == D_SEP2 ~ "September 2",
      TRUE ~ NA_character_
    )
  ) %>%
  group_by(Date) %>%
  summarise(
    n = n(),
    n_zero = sum(Green_Algae == 0, na.rm = TRUE),
    pct_zero = mean(Green_Algae == 0, na.rm = TRUE) * 100,
    min = min(Green_Algae, na.rm = TRUE),
    median = median(Green_Algae, na.rm = TRUE),
    max = max(Green_Algae, na.rm = TRUE),
    .groups = "drop"
  )

print(zero_tbl)

# Several subsequent dates have completely 0 values (or near 0)
# Gamma models require strictly positive responses, meaning Gamma would be invalid without adding a psuedocount (and thus changing inferences)
# Using alternative models (e.g. Gaussian) would lose interpretability and consistency of modelling approach throughout summer, and would provide little practical inference due to such low values regardless.
# Therefore, additional GLMs and any whole-season GAM are not valuable for statistical inference.

########## 3) June-only: Gamma(identity) vs Gaussian(identity) ##########

JUN <- dplyr::filter(brighton_data, Sampling_Date == D_JUN)
JUN$Pad_Colour <- relevel(JUN$Pad_Colour, ref = "W")
JUN$Pad_Region <- relevel(JUN$Pad_Region, ref = "N")

# Gaussian(identity)
glm_JUN_gaus <- glm(Green_Algae ~ Pad_Colour * Pad_Region,
                    data = JUN, family = gaussian(link = "identity"))
print(summary(glm_JUN_gaus))

par(mfrow = c(1, 2))
hist(residuals(glm_JUN_gaus, type = "pearson"),
     main = "June Gaussian: residual hist", xlab = "Pearson residuals",
     col = "grey80", border = "white")
qqnorm(residuals(glm_JUN_gaus, type = "pearson"),
       main = "June Gaussian: QQ plot")
qqline(residuals(glm_JUN_gaus, type = "pearson"))

# Gamma(identity)
glm_JUN_gamma <- glm(Green_Algae ~ Pad_Colour * Pad_Region,
                    data = JUN, family = Gamma(link = "identity"))
print(summary(glm_JUN_gamma))

par(mfrow = c(1, 2))
hist(residuals(glm_JUN_gamma, type = "pearson"),
     main = "June Gamma: residual hist", xlab = "Pearson residuals",
     col = "grey80", border = "white")
qqnorm(residuals(glm_JUN_gamma, type = "pearson"),
       main = "June Gamma: QQ plot")
qqline(residuals(glm_JUN_gaus, type = "pearson"))

print(AIC(glm_JUN_gamma, glm_JUN_gaus))
# Comparable visual diagnostics, and Gamma has a lower AIC; proceed with Gamma for June