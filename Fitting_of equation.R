library(dplyr)

# -------------------------------
# 1. INPUT DATA
# -------------------------------
df <- tribble(
  ~RH, ~Temp, ~FL, ~NC, ~CA,
  92, 30, 1, 1, 1,
  92, 27, 1, 1, 1,
  92, 23, 1, 1, 1,
  92, 20, 1, 1, 1,
  75, 30, 1, 1, 1,
  75, 27, 1, 1, 1,
  75, 23, 1, 1, 1,
  75, 20, 1, 1, 1,
  52, 30, 1, 1, 1,
  52, 27, 1, 1, 1,
  52, 23, 1, 1, 1,
  52, 20, 0.666666667, 1, 1,
  33, 30, 1, 1, 1,
  33, 27, 1, 1, 1,
  33, 23, 1, 1, 1,
  33, 20, 0.666666667, 1, 1
)

# --------------------------------------
# 2. FUNCTION TO FIT MODELS + REPORT R²
# --------------------------------------
fit_model <- function(column) {
  
  form <- as.formula(paste(column, "~ Temp + RH"))
  
  model <- lm(form, data = df)
  
  r2   <- summary(model)$r.squared
  rmse <- sqrt(mean(residuals(model)^2))
  
  list(model = model, R2 = r2, RMSE = rmse)
}

# --------------------------------------
# 3. FIT MODELS FOR EACH STRAIN
# --------------------------------------
fit_FL <- fit_model("FL")
fit_NC <- fit_model("NC")
fit_CA <- fit_model("CA")

# --------------------------------------
# 4. PRINT RESULTS
# --------------------------------------
cat("\n=== FL STRAIN ===\n")
print(summary(fit_FL$model))
cat("R² =", fit_FL$R2, "\nRMSE =", fit_FL$RMSE, "\n")

cat("\n=== NC STRAIN ===\n")
print(summary(fit_NC$model))
cat("R² =", fit_NC$R2, "\nRMSE =", fit_NC$RMSE, "\n")

cat("\n=== CA STRAIN ===\n")
print(summary(fit_CA$model))
cat("R² =", fit_CA$R2, "\nRMSE =", fit_CA$RMSE, "\n")
library(nls2)
install.packages("nls2")
library(nls2)
# -----------------------
# FL strain data
# -----------------------
df_FL <- data.frame(
  T = c(20, 23, 27, 30),
  S = c(0.67, 1, 1, 1)
)

# -----------------------
# Logistic function
# -----------------------
logistic_fun <- function(T, a, b) {
  1 / (1 + exp(-(a + b*T)))
}

# -----------------------
# Fit logistic using nls2
# -----------------------
fit_FL <- nls2(
  S ~ 1 / (1 + exp(-(a + b*T))),
  data = df_FL,
  start = list(a = -5, b = 0.2),
  algorithm = "brute-force"
)

summary(fit_FL)
coef(fit_FL)
########################
library(dplyr)
library(Metrics)

# -----------------------------
# 1. INPUT DATA
# -----------------------------
df <- tribble(
  ~Temp, ~RH, ~mu_CA, ~mu_NC, ~mu_FL,
  20, 33, 0.2857143, 0.16095238, 0.24952381,
  20, 52, 0.1900000, 0.07809524, 0.15619048,
  20, 75, 0.03095238, 0.01666667, 0.02904762,
  20, 92, 0.01095238, 0.01000000, 0.01285714,
  23, 33, 0.2857143, 0.15000000, 0.28571429,
  23, 52, 0.19190476, 0.06380952, 0.28571429,
  23, 75, 0.03000000, 0.02428571, 0.02619048,
  23, 92, 0.01809524, 0.01380952, 0.01666667,
  27, 33, 0.2857143, 0.26523810, 0.28571429,
  27, 52, 0.2857143, 0.20285714, 0.28571429,
  27, 75, 0.06000000, 0.02761905, 0.03476190,
  27, 92, 0.03285714, 0.01619048, 0.01619048,
  30, 33, 0.2857143, 0.27904762, 0.28571429,
  30, 52, 0.2857143, 0.17571429, 0.28571429,
  30, 75, 0.06380952, 0.06047619, 0.03476190,
  30, 92, 0.04523810, 0.01714286, 0.02000000,
  35, 33, 0.2857143, 0.28571429, 0.28571429,
  35, 52, 0.2857143, 0.23428571, 0.28571429,
  35, 75, 0.14714286, 0.04619048, 0.14333333,
  35, 92, 0.04857143, 0.03095238, 0.03095238
)

# -----------------------------
# 2. FITTING FUNCTION
# -----------------------------
fit_model <- function(column) {
  
  form <- as.formula(paste(column, "~ Temp + RH"))
  model <- lm(form, data = df)
  
  r2   <- summary(model)$r.squared
  rmse <- sqrt(mean(residuals(model)^2))
  
  list(model = model, R2 = r2, RMSE = rmse)
}

# -----------------------------
# 3. FIT MODELS
# -----------------------------
fit_CA <- fit_model("mu_CA")
fit_NC <- fit_model("mu_NC")
fit_FL <- fit_model("mu_FL")

# -----------------------------
# 4. PRINT RESULTS
# -----------------------------
cat("\n=== CA Strain ===\n")
print(summary(fit_CA$model))
cat("R² =", fit_CA$R2, "\nRMSE =", fit_CA$RMSE, "\n")

cat("\n=== NC Strain ===\n")
print(summary(fit_NC$model))
cat("R² =", fit_NC$R2, "\nRMSE =", fit_NC$RMSE, "\n")

cat("\n=== FL Strain ===\n")
print(summary(fit_FL$model))
cat("R² =", fit_FL$R2, "\nRMSE =", fit_FL$RMSE, "\n")
###########################################################
library(dplyr)
library(tibble)

# ------------------ Build the Nymph dataset ------------------
nymph <- tribble(
  ~Temp, ~RH, ~mu_CA, ~mu_NC, ~mu_FL,
  20,33,0.024285714,0.024285714,0.025238095,
  20,52,0.019047619,0.019047619,0.013333333,
  20,75,0.011904762,0.011904762,0.009047619,
  20,92,0.007619048,0.007619048,0.009523810,
  23,33,0.033809524,0.033809524,0.023809524,
  23,52,0.020476190,0.020476190,0.015714286,
  23,75,0.015238095,0.015238095,0.008095238,
  23,92,0.013809524,0.013809524,0.008095238,
  27,33,0.064285714,0.064285714,0.027142857,
  27,52,0.040476190,0.040476190,0.018095238,
  27,75,0.033333333,0.033333333,0.010000000,
  27,92,0.016666667,0.016666667,0.009523810,
  30,33,0.038571429,0.038571429,0.044761905,
  30,52,0.044285714,0.044285714,0.026666667,
  30,75,0.028095238,0.028095238,0.013333333,
  30,92,0.022857143,0.022857143,0.011904762,
  35,33,0.168095238,0.168095238,0.253333333,
  35,52,0.123809524,0.123809524,0.049523810,
  35,75,0.057619048,0.057619048,0.020000000,
  35,92,0.053333333,0.053333333,0.013333333
)

# ------------------ Fit Function ------------------
fit_fun <- function(y, df){
  model <- lm(y ~ Temp + RH, data = df)
  pred <- predict(model)
  R2 <- cor(pred, y)^2
  RMSE <- sqrt(mean((pred - y)^2))
  list(model = model, R2 = R2, RMSE = RMSE)
}

# ------------------ Run the Models ------------------
fit_CA_n <- fit_fun(nymph$mu_CA, nymph)
fit_NC_n <- fit_fun(nymph$mu_NC, nymph)
fit_FL_n <- fit_fun(nymph$mu_FL, nymph)

# ------------------ Print Results ------------------
cat("\n=== CA STRAIN ===\n")
print(summary(fit_CA_n$model))
cat("R² =", fit_CA_n$R2, "\nRMSE =", fit_CA_n$RMSE, "\n")

cat("\n=== NC STRAIN ===\n")
print(summary(fit_NC_n$model))
cat("R² =", fit_NC_n$R2, "\nRMSE =", fit_NC_n$RMSE, "\n")

cat("\n=== FL STRAIN ===\n")
print(summary(fit_FL_n$model))
cat("R² =", fit_FL_n$R2, "\nRMSE =", fit_FL_n$RMSE, "\n")
###############
# -----------------------------
# 1. Load packages
# -----------------------------
library(dplyr)
library(Metrics)

# -----------------------------
# 2. Adult mortality dataset
# -----------------------------
df <- tribble(
  ~Temp, ~RH, ~mu_CA, ~mu_NC, ~mu_FL,
  20, 33, 0.007619048, 0.006666667, 0.005238095,
  20, 52, 0.006666667, 0.006666667, 0.003809524,
  20, 75, 0.006666667, 0.002857143, 0.002380952,
  20, 92, 0.006190476, 0.002857143, 0.003333333,
  23, 33, 0.005714286, 0.005714286, 0.005238095,
  23, 52, 0.004285714, 0.003809524, 0.007142857,
  23, 75, 0.004285714, 0.002857143, 0.003809524,
  23, 92, 0.005714286, 0.002857143, 0.003333333,
  30, 33, 0.019047619, 0.010476190, 0.006666667,
  30, 52, 0.013333333, 0.007142857, 0.006190476,
  30, 75, 0.010476190, 0.002857143, 0.004285714,
  30, 92, 0.012857143, 0.003333333, 0.004285714,
  35, 33, 0.031904762, 0.013333333, 0.017619048,
  35, 52, 0.025714286, 0.011904762, 0.012380952,
  35, 75, 0.017142857, 0.005238095, 0.010476190,
  35, 92, 0.012857143, 0.004285714, 0.007619048
)

# -----------------------------
# 3. Function to fit linear model
# -----------------------------
fit_model <- function(y_col) {
  form <- as.formula(paste(y_col, "~ Temp + RH"))
  model <- lm(form, data = df)
  
  # Calculate metrics
  r2   <- summary(model)$r.squared
  rmse <- rmse(df[[y_col]], predict(model))
  
  list(model = model, R2 = r2, RMSE = rmse)
}

# -----------------------------
# 4. Fit models for CA, NC, FL
# -----------------------------
fit_CA <- fit_model("mu_CA")
fit_NC <- fit_model("mu_NC")
fit_FL <- fit_model("mu_FL")

# -----------------------------
# 5. Print model summaries
# -----------------------------
cat("\n=== CA STRAIN ===\n")
print(summary(fit_CA$model))
cat("R² =", fit_CA$R2, " RMSE =", fit_CA$RMSE, "\n")

cat("\n=== NC STRAIN ===\n")
print(summary(fit_NC$model))
cat("R² =", fit_NC$R2, " RMSE =", fit_NC$RMSE, "\n")

cat("\n=== FL STRAIN ===\n")
print(summary(fit_FL$model))
cat("R² =", fit_FL$R2, " RMSE =", fit_FL$RMSE, "\n")

# -----------------------------
# 6. Create results table
# -----------------------------
results <- tibble(
  Strain = c("CA", "NC", "FL"),
  Equation = c(
    paste0("mu_CA = ",
           round(coef(fit_CA$model)[1],6)," + ",
           round(coef(fit_CA$model)[2],6),"*Temp - ",
           abs(round(coef(fit_CA$model)[3],6)),"*RH"),
    
    paste0("mu_NC = ",
           round(coef(fit_NC$model)[1],6)," + ",
           round(coef(fit_NC$model)[2],6),"*Temp - ",
           abs(round(coef(fit_NC$model)[3],6)),"*RH"),
    
    paste0("mu_FL = ",
           round(coef(fit_FL$model)[1],6)," + ",
           round(coef(fit_FL$model)[2],6),"*Temp - ",
           abs(round(coef(fit_FL$model)[3],6)),"*RH")
  ),
  R2   = c(fit_CA$R2, fit_NC$R2, fit_FL$R2),
  RMSE = c(fit_CA$RMSE, fit_NC$RMSE, fit_FL$RMSE)
)

print(results)
########################################
library(dplyr)
library(Metrics)

# -----------------------------
# 1. Larval mortality dataset
# -----------------------------
df <- tribble(
  ~Temp, ~RH, ~mu_CA, ~mu_NC, ~mu_FL,
  20, 33, 0.285714286, 0.160952381, 0.249523810,
  20, 52, 0.190000000, 0.078095238, 0.156190476,
  20, 75, 0.030952381, 0.016666667, 0.029047619,
  20, 92, 0.010952381, 0.010000000, 0.012857143,
  
  23, 33, 0.285714286, 0.150000000, 0.285714286,
  23, 52, 0.191904762, 0.063809524, 0.285714286,
  23, 75, 0.030000000, 0.024285714, 0.026190476,
  23, 92, 0.018095238, 0.013809524, 0.016666667,
  
  27, 33, 0.285714286, 0.265238095, 0.285714286,
  27, 52, 0.285714286, 0.202857143, 0.285714286,
  27, 75, 0.060000000, 0.027619048, 0.034761905,
  27, 92, 0.032857143, 0.016190476, 0.016190476,
  
  30, 33, 0.285714286, 0.279047619, 0.285714286,
  30, 52, 0.285714286, 0.175714286, 0.285714286,
  30, 75, 0.063809524, 0.060476190, 0.034761905,
  30, 92, 0.045238095, 0.017142857, 0.020000000,
  
  35, 33, 0.285714286, 0.285714286, 0.285714286,
  35, 52, 0.285714286, 0.234285714, 0.285714286,
  35, 75, 0.147142857, 0.046190476, 0.143333333,
  35, 92, 0.048571429, 0.030952381, 0.030952381
)

# -----------------------------
# 2. Function to fit model
# -----------------------------
fit_model <- function(y_col){
  form <- as.formula(paste(y_col, "~ Temp + RH"))
  model <- lm(form, data=df)
  
  r2   <- summary(model)$r.squared
  rmse <- rmse(df[[y_col]], predict(model))
  
  list(model=model, R2=r2, RMSE=rmse)
}

# -----------------------------
# 3. Fit the models
# -----------------------------
fit_CA <- fit_model("mu_CA")
fit_NC <- fit_model("mu_NC")
fit_FL <- fit_model("mu_FL")

# -----------------------------
# 4. Print summaries
# -----------------------------
cat("\n=== CA STRAIN ===\n")
print(summary(fit_CA$model))
cat("R² =", fit_CA$R2, " RMSE =", fit_CA$RMSE, "\n")

cat("\n=== NC STRAIN ===\n")
print(summary(fit_NC$model))
cat("R² =", fit_NC$R2, " RMSE =", fit_NC$RMSE, "\n")

cat("\n=== FL STRAIN ===\n")
print(summary(fit_FL$model))
cat("R² =", fit_FL$R2, " RMSE =", fit_FL$RMSE, "\n")

# -----------------------------
# 5. Word-ready results table
# -----------------------------
results <- tibble(
  Strain = c("CA", "NC", "FL"),
  Equation = c(
    paste0("mu_CA = ", round(coef(fit_CA$model)[1],6),
           " + ", round(coef(fit_CA$model)[2],6),"*Temp - ",
           abs(round(coef(fit_CA$model)[3],6)),"*RH"),
    
    paste0("mu_NC = ", round(coef(fit_NC$model)[1],6),
           " + ", round(coef(fit_NC$model)[2],6),"*Temp - ",
           abs(round(coef(fit_NC$model)[3],6)),"*RH"),
    
    paste0("mu_FL = ", round(coef(fit_FL$model)[1],6),
           " + ", round(coef(fit_FL$model)[2],6),"*Temp - ",
           abs(round(coef(fit_FL$model)[3],6)),"*RH")
  ),
  R2 = c(fit_CA$R2, fit_NC$R2, fit_FL$R2),
  RMSE = c(fit_CA$RMSE, fit_NC$RMSE, fit_FL$RMSE)
)

print(results)
########################
library(dplyr)
library(Metrics)

# ================================
# 1. INPUT DATA (Egg development time)
# ================================
df <- tribble(
  ~RH, ~Temp, ~FL, ~NC, ~CA,
  75, 30, 22, 20, 19.5,
  75, 27, 25, 22, 22.5,
  75, 23, 42, 40, 41.5,
  75, 20, 80, 55, 47,
  
  52, 30, 21, 20, 20,
  52, 27, 23, 21, 22,
  52, 23, 42, 38, 39.5,
  52, 20, 80, 60, 45,
  
  33, 30, 20, 19.5, 19.5,
  33, 27, 25, 22, 24,
  33, 23, 42, 38, 40,
  33, 20, 84, 60, 48
)

# ======================================
# 2. CONVERT DEVELOPMENTAL TIME → RATE
# ======================================
df <- df %>%
  mutate(
    FL_rate = 1 / FL,
    NC_rate = 1 / NC,
    CA_rate = 1 / CA
  )

# ======================================
# 3. FUNCTION TO FIT MODEL + METRICS
# ======================================
fit_model <- function(rate_col) {
  
  form <- as.formula(paste(rate_col, "~ Temp + RH"))
  model <- lm(form, data = df)
  
  r2 <- summary(model)$r.squared
  rmse <- sqrt(mean(residuals(model)^2))
  
  list(model = model, R2 = r2, RMSE = rmse)
}

# ======================================
# 4. RUN MODELS FOR EACH STRAIN
# ======================================
fit_FL <- fit_model("FL_rate")
fit_NC <- fit_model("NC_rate")
fit_CA <- fit_model("CA_rate")

# ======================================
# 5. PRINT RESULTS
# ======================================
cat("\n=== FL STRAIN ===\n")
print(summary(fit_FL$model))
cat("R² =", fit_FL$R2, " RMSE =", fit_FL$RMSE, "\n")

cat("\n=== NC STRAIN ===\n")
print(summary(fit_NC$model))
cat("R² =", fit_NC$R2, " RMSE =", fit_NC$RMSE, "\n")

cat("\n=== CA STRAIN ===\n")
print(summary(fit_CA$model))
cat("R² =", fit_CA$R2, " RMSE =", fit_CA$RMSE, "\n")
###########################################
# =============================================================
# LARVAL DEVELOPMENTAL RATE MODEL (FL, NC, CA)
# =============================================================

library(dplyr)

# 1. Raw dataset
df <- tribble(
  ~RH, ~Temp, ~FL_time, ~NC_time, ~CA_time,
  75, 35, 7.5, 8, 7.5,
  75, 30, 9, 8, 9,
  75, 27, 10, 9, 10,
  75, 23, 18, 17.5, 17,
  75, 20, 27, 18, 24,
  52, 35, 7, 7.5, 7,
  52, 30, 7.5, 7, 7,
  52, 27, 10, 9, 10,
  52, 23, 18, 19, 17.5,
  52, 20, 27, 19, 24,
  33, 35, 7, 8, 7,
  33, 30, 8, 7, 8,
  33, 27, 10, 9, 10,
  33, 23, 17.5, 17.5, 13,
  33, 20, 27, 17, 26
)

# 2. Convert time (days) → rate (1/day)
df <- df %>%
  mutate(
    FL = 1 / FL_time,
    NC = 1 / NC_time,
    CA = 1 / CA_time
  )

# 3. Fit models for each strain
form <- as.formula("rate ~ Temp + RH")

fit_FL <- list(
  model = lm(rate ~ Temp + RH, data = df %>% rename(rate = FL)),
  R2    = summary(lm(rate ~ Temp + RH, data = df %>% rename(rate = FL)))$r.squared,
  RMSE  = sigma(lm(rate ~ Temp + RH, data = df %>% rename(rate = FL)))
)

fit_NC <- list(
  model = lm(rate ~ Temp + RH, data = df %>% rename(rate = NC)),
  R2    = summary(lm(rate ~ Temp + RH, data = df %>% rename(rate = NC)))$r.squared,
  RMSE  = sigma(lm(rate ~ Temp + RH, data = df %>% rename(rate = NC)))
)

fit_CA <- list(
  model = lm(rate ~ Temp + RH, data = df %>% rename(rate = CA)),
  R2    = summary(lm(rate ~ Temp + RH, data = df %>% rename(rate = CA)))$r.squared,
  RMSE  = sigma(lm(rate ~ Temp + RH, data = df %>% rename(rate = CA)))
)

# 4. Print summaries
cat("\n=== FL STRAIN ===\n")
print(summary(fit_FL$model))
cat("R² =", fit_FL$R2, " RMSE =", fit_FL$RMSE, "\n")

cat("\n=== NC STRAIN ===\n")
print(summary(fit_NC$model))
cat("R² =", fit_NC$R2, " RMSE =", fit_NC$RMSE, "\n")

cat("\n=== CA STRAIN ===\n")
print(summary(fit_CA$model))
cat("R² =", fit_CA$R2, " RMSE =", fit_CA$RMSE, "\n")
####################################
library(dplyr)
#Nymph

# ===============================
# 1. Create the dataset
# ===============================
df <- tribble(
  ~RH, ~Temp, ~FL_days, ~NC_days, ~CA_days,
  75, 35, 8, 11, 10,
  75, 30, 8, 14, 13,
  75, 27, 17, 18, 17,
  75, 23, 30, 27, 27,
  75, 20, 50, 37, 40,
  52, 35, 8, 11, 10,
  52, 30, 8, 14, 13,
  52, 27, 18, 19, 18,
  52, 23, 30, 27, 27.5,
  52, 20, 50, 37, 44,
  33, 35, 8, 12, 10,
  33, 30, 8, 14, 13,
  33, 27, 17, 19, 17,
  33, 23, 30, 25, 26,
  33, 20, 50, 37, 40
)

# Convert to rates
df <- df %>%
  mutate(
    FL = 1 / FL_days,
    NC = 1 / NC_days,
    CA = 1 / CA_days
  )

# ===============================
# 2. Fit models for each strain
# ===============================
fit_model <- function(yvar) {
  f <- as.formula(paste(yvar, "~ Temp + RH"))
  m <- lm(f, data = df)
  preds <- predict(m)
  RMSE <- sqrt(mean((df[[yvar]] - preds)^2))
  list(model = m, R2 = summary(m)$r.squared, RMSE = RMSE)
}

fit_FL <- fit_model("FL")
fit_NC <- fit_model("NC")
fit_CA <- fit_model("CA")

# ===============================
# 3. Print results
# ===============================
cat("\n=== FL STRAIN ===\n")
print(summary(fit_FL$model))
cat("R² =", fit_FL$R2, " RMSE =", fit_FL$RMSE, "\n")

cat("\n=== NC STRAIN ===\n")
print(summary(fit_NC$model))
cat("R² =", fit_NC$R2, " RMSE =", fit_NC$RMSE, "\n")

cat("\n=== CA STRAIN ===\n")
print(summary(fit_CA$model))
cat("R² =", fit_CA$R2, " RMSE =", fit_CA$RMSE, "\n")
