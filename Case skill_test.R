library(terra)
library(readr)
library(dplyr)

# ------------------------------------------------------------
# 0. Clear workspace
# ------------------------------------------------------------
rm(list = ls())
cat("\014")

# ------------------------------------------------------------
# 1. Define development / mortality / oviposition functions
# ------------------------------------------------------------
egg_rate_FL <- function(T, RH){
  pmax(0, -0.05866 + 0.003649 * T - 0.00002482 * RH)
}
larva_rate_FL <- function(T, RH){
  pmax(0, -0.09847 + 0.007270 * T - 0.0001245 * RH)
}
nymph_rate_FL <- function(T, RH){
  pmax(0, -0.1434 + 0.007983 * T + 0.000000985 * RH)
}
muL_FL <- function(T, RH){
  pmax(0, 0.35688 + 0.00408 * T - 0.00499 * RH)
}
muN_FL <- function(T, RH){
  pmax(0, -0.01980 + 0.004269 * T - 0.001031 * RH)
}
muA_FL <- function(T, RH){
  pmax(0, -0.002067 + 0.0004848 * T - 0.00007199 * RH)
}
oviposition_FL <- function(T, RH = NULL){
  ifelse(
    T >= 15 & T <= 40,
    1,
    ifelse((T >= 10 & T < 15) | (T > 40 & T <= 45),
           0.5,
           0
    )
  )
}

R_tick <- function(H,
                   ratio_free = 828,
                   ratio_owned = 3.3,
                   roaming_fraction = 0.8,
                   beta = 0.02){
  dogs_free  <- H / ratio_free
  dogs_owned <- H / ratio_owned
  dogs_roaming <- dogs_free + roaming_fraction * dogs_owned
  R <- 1 - exp(-beta * dogs_roaming)
  pmax(0, pmin(R, 1))
}

R0_strain <- function(T, RH, H,
                      dev_fun_egg, dev_fun_larva, dev_fun_nymph,
                      mort_fun_L, mort_fun_N, mort_fun_A,
                      ovip_fun,
                      Tp = 12.9,
                      F0 = 805.4,
                      p_female = 0.5){
  sigma1 <- dev_fun_egg(T, RH)
  sigma2 <- dev_fun_larva(T, RH)
  sigma3 <- dev_fun_nymph(T, RH)
  muL <- mort_fun_L(T, RH)
  muN <- mort_fun_N(T, RH)
  muA <- mort_fun_A(T, RH)
  muA[muA < 0.02] <- 0.02
  Rhost <- R_tick(H)
  hatch <- rep(1, length(T))
  muE_eff <- 0.07
  F_total <- F0 * ovip_fun(T, RH) * hatch * Rhost * p_female
  F_daily <- F_total * muA
  R0 <- exp(-muA * Tp) *
    (F_daily / muA) *
    (sigma1 / (sigma1 + muE_eff)) *
    (sigma2 / (sigma2 + muL)) *
    (sigma3 / (sigma3 + muN))
  R0[is.nan(R0) | R0 < 0 | is.infinite(R0)] <- 0
  return(R0)
}

R0_FL_fun <- function(T, RH, H){
  R0_strain(T = T, RH = RH, H = H,
            dev_fun_egg   = egg_rate_FL,
            dev_fun_larva = larva_rate_FL,
            dev_fun_nymph = nymph_rate_FL,
            mort_fun_L    = muL_FL,
            mort_fun_N    = muN_FL,
            mort_fun_A    = muA_FL,
            ovip_fun      = oviposition_FL)
}

# ------------------------------------------------------------
# 2. Read occurrence / case CSV
# ------------------------------------------------------------
df_pts <- read_csv("/Users/kagboka/Desktop/Anaplasmosis_case_Africa.csv", na = c("", "NA"))

# ------------------------------------------------------------
# 3. Load and prepare environmental rasters (temperature, humidity, human density → dog density)
# ------------------------------------------------------------
temp_folder <- '/Users/kagboka/Desktop/Ritter_work/LTM_tif/factual/tas/'
temp_stack  <- rast(list.files(temp_folder, pattern = "\\.tif$", full.names = TRUE))
r_T <- mean(temp_stack, na.rm = TRUE)

rh_folder <- '/Users/kagboka/Desktop/Ritter_work/LTM_tif/factual/hurs/'
r_RH <- mean(rast(list.files(rh_folder, pattern = "\\.tif$", full.names = TRUE)), na.rm = TRUE)

r_H <- rast('/Users/kagboka/Desktop/Lesh2.0/gpw_v4_population_density_adjusted_to_2015_unwpp_country_totals_rev11_2000_2pt5_min.tif')

r_RH <- project(r_RH, r_T, method = "bilinear")
r_H  <- project(r_H,  r_T, method = "bilinear")

ext_common <- intersect(ext(r_T), intersect(ext(r_RH), ext(r_H)))
r_T   <- crop(r_T,  ext_common)
r_RH  <- crop(r_RH, ext_common)
r_H   <- crop(r_H,  ext_common)

dog_density_fun <- function(h, ratio_free = 828, ratio_owned = 3.3, roaming_fraction = 0.8){
  dogs_free  <- h / ratio_free
  dogs_owned <- h / ratio_owned
  dogs_free + roaming_fraction * dogs_owned
}

r_H_dogs <- app(r_H, fun = dog_density_fun)

env_stack <- c(r_T, r_RH, r_H_dogs)
names(env_stack) <- c("mean_temp", "mean_RH", "H")

# ------------------------------------------------------------
# 4. Extract env values at occurrence coordinates
# ------------------------------------------------------------
pts_vect <- vect(df_pts,
                 geom = c("Sampling_longitude", "Sampling_latitude"),
                 crs = crs(env_stack))

extracted <- terra::extract(env_stack, pts_vect)
df_out <- bind_cols(df_pts, extracted)

# ------------------------------------------------------------
# 5. Compute predicted R0 and correlations
# ------------------------------------------------------------
df_out$R0_pred <- R0_FL_fun(
  T  = df_out$mean_temp,
  RH = df_out$mean_RH,
  H  = df_out$H
)

df_valid <- df_out[ !is.na(df_out$R0_pred) & !is.na(df_out$N_positive), ]

cat("Number of valid observations:", nrow(df_valid), "\n")

pearson_res  <- cor.test(df_valid$R0_pred, df_valid$N_positive, method = "pearson")
spearman_res <- cor.test(df_valid$R0_pred, df_valid$N_positive, method = "spearman")

cat("=== Correlation results ===\n")
cat("Pearson r:", pearson_res$estimate, " (p =", pearson_res$p.value, ")\n")
cat("Spearman rho:", spearman_res$estimate, " (p =", spearman_res$p.value, ")\n")
# --- After your correlation block ---

## Linear regression (least-squares) — simple baseline
lm_fit <- lm(N_positive ~ R0_pred, data = df_valid)
summary(lm_fit)

cat("=== Linear regression (OLS) ===\n")
cat("R-squared:", summary(lm_fit)$r.squared, "\n")
cat("Coefficient (R0_pred):", coef(lm_fit)["R0_pred"], " (p =", summary(lm_fit)$coefficients["R0_pred","Pr(>|t|)"], ")\n\n")

## Poisson (or Negative-Binomial) regression — for count data
# If you expect overdispersion, you might switch to negative binomial (e.g. using MASS::glm.nb)
glm_pois <- glm(N_positive ~ R0_pred,
                data = df_valid,
                family = poisson(link = "log"))

cat("=== Poisson regression ===\n")
print(summary(glm_pois))

# Check over-dispersion
disp_ratio <- sum(residuals(glm_pois, type = "pearson")^2) / df.residual(glm_pois)
cat("Dispersion ratio (pearson):", disp_ratio, "\n\n")

if(disp_ratio > 1.5) {
  # Over-dispersion detected — try negative binomial
  if(requireNamespace("MASS", quietly = TRUE)) {
    nb_fit <- MASS::glm.nb(N_positive ~ R0_pred, data = df_valid)
    cat("=== Negative-Binomial regression ===\n")
    print(summary(nb_fit))
  } else {
    warning("MASS package needed for negative binomial regression; please install.")
  }
}

## Generalized Additive Model (GAM) — allows non-linear relationship
if(requireNamespace("mgcv", quietly = TRUE)) {
  library(mgcv)
  gam_fit <- gam(N_positive ~ s(R0_pred), data = df_valid, family = poisson(link = "log"))
  cat("=== GAM (Poisson) ===\n")
  print(summary(gam_fit))
  
  # Optionally plot the smooth effect
  plot(gam_fit, shade = TRUE, seWithMean = TRUE,
       main = "Smooth effect of R0_pred on N_positive")
} else {
  warning("mgcv package not installed — cannot fit GAM.")
}
###################################
# install required packages if not already
if (!requireNamespace("MASS", quietly = TRUE)) install.packages("MASS")
if (!requireNamespace("mgcv", quietly = TRUE)) install.packages("mgcv")
if (!requireNamespace("pscl", quietly = TRUE)) install.packages("pscl")
if (!requireNamespace("DHARMa", quietly = TRUE)) install.packages("DHARMa")

library(MASS)
library(mgcv)
library(pscl)
library(DHARMa)

# Subset to valid data
df_valid <- df_out %>%
  filter(!is.na(R0_pred), !is.na(N_positive))

cat("Number of valid observations:", nrow(df_valid), "\n")

# ----- 1) Negative-Binomial GLM -----
nb_glm <- glm.nb(N_positive ~ R0_pred, data = df_valid)
summary(nb_glm)

# Check dispersion / goodness-of-fit
res <- simulateResiduals(nb_glm)
plot(res)  # diagnostic plots
cat("Dispersion (Pearson):", sum(residuals(nb_glm, type = "pearson")^2) / nb_glm$df.residual, "\n")

# ----- 2) Zero-Inflated negative binomial (if many zeros) -----
zinb <- zeroinfl(N_positive ~ R0_pred | 1,
                 data = df_valid,
                 dist = "negbin")
summary(zinb)

# Compare AIC / BIC of NB vs ZINB
cat("AIC NB:", AIC(nb_glm), "\n")
cat("AIC ZINB:", AIC(zinb), "\n")

# ----- 3) GAM with Negative-Binomial family (allows non-linear effect) -----
gam_nb <- gam(N_positive ~ s(R0_pred),
              data = df_valid,
              family = nb(link = "log"),
              method = "REML")
summary(gam_nb)
plot(gam_nb, shade = TRUE, main = "Smooth effect of predicted R0 on case counts")

# Diagnostics using DHARMa
res_gam <- simulateResiduals(gam_nb)
plot(res_gam)

# ----- 4) (Optional) Cross-validation (k-fold) to assess predictive stability -----
set.seed(123)
k <- 5
folds <- sample(rep(1:k, length.out = nrow(df_valid)))
cv_res <- lapply(1:k, function(i) {
  train <- df_valid[folds != i, ]
  test  <- df_valid[folds == i, ]
  
  m <- glm.nb(N_positive ~ R0_pred, data = train)
  preds <- predict(m, newdata = test, type = "response")
  cor_obs_pred <- cor(test$N_positive, preds, method = "spearman")
  data.frame(fold = i,
             spearman = cor_obs_pred,
             rmse = sqrt(mean((test$N_positive - preds)^2)))
})
cv_df <- do.call(rbind, cv_res)
print(cv_df)
cat("Mean CV Spearman:", mean(cv_df$spearman), "\n")
cat("Mean CV RMSE:", mean(cv_df$rmse), "\n")

# ----- 5) (Optional) Add covariates if available and refit -----
# e.g. suppose you have sampling_year or country or population_density or other variables in df_out
#  nb_glm2 <- glm.nb(N_positive ~ R0_pred + pop_density + factor(country), data = df_valid)

# Then inspect coefficients, diagnostics, etc.
#######################################################
# -----------------------------------------
# LOG-TRANSFORMED CASES
# -----------------------------------------

# Add log(N_positive + 1) to avoid log(0)
df_valid$log_cases <- log(df_valid$N_positive + 1)

cat("=== Log-transformed cases summary ===\n")
summary(df_valid$log_cases)

# Correlations with log cases
pearson_log <- cor.test(df_valid$R0_pred, df_valid$log_cases, method = "pearson")
spearman_log <- cor.test(df_valid$R0_pred, df_valid$log_cases, method = "spearman")

cat("=== Correlation with log(N_positive + 1) ===\n")
cat("Pearson r (log):", pearson_log$estimate, " (p =", pearson_log$p.value, ")\n")
cat("Spearman rho (log):", spearman_log$estimate, " (p =", spearman_log$p.value, ")\n")

# -----------------------------------------
# Linear regression on log cases
# -----------------------------------------
lm_log <- lm(log_cases ~ R0_pred, data = df_valid)
cat("=== Linear regression (log cases) ===\n")
print(summary(lm_log))

cat("R-squared (log model):", summary(lm_log)$r.squared, "\n")
cat("Coefficient (R0_pred):", coef(lm_log)["R0_pred"],
    " (p =", summary(lm_log)$coefficients["R0_pred", "Pr(>|t|)"], ")\n\n")

# -----------------------------------------
# GAM on log cases
# -----------------------------------------
if(requireNamespace("mgcv", quietly = TRUE)) {
  library(mgcv)
  gam_log <- gam(log_cases ~ s(R0_pred),
                 data = df_valid,
                 method = "REML")
  
  cat("=== GAM (log cases) ===\n")
  print(summary(gam_log))
  
  # Plot effect
  plot(gam_log, shade = TRUE,
       main = "Smooth effect of R0_pred on log(N_positive + 1)")
  
} else {
  warning("mgcv package not installed; cannot fit GAM on log cases.")
}

# -----------------------------------------
# Diagnostics for transformed models
# -----------------------------------------
if(requireNamespace("DHARMa", quietly = TRUE)) {
  library(DHARMa)
  
  cat("=== Linear model residual diagnostics ===\n")
  sim_lm <- simulateResiduals(lm_log)
  plot(sim_lm)
  
  if(exists("gam_log")) {
    cat("=== GAM (log cases) residual diagnostics ===\n")
    sim_gam_log <- simulateResiduals(gam_log)
    plot(sim_gam_log)
  }
}
