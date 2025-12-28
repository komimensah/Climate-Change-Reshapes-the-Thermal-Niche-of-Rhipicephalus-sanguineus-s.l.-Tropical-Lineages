############################################################
## 5. STABILITY ANALYSIS: 4-STAGE ODE + EIGENVALUE λ_max
############################################################

# This uses your existing functions:
# egg_rate_FL(), larva_rate_FL(), nymph_rate_FL()
# muL_FL(), muN_FL(), muA_FL()
# oviposition_FL(), R_tick()

# Biological constants (same as in your R0 formulation)
Tp_default    <- 12.9    # pre-infectious period / delay [days]
muE_default   <- 0.4605  # non-development egg loss
F0_default    <- 805.4   # baseline fecundity
p_female_def  <- 0.5     # female proportion

# 5.1 Function that returns λ_max for given T, RH, H
lambda_max_FL <- function(T, RH, H,
                          Tp = Tp_default,
                          muE = muE_default,
                          F0 = F0_default,
                          p_female = p_female_def) {
  
  # --- Development rates (per day) ---
  sigma1 <- egg_rate_FL(T, RH)    # E -> L
  sigma2 <- larva_rate_FL(T, RH)  # L -> N
  sigma3 <- nymph_rate_FL(T, RH)  # N -> A
  
  # --- Mortality rates (per day) ---
  muL <- muL_FL(T, RH)
  muN <- muN_FL(T, RH)
  muA <- muA_FL(T, RH)
  muA[muA < 0.01] <- 0.01  # protect against 0
  
  # --- Host contact / bite success (0–1) ---
  Rhost <- R_tick(H)
  
  # --- Oviposition success (0–1) ---
  Ovi <- oviposition_FL(T, RH)
  
  # --- Total eggs per female over lifetime ---
  F_total <- F0 * Ovi * Rhost * p_female
  
  # --- Egg production rate per adult per day ---
  # (eggs per female lifetime) * (1 / life expectancy)
  F_perA <- F_total * muA
  
  # ---- Build linear ODE system: dX/dt = J %*% X ----
  # States: X = (E, L, N, A)^T
  # dE/dt = F_perA * A        - (muE + sigma1) * E
  # dL/dt = sigma1 * E        - (muL + sigma2) * L
  # dN/dt = sigma2 * L        - (muN + sigma3) * N
  # dA/dt = sigma3 * N        - muA           * A
  
  J <- matrix(0, nrow = 4, ncol = 4)
  rownames(J) <- colnames(J) <- c("E", "L", "N", "A")
  
  J["E", "E"] <- -(muE + sigma1)
  J["E", "A"] <-  F_perA
  
  J["L", "E"] <-  sigma1
  J["L", "L"] <- -(muL + sigma2)
  
  J["N", "L"] <-  sigma2
  J["N", "N"] <- -(muN + sigma3)
  
  J["A", "N"] <-  sigma3
  J["A", "A"] <- -muA
  
  # --- Eigenvalues of the Jacobian ---
  ev <- eigen(J, only.values = TRUE)$values
  # dominant real part
  lambda_max <- max(Re(ev), na.rm = TRUE)
  
  return(lambda_max)
}

############################################################
## 6. QUICK CHECK AT A SINGLE CLIMATE POINT
############################################################

T_test  <- 28   # °C
RH_test <- 75   # %
H_test  <- 500  # dogs/km²

lambda_max_FL(T_test, RH_test, H_test)
# >  λ_max > 0  -> population tends to grow (persistence)
# >  λ_max < 0  -> population tends to decline (extinction)

############################################################
## 7. STABILITY SURFACE OVER T–RH FOR FIXED DOG DENSITY
############################################################

library(dplyr)
library(tidyr)
library(ggplot2)

T_vec  <- seq(10, 40, by = 1)
RH_vec <- seq(20, 95, by = 5)
H_fixed <- 300   # choose a representative dog density

grid <- expand.grid(T = T_vec, RH = RH_vec)

grid$lambda_max <- mapply(
  lambda_max_FL,
  T  = grid$T,
  RH = grid$RH,
  H  = H_fixed
)

# Optional: classify stable vs unstable around 0
grid <- grid %>%
  mutate(
    stability = ifelse(lambda_max > 0, "Persistence (λmax > 0)", 
                       "Decline (λmax ≤ 0)")
  )

# Nice stability plot (heat + contour + boundary line)
p_stab <- ggplot(grid, aes(x = T, y = RH, fill = lambda_max)) +
  geom_tile() +
  geom_contour(aes(z = lambda_max),
               breaks = 0, colour = "black", linewidth = 0.7) +
  scale_fill_viridis_c(
    option = "plasma",
    name   = expression(lambda[max]),
    limits = range(grid$lambda_max, na.rm = TRUE)
  ) +
  labs(
    x = "Temperature (°C)",
    y = "Relative humidity (%)",
    title = expression("Climate-dependent stability of " * R[0]),
    subtitle = bquote("FL strain; dog density H = " * .(H_fixed) * " dogs/km"^2)
  ) +
  
  
  theme_bw(base_size = 13) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title       = element_text(face = "bold"),
    plot.subtitle    = element_text(size = 11)
  )

print(p_stab)

p_class <- ggplot(grid, aes(x = T, y = RH, fill = stability)) +
  geom_tile() +
  scale_fill_manual(values = c(
    "Persistence (λmax > 0)" = "#1a9850",
    "Decline (λmax ≤ 0)"     = "#d73027"
  )) +
  labs(
    x = "Temperature (°C)",
    y = "Relative humidity (%)",
    fill = "Local dynamics",
    title = "Tick population persistence boundary",
    subtitle = bquote("FL strain; dog density H = " * .(H_fixed) * " dogs/km"^2)
  ) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title       = element_text(face = "bold"),
    plot.subtitle    = element_text(size = 11)
  )
print(p_class)
ggsave("p_class.png", p_class, width = 8, height = 6, units = "in", dpi = 300)
###############################
library(dplyr)
library(tidyr)
library(ggplot2)

set.seed(2025)

# Number of Monte Carlo draws per grid cell
n_draws <- 200  # you can increase to 500, 1000, depending on computing time

# Define distributions/ranges for H (dog density) — adjust to plausible range
H_min <- 50      # e.g. 50 dogs/km²
H_max <- 1000    # e.g. 1000 dogs/km²

T_vec  <- seq(10, 40, by = 2)
RH_vec <- seq(20, 95, by = 5)

grid0 <- expand.grid(T = T_vec, RH = RH_vec)

# For each grid cell, draw n_draws dog-density values, compute lambda_max
results <- lapply(seq_len(nrow(grid0)), function(i) {
  T0  <- grid0$T[i]
  RH0 <- grid0$RH[i]
  
  H_samples <- runif(n_draws, min = H_min, max = H_max)
  
  lambdas <- sapply(H_samples, function(H0) {
    lambda_max_FL(T = T0, RH = RH0, H = H0)
  })
  
  data.frame(
    T = T0,
    RH = RH0,
    lambda_median = median(lambdas, na.rm = TRUE),
    lambda_p05 = quantile(lambdas, probs = 0.05, na.rm = TRUE),
    lambda_p95 = quantile(lambdas, probs = 0.95, na.rm = TRUE)
  )
})

df_uncert <- do.call(rbind, results)

# Example plot: median + uncertainty ribbon contours
p_med <- ggplot(df_uncert, aes(x = T, y = RH, fill = lambda_median)) +
  geom_tile() +
  scale_fill_viridis_c(option = "plasma") +
  labs(
    x = "Temperature (°C)",
    y = "Relative humidity (%)",
    fill = expression(lambda[max]),
    title = "Median λ[max] over dog-density uncertainty",
    subtitle = paste0(n_draws, " draws; H ~ Uniform(", H_min,",",H_max,") dogs/km²")
  ) +
  theme_bw()

print(p_med)
ggsave("p_med.png", p_med, width = 8, height = 6, units = "in", dpi = 300)
p_bound <- ggplot(df_uncert, aes(x = T, y = RH)) +
  geom_tile(aes(fill = lambda_p05)) +
  scale_fill_viridis_c(option = "magma") +
  labs(
    x = "Temperature (°C)",
    y = "Relative humidity (%)",
    fill = "λmax (5th percentile)",
    title = "Lower 5% quantile of λmax under dog-density uncertainty"
  ) +
  theme_bw()
print(p_bound)
ggsave("p_bound.png", p_bound, width = 8, height = 6, units = "in", dpi = 300)
