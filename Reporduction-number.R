
############################################################
## 1. Clear workspace
############################################################
rm(list = ls())
cat("\014")

############################################################
## 2. DEVELOPMENT & MORTALITY FUNCTIONS (unchanged)
############################################################

egg_rate_FL <- function(T, RH){
  pmax(0, -0.05866 + 0.003649 * T - 0.00002482 * RH)
}
egg_rate_NC <- function(T, RH){
  pmax(0, -0.05389 + 0.003567 * T - 0.00000696 * RH)
}
egg_rate_CA <- function(T, RH){
  pmax(0, -0.04524 + 0.003192 * T + 0.00001260 * RH)
}

larva_rate_FL <- function(T, RH){
  pmax(0, -0.09847 + 0.007270 * T - 0.0001245 * RH)
}
larva_rate_NC <- function(T, RH){
  pmax(0, -0.05767 + 0.005947 * T - 0.0001028 * RH)
}
larva_rate_CA <- function(T, RH){
  pmax(0, -0.08265 + 0.006921 * T - 0.0001880 * RH)
}

nymph_rate_FL <- function(T, RH){
  pmax(0, -0.1434 + 0.007983 * T + 0.000000985 * RH)
}
nymph_rate_NC <- function(T, RH){
  pmax(0, -0.05969 + 0.004203 * T + 0.00003553 * RH)
}
nymph_rate_CA <- function(T, RH){
  pmax(0, -0.07980 + 0.005159 * T - 0.000004675 * RH)
}

muL_FL <- function(T, RH){
  pmax(0, 0.35688 + 0.00408 * T - 0.00499 * RH)
}
muL_NC <- function(T, RH){
  pmax(0, 0.17432 + 0.00635 * T - 0.00378 * RH)
}
muL_CA <- function(T, RH){
  pmax(0, 0.33923 + 0.00444 * T - 0.00478 * RH)
}

muN_FL <- function(T, RH){
  pmax(0, -0.01980 + 0.004269 * T - 0.001031 * RH)
}
muN_NC <- function(T, RH){
  pmax(0, -0.05026 + 0.005165 * T - 0.0007509 * RH)
}
muN_CA <- function(T, RH){
  pmax(0, -0.05026 + 0.005165 * T - 0.0007509 * RH)
}

muA_FL <- function(T, RH){
  pmax(0, -0.002067 + 0.0004848 * T - 0.00007199 * RH)
}
muA_NC <- function(T, RH){
  pmax(0, 0.004948 + 0.0002812 * T - 0.0001069 * RH)
}
muA_CA <- function(T, RH){
  pmax(0, -0.01009 + 0.001084 * T - 0.0001152 * RH)
}

############################################################
## 3. OVIPOSITION SUCCESS (unchanged)
############################################################

oviposition_FL <- function(T, RH = NULL) {
  ifelse(
    T >= 15 & T <= 40, 
    1,
    ifelse(
      (T >= 10 & T < 15) | (T > 40 & T <= 45),
      0.5,
      0
    )
  )
}
oviposition_NC <- oviposition_FL
oviposition_CA <- oviposition_FL

############################################################
## 4. HOST CONTACT RATE (unchanged)
############################################################

R_tick <- function(H,
                   ratio_free = 828,
                   ratio_owned = 3.3,
                   roaming_fraction = 0.8,
                   beta = 0.02) {
  dogs_free  <- H / ratio_free
  dogs_owned <- H / ratio_owned
  dogs_roaming <- dogs_free + roaming_fraction * dogs_owned
  R <- 1 - exp(-beta * dogs_roaming)
  return(pmax(0, pmin(R, 1)))
}

############################################################
## 5. Egg-hatch probability (vectorized)
############################################################


############################################################
## 6. R0 formula (no w_lc)
############################################################

R0_strain <- function(T, RH, H,
                      dev_fun_egg, dev_fun_larva, dev_fun_nymph,
                      mort_fun_L, mort_fun_N, mort_fun_A,
                      ovip_fun,
                      Tp     = 12.9,
                      F0     = 805.4,
                      p_female = 0.5) {
  
  sigma1 <- dev_fun_egg(T, RH)
  sigma2 <- dev_fun_larva(T, RH)
  sigma3 <- dev_fun_nymph(T, RH)
  
  muL <- mort_fun_L(T, RH)
  muN <- mort_fun_N(T, RH)
  muA <- mort_fun_A(T, RH)
  muA[muA < 0.02] <- 0.02
  
  Rhost <- R_tick(H)
  
  # Egg stage
  # hatch <- egg_hatch_prob(T)   # REMOVE THIS (causes raster crash)
  hatch <- rep(1, length(T))     # safe dummy
  muE_eff <- 0.07                # constant egg mortality
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

############################################################
## 7. Build Grid (T, RH, H → dog density)
############################################################

T_vec  <- seq(-50, 40, by = 0.5)
RH_vec <- seq(17, 92, by = 10)
H_humans <- c(100, 300, 600, 1000, 7000)

dog_density <- function(H,
                        ratio_free = 828,
                        ratio_owned = 3.3,
                        roaming_fraction = 0.8) {
  dogs_free  <- H / ratio_free
  dogs_owned <- H / ratio_owned
  dogs_roaming <- dogs_free + roaming_fraction * dogs_owned
  return(dogs_roaming)
}

H_vec <- round(dog_density(H_humans))

grid <- expand.grid(
  T  = T_vec,
  RH = RH_vec,
  H  = H_vec
)

############################################################
## 8. Compute R0 for each strain (FL, NC, CA)
############################################################

grid$R0_FL <- R0_strain(grid$T, grid$RH, grid$H,
                        egg_rate_FL, larva_rate_FL, nymph_rate_FL,
                        muL_FL, muN_FL, muA_FL,
                        oviposition_FL)
#grid$R0_NC <- R0_strain(grid$T, grid$RH, grid$H,
#   egg_rate_NC, larva_rate_NC, nymph_rate_NC,
#   muL_NC, muN_NC, muA_NC,
#  oviposition_NC)
grid$R0_CA <- R0_strain(grid$T, grid$RH, grid$H,
                        egg_rate_CA, larva_rate_CA, nymph_rate_CA,
                        muL_CA, muN_CA, muA_CA,
                        oviposition_CA)

############################################################
## 9. Format & Plot (long + ggplot)
############################################################

library(dplyr)
library(tidyr)
library(ggplot2)

grid_long <- grid %>%
  pivot_longer(cols = starts_with("R0_"),
               names_to = "Strain",
               values_to = "R0") %>%
  mutate(Strain = recode(Strain,
                         R0_FL = "FL",
                         R0_NC = "NC",
                         R0_CA = "CA")) %>%
  
  # Create R0 classes (no 0 class)
  mutate(R0_class = cut(
    R0,
    breaks = c(-Inf, 1, 5, 10, 20, 40, Inf),
    labels = c(
      "<1",
      "1–5",
      "5–10",
      "10–20",
      "20–40",
      ">40"
    ),
    right = FALSE
  ))

# Colors for each class


# Colour palette for discrete VC classes (red-based, high risk = dark)
cols <- c(
  "#fff5eb",  # very light peach (lowest class, <1)
  "#fdd0a2",  # light orange
  "#fdae6b",  # mid orange
  "#fd8d3c",  # strong orange
  "#e6550d",  # red-orange
  "#a63603"   # dark red (highest risk, >40)
)
# Plot
p_heat <- ggplot(grid_long, aes(x = T, y = RH, fill = R0_class)) +
  geom_tile() +
  scale_fill_manual(values = cols, name = "Ro class") +
  facet_grid(Strain ~ H,
             labeller = labeller(
               H = function(x) paste0(x, " dogs/km²")
             )) +
  labs(
    x = "Temperature (°C)",
    y = "Relative Humidity (%)",
    title = expression()
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "grey90", colour = NA),
    legend.position = "right"
  )

print(p_heat)
# Save as high-resolution PNG
ggsave(
  filename = "R0_heatmap.png",
  plot = p_heat,
  width = 14,       # in inches
  height = 10,
  dpi = 300,        # high-resolution
  bg = "white"
)

# Save as high-resolution TIFF (journal-friendly)
ggsave(
  filename = "R0_heatmap.tiff",
  plot = p_heat,
  width = 14,
  height = 10,
  dpi = 300,
  compression = "lzw",
  bg = "white"
)

# Save as vector PDF (infinite resolution)
ggsave(
  filename = "R0_heatmap.pdf",
  plot = p_heat,
  width = 14,
  height = 10,
  bg = "white"
)
