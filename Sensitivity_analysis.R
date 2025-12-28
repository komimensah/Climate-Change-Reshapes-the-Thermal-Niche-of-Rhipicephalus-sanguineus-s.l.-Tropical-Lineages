library(lhs)
library(sensitivity)
library(ggplot2)

############################################################
## 1. DEVELOPMENT & MORTALITY FUNCTIONS (unchanged)
############################################################

egg_rate_FL <- function(T, RH){
  pmax(0, -0.05866 + 0.003649 * T - 0.00002482 * RH)
}

larva_rate_FL <- function(T, RH){
  pmax(0, -0.09847 + 0.007270 * T - 0.0001245 * RH)
}

nymph_rate_FL <- function(T, RH){
  pmax(0, -0.1434 + 0.007983 * T + 0.000000985 * RH)
}

## Adult mortality now EXPOSES β_A as a parameter
muA_FL <- function(T, RH, betaA){
  pmax(0, -0.002067 + betaA * T - 0.00007199 * RH)
}

muL_FL <- function(T, RH){
  pmax(0, 0.35688 + 0.00408 * T - 0.00499 * RH)
}

muN_FL <- function(T, RH){
  pmax(0, -0.01980 + 0.004269 * T - 0.001031 * RH)
}

############################################################
## 2. OVIPOSITION SUCCESS (unchanged)
############################################################

oviposition_FL <- function(T, RH){
  base <- ifelse(T < 20, 0,
                 ifelse(T == 20, 0.67,
                        ifelse(T > 20 & T <= 35, 1, 0)))
  base
}

############################################################
## 3. HOST CONTACT FUNCTION — β NOW A PARAMETER
############################################################

R_tick <- function(H, beta,
                   ratio_free = 828,
                   ratio_owned = 3.3,
                   roaming_fraction = 0.8){
  
  dogs_free  <- H / ratio_free
  dogs_owned <- H / ratio_owned
  dogs_roaming <- dogs_free + roaming_fraction * dogs_owned
  
  R <- 1 - exp(-beta * dogs_roaming)
  return(pmax(0, pmin(R, 1)))
}

############################################################
## 4. UPDATED R0 MODEL INCLUDING sampled: β, F0, Tp, βA
############################################################

R0_model <- function(T, RH, H, beta, F0, Tp, betaA,
                     muE = 0.4605, p_female = 0.5){
  
  # Development
  sigma1 <- egg_rate_FL(T, RH)
  sigma2 <- larva_rate_FL(T, RH)
  sigma3 <- nymph_rate_FL(T, RH)
  
  # Mortality (adult uses betaA as slope)
  muL <- muL_FL(T, RH)
  muN <- muN_FL(T, RH)
  muA <- muA_FL(T, RH, betaA)
  muA[muA < 0.01] <- 0.01
  
  # Contact rate
  Rhost <- R_tick(H, beta)
  
  # Oviposition
  Ovi <- oviposition_FL(T, RH)
  
  # Fecundity (now F0 varies)
  F_total <- F0 * Ovi * Rhost * p_female
  
  F <- F_total * muA
  
  # R0 formulation
  R0 <- exp(-muA * Tp) *
    (F / muA) *
    (sigma1 / (sigma1 + muE)) *
    (sigma2 / (sigma2 + muL)) *
    (sigma3 / (sigma3 + muN))
  
  R0[is.nan(R0) | is.infinite(R0) | R0 < 0] <- 0
  return(R0)
}

############################################################
## 5. LHS SAMPLING OF THE 7 PARAMETERS
############################################################

set.seed(2025)
n_sims <- 1000

param_names <- c("T","RH","H","beta","F0","Tp","betaA")

param_min <- c(T=10, RH=20, H=10, beta=0.005, F0=500, Tp=8, betaA=0.0001)
param_max <- c(T=35, RH=90, H=1000, beta=0.05,  F0=1100, Tp=20, betaA=0.002)

X_u <- randomLHS(n_sims, length(param_names))
colnames(X_u) <- param_names
X_df <- as.data.frame(X_u)

# Scale to parameter ranges
for (v in param_names) {
  X_df[[v]] <- param_min[v] + X_df[[v]] * (param_max[v] - param_min[v])
}

############################################################
## 6. EVALUATE THE MODEL
############################################################

Y <- numeric(n_sims)
for(i in seq_len(n_sims)){
  Y[i] <- R0_model(
    T = X_df$T[i],
    RH = X_df$RH[i],
    H = X_df$H[i],
    beta = X_df$beta[i],
    F0 = X_df$F0[i],
    Tp = X_df$Tp[i],
    betaA = X_df$betaA[i]
  )
}

summary(Y)

############################################################
## 7. PRCC
############################################################

prcc_res <- pcc(X = X_df, y = Y, rank = TRUE, nboot = 200, conf = 0.95)
print(prcc_res$PRCC)

############################################################
## 8. HIGH-QUALITY PRCC PLOT
############################################################

prcc_df <- prcc_res$PRCC

df_plot <- data.frame(
  Parameter = rownames(prcc_df),
  PRCC = prcc_df$original
)

df_plot <- df_plot[order(abs(df_plot$PRCC), decreasing = TRUE), ]

library(ggplot2)

# Extract PRCC + 95% CI
df_plot <- data.frame(
  Parameter = rownames(prcc_res$PRCC),
  PRCC      = prcc_res$PRCC$original,
  lowerCI   = prcc_res$PRCC$`min. c.i.`,
  upperCI   = prcc_res$PRCC$`max. c.i.`
)

# Order by effect size
df_plot <- df_plot[order(abs(df_plot$PRCC), decreasing = TRUE), ]

# Extract PRCC table safely
prcc_tab <- as.data.frame(prcc_res$PRCC)

# Identify possible CI columns
ci_cols <- grep("c.i", colnames(prcc_tab), value = TRUE)

# If CI exists, extract it; otherwise use NA
if (length(ci_cols) >= 2) {
  lowerCI <- prcc_tab[[ci_cols[1]]]
  upperCI <- prcc_tab[[ci_cols[2]]]
} else {
  lowerCI <- upperCI <- rep(NA, nrow(prcc_tab))
}

# Build plotting dataframe
df_plot <- data.frame(
  Parameter = rownames(prcc_tab),
  PRCC      = prcc_tab$original,
  lowerCI   = lowerCI,
  upperCI   = upperCI
)

# Order by absolute PRCC magnitude
df_plot <- df_plot[order(abs(df_plot$PRCC), decreasing = TRUE), ]
library(ggplot2)

p <- ggplot(df_plot, aes(x = reorder(Parameter, PRCC), y = PRCC)) +
  
  # Main bars
  geom_col(fill = "#2b8cbe", width = 0.65) +
  
  # CI bars only if they exist
  geom_errorbar(
    aes(ymin = lowerCI, ymax = upperCI),
    width = 0.15, color = "black", linewidth = 0.6,
    na.rm = TRUE
  ) +
  
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.7) +
  
  coord_flip() +
  
  theme_bw(base_size = 15) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.title = element_text(face = "bold"),
    axis.text  = element_text(size = 12),
    plot.title = element_text(face = "bold", size = 18)
  ) +
  
  labs(
    title = "",
    subtitle = "",
    x = "Parameters",
    y = "Partial Rank Correlation Coefficient (PRCC)"
  )

print(p)


ggsave("PRCC_R0_FL.pdf", p, width = 8, height = 6, units = "in")
ggsave("PRCC_R0_FL.png", p, width = 8, height = 6, units = "in", dpi = 300)
#######################
