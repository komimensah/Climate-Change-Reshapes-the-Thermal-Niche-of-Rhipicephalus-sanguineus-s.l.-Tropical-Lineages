############################################################
## Host–Tick Encounter Rate (Hofmeester Gas Model Adapted)
## R_tick = D_total * v * (2 * r) / pi
##
## Where:
##   H = human density (humans per km^2)
##   D_free = H / 828         (free-roaming dogs)
##   D_dom  = H / 3.3         (owned dogs)
##   D_total = D_free + D_dom (dogs per km^2)
##   v = dog daily movement (km/day) [default = 8]
##   r = 0.5 m = 0.0005 km (host–tick interaction radius)
############################################################

R_tick <- function(H,
                   v = 8,             # dog movement speed (km/day)
                   ratio_free = 828,  # humans per free dog
                   ratio_dom  = 3.3,  # humans per domestic dog
                   r = 0.0005) {      # interaction radius (km)
  
  # Dog densities
  D_free <- H / ratio_free
  D_dom  <- H / ratio_dom
  
  # Total dog density
  D_total <- D_free + D_dom
  
  # Host–tick encounter rate
  R <- D_total * v * (2 * r) / pi
  
  return(R)
}

############################################################
## Example usage
############################################################

# Example human densities (humans per km^2)
H_vals <- c(10, 50, 100, 500, 1000)

# Compute encounter rates
R_tick(H_vals)

# If you want low, mid, and high movement scenarios:
R_low  <- R_tick(H_vals, v = 5)     # lower bound from literature
R_mid  <- R_tick(H_vals, v = 8)     # central value
R_high <- R_tick(H_vals, v = 13.5)  # upper bound

# Print results
data.frame(
  Human_Density = H_vals,
  R_low  = R_low,
  R_mid  = R_mid,
  R_high = R_high
)
