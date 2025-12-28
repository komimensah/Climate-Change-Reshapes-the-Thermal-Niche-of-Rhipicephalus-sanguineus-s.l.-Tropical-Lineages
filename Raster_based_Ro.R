
setwd('/Users/kagboka/Desktop/Anaplasma/')

library(terra)
library(ggplot2)
library(viridis)
library(tidyverse)

############################################################
## 1. LOAD RASTERS (Temperature, Humidity, Human Density)
############################################################

# ---- Temperature (12 months → mean) ----
#temp_folder <- '/Users/kagboka/Desktop/Ritter_work/LTM_tif/counterfactual/tas/'
#temp_files  <- list.files(temp_folder, pattern = "\\.tif$", full.names = TRUE)
temp_folder <- '/Users/kagboka/Desktop/Ritter_work/LTM_tif/factual/tas/'
temp_files  <- list.files(temp_folder, pattern = "\\.tif$", full.names = TRUE)
temp_stack  <- rast(temp_files)
r_T         <- mean(temp_stack, na.rm = TRUE)

# ---- Relative Humidity (12 months → mean) ----
#rh_folder <- '/Users/kagboka/Desktop/Ritter_work/LTM_tif/counterfactual/hurs/'
#rh_files  <- list.files(rh_folder, pattern = "\\.tif$", full.names = TRUE)
rh_folder <- '/Users/kagboka/Desktop/Ritter_work/LTM_tif/factual/hurs/'
rh_files  <- list.files(rh_folder, pattern = "\\.tif$", full.names = TRUE)
rh_stack  <- rast(rh_files)
r_RH      <- mean(rh_stack, na.rm = TRUE)

# ---- Human population density ----
r_H <- rast('/Users/kagboka/Desktop/Lesh2.0/gpw_v4_population_density_adjusted_to_2015_unwpp_country_totals_rev11_2000_2pt5_min.tif')
names(r_H) <- "H"

############################################################
## 2. HARMONISE RESOLUTION / EXTENT / CRS
############################################################

r_RH <- project(r_RH, r_T, method = "bilinear")
r_H  <- project(r_H,  r_T, method = "bilinear")

# intersection of rasters
ext_common <- intersect(ext(r_T), intersect(ext(r_RH), ext(r_H)))

r_T  <- crop(r_T, ext_common)
r_RH <- crop(r_RH, ext_common)
r_H  <- crop(r_H,  ext_common)

############################################################
## 3. DOG-DENSITY FROM HUMAN-DENSITY
############################################################

dog_density_fun <- function(H,
                            ratio_free = 828,
                            ratio_owned = 3.3,
                            roaming_fraction = 0.8) {
  dogs_free  <- H / ratio_free
  dogs_owned <- H / ratio_owned
  dogs_free + roaming_fraction * dogs_owned
}

# Convert humans/km² → dogs/km² pixel-wise
r_H_dogs <- app(r_H, fun = function(h) dog_density_fun(h))
names(r_H_dogs) <- "H_dogs"

############################################################
## 4. PREPARE R0 INPUT RASTER STACK (T, RH, H_dogs)
############################################################

r_in <- c(r_T, r_RH, r_H_dogs)
names(r_in) <- c("T", "RH", "H")

############################################################
## 5. R0 FUNCTIONS (NO LAND COVER ARGUMENT)
############################################################

R0_FL_fun <- function(T, RH, H){
  R0_strain(
    T = T, RH = RH, H = H,
    dev_fun_egg   = egg_rate_FL,
    dev_fun_larva = larva_rate_FL,
    dev_fun_nymph = nymph_rate_FL,
    mort_fun_L    = muL_FL,
    mort_fun_N    = muN_FL,
    mort_fun_A    = muA_FL,
    ovip_fun      = oviposition_FL
  )
}

#R0_NC_fun <- function(T, RH, H){
# R0_strain(
#  T = T, RH = RH, H = H,
#  dev_fun_egg   = egg_rate_NC,
#  dev_fun_larva = larva_rate_NC,
# dev_fun_nymph = nymph_rate_NC,
# mort_fun_L    = muL_NC,
#  mort_fun_N    = muN_NC,
#  mort_fun_A    = muA_NC,
#  ovip_fun      = oviposition_NC
#  )
#}

R0_CA_fun <- function(T, RH, H){
  R0_strain(
    T = T, RH = RH, H = H,
    dev_fun_egg   = egg_rate_CA,
    dev_fun_larva = larva_rate_CA,
    dev_fun_nymph = nymph_rate_CA,
    mort_fun_L    = muL_CA,
    mort_fun_N    = muN_CA,
    mort_fun_A    = muA_CA,
    ovip_fun      = oviposition_CA
  )
}

############################################################
## 6. APPLY R0 MODEL TO PIXELS (terra::lapp)
############################################################

R0_FL_raster <- lapp(r_in, fun = function(T, RH, H){ R0_FL_fun(T, RH, H) })
#R0_NC_raster <- lapp(r_in, fun = function(T, RH, H){ R0_NC_fun(T, RH, H) })
R0_CA_raster <- lapp(r_in, fun = function(T, RH, H){ R0_CA_fun(T, RH, H) })

names(R0_FL_raster) <- "R0_FL"
#names(R0_NC_raster) <- "R0_NC"
names(R0_CA_raster) <- "R0_CA"

R0_stack <- c(R0_FL_raster, R0_CA_raster)
R0_stack_bin <- R0_stack > 1   # Threshold map

plot(R0_stack_bin)

############################################################
## 7. GGPlot Example
############################################################

df_FL <- as.data.frame(R0_FL_raster, xy = TRUE, na.rm = TRUE)

ggplot(df_FL, aes(x = x, y = y, fill = R0_FL)) +
  geom_raster() +
  scale_fill_viridis(option = "magma") +
  coord_equal() +
  labs(
    x = "Longitude",
    y = "Latitude",
    title = expression(R[0] ~ "for" ~ italic(R.~sanguineus) ~ "FL strain")
  ) +
  theme_bw(base_size = 12)

############################################################
## 8. EXPORT RESULTS
############################################################

writeRaster(R0_FL_raster, "R0_FL_noLandcover_factual.tif", overwrite = TRUE)
#writeRaster(R0_NC_raster, "R0_NC_noLandcover.tif", overwrite = TRUE)
writeRaster(R0_CA_raster, "R0_CA_noLandcover_factual.tif", overwrite = TRUE)
