
##########################
library(terra)
library(dplyr)
library(pROC)
library(ecospat)
library(ggplot2)

############################################################
# 1. Load your R0 raster (replace with your object if needed)
############################################################
# Example:
# r_R0 <- rast("R0_FL_strain.tif")

# Make sure r_R0 is loaded before continuing


############################################################
# 2. Load occurrence data
############################################################
# occ_df must contain: decimalLongitude, decimalLatitude
# occ_df <- read.csv("occurrences.csv")

library(readr)
# Make sure both rasters match in CRS, extent, resolution
vc_ssp126_hi <- project(vc_ssp126_hi, vc_ssp245_hi)
vc_ssp126_hi <- resample(vc_ssp126_hi, vc_ssp245_hi)

# Combine only the two desired layers
r_R0 <- (vc_ssp245_hi + vc_ssp126_hi) / 2

# Plot to confirm
plot(r_R0)
# --- Tick occurrence data (lon / lat in WGS84) ---
occ <- read_csv('/Users/kagboka/Desktop/Anaplasma/Rhipicephalus Sanguineus.csv')  # columns: lon, lat, (optionally other fields)

# Create coordinate matrix for terra
occ_xy <- as.matrix(occ[, c("decimalLongitude", "decimalLatitude")])


############################################################
# 3. Extract R0 values at occurrence points (CORRECT VERSION)
############################################################
occ_extract <- terra::extract(r_R0, occ_xy)

# Extract values
occ_extract <- terra::extract(r_R0, occ_xy)

# Automatically detect the raster value column
val_col <- setdiff(names(occ_extract), "ID")

# Extract R0 at occurrences
occ_vals <- occ_extract[[val_col]]

# Remove NA
occ_vals <- occ_vals[!is.na(occ_vals)]


# Add presence column
occ_pres_df <- data.frame(
  R0 = occ_vals,
  presence = 1
)


############################################################
# 4. Generate background points
############################################################
bg_pts <- spatSample(
  r_R0,
  size      = 10000,
  method    = "random",
  na.rm     = TRUE,
  as.points = TRUE,
  values    = TRUE
)

bg_df <- as.data.frame(bg_pts)

# Identify the raster column
valcol <- setdiff(names(bg_df), c("x", "y"))
names(bg_df)[names(bg_df) == valcol] <- "R0"

# Remove NA
bg_df <- bg_df %>% filter(!is.na(R0))

# Add presence = 0
bg_df$presence <- 0

# Background only needs R0 + presence
bg_pres_df <- bg_df[, c("R0", "presence")]


############################################################
# 5. Combine presence + background data
############################################################
eval_df <- bind_rows(occ_pres_df, bg_pres_df)


############################################################
# 6. AUC calculation
############################################################
roc_obj <- roc(eval_df$presence, eval_df$R0, quiet = TRUE)
auc_val <- auc(roc_obj)
cat("AUC: ", auc_val, "\n")


############################################################
# 7. Compute sensitivity-maximizing threshold (for TSS)
############################################################
opt <- coords(roc_obj, x = "best", best.method = "youden", transpose = TRUE)
threshold <- opt["threshold"]

# Sensitivity at threshold
TSS <- sens <- sensitivity(roc_obj, threshold = threshold)

cat("TSS (sensitivity at threshold): ", TSS, "\n")
cat("Threshold: ", threshold, "\n")


############################################################
# 8. Boyce Index (continuous)
############################################################
# Full raster values
bg_vals <- values(r_R0)
bg_vals <- bg_vals[!is.na(bg_vals)]

boyce_res <- ecospat.boyce(
  fit = bg_vals,
  obs = occ_vals,
  nclass = 0,
  window.w = "default",
  res = 100
)

cat("Boyce Index: ", boyce_res$cor, "\n")


############################################################
# 9. Plot Boyce curve for publication
############################################################
plot(boyce_res, type = "l", lwd = 3, col = "#1d91c0",
     main = "Continuous Boyce Curve",
     xlab = "Suitability bins", 
     ylab = "Predicted-to-Expected Ratio (P/E)")
abline(h = 1, lty = 2)
#############################
############################################################
library(terra)
library(sf)
library(blockCV)
library(dplyr)

# --- 0. Setup -----------------------------------------------------------------

# Your data: 
# r_R0 : SpatRaster
# occ_df : data.frame with decimalLongitude, decimalLatitude, presence = 1

occ_sf <- st_as_sf(
  occ_df,
  coords = c("decimalLongitude","decimalLatitude"),
  crs = 4326
)
r_R0 <- terra::project(r_R0, "EPSG:4326")

# --- 1. Define CV configurations to test -------------------------------------

configs <- expand.grid(
  k    = c(2, 3, 4, 5),         # number of folds to test
  size = c(100000, 200000, 500000, 1000000)  # block size in metres
)

valid_configs <- list()

# --- 2. Loop over each configuration, test for empty test folds ----------------

for (i in seq_len(nrow(configs))) {
  k_i    <- configs$k[i]
  size_i <- configs$size[i]
  
  set.seed(123 + i)  # vary seed so folds differ
  cv_obj <- cv_spatial(
    x      = occ_sf,
    column = "presence",
    r      = r_R0,
    k      = k_i,
    size   = size_i,
    hexagon = FALSE
  )
  
  # count presences in each test fold
  counts <- sapply(cv_obj$test_ids, length)
  
  if (all(counts > 0)) {
    valid_configs[[length(valid_configs)+1]] <-
      list(k = k_i, size = size_i, counts = counts)
  }
}

# --- 3. Print working configurations ------------------------------------------

if (length(valid_configs) == 0) {
  cat("No configuration found that assigns ≥1 presence to every fold.\n")
} else {
  cat("Valid block-CV configurations (no empty test folds):\n")
  for (cfg in valid_configs) {
    cat("  k =", cfg$k,
        "; block size =", cfg$size,
        "; test-fold counts =", paste(cfg$counts, collapse = ","),
        "\n")
  }
}
##############################
library(terra)
library(dplyr)
library(pROC)

# Prepare occurrence values
occ_xy   <- as.matrix(occ_df[, c("decimalLongitude","decimalLatitude")])
occ_vals <- terra::extract(r_R0, occ_xy)[,1]
occ_vals <- occ_vals[!is.na(occ_vals)]
occ_df2  <- data.frame(R0 = occ_vals, presence = 1)

# Background sampling
bg_pts <- spatSample(
  r_R0,
  size      = 10000,
  method    = "random",
  na.rm     = TRUE,
  as.points = TRUE,
  values    = TRUE
)
bg_df <- as.data.frame(bg_pts)
valcol <- setdiff(names(bg_df), c("x","y"))
names(bg_df)[names(bg_df) == valcol] <- "R0"
bg_df <- bg_df %>% filter(!is.na(R0))
bg_df$presence <- 0

# Combine
eval_df <- bind_rows(occ_df2, bg_df)

# AUC
roc_obj <- roc(eval_df$presence, eval_df$R0, quiet = TRUE)
auc_val <- auc(roc_obj)

# Threshold & TSS
cb <- coords(roc_obj, x = "best", ret = c("threshold", "sensitivity", "specificity"),
             best.method = "youden")
thr <- as.numeric(cb["threshold"])
pred   <- ifelse(eval_df$R0 > thr, 1, 0)
sens   <- mean(pred[eval_df$presence == 1] == 1)
spec   <- mean(pred[eval_df$presence == 0] == 0)
tss    <- sens + spec - 1

cat("AUC =", round(auc_val,3),
    "; threshold =", round(thr,3),
    "; TSS =", round(tss,3), "\n")
