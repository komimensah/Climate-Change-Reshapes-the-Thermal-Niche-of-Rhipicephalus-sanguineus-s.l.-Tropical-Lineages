#################
library(terra)
library(ggplot2)
library(patchwork)
library(dplyr)
library(sf)
library(rnaturalearth)
setwd('/Users/kagboka/Desktop/Anaplasma/')
# ------------------------------------------------------------
# 1. Load world boundaries and filter continents of interest
# ------------------------------------------------------------

world_sf <- dplyr::select(world_sf, continent)
world_sf <- dplyr::filter(world_sf, continent %in% target_continents)
world_sf <- ne_countries(scale = "medium", returnclass = "sf") %>%
  select(continent)

target_continents <- c("Africa", "Asia", "South America")
world_sf <- world_sf %>% filter(continent %in% target_continents)

world_vect <- vect(world_sf)

plot(world_sf)

# Update with your actual folder paths
factual_folder <- "/Users/kagboka/Desktop/Ritter_work/LTM_tif/factual/tas/"
counterfactual_folder <- "/Users/kagboka/Desktop/Ritter_work/LTM_tif/counterfactual/tas/"

shapefile_path <- "/Users/kagboka/Desktop/Ritter_work/MAin_new_project/Results/sahara_countries.shp"  # <- update this path

# ========== 2. Load raster files ==========
factual_files <- list.files(factual_folder, pattern = "\\.tif$", full.names = TRUE)
counter_files <- list.files(counterfactual_folder, pattern = "\\.tif$", full.names = TRUE)

factual_files <- sort(factual_files)
counter_files <- sort(counter_files)

# ========== 3. Load shapefile and convert to same CRS ==========
shape <- world_vect

# ========== 4. Stack and crop/mask ==========
factual_stack <- rast(factual_files)
counter_stack <- rast(counter_files)

# Reproject shapefile to match raster CRS
shape <- project(shape, crs(factual_stack))

# Clip rasters using shapefile
factual_stack <- mask(crop(factual_stack, shape), shape)
counter_stack <- mask(crop(counter_stack, shape), shape)

# ========== 5. Prepare for plotting ==========
month_labels <- month.name[1:12]   # <-- FULL month names
plot_list <- list()

for (i in 1:12) {
  factual_vals <- values(factual_temp_stack[[i]])
  counter_vals <- values(counterfactual_temp_stack[[i]])
  
  df <- data.frame(
    Temp = c(factual_vals, counter_vals),
    Scenario = rep(c("Factual", "Counterfactual"), each = length(factual_vals))
  )
  
  df <- df[complete.cases(df) & df$Temp > -10 & df$Temp < 60, ]
  
  p <- ggplot(df, aes(x = Temp, fill = Scenario)) +
    geom_density(alpha = 0.5) +
    scale_fill_manual(values = c("Factual" = "#377eb8",
                                 "Counterfactual" = "#e41a1c")) +
    labs(
      title = month_labels[i],      # <-- NO "Month:"
      x = "Temperature (°C)",
      y = "Density"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      panel.grid = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.6),
      axis.ticks = element_line(color = "black"),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 12)
    )
  
  plot_list[[i]] <- p
}

final_plotA <- wrap_plots(plot_list, ncol = 4, guides = "collect") +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

# Add ONE big panel label "A" (figure-level)
final_plotA <- final_plotA +
  plot_annotation(
    title = "A"
  ) &
  theme(
    plot.title = element_text(
      size = 22,
      face = "bold",
      hjust = 0   # left-aligned like journal figures
    )
  )

ggsave(
  "Monthly_Temperature_Factual_vs_Counterfactual_NEW_COLORS.png",
  final_plotA,
  width = 14,
  height = 10,
  dpi = 300
)

print(final_plotA)
#####################################################Maping climate varibales humidity

# Update with your actual folder paths
factual_folder <- "/Users/kagboka/Desktop/Ritter_work/LTM_tif/factual/hurs/"
counterfactual_folder <- "/Users/kagboka/Desktop/Ritter_work/LTM_tif/counterfactual/hurs/"

shapefile_path <- "/Users/kagboka/Desktop/Ritter_work/MAin_new_project/Results/sahara_countries.shp"  # <- update this path

# ========== 2. Load raster files ==========
factual_files <- list.files(factual_folder, pattern = "\\.tif$", full.names = TRUE)
counter_files <- list.files(counterfactual_folder, pattern = "\\.tif$", full.names = TRUE)

factual_files <- sort(factual_files)
counter_files <- sort(counter_files)

# ========== 3. Load shapefile and convert to same CRS ==========
shape <- world_vect

# ========== 4. Stack and crop/mask ==========
factual_stack <- rast(factual_files)
counter_stack <- rast(counter_files)

# Reproject shapefile to match raster CRS
shape <- project(shape, crs(factual_stack))

# Clip rasters using shapefile
factual_stack <- mask(crop(factual_stack, shape), shape)
counter_stack <- mask(crop(counter_stack, shape), shape)

# ========== 5. Prepare for plotting ==========
month_labels <- month.name[1:12]   # FULL month names
plot_list <- list()

factual_rh_stack <- factual_stack
counterfactual_rh_stack <- counter_stack

for (i in 1:12) {
  factual_vals <- values(factual_rh_stack[[i]])
  counter_vals <- values(counterfactual_rh_stack[[i]])
  
  df <- data.frame(
    RH = c(factual_vals, counter_vals),
    Scenario = rep(c("Factual", "Counterfactual"), each = length(factual_vals))
  )
  
  df <- df[complete.cases(df) & df$RH >= 0 & df$RH <= 100, ]
  
  p <- ggplot(df, aes(x = RH, fill = Scenario)) +
    geom_density(alpha = 0.5) +
    scale_fill_manual(values = c("Factual" = "#984ea3",
                                 "Counterfactual" = "#4daf4a")) +
    labs(
      title = month_labels[i],      # NO "Month:"
      x = "Relative Humidity (%)",
      y = "Density"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      panel.grid = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.6),
      axis.ticks = element_line(color = "black"),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 12)
    )
  
  plot_list[[i]] <- p
}

final_plotB <- wrap_plots(plot_list, ncol = 4, guides = "collect") +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

# Add big panel label "B"
final_plotB <- final_plotB +
  plot_annotation(tag_levels = list("B")) &
  theme(plot.tag = element_text(size = 22, face = "bold"))

ggsave("Monthly_Relative_Humidity_Factual_vs_Counterfactual.png",
       final_plotB, width = 14, height = 10, dpi = 300)

print(final_plotB)


library(magick)
library(cowplot)
getwd()
# ---------------------------------------------
# Step 1: Read the two images
# ---------------------------------------------
library(magick)

# Step 1: Read and flatten
img1 <- image_read('/Users/kagboka/Desktop/Anaplasma/Monthly_Temperature_Factual_vs_Counterfactual_NEW_COLORS.png') %>% image_flatten()
img2 <- image_read('/Users/kagboka/Desktop/Anaplasma/Monthly_Relative_Humidity_Factual_vs_Counterfactual.png') %>% image_flatten()




# Step 2: Resize to same width
common_width <- min(image_info(img1)$width, image_info(img2)$width)
img1 <- image_resize(img1, geometry = paste0(common_width))
img2 <- image_resize(img2, geometry = paste0(common_width))

# Step 3: Stack vertically (a on top of b)
combined <- image_append(c(img1, img2), stack = TRUE)

# Step 4: Add labels
labeled_plot <- image_annotate(combined, "", size = 80, gravity = "NorthWest", location = "+30+30", color = "black")
# Estimate label position for second image based on img1's height
label_b_y <- image_info(img1)$height + 30
labeled_plot <- image_annotate(labeled_plot, "", size = 80, gravity = "NorthWest",
                               location = paste0("+30+", label_b_y), color = "black")

# Step 5: Save the combined image
image_write(labeled_plot, path = "combined_facet_climate_trends.png", format = "png")



# ---------------------------------------------------------
# VC Delta & Diffusion under SSP1-2.6 and SSP2-4.5
# 4-panel figure with separate legends, using patchwork
# ---------------------------------------------------------

library(terra)
library(ggplot2)
library(tidyverse)
library(sf)
library(rnaturalearth)
library(cmocean)
library(patchwork)

# 1. Load shapefiles for clipping and borders
world_sf <- ne_countries(scale = "medium", returnclass = "sf")
world_vect <- vect(world_sf)

# 2. Load and clip rasters
load_and_prepare <- function(path) {
  r_hi <- rast(path)
  r_hi <- mask(crop(r_hi, ext(world_vect)), world_vect)
  r_lo <- r_hi   # keep original resolution
  return(as.data.frame(r_lo, xy = TRUE, na.rm = TRUE) %>%
           rename(vc = 3))
}

R0_CA_CF <- load_and_prepare('/Users/kagboka/Desktop/Anaplasma/R0_CA_noLandcover_counterfactual.tif')
R0_CA_F <- load_and_prepare('/Users/kagboka/Desktop/Anaplasma/R0_CA_noLandcover_factual.tif')
R0_FL_CF  <- load_and_prepare('/Users/kagboka/Desktop/Anaplasma/R0_FL_noLandcover_counterfactual.tif')
R0_FL_F  <- load_and_prepare('/Users/kagboka/Desktop/Anaplasma/R0_FL_noLandcover_factual.tif')

############################################################
## VC CLASS BREAKS, LABELS, AND COLOURS
############################################################

# Define class intervals
vc_breaks <- c(-Inf, 1, 5, 10, 20, 40, Inf)

# Labels starting explicitly at "<1"
vc_labels <- c(
  "<1",
  "1–5",
  "5–10",
  "10–20",
  "20–40",
  ">40"
)

# Colour palette for discrete VC classes (red-based, high risk = dark)
vc_colors <- c(
  "#fff5eb",  # very light peach (lowest class, <1)
  "#fdd0a2",  # light orange
  "#fdae6b",  # mid orange
  "#fd8d3c",  # strong orange
  "#e6550d",  # red-orange
  "#a63603"   # dark red (highest risk, >40)
)

############################################################
## FUNCTION TO PLOT A PANEL WITH FIXED CLASSES
############################################################
# Bold outline only for Africa, Asia, Oceania
world_bold <- world_sf %>% dplyr::filter(continent %in% c("Africa", "Asia", "Oceania", "South America"))
make_vc_plot <- function(df, title, legend_title = "VC class") {
  
  # Classify VC into fixed bins starting at "<1"
  df$vc_class <- cut(
    df$vc,
    breaks = vc_breaks,
    labels = vc_labels,
    include.lowest = TRUE,
    ordered_result = TRUE
  )
  
  ggplot() +
    geom_tile(data = df, aes(x = x, y = y, fill = vc_class)) +
   # geom_sf(data = world_sf, fill = NA, color = "gray10", linewidth = 0.2) +
    geom_sf(data = world_bold, fill = NA, color = "black", linewidth = 0.6) +
    coord_sf(xlim = c(-180, 180), ylim = c(-62, 50), expand = FALSE) +
    
    # Fixed discrete colour scale
    scale_fill_manual(
      values = vc_colors,
      drop = FALSE,
      name = legend_title
    ) +
    
    labs(title = title) +
    theme_void() +
    theme(
      plot.title = element_text(size = 13, face = "bold", hjust = 0),
      legend.position = c(0.75, 0.17),
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 8),
      legend.key.height = unit(0.3, "cm"),
      legend.key.width  = unit(0.4, "cm")
    )
}

# 4. Create all four panels
p1 <- make_vc_plot(R0_CA_CF, "a) Florida (tropical)", "counterfactual")
p2 <- make_vc_plot(R0_CA_F, "b) Florida (tropical)", "factual")
p3 <- make_vc_plot(R0_FL_CF,  "c) California (tropical)", "counterfactual")
p4 <- make_vc_plot(R0_FL_F,  "d) California (tropical)", "factual")
# 5. Combine into 4-panel layout
final_plot <- p1 / p2 / p3 / p4  # vertical stack
print(final_plot)
# Alternative 2x2 layout:
# final_plot <- (p1 | p2) / (p3 | p4)

# 6. Save output
ggsave("figureS_4no_comp.png",
       plot = final_plot, width = 8, height = 10, dpi = 300, bg = "white")
###########################################
#Figure 2
###########################################
library(terra)
library(sf)
library(tidyverse)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggplot2)
library(patchwork)

# ------------------------------------------------------------
# 1. Load world borders
# ------------------------------------------------------------
world_sf   <- ne_countries(scale = "medium", returnclass = "sf")
world_vect <- vect(world_sf)

# ------------------------------------------------------------
# 2. Function to load + crop + mask + downsample raster
# ------------------------------------------------------------
load_and_prepare <- function(path) {
  r <- rast(path)
  
  # Clip to world extent
  r <- crop(r, ext(world_vect))
  r <- mask(r, world_vect)
  
  # Reduce resolution for mapping (important)
  r_low <- aggregate(r, fact = 10, fun = mean)
  
  return(r_low)   # return raster, NOT dataframe
}

# ------------------------------------------------------------
# 3. Load all R0 rasters (Florida & California, factual & counterfactual)
# ------------------------------------------------------------
R0_CA_CF <- load_and_prepare('/Users/kagboka/Desktop/Anaplasma/R0_CA_noLandcover_counterfactual.tif')
R0_CA_F  <- load_and_prepare('/Users/kagboka/Desktop/Anaplasma/R0_CA_noLandcover_factual.tif')

R0_FL_CF <- load_and_prepare('/Users/kagboka/Desktop/Anaplasma/R0_FL_noLandcover_counterfactual.tif')
R0_FL_F  <- load_and_prepare('/Users/kagboka/Desktop/Anaplasma/R0_FL_noLandcover_factual.tif')

# ------------------------------------------------------------
# 4. Compute ΔR0 maps
# ------------------------------------------------------------
dCA <- R0_CA_F  - R0_CA_CF   # factual - counterfactual
dFL <- R0_FL_F  - R0_FL_CF   # factual - counterfactual

names(dCA) <- "California (temperate)"
names(dFL) <- "Florida (tropical)"

# ------------------------------------------------------------
# 5. Stack ΔR0 rasters and convert to dataframe
# ------------------------------------------------------------
stack_all <- c(dCA, dFL)

df <- as.data.frame(stack_all, xy = TRUE, na.rm = TRUE)
df$id_row <- seq_len(nrow(df))  # unique ID for joining

# ------------------------------------------------------------
# 6. Assign continent to each point
# ------------------------------------------------------------
points_sf <- st_as_sf(df, coords = c("x", "y"), crs = st_crs(world_sf))
points_sf$id_row <- df$id_row

world_cont <- ne_countries(scale = "medium", returnclass = "sf")[, c("continent")]

joined_sf <- suppressWarnings(st_join(points_sf, world_cont, join = st_within))

continent_df <- joined_sf %>% st_drop_geometry() %>% select(id_row, continent)

df <- left_join(df, continent_df, by = "id_row")

# Filter continents of interest
df <- df %>% filter(continent %in% c("Africa", "Asia", "South America"))

# ------------------------------------------------------------
# 7. Long format for plotting
# ------------------------------------------------------------
df_long <- df %>%
  pivot_longer(
    cols = c(`California (temperate)`, `Florida (tropical)`),
    names_to = "strain",
    values_to = "delta_R0"
  )

# ------------------------------------------------------------
# 8. Panel plot comparing changes
# ------------------------------------------------------------
p <- ggplot(df_long, aes(x = delta_R0, y = delta_R0, color = continent)) +
  geom_point(alpha = 0.22, size = 1.4) +
  geom_abline(slope = 1, intercept = 0, color = "black", linewidth = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
  xlab("ΔR₀ (change from counterfactual → factual)") +
  ylab("ΔR₀ (self-comparison)") +
  facet_wrap(~strain, scales = "free") +
  scale_color_brewer(palette = "Set2") +
  theme_classic(base_size = 11) +
  theme(
    strip.text     = element_text(face = "bold", size = 12),
    legend.position = "bottom",
    panel.grid = element_blank()
  )

print(p)
# ------------------------------------------------------------
# 8. Panel plot comparing Florida vs California ΔR0
# ------------------------------------------------------------

# Pivot wider to get two columns: delta_CA and delta_FL
df_wide <- df_long %>%
  select(x, y, continent, strain, delta_R0) %>%
  pivot_wider(names_from = strain, values_from = delta_R0)

# Rename for clarity
df_wide <- df_wide %>%
  rename(
    delta_CA = `California (temperate)`,
    delta_FL = `Florida (tropical)`
  )

# Plot FL vs CA
p <- ggplot(df_wide, aes(x = delta_CA, y = delta_FL, color = continent)) +
  geom_point(alpha = 0.45, size = 2) +
  geom_abline(slope = 1, intercept = 0, color = "black", linewidth = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  xlab("ΔRo — California (tropical strain)") +
  ylab("ΔRo — Florida (tropical strain)") +
  scale_color_brewer(palette = "Set2", name = "Continent") +
  theme_classic(base_size = 11) +
  theme(
    strip.text     = element_text(face = "bold", size = 12),
    legend.position = "bottom",
    panel.grid      = element_blank()
  )

print(p)
# ------------------------------------------------------------
# 10. Save high-res output
# ------------------------------------------------------------
ggsave("figureS13_style_main_vs_continent_from_rasters.png",
       p, width = 9, height = 6.5, dpi = 300, bg = "white")
###############################################



library(terra)
library(sf)
library(tidyverse)
library(rnaturalearth)
library(ggplot2)

# ------------------------------------------------------------
# 1. Load world borders
# ------------------------------------------------------------
world_sf <- ne_countries(scale = "medium", returnclass = "sf")
world_vect <- vect(world_sf)

# ------------------------------------------------------------
# 2. Function to load + crop + mask + downsample raster
# ------------------------------------------------------------
load_and_prepare <- function(path) {
  r <- rast(path)
  r <- crop(r, ext(world_vect))
  r <- mask(r, world_vect)
  r <- aggregate(r, fact = 10, fun = mean)
  return(r)
}

# ------------------------------------------------------------
# 3. Load R0 rasters
# ------------------------------------------------------------
R0_CA_CF <- load_and_prepare("R0_CA_noLandcover_counterfactual.tif")
R0_CA_F  <- load_and_prepare("R0_CA_noLandcover_factual.tif")

R0_FL_CF <- load_and_prepare("R0_FL_noLandcover_counterfactual.tif")
R0_FL_F  <- load_and_prepare("R0_FL_noLandcover_factual.tif")

# ------------------------------------------------------------
# 4. Compute ΔR0
# ------------------------------------------------------------
dCA <- R0_CA_F - R0_CA_CF
dFL <- R0_FL_F - R0_FL_CF

names(dCA) <- "delta_CA"
names(dFL) <- "delta_FL"

# ------------------------------------------------------------
# 5. Convert to dataframe
# ------------------------------------------------------------
stack_all <- c(dCA, dFL)
df <- as.data.frame(stack_all, xy = TRUE, na.rm = TRUE)

# ------------------------------------------------------------
# 6. Assign continent
# ------------------------------------------------------------
points_sf <- st_as_sf(df, coords = c("x", "y"), crs = st_crs(world_sf))
world_cont <- world_sf[, "continent"]

joined <- st_join(points_sf, world_cont, join = st_within)
df$continent <- joined$continent

df <- df %>%
  filter(continent %in% c("Africa", "Asia", "South America"))

# ------------------------------------------------------------
# 7. LINEAR REGRESSION (THIS IS WHAT THE REVIEWER WANTS)
# ------------------------------------------------------------
lm_fit <- lm(delta_FL ~ delta_CA, data = df)
lm_sum <- summary(lm_fit)

beta0 <- coef(lm_fit)[1]
beta1 <- coef(lm_fit)[2]
R2    <- lm_sum$r.squared

eq_label <- paste0(
  "ΔR₀_FL = ",
  round(beta0, 3), " + ",
  round(beta1, 3), "·ΔR₀_CA\n",
  "R² = ", round(R2, 3)
)

# ------------------------------------------------------------
# 8. Scatter plot WITH regression
# ------------------------------------------------------------
p <- ggplot(df, aes(x = delta_CA, y = delta_FL, color = continent)) +
  geom_point(alpha = 0.35, size = 1.6) +
  
  # 1:1 reference line
  #geom_abline(slope = 1, intercept = 0, color = "black", linewidth = 0.5) +
  
  # Fitted regression
  geom_smooth(method = "lm", se = FALSE,
              color = "black", linewidth = 0.8) +
  
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  
  annotate(
    "text",
    x = quantile(df$delta_CA, 0.05),
    y = quantile(df$delta_FL, 0.95),
    label = eq_label,
    hjust = 0,
    vjust = 1,
    size = 3.6
  ) +
  
  xlab("ΔR₀ — California (tropical strain)") +
  ylab("ΔR₀ — Florida (tropical strain)") +
  scale_color_brewer(palette = "Set2", name = "Continent") +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "bottom",
    panel.grid = element_blank()
  )

print(p)

# ------------------------------------------------------------
# 9. Save
# ------------------------------------------------------------
ggsave(
  "figureS13_deltaR0_regression.png",
  p, width = 9, height = 6.5, dpi = 300, bg = "white"
)


#Figure 3
#############################################
# Load packages
library(terra)
library(ggplot2)
library(tidyverse)
library(sf)
library(rnaturalearth)
library(cmocean)

# ---------------------------------------------------------
# 1. Load global shapefile (for land mask and borders)
# ---------------------------------------------------------

world_sf <- ne_countries(scale = "medium", returnclass = "sf")  # ~1:50m scale
world_vect <- vect(world_sf)  # convert to terra format for masking

# ---------------------------------------------------------
# 2. Load rasters and clip to land
# ---------------------------------------------------------
# ------------------------------------------------------------
# 3. Load all R0 rasters (Florida & California, factual & counterfactual)
# ------------------------------------------------------------
R0_CA_CF <- rast('/Users/kagboka/Desktop/Anaplasma/R0_CA_noLandcover_counterfactual.tif')
R0_CA_F  <- rast('/Users/kagboka/Desktop/Anaplasma/R0_CA_noLandcover_factual.tif')

R0_FL_CF <- rast('/Users/kagboka/Desktop/Anaplasma/R0_FL_noLandcover_counterfactual.tif')
R0_FL_F  <- rast('/Users/kagboka/Desktop/Anaplasma/R0_FL_noLandcover_factual.tif')

# ------------------------------------------------------------
# 4. Compute ΔR0 maps
# ------------------------------------------------------------
dCA <- R0_CA_F  - R0_CA_CF   # factual - counterfactual
dFL <- R0_FL_F  - R0_FL_CF   # factual - counterfactual

dCA <- R0_CA_F  - R0_CA_CF   # factual - counterfactual
dFL <- R0_FL_F  - R0_FL_CF   # factual - counterfactual
#r126 <- dCA
#r245 <- dFL
## ---- Symmetric percentage change Δs for California strain ----
# Denominator = mean(F, CF); add tiny epsilon to avoid /0
den_CA <- (R0_CA_F + R0_CA_CF) / 2

# Optionally drop cells where both are (almost) zero – change is meaningless there
den_CA[abs(den_CA) < 1e-6] <- NA

r126 <- (R0_CA_F - R0_CA_CF) / (den_CA + 1e-6)
names(r126) <- "Δs — California (tropical strain)"   # dimensionless

## ---- Symmetric percentage change Δs for Florida strain ----
den_FL <- (R0_FL_F + R0_FL_CF) / 2
den_FL[abs(den_FL) < 1e-6] <- NA

r245 <- (R0_FL_F - R0_FL_CF) / (den_FL + 1e-6)
names(r245) <- "Δs — Florida (tropical strain)"




vc_ssp126_hi <- r126
vc_ssp245_hi <- r245

# Crop and mask to land area only
vc_ssp126_hi <- mask(crop(vc_ssp126_hi, ext(world_vect)), world_vect)
vc_ssp245_hi <- mask(crop(vc_ssp245_hi, ext(world_vect)), world_vect)

# ---------------------------------------------------------
# 3. Downsample to ~1° resolution (change 'fact' as needed)
# ---------------------------------------------------------

vc_ssp126 <- vc_ssp126_hi
vc_ssp245 <- vc_ssp245_hi

# ---------------------------------------------------------
# 4. Convert to data frame for ggplot
# ---------------------------------------------------------

df_126 <- as.data.frame(vc_ssp126, xy = TRUE, na.rm = TRUE)
df_245 <- as.data.frame(vc_ssp245, xy = TRUE, na.rm = TRUE)

names(df_126)[3] <- "vc"
names(df_245)[3] <- "vc"
df_245$scenario <- "a) ΔRo — Florida (tropical strain)"
#df_126$scenario <- "a) ΔRo — California (tropical strain)"
#df_245$scenario <- "b) ΔRo — Florida (tropical strain)"
df_126$scenario <- "b) ΔRo — California (tropical strain)"
df_all <- bind_rows(df_126, df_245)

library(ggplot2)
library(dplyr)

df_all2 <- df_all %>%
  mutate(
    vc_class = cut(
      vc,
      breaks = c(-Inf, -2, -1, -0.5, 0, 0.5, 1, 2, Inf),
      labels = c("< -2", "-2 – -1", "-1 – -0.5", "-0.5 – 0",
                 "0 – 0.5", "0.5 – 1", "1 – 2", "> 2"),
      right = TRUE,
      include.lowest = TRUE
    )
  )

# Define your custom colour palette (named vector)
my_colors <- c(
  "#d73027",  # class 1: < -2 (dark red — strong negative change)
  "#fc8d59",  # class 2: -2 – -1
  "#fee08b",  # class 3: -1 – -0.5
  "#ffffbf",  # class 4: -0.5 – 0
  "#d9ef8b",  # class 5: 0 – 0.5
  "#91cf60",  # class 6: 0.5 – 1
  "#1a9850",  # class 7: 1 – 2
  "#006837"   # class 8: > 2 (dark green — strong positive change)
)
p <- ggplot() +
  geom_tile(data = df_all2, aes(x = x, y = y, fill = vc_class)) +
  geom_sf(data = world_sf, fill = NA, color = "gray10", linewidth = 0.2) +
  coord_sf(xlim = c(-180, 180), ylim = c(-62, 55), expand = FALSE) +
  facet_wrap(~ scenario, ncol = 1) +
  scale_fill_manual(
    values = my_colors,
    na.value = "grey90",
    name = "Ro (class)"
  ) +
  theme_void() +
  theme(
    strip.text = element_text(size = 13, face = "bold", hjust = 0),
    legend.position = c(0.70, 0.15),
    legend.title = element_text(size = 9),
    legend.text  = element_text(size = 8),
    legend.key.height = unit(0.3, "cm"),
    legend.key.width  = unit(0.4, "cm"),
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin = unit(c(10.5, 5.5, 5.5, 5.5), "points")
  )

print(p)
# ---------------------------------------------------------
# 6. Save the map
# ---------------------------------------------------------

ggsave("figureS_comp_nocomp67.png",
       plot = p, width = 8, height = 5.5, bg = "white", dpi = 300)
###################################################
#Figure 4
# Load packages
library(terra)
library(tidyverse)
setwd('/Users/kagboka/Desktop/Anophelese stephensis/')

# Load delta T (future - present) rasters per scenario

# ---------------------------------------------------------
# 2. Load rasters and clip to land
# ---------------------------------------------------------
# ------------------------------------------------------------
# 3. Load all R0 rasters (Florida & California, factual & counterfactual)
# ------------------------------------------------------------
R0_CA_CF <- rast('/Users/kagboka/Desktop/Anaplasma/R0_CA_noLandcover_counterfactual.tif')
R0_CA_F  <- rast('/Users/kagboka/Desktop/Anaplasma/R0_CA_noLandcover_factual.tif')

R0_FL_CF <- rast('/Users/kagboka/Desktop/Anaplasma/R0_FL_noLandcover_counterfactual.tif')
R0_FL_F  <- rast('/Users/kagboka/Desktop/Anaplasma/R0_FL_noLandcover_factual.tif')

# ------------------------------------------------------------
# 4. Compute ΔR0 maps
# ------------------------------------------------------------
dCA <- R0_CA_F  - R0_CA_CF   # factual - counterfactual
dFL <- R0_FL_F  - R0_FL_CF   # factual - counterfactual
#r126 <- dCA
#r245 <- dFL
## ---- Symmetric percentage change Δs for California strain ----
# Denominator = mean(F, CF); add tiny epsilon to avoid /0
den_CA <- (R0_CA_F + R0_CA_CF) / 2

# Optionally drop cells where both are (almost) zero – change is meaningless there
den_CA[abs(den_CA) < 1e-6] <- NA

r126 <- (R0_CA_F - R0_CA_CF) / (den_CA + 1e-6)
names(r126) <- "Δs — California (tropical strain)"   # dimensionless

## ---- Symmetric percentage change Δs for Florida strain ----
den_FL <- (R0_FL_F + R0_FL_CF) / 2
den_FL[abs(den_FL) < 1e-6] <- NA

r245 <- (R0_FL_F - R0_FL_CF) / (den_FL + 1e-6)
names(r245) <- "Δs — Florida (tropical strain)"

## Quick check of ranges
summary(r126)
summary(r245)



obs_temp <- rast('/Users/kagboka/Desktop/Lesh2.0/gpw_v4_population_density_adjusted_to_2015_unwpp_country_totals_rev11_2000_2pt5_min.tif')

obs_temp = obs_temp/828  + 0.8 * (obs_temp/3.3)
plot(obs_temp)
# --- 1. Load world shapefile ---
world_sf <- ne_countries(scale = "medium", returnclass = "sf")
world_vect <- vect(world_sf)


r <- r126
library(terra)



# Compute global summary statistics
stats <- global(r, c("min", "mean", "max", "std"), na.rm = TRUE)

print(stats)
r126 <- mask(crop(r126, ext(world_vect)), world_vect)
r245  <- mask(crop(r245,  ext(world_vect)), world_vect)

obs_temp  <- mask(crop(obs_temp,  ext(world_vect)), world_vect)



obs_temp <- project(obs_temp, r245)
r126     <- project(r126, r245)


library(terra)
library(tidyverse)
library(sf)
library(rnaturalearth)

# ------------------------------------------------------------
# 1. Load world boundaries and filter continents of interest
# ------------------------------------------------------------
world_sf <- ne_countries(scale = "medium", returnclass = "sf") %>%
  select(continent)

target_continents <- c("Africa", "Asia", "South America")
world_sf <- world_sf %>% filter(continent %in% target_continents)

world_vect <- vect(world_sf)

# ------------------------------------------------------------
# 2. Harmonise extents + resample
# ------------------------------------------------------------
ext_common <- intersect(ext(obs_temp), intersect(ext(r126), ext(r245)))

obs_temp <- crop(obs_temp, ext_common)
r126     <- crop(r126, ext_common)
r245     <- crop(r245, ext_common)

obs_temp <- resample(obs_temp, r245)
r126     <- resample(r126, r245)

names(r126) <- "ΔRo — California (tropical strain)"
names(r245) <- "ΔRo — Florida (tropical strain)"
names(obs_temp) <- "dog_density"

# ------------------------------------------------------------
# 3. Clip rasters to selected continents ONLY
# ------------------------------------------------------------
obs_temp_c <- mask(crop(obs_temp, world_vect), world_vect)
r126_c     <- mask(crop(r126, world_vect), world_vect)
r245_c     <- mask(crop(r245, world_vect), world_vect)

# ------------------------------------------------------------
# 4. Stack and convert to dataframe
# ------------------------------------------------------------
stack_all <- c(obs_temp_c, r126_c, r245_c)

df <- as.data.frame(stack_all, xy = TRUE, na.rm = TRUE)

# ------------------------------------------------------------
# 5. Assign continent to each pixel
# ------------------------------------------------------------
pts_sf <- st_as_sf(df, coords = c("x", "y"),
                   crs = st_crs(world_sf))
pts_sf$id_row <- seq_len(nrow(df))

cont_join <- st_join(pts_sf, world_sf)
df$continent <- cont_join$continent

df <- df %>% drop_na(continent, dog_density)

# ------------------------------------------------------------
# 6. Prepare long format for plotting
# ------------------------------------------------------------
df_long <- df %>%
  pivot_longer(
    cols = starts_with("ΔRo"),
    names_to = "scenario",
    values_to = "dRo"
  ) %>%
  drop_na(dRo)

# Mean ΔRo per scenario (global)
scenario_stats <- df_long %>%
  group_by(scenario) %>%
  summarise(scenario_mean = mean(dRo, na.rm = TRUE))

df_long <- df_long %>%
  left_join(scenario_stats, by = "scenario")

df_long <- df_long %>%
  mutate(
    scenario_label = paste0(scenario,
                            "\nMean ΔRo = ",
                            round(scenario_mean, 3))
  )

# ------------------------------------------------------------
# 7. Plot with continents colored
# ------------------------------------------------------------
continent_colors <- c(
  "Africa"        = "#d73027",  # red
  "Asia"          = "#4575b4",  # blue
  "South America" = "#7fc97f",  # green
  "Oceania"       = "#984ea3"   # purple
)

p <- ggplot(df_long, aes(x = dog_density, 
                         y = dRo, 
                         color = continent)) +
  
  geom_point(alpha = 0.45, shape = 16, size = 2) +   # <--- slightly stronger alpha
  
  geom_hline(aes(yintercept = scenario_mean),
             color = "black", linewidth = 0.8, linetype = "dashed") +
  
  scale_color_manual(values = continent_colors) +
  
  facet_wrap(~scenario_label, scales = "free") +
  
  labs(
    x = "Dog density (dogs per km²)",
    y = "ΔRo",
    color = "Continent"
  ) +
  
  theme_classic(base_size = 13) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom"
  )

print(p)
setwd('/Users/kagboka/Desktop/Anaplasma/')
# Save high-res version
ggsave("figure_dT_scenario_from_rasters.png", p, width = 7, height = 5, dpi = 300)
################################################
