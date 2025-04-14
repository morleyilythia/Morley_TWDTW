################################################################################
# Time series modelling spatiotemporal changes in Biogeoclimatic Ecosystem
# Classification (BEC) zones between 1997 and 2019 in west-central 
# British Columbia, Canada.
# Author: Ilythia Morley
# Date: March 1st 2025
#
# Part 2: Raster Data Processing
################################################################################

################################################################################
# 📦 1. LOAD LIBRARIES
################################################################################
library(sf)
library(dplyr)
library(raster)
library(terra)
library(tools)
library(foreign)

################################################################################
# 📂 2. DEFINE INPUTS AND PARAMETERS
################################################################################
base_dir <- "C:\Users\User\Desktop\GLAD_Output"
mosaic_dir <- file.path(base_dir, "Mosaic")
dir.create(mosaic_dir, showWarnings = FALSE)

proj <- CRS("+proj=longlat +datum=NAD83 +no_defs")

valid_intervals <- as.character(c(
  402, 404, 405, 408, 424, 422, 426, 428, 431, 447, 449, 451, 453, 470, 472, 474,
  476, 493, 495, 497, 499, 516, 518, 520, 522, 539, 541, 543, 545, 562, 564, 566,
  568, 585, 587, 589, 591, 608, 610, 612, 614, 631, 633, 635, 637, 654, 656, 658,
  660, 677, 679, 681, 683, 700, 702, 704, 706, 723, 725, 727, 729, 747, 748, 750,
  752, 769, 771, 773, 775, 792, 794, 796, 798, 815, 817, 819, 821, 838, 840, 842,
  844, 861, 863, 865, 867, 884, 886, 888, 890, 907, 909, 911, 913, 423, 450, 452,
  454, 455, 471, 477, 492, 498, 501, 514, 519, 523, 542, 544, 547, 563, 570, 588,
  590, 593, 611, 613, 616, 638, 639, 652, 653, 657, 661, 662, 676, 680, 685, 703,
  705, 708, 721, 730, 731, 749, 754, 767, 772, 790, 795, 797, 814, 822, 836, 837,
  842, 843, 845, 868, 882, 883, 887, 889, 912
))

################################################################################
# 🧹 3. CLEAN FILE LIST TO KEEP VALID INTERVALS
################################################################################
tif_files <- list.files(base_dir, pattern = "\\.tif$", recursive = TRUE, full.names = TRUE)

for (file in tif_files) {
  filename <- file_path_sans_ext(basename(file))
  if (!(filename %in% valid_intervals)) {
    cat("Deleting:", file, "\n")
    file.remove(file)
  } else {
    cat("Keeping:", file, "\n")
  }
}

################################################################################
# ✅ 4. VERIFY FILE PRESENCE ACROSS YEARS
################################################################################
required_files <- paste0(valid_intervals, ".tif")
subdirs <- list.dirs(base_dir, full.names = TRUE, recursive = FALSE)

for (subdir in subdirs) {
  existing_files <- list.files(subdir, pattern = "\\.tif$", full.names = FALSE)
  missing_files <- setdiff(required_files, existing_files)
  
  if (length(missing_files) > 0) {
    cat("\n❌ Missing files in:", subdir, "\n")
    print(missing_files)
  } else {
    cat("\n✅ All files present in:", subdir, "\n")
  }
}

################################################################################
# 🗂️ 5. ORGANIZE FILES INTO GROUPS BY INTERVAL
################################################################################
all_tifs <- list.files(base_dir, pattern = "\\.tif$", full.names = TRUE, recursive = TRUE)
valid_tifs <- all_tifs[basename(all_tifs) %in% paste0(valid_intervals, ".tif")]
final_list <- split(valid_tifs, tools::file_path_sans_ext(basename(valid_tifs)))

################################################################################
# 🧩 6. MOSAIC FILES BY INTERVAL
################################################################################
for (interval in names(final_list)) {
  paths <- final_list[[interval]]
  if (length(paths) > 0) {
    cat("Mosaicking", interval, "with", length(paths), "rasters...\n")
    rasters <- lapply(paths, rast)
    mos <- do.call(mosaic, rasters)
    out_file <- file.path(mosaic_dir, paste0(interval, ".tif"))
    writeRaster(mos, out_file, overwrite = TRUE)
  }
}

################################################################################
# ⏳ 7. CREATE TIME INDEX
################################################################################
setwd(mosaic_dir)
idx_vec <- c()
for (file in dir(no.. = TRUE)) {
  idx_vec <- c(idx_vec, strtoi(file_path_sans_ext(file)))
}
timeline <- all_dates[idx_vec]

################################################################################
# 📏 8. DEFINE STUDY AREA AND CLIP RASTERS
################################################################################
studyarea_extent <- c(-124.37, -123.55, 50.9995, 52.0005)
gridRaster <- raster(extent(studyarea_extent), crs = crs(proj))
gridPoly <- rasterToPolygons(gridRaster)
gridPoly$ID <- 1:nrow(gridPoly)

################################################################################
# ⚙️ 9. PROCESS RASTER STACKS BY BAND
################################################################################
process_band <- function(band_num, mask_poly) {
  files <- list.files(mosaic_dir, full.names = TRUE)
  stack_band <- stack()
  for (file in files) {
    r <- raster(file, band = band_num)
    r <- crop(r, mask_poly)
    r <- mask(r, mask_poly)
    stack_band <- stack(stack_band, r)
  }
  return(stack_band)
}

bands <- lapply(1:7, function(i) approxNA(process_band(i, gridPoly), method = "linear", rule = 2, f = 0, ties = mean, NArule = 1))
names(bands) <- paste0("stack_mosaic_b", 1:7)
list2env(bands, .GlobalEnv)

################################################################################
# 🌱 10. CALCULATE VEGETATION INDICES: NDVI, NDPI, NBR
################################################################################
calc_vi_stack <- function(band_func, gridPoly) {
  vi_stack <- stack()
  for (file in dir(mosaic_dir, pattern = "\\.tif$", full.names = TRUE)) {
    vi <- band_func(file)
    vi <- crop(vi, gridPoly[])
    vi <- mask(vi, gridPoly[])
    vi_stack <- stack(vi_stack, vi)
  }
  approxNA(vi_stack, method = "linear", rule = 2, f = 0, ties = mean, NArule = 1)
}

ndvi_stack <- calc_vi_stack(function(file) {
  NIR <- raster(file, band = 4)
  red <- raster(file, band = 3)
  (NIR - red) / (NIR + red)
}, gridPoly)

ndpi_stack <- calc_vi_stack(function(file) {
  NIR <- raster(file, band = 4)
  red <- raster(file, band = 3)
  swir1 <- raster(file, band = 5)
  (NIR - (0.74 * red + (1 - 0.74) * swir1)) / (NIR + (0.74 * red + (1 - 0.74) * swir1))
}, gridPoly)

nbr_stack <- calc_vi_stack(function(file) {
  NIR <- raster(file, band = 4)
  swir1 <- raster(file, band = 5)
  (NIR - swir1) / (NIR + swir1)
}, gridPoly)

################################################################################
# ✅ 11. VALIDATE STACKS
################################################################################
stopifnot(nlayers(ndvi_stack) == nlayers(ndpi_stack), nlayers(ndpi_stack) == nlayers(nbr_stack))
stopifnot(is.na(freq(ndvi_stack, value = NA)) == FALSE)
stopifnot(is.na(freq(ndpi_stack, value = NA)) == FALSE)
stopifnot(is.na(freq(nbr_stack, value = NA)) == FALSE)

################################################################################
# 📌 12. LOAD VALIDATION DATA
################################################################################
val_df <- read.dbf("C:/Users/User/Desktop/Morley_Scripts_Code/03_ValidationData/03_VectorScipt_Output/ValidationData_V3_11_Combined.dbf")
xy <- val_df[, c(2, 3)]
val_spdf <- SpatialPointsDataFrame(coords = xy, data = val_df, proj4string = proj)

################################################################################
# 🧭 13. ALIGN VALIDATION DATA TO NDVI EXTENT
################################################################################
r <- ndvi_stack[[1]]
rc <- mask(r, gridPoly[1,])
val_spdf@bbox <- as.matrix(extent(rc))
stopifnot(extent(val_spdf) == extent(rc))

################################################################################
# 📊 14. FINALIZE FIELD SAMPLES
################################################################################
field_samples <- as.data.frame(val_spdf)
field_samples <- subset(field_samples, select = -c(longitude.1, latitude.1))

################################################################################
# ⏱️ 15. CREATE RASTER TIME SERIES OBJECTS
################################################################################
rts <- twdtwRaster(stack_mosaic_b1, stack_mosaic_b2, stack_mosaic_b3, stack_mosaic_b4, stack_mosaic_b5, stack_mosaic_b6, timeline = timeline)
rts_VI <- twdtwRaster(ndvi_stack, ndpi_stack, nbr_stack, timeline = timeline)

################################################################################
# 🔍 16. SPLIT FIELD SAMPLES INTO TRAINING AND VALIDATION SETS
################################################################################
set.seed(1)
I <- unlist(createDataPartition(field_samples$label, p = 0.5))
training_samples <- field_samples[I, ]
validation_samples <- field_samples[-I, ]

table(training_samples[["label"]])

################################################################################
# 📉 17. EXTRACT TIME SERIES FOR TRAINING AND VALIDATION SAMPLES
################################################################################
proj <- crs("+proj=aea +lat_0=45 +lon_0=-126 +lat_1=50 +lat_2=58.5 +x_0=1000000 +y_0=0 +datum=NAD83 +units=m +no_defs ")

training_ts <- getTimeSeries(rts_VI, y = training_samples, proj4string = proj)
validation_ts <- getTimeSeries(rts_VI, y = validation_samples, proj4string = proj)

################################################################################
# 📈 18. CREATE TEMPORAL PATTERNS
################################################################################
temporal_patterns <- createPatterns(training_ts, from = "1998-05-25", to = "1998-09-14", freq = 16, formula = y ~ s(x, k = 3))
temporal_patterns_VI <- createPatterns(training_ts, from = "1997-05-25", to = "1997-09-14", freq = 16, formula = y ~ s(x, k = 3))

# Plot temporal patterns - All Bands
temp_pat_1997 <- plot(temporal_patterns, type = "patterns") +
  theme_bw() +
  scale_color_discrete(labels = c("Blue", "Green", "Red", "NIR", "SWIR1", "SWIR2"), name = "Bands") +
  theme(legend.title = element_text(colour = "black", size = 14), 
        legend.text = element_text(colour = "black", size = 14), 
        axis.text = element_text(colour = "black", size = 14),
        axis.title = element_text(colour = "black", size = 14), 
        plot.background = element_rect(color = "white"), 
        legend.position = "bottom") +
  labs(y = "Reflectance value")

# Plot temporal patterns - VI Only
temp_pat_1997_VI <- plot(temporal_patterns_VI, type = "patterns") +
  theme_bw() +
  scale_color_discrete(labels = c("NDVI", "NDPI", "NBR")) +
  theme(legend.title = element_text(colour = "black", size = 14), 
        legend.text = element_text(colour = "black", size = 14), 
        axis.text = element_text(colour = "black", size = 14),
        axis.title = element_text(colour = "black", size = 14), 
        plot.background = element_rect(color = "white"), 
        legend.position = "bottom") +
  labs(y = "Reflectance value") +
  guides(color = guide_legend(title = "Derived indices"))

# Plot temporal patterns - AT Class Only
temp_pat_AT <- plot(temporal_patterns, type = "patterns", labels = c(1)) +
  theme_bw() +
  scale_color_discrete(labels = c("Blue", "Green", "Red", "NIR", "SWIR1", "SWIR2"), name = "Bands") +
  theme(legend.title = element_text(colour = "black", size = 14), 
        legend.text = element_text(colour = "black", size = 14), 
        axis.text = element_text(colour = "black", size = 14),
        axis.title = element_text(colour = "black", size = 14), 
        plot.background = element_rect(color = "white"), 
        legend.position = "bottom")

################################################################################
# 🧠 19. CLASSIFY RASTER TIME SERIES USING TWDTW
################################################################################
log_fun <- logisticWeight(-0.1, 50)

beginCluster()
r_twdtw_final <- twdtwApply(
  x = rts, y = temporal_patterns, dates = timeline,
  weight.fun = log_fun, progress = 'text'
)
endCluster()

r_lucc <- twdtwClassify(r_twdtw_final, progress = 'text')

plotChanges(
  r_lucc, time.levels = c(1),
  class.labels = c("AT ", "CWH ", "CWH Dry", "ESSF ", "ESSF Dry", "IDF Dry", "MH ", "MS Dry", "SBPS Dry"),
  class.levels = 1:9,
  ylim = c("AT ", "CWH ", "CWH Dry", "ESSF ", "ESSF Dry", "IDF Dry", "MH ", "MS Dry", "SBPS Dry")
)

################################################################################
# 📏 20. ASSESS CLASSIFICATION ACCURACY
################################################################################
twdtw_assess <- twdtwAssess(
  object = r_lucc, y = validation_samples,
  proj4string = proj, conf.int = 0.95
)

plot(twdtw_assess, type = "accuracy") +
  xlim(c("AT ", "CWH ", "CWH Dry", "ESSF ", "ESSF Dry", "IDF Dry", "MH ", "MS Dry", "SBPS Dry"))

plotAccuracy(
  twdtw_assess, perc = TRUE, conf.int = 0.95,
  time.labels = NULL, category.name = NULL, category.type = NULL
)
################################ END OF SCRIPT #################################
################################################################################