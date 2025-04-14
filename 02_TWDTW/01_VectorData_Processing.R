################################################################################
# Time series modelling spatiotemporal changes in Biogeoclimatic Ecosystem Classification (BEC) zones between 1997 and 2019 in west-central British Columbia, Canada.
# Ilythia Morley
# March 1st 2025
#
# Part 1: Vector Data Processing
################################################################################
################################################################################
# 📦 1. LOAD LIBRARIES
################################################################################
library(sf)
library(dplyr)
library(raster)
library(terra)

################################ WATER POINTS ##################################
################################################################################
# 💧 2. DEFINE WATER POINT PROCESSING FUNCTION
################################################################################
# Function to generate water points with dates
generate_water_points <- function(shapefile_path, proj_string, from_dates, to_dates, counts, zone_label = "Water") {
  water_pts <- read_sf(shapefile_path) %>%
    st_transform(crs = proj_string)
  
  water_pts_all <- cbind(as.data.frame(water_pts),
                         st_coordinates(water_pts) %>%
                           as.data.frame() %>%
                           rename(longitude = X, latitude = Y))
  
  water_pts_all$from <- as.Date(unlist(Map(rep, from_dates, counts)))
  water_pts_all$to   <- as.Date(unlist(Map(rep, to_dates, counts)))
  water_pts_all$label <- zone_label
  
  return(water_pts_all)
}

################################################################################
# 💧 3. GENERATE WATER POINT VERSIONS V3–V11
################################################################################
# Define projection and shapefile
proj <- "+proj=longlat +datum=WGS84 +no_defs"
shapefile_path <- "C:/Users/User/Desktop/Morley_Scripts_Code/03_ValidationData/01_WaterPoints/01_Sampled_WaterPoints.shp"

# ----- V3 -----
from <- c("1997-05-25", "1998-05-25", "1999-05-25", "2000-05-25")
to   <- c("1997-09-14", "1998-09-14", "1999-09-14", "2000-09-14")
count <- c(1992, 1993, 1993, 1993)
water_V3 <- generate_water_points(shapefile_path, proj, from, to, count)

# ----- V4 -----
from <- c("2000-05-25", "2001-05-25", "2002-05-25", "2003-05-25")
to   <- c("2000-09-14", "2001-09-14", "2002-09-14", "2003-09-14")
count <- c(1992, 1993, 1993, 1993)
water_V4 <- generate_water_points(shapefile_path, proj, from, to, count)

# ----- V5 -----
from <- c("2003-05-25", "2004-05-25", "2005-05-25", "2006-05-25")
to   <- c("2003-09-14", "2004-09-14", "2005-09-14", "2006-09-14")
count <- c(1992, 1993, 1993, 1993)
water_V5 <- generate_water_points(shapefile_path, proj, from, to, count)

# ----- V6 -----
from <- c("2006-05-25", "2007-05-25", "2008-05-25")
to   <- c("2006-09-14", "2007-09-14", "2008-09-14")
count <- c(2657, 2657, 2657)
water_V6 <- generate_water_points(shapefile_path, proj, from, to, count)

# ----- V7 -----
from <- c("2003-05-25", "2004-05-25", "2005-05-25", "2006-05-25")
to   <- c("2003-09-14", "2004-09-14", "2005-09-14", "2006-09-14")
count <- c(1992, 1993, 1993, 1993)
water_V7 <- generate_water_points(shapefile_path, proj, from, to, count)

# ----- V8 -----
from <- c("2006-05-25", "2007-05-25", "2008-05-25")
to   <- c("2006-09-14", "2007-09-14", "2008-09-14")
count <- c(2657, 2657, 2657)  # fixed typo
water_V8 <- generate_water_points(shapefile_path, proj, from, to, count)

# ----- V9 -----
from <- c("2014-05-25", "2015-05-25", "2016-05-25")
to   <- c("2014-09-14", "2015-09-14", "2016-09-14")
count <- c(2657, 2657, 2657)
water_V9 <- generate_water_points(shapefile_path, proj, from, to, count)

# ----- V10 -----
from <- c("2016-05-25", "2017-05-25", "2018-05-25")
to   <- c("2016-09-14", "2017-09-14", "2018-09-14")
count <- c(2657, 2657, 2657)
water_V10 <- generate_water_points(shapefile_path, proj, from, to, count)

# ----- V11 -----
from <- c("2018-05-25", "2019-05-25")
to   <- c("2018-09-14", "2019-09-14")
count <- c(3985, 3986)
water_V11 <- generate_water_points(shapefile_path, proj, from, to, count)

# ----- Combine All Water Data -----

V3_11_Water_Combined <-  do.call("rbind", list(water_V3, water_V4, water_V5, water_V6, water_V7, water_V8, water_V9, water_V10, water_V11))

################################# BEC POINTS ###################################
################################################################################
# 🌲 4. DEFINE BEC POINT PROCESSING FUNCTION
################################################################################
process_bec_points <- function(shapefile_path, proj_string, from_dates, to_dates, zone_col = "ZONE") {
  # Read and transform shapefile
  sf_data <- read_sf(shapefile_path) %>%
    st_transform(crs = proj_string)
  
  # Extract coordinate columns
  coords <- st_coordinates(sf_data) %>%
    as.data.frame() %>%
    rename(longitude = X, latitude = Y)
  
  # Convert to data frame and bind coordinates
  df <- as.data.frame(sf_data) %>%
    bind_cols(coords)

  # Remove rows with NA in zone column
  df <- df %>% filter(!is.na(.data[[zone_col]]))
  
  # Generate even count split
  total_n <- nrow(df)
  n_bins <- length(from_dates)
  base_count <- floor(total_n / n_bins)
  remainder <- total_n %% n_bins
  counts <- rep(base_count, n_bins)
  if (remainder > 0) counts[1:remainder] <- counts[1:remainder] + 1
  
  # Add from/to dates
  df$from <- as.Date(unlist(Map(rep, from_dates, counts)))
  df$to   <- as.Date(unlist(Map(rep, to_dates, counts)))
  
  return(df)
}

################################################################################
# 🌲 5. GENERATE BEC POINT VERSIONS V3–V11
################################################################################
proj <- "+proj=longlat +datum=WGS84 +no_defs"

# ----- V3 -----
V3_POINTS <- process_bec_points(
  shapefile_path = "C:/Users/User/Desktop/Morley_Scripts_Code/03_ValidationData/02_BEC/BEC_V3_Points.shp",
  proj_string = proj,
  from_dates = c("1997-05-25", "1998-05-25", "1999-05-25", "2000-05-25"),
  to_dates = c("1997-09-14", "1998-09-14", "1999-09-14", "2000-09-14")
)

# ----- V4 -----
V4_POINTS <- process_bec_points(
  shapefile_path = "C:/Users/User/Desktop/Morley_Scripts_Code/03_ValidationData/02_BEC/BEC_V4_Points.shp",
  proj_string = proj,
  from_dates = c("2000-05-25", "2001-05-25", "2002-05-25", "2003-05-25"),
  to_dates = c("2000-09-14", "2001-09-14", "2002-09-14", "2003-09-14")
)

# ----- V5 -----
V5_POINTS <- process_bec_points(
  shapefile_path = "C:/Users/User/Desktop/Morley_Scripts_Code/03_ValidationData/02_BEC/BEC_V5_Points.shp",
  proj_string = proj,
  from_dates = c("2003-05-25", "2004-05-25", "2005-05-25", "2006-05-25"),
  to_dates = c("2003-09-14", "2004-09-14", "2005-09-14", "2006-09-14")
)

# ----- V6 -----
V6_POINTS <- process_bec_points(
  shapefile_path = "C:/Users/User/Desktop/Morley_Scripts_Code/03_ValidationData/02_BEC/BEC_V6_Points.shp",
  proj_string = proj,
  from_dates = c("2006-05-25", "2007-05-25", "2008-05-25"),
  to_dates = c("2006-09-14", "2007-09-14", "2008-09-14")
)

# ----- V7 -----
V7_POINTS <- process_bec_points(
  shapefile_path = "C:/Users/User/Desktop/Morley_Scripts_Code/03_ValidationData/02_BEC/BEC_V7_Points.shp",
  proj_string = proj,
  from_dates = c("2008-05-25", "2009-05-25", "2010-05-25"),
  to_dates = c("2008-09-14", "2009-09-14", "2010-09-14")
)

# ----- V8 -----
V8_POINTS <- process_bec_points(
  shapefile_path = "C:/Users/User/Desktop/Morley_Scripts_Code/03_ValidationData/02_BEC/BEC_V8_Points.shp",
  proj_string = proj,
  from_dates = c("2010-05-25", "2011-05-25", "2012-05-25", "2013-05-25", "2014-05-25"),
  to_dates = c("2010-09-14", "2011-09-14", "2012-09-14", "2013-09-14", "2014-09-14")
)

# ----- V9 -----
V9_POINTS <- process_bec_points(
  shapefile_path = "C:/Users/User/Desktop/Morley_Scripts_Code/03_ValidationData/02_BEC/BEC_V9_Points.shp",
  proj_string = proj,
  from_dates = c("2014-05-25", "2015-05-25", "2016-05-25"),
  to_dates = c("2014-09-14", "2015-09-14", "2016-09-14")
)

# ----- V10 -----
V10_POINTS <- process_bec_points(
  shapefile_path = "C:/Users/User/Desktop/Morley_Scripts_Code/03_ValidationData/02_BEC/BEC_V10_Points.shp",
  proj_string = proj,
  from_dates = c("2016-05-25", "2017-05-25", "2018-05-25"),
  to_dates = c("2016-09-14", "2017-09-14", "2018-09-14")
)

# ----- V11 -----
V11_POINTS <- process_bec_points(
  shapefile_path = "C:/Users/User/Desktop/Morley_Scripts_Code/03_ValidationData/02_BEC/BEC_V11_Points.shp",
  proj_string = proj,
  from_dates = c("2018-05-25", "2019-05-25"),
  to_dates = c("2018-09-14", "2019-09-14")
)

# ----- Combine All BEC Data -----

# Define the desired columns
keep_cols <- c("ZONE", "SUBZONE", "longitude", "latitude", "from", "to")

# Apply to each version
V3_POINTS       <- V3_POINTS[, keep_cols]
V4_POINTS       <- V4_POINTS[, keep_cols]
V5_POINTS       <- V5_POINTS[, keep_cols]
V6_POINTS       <- V6_POINTS[, keep_cols]
V7_POINTS       <- V7_POINTS[, keep_cols]
V8_POINTS       <- V8_POINTS[, keep_cols]
V9_POINTS       <- V9_POINTS[, keep_cols]
V10_POINTS      <- V10_POINTS[, keep_cols]
V11_POINTS      <- V11_POINTS[, keep_cols]

#Combine all Versions
V3_11_Combined <-  do.call("rbind", list(V3_POINTS, V4_POINTS, V5_POINTS, V6_POINTS, V7_POINTS, V8_POINTS, V9_POINTS, V10_POINTS, V11_POINTS))

################################################################################
# 🧼 6. CLEAN, COMBINE, AND SIMPLIFY BEC DATA
################################################################################
V3_11_Combined <- V3_11_Combined %>%
  mutate(
    # Merge IMA/CMA/BAFA into AT
    ZONE = case_when(
      ZONE %in% c("IMA", "CMA", "BAFA") ~ "AT",
      TRUE ~ ZONE
    ),
    
    # Recode SUBZONE to simplified categories
    SUBZONE = case_when(
      SUBZONE %in% c("dv", "dcp", "dcw", "dw", "dc", "ds", "xv", "dvp", "dm", "dvw",
                     "dk", "xc", "xm", "xw", "xh", "xk", "xcp", "xvw", "xvp", "xcw") ~ "Dry",
      SUBZONE %in% c("wk", "wc", "ws", "vm", "vk", "vh", "wcp", "ww", "wh", "wcw") ~ "Wet",
      SUBZONE %in% c("mm", "mc", "mv", "mk", "ms", "mw", "mh", "mwp", "mmp", "mcp", "mww",
                     "p", "unv", "unk", "unc", "unp", "un") ~ "",  # Moist & Undiff → blank
      TRUE ~ SUBZONE
    ),
    
    # Create final label
    label = paste(ZONE, SUBZONE)
  ) %>%
  dplyr::select(-ZONE, -SUBZONE)  # Drop the originals


################################################################################
# 🔗 7. COMBINE BEC + WATER DATASETS
################################################################################
keep_cols <- c("label", "longitude", "latitude", "from", "to")
V3_11_Water_Combined       <- V3_11_Water_Combined[, keep_cols]
V3_11_Combined       <- V3_11_Combined[, keep_cols]
V3_11_Combined <-  do.call("rbind", list(V3_11_Water_Combined, V3_11_Combined))


st_write(V3_11_Combined, "C:/Users/User/Desktop/Morley_Scripts_Code/03_ValidationData/03_VectorScipt_Output/ValidationData_V3_11_Combined.dbf")

rm("V3_POINTS", "V4_POINTS", "V5_POINTS", "V6_POINTS", "V7_POINTS", "V8_POINTS",
   "V9_POINTS", "V10_POINTS", "V11_POINTS", "water_pts", "Water_pts", "water_pts_all",
   "water_COORDS", "water_data_frame", "water_V3", "water_V4", "water_V5", "water_V6",
   "water_V7", "water_V8", "water_V9", "water_V10", "water_V11", "V3_11_Water_Combined")
################################ END OF SCRIPT #################################
################################################################################
