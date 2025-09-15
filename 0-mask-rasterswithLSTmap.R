library(sf)
library(terra)
library(exactextractr)
library(dplyr)

# --- Inputs ---
shp_path <- "/Users/yuchang/PhD/PhD_research/swiss_heat_exposure/data/revision_r1/muni_.shp"

raster_paths <- list(
  albedo  = "/Users/yuchang/PhD/PhD_research/swiss_heat_exposure/data/revision_r1/vars_raster/re-masked-rasters/Swiss_L8_Albedo_2018_100m_noUImask.tif",
  evi     = "/Users/yuchang/PhD/PhD_research/swiss_heat_exposure/data/revision_r1/vars_raster/re-masked-rasters/Swiss_L8_EVI_2018_100m_noUImask_ge0_le1_masked.tif",
  imperv  = "/Users/yuchang/PhD/PhD_research/swiss_heat_exposure/data/revision_r1/vars_raster/re-masked-rasters/impervious_Switzerland_inverted_masked.tif",
  bheight = "/Users/yuchang/PhD/PhD_research/swiss_heat_exposure/data/revision_r1/vars_raster/re-masked-rasters/Swiss_Building_Height_2018_full.tif",
  dem     = "/Users/yuchang/PhD/PhD_research/swiss_heat_exposure/data/revision_r1/vars_raster/re-masked-rasters/Swiss_DEM_full_100m.tif"
)

out_shp <- "/Users/yuchang/PhD/PhD_research/swiss_heat_exposure/data/revision_r1/vars_raster/re-masked-rasters/muni_with_zonal_means.shp"

# --- Load polygons ---
muni <- st_read(shp_path, quiet = TRUE)

# Use 'id' as the join key
id_field <- "id"
stopifnot(id_field %in% names(muni))

# Helper to calculate unweighted mean for one raster
zonal_mean_one <- function(path, poly_sf, id_field) {
  r  <- rast(path)
  rr <- raster::raster(r)
  poly_r <- st_transform(poly_sf, crs(rr))
  vals <- exact_extract(rr, poly_r, fun = "mean", progress = FALSE)
  tibble(!!id_field := poly_r[[id_field]], mean = vals)
}

# Loop over rasters and merge into muni
for (nm in names(raster_paths)) {
  message("Processing: ", nm)
  res_tbl <- zonal_mean_one(raster_paths[[nm]], muni, id_field) |>
    rename(!!paste0(nm, "_mean") := mean)
  muni <- muni |> left_join(res_tbl, by = id_field)
}

# --- Save as new shapefile ---
st_write(muni, out_shp, delete_layer = TRUE)

message("Done. Wrote shapefile: ", out_shp)



library(sf)
library(terra)
library(exactextractr)
library(dplyr)

# ---- Inputs ----
# Use the shapefile that already contains your former five fields:
shp_path <- "/Users/yuchang/PhD/PhD_research/swiss_heat_exposure/data/revision_r1/muni_with_zonal_means.shp"
# (If you prefer to start from the base layer instead, set shp_path to muni_.shp)

heat_rasters <- list(
  heatwave_pr = "/Users/yuchang/PhD/PhD_research/swiss_heat_exposure/data/heat_stress_tif/heatwave_probability_max_95th_30_summer.tif",
  avgtemp25   = "/Users/yuchang/PhD/PhD_research/swiss_heat_exposure/data/heat_stress_tif/average_temp_days_2018_25.tif"
)

# Output (overwrite same file)
out_shp <- shp_path

# ---- Load polygons ----
muni <- st_read(shp_path, quiet = TRUE)

# Join key
id_field <- "id"
stopifnot(id_field %in% names(muni))
if (any(duplicated(muni[[id_field]]))) stop("Duplicates found in 'id' field.")

# Helper: unweighted mean (reproject polygons to raster CRS; raster unchanged)
zonal_mean_one <- function(path, poly_sf, id_field) {
  r  <- rast(path)
  rr <- raster::raster(r)
  poly_r <- st_transform(poly_sf, crs(rr))
  vals <- exactextractr::exact_extract(rr, poly_r, fun = "mean", progress = FALSE)
  tibble::tibble(!!id_field := poly_r[[id_field]], mean = vals)
}

# ---- Compute & append the two fields ----
for (nm in names(heat_rasters)) {
  message("Processing: ", nm)
  res_tbl <- zonal_mean_one(heat_rasters[[nm]], muni, id_field) |>
    dplyr::rename(!!paste0(nm, "_mean") := mean)
  muni <- dplyr::left_join(muni, res_tbl, by = id_field)
}

# ---- Save (overwrite layer) ----
sf::st_write(muni, out_shp, delete_layer = TRUE, quiet = TRUE)
message("Done. Updated shapefile: ", out_shp)
