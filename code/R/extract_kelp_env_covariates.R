# =============================================================================
# extract_kelp_env_covariates.R
# =============================================================================
#
# PURPOSE:
#   Derive per-MPA environmental covariates from Tom Bell's island-inclusive
#   kelp-canopy + environment product, for the resilience moderator analyses
#   (16_environmental_moderators.R, 18_mpa_effectiveness_predictors.R). Covers
#   the Channel Islands (the per-pixel product Wanner et al. 2024 used there).
#   Supersedes the wave-only extract_wave_exposure.R.
#
# SOURCE (open-access, EDI):
#   Bell, T. W. 2023. "SBC LTER: Reef: California kelp canopy and environmental
#   variable dynamics ver 1." Environmental Data Initiative.
#   Package knb-lter-sbc.162.1   DOI: 10.6073/pasta/c40db2c8629cfa3fbe80fdc9e086a9aa
#   NetCDF "Kelp canopy and environmental variables" (~201 MB): 332,640 Landsat-
#   pixel "stations" x 152 quarters (1984-2021). Same product behind Kumagai et al.
#   (2024)'s hsmax. Variables extracted here:
#     hsmax       -> wave_hs    : max significant wave height (m), CDIP MOP v1.1
#     nitrate     -> nitrate    : seawater nitrate (uM), from SST empirical relationship
#     temperature -> temperature: sea surface temperature (deg C)
#     npp         -> npp        : net primary productivity proxy
#     depth       -> depth_m    : station depth (m)
#
#   Plus HUMAN GRAVITY (human-pressure / accessibility index, Cinner et al. 2018
#   framework: population within 50 km weighted by inverse-square travel distance),
#   per kelp patch, from the vendored Kumagai comparator file in
#   data/external/kumagai2024/ (or KUMAGAI_HUMAN_GRAVITY); 0 = no population within
#   50 km. Joined by nearest grid point; NA if the file is absent. This completes
#   the set of environmental covariates Kumagai mapped.
#
# METHOD:
#   Per covariate, take the climatological (time-mean) value at each kelp-pixel
#   station, then average over all stations within 3 km of the MPA (nearest station
#   if none within 3 km). Haversine distance. One row per MPA.
#
# INPUT  (external, not tracked):  SBC_KELP_ENV_NETCDF or the default local
#        Bell 2023 NetCDF development path.
# OUTPUT (tracked):                data/per_mpa_kelp_env.csv
#        columns: MPA, Lat, Lon, nearest_km, n_stations, wave_hs, nitrate,
#                 temperature, npp, depth_m, gravity, island, coverage
#
# Standalone; only re-run if the source NetCDF changes. The derived CSV is tracked
# so the analyses run without it. Skips gracefully if the NetCDF is absent.
#
# AUTHORS: Emily Donham & Adrian Stier
# =============================================================================

cat("Extracting per-MPA kelp-forest environmental covariates from Bell 2023 NetCDF...\n")
suppressMessages({library(ncdf4); library(data.table)})

NC <- path.expand(Sys.getenv(
  "SBC_KELP_ENV_NETCDF",
  unset = "~/sbc-kelp-env/CAkelpCanopyEnv_sbc162.nc"
))
if (!file.exists(NC)) {
  message("  [env] Source NetCDF not found at ", NC,
          "\n        Download EDI knb-lter-sbc.162.1 entity 5fbcb7b9780ad84157e3d4bbb0ab0947",
          "\n        to that path or set SBC_KELP_ENV_NETCDF. Skipping; the tracked",
          "\n        data/per_mpa_kelp_env.csv is used by the analyses as-is.")
} else {
  RADIUS_KM <- 3
  vars <- c(wave_hs = "hsmax", nitrate = "nitrate", temperature = "temperature",
            npp = "npp", depth_m = "depth")
  hav <- function(la1, lo1, la2, lo2) { R <- 6371
    dla <- (la2 - la1) * pi / 180; dlo <- (lo2 - lo1) * pi / 180
    a <- sin(dla / 2)^2 + cos(la1 * pi / 180) * cos(la2 * pi / 180) * sin(dlo / 2)^2
    2 * R * asin(sqrt(a)) }

  nc  <- nc_open(NC)
  lat <- ncvar_get(nc, "lat"); lon <- ncvar_get(nc, "lon")
  st <- data.table(lat = lat, lon = lon)            # per-station climatological means
  for (out in names(vars)) {                        # read one variable at a time (memory)
    v <- ncvar_get(nc, vars[out])
    st[[out]] <- if (is.matrix(v)) rowMeans(v, na.rm = TRUE) else as.numeric(v)
    rm(v); gc()
  }
  nc_close(nc)

  meta <- as.data.table(read.csv(here::here("data", "harmonized", "harmonized_site_metadata.csv"),
                                 stringsAsFactors = FALSE))
  res <- meta[, {
    d <- hav(Lat, Lon, st$lat, st$lon)
    near <- which(d <= RADIUS_KM); if (length(near) == 0) near <- which.min(d)
    sub <- st[near]
    c(list(nearest_km = round(min(d), 2), n_stations = length(near)),
      lapply(names(vars), function(x) round(mean(sub[[x]], na.rm = TRUE), 3)) |> setNames(names(vars)),
      list(island = !is.na(ChannelIsland)))
  }, by = .(MPA = CA_MPA_Name_Short, Lat, Lon)]
  res[, coverage := ifelse(nearest_km <= 5, "valid", "far")]

  # human gravity (human-pressure index) from Kumagai's per-kelp-patch grid (nearest point)
  GRAV <- {
    override <- Sys.getenv("KUMAGAI_HUMAN_GRAVITY", unset = "")
    local <- here::here("data", "external", "kumagai2024", "human_gravity_for_kelp_patches.csv")
    if (nzchar(override)) path.expand(override) else if (file.exists(local)) local else
      path.expand("~/kumagai2024-comparison/repo/Data/Population/human_gravity_for_kelp_patches.csv")
  }
  if (file.exists(GRAV)) {
    g <- as.data.table(read.csv(GRAV, stringsAsFactors = FALSE))
    res[, gravity := vapply(seq_len(.N), function(i) {
      d <- hav(Lat[i], Lon[i], g$lat, g$lon); g$gravity[which.min(d)] }, numeric(1))]
    res[, gravity := round(gravity, 1)]
  } else { message("  [env] Kumagai gravity file not found; gravity = NA."); res[, gravity := NA_real_] }
  setcolorder(res, c("MPA", "Lat", "Lon", "nearest_km", "n_stations", "wave_hs", "nitrate",
                     "temperature", "npp", "depth_m", "gravity", "island", "coverage"))

  fwrite(res, here::here("data", "per_mpa_kelp_env.csv"))
  cat("  Wrote data/per_mpa_kelp_env.csv:", nrow(res), "MPAs (",
      sum(res$island), "island,", sum(!res$island), "mainland);",
      sum(res$coverage == "valid"), "within 5 km.\n")
  cat("  wave_hs:", paste(round(range(res$wave_hs, na.rm = TRUE), 2), collapse = "-"), "m;",
      " nitrate:", paste(round(range(res$nitrate, na.rm = TRUE), 2), collapse = "-"), "uM;",
      " temp:", paste(round(range(res$temperature, na.rm = TRUE), 1), collapse = "-"), "C\n")
}
