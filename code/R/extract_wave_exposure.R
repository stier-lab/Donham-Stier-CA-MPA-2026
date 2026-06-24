# =============================================================================
# extract_wave_exposure.R
# =============================================================================
#
# PURPOSE:
#   Derive a per-MPA wave-exposure covariate from Tom Bell's island-inclusive
#   kelp-canopy + environment product, for use in the resilience moderator
#   analyses (16_environmental_moderators.R, 18_mpa_effectiveness_predictors.R).
#
# SOURCE (open-access, EDI):
#   Bell, T. W. 2023. "SBC LTER: Reef: California kelp canopy and environmental
#   variable dynamics ver 1." Environmental Data Initiative.
#   Package knb-lter-sbc.162.1   DOI: 10.6073/pasta/c40db2c8629cfa3fbe80fdc9e086a9aa
#   Single NetCDF entity "Kelp canopy and environmental variables" (~201 MB):
#     332,640 Landsat-pixel "stations" x 152 quarters (1984-2021), variables incl.
#     lat, lon, depth, area, biomass, temperature, npp, nitrate, and HSMAX = the
#     maximum significant wave height (m) from the CDIP MOP v1.1 swell-propagation
#     model. This is the SAME product Kumagai et al. (2024) used for `hsmax` and
#     that Wanner et al. (2024, Ecology) used for the CHANNEL ISLANDS -- it covers
#     the offshore islands (94,956 of the stations fall in the Channel-Islands box),
#     unlike the mainland-only 500-m coastline-segment product (knb-lter-sbc.144).
#
# METHOD:
#   wave exposure per MPA = mean of the climatological (time-mean) HSMAX over all
#   kelp-pixel stations within 3 km of the MPA (nearest station if none within 3 km).
#   Distance is haversine. Output is one row per MPA.
#
# INPUT  (external, not tracked):  ~/sbc-kelp-env/CAkelpCanopyEnv_sbc162.nc
#        (download: the EDI entity above; see ~/sbc-kelp-env/PROVENANCE.md)
# OUTPUT (tracked):                data/per_mpa_wave_exposure.csv
#        columns: MPA, Lat, Lon, nearest_km (km to nearest station),
#                 n_stations (within 3 km), wave_hs (mean HSMAX, m), island, coverage
#
# This script is STANDALONE and only needs re-running if the source NetCDF changes;
# the derived CSV is tracked so the analyses run without it. Skips gracefully if the
# NetCDF is absent.
#
# AUTHORS: Emily Donham & Adrian Stier
# =============================================================================

cat("Extracting per-MPA wave exposure from Bell 2023 kelp-canopy-env NetCDF...\n")
suppressMessages({library(ncdf4); library(data.table)})

NC <- path.expand("~/sbc-kelp-env/CAkelpCanopyEnv_sbc162.nc")
if (!file.exists(NC)) {
  message("  [wave] Source NetCDF not found at ", NC,
          "\n         Download EDI knb-lter-sbc.162.1 entity 5fbcb7b9780ad84157e3d4bbb0ab0947",
          "\n         to that path (see ~/sbc-kelp-env/PROVENANCE.md). Skipping; the tracked",
          "\n         data/per_mpa_wave_exposure.csv is used by the analyses as-is.")
} else {
  RADIUS_KM <- 3
  hav <- function(la1, lo1, la2, lo2) { R <- 6371
    dla <- (la2 - la1) * pi / 180; dlo <- (lo2 - lo1) * pi / 180
    a <- sin(dla / 2)^2 + cos(la1 * pi / 180) * cos(la2 * pi / 180) * sin(dlo / 2)^2
    2 * R * asin(sqrt(a)) }

  nc  <- nc_open(NC)
  lat <- ncvar_get(nc, "lat"); lon <- ncvar_get(nc, "lon")
  hs  <- ncvar_get(nc, "hsmax")                  # [station x time]
  nc_close(nc)
  hs_st <- rowMeans(hs, na.rm = TRUE); rm(hs); gc()   # climatological mean HSMAX per station
  keep <- is.finite(hs_st) & is.finite(lat) & is.finite(lon)
  lat <- lat[keep]; lon <- lon[keep]; hs_st <- hs_st[keep]

  meta <- as.data.table(read.csv(here::here("data", "harmonized", "harmonized_site_metadata.csv"),
                                 stringsAsFactors = FALSE))
  res <- meta[, {
    d <- hav(Lat, Lon, lat, lon)
    near <- which(d <= RADIUS_KM); if (length(near) == 0) near <- which.min(d)
    .(nearest_km = round(min(d), 2), n_stations = length(near),
      wave_hs = round(mean(hs_st[near]), 3), island = !is.na(ChannelIsland))
  }, by = .(MPA = CA_MPA_Name_Short, Lat, Lon)]
  res[, coverage := ifelse(nearest_km <= 5, "valid", "far")]

  fwrite(res, here::here("data", "per_mpa_wave_exposure.csv"))
  cat("  Wrote data/per_mpa_wave_exposure.csv:", nrow(res), "MPAs (",
      sum(res$island), "island,", sum(!res$island), "mainland);",
      sum(res$coverage == "valid"), "within 5 km.\n")
  cat("  HSMAX range:", paste(round(range(res$wave_hs), 2), collapse = " - "), "m",
      " (most exposed:", res$MPA[which.max(res$wave_hs)],
      "; most sheltered:", res$MPA[which.min(res$wave_hs)], ")\n")
}
