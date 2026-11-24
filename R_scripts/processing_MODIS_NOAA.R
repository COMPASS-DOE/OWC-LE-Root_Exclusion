# Roberta Peixoto 2026-11-24
# processing_external.R  -  MODIS_NOAA data downloaded from the web (run once; the analysis reads the CSVs)
#   MODIS GPP (MOD17A2HGF, 8-day composites, 500 m pixel at the OWC site)
#   Lake Erie water level (NOAA station 9063079, Marblehead OH)
# Saves:
#   data/MODIS_NOAA/processing_outputs/modis_gpp_8day.csv   DATE (start of the 8-day composite), gpp (g C m-2 d-1)
#   data/MODIS_NOAA/processing_outputs/lake_erie_daily.csv  DATE, lake (m, IGLD; daily mean of hourly heights)

require(pacman)
p_load(tidyverse, lubridate, janitor, MODISTools, here)

out_dir <- here::here("data", "MODIS_NOAA", "processing_outputs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

site_lat <- 41.378165
site_lon <- -82.508944
period   <- as.Date(c("2022-01-01", "2025-12-31"))
noaa_station <- "9063079"

# 1. MODIS GPP 
# Value = kg C m-2 per 8 days x 0.0001 (scale); > 32760 = fill. Converted to g C m-2 d-1.
# If the download fails (SSL/server), a previously downloaded copy is used when present.
gpp_cache <- here::here("data", "MODIS_NOAA", "modis_gpp_OWC_8day_2024_2025.rds")
gpp <- tryCatch(
  mt_subset(product = "MOD17A2HGF", band = "Gpp_500m", lat = site_lat, lon = site_lon,
            start = as.character(period[1]), end = as.character(period[2]),
            km_lr = 0, km_ab = 0, site_name = "OWC", internal = TRUE, progress = FALSE) %>%
    transmute(DATE = as_date(calendar_date),
              gpp  = if_else(value > 32760, NA_real_, value * 0.0001 * 1000 / 8)),
  error = function(e) {
    if (!file.exists(gpp_cache)) stop("MODIS download failed and no local copy at ", gpp_cache)
    message("MODIS download failed (", conditionMessage(e), "); using local copy ", basename(gpp_cache))
    readRDS(gpp_cache)
  })
message(nrow(gpp), " GPP composites, ", sum(is.na(gpp$gpp)), " fill values")

# 2. Lake Erie level (NOAA hourly heights, one request per month, local standard time) 
get_month <- function(m) {
  url <- paste0("https://api.tidesandcurrents.noaa.gov/api/prod/datagetter?",
                "product=hourly_height&station=", noaa_station,
                "&begin_date=", format(m, "%Y%m%d"),
                "&end_date=", format(ceiling_date(m, "month") - days(1), "%Y%m%d"),
                "&datum=IGLD&units=metric&time_zone=lst&format=csv&application=OWC_root_exclusion")
  out <- tryCatch(clean_names(read_csv(url, show_col_types = FALSE)), error = function(e) NULL)
  if (is.null(out) || !"water_level" %in% names(out)) { message("No data: ", format(m, "%Y-%m")); return(NULL) }
  tibble(datetime = if (inherits(out$date_time, "POSIXct")) out$date_time else ymd_hm(out$date_time),
         lake_level = as.numeric(out$water_level))
}
lake_hourly <- purrr::map_dfr(seq(floor_date(period[1], "month"), floor_date(period[2], "month"), by = "month"),
                              get_month)

lake_daily <- lake_hourly %>%
  group_by(DATE = as_date(datetime)) %>%
  summarise(lake = mean(lake_level, na.rm = TRUE), n_hours = sum(!is.na(lake_level)), .groups = "drop")
message(nrow(lake_daily), " lake days")

# 3. Check and save 
ggplot(gpp, aes(DATE, gpp)) + geom_line() + labs(y = "GPP (g C m-2 d-1)", title = "MODIS GPP") + theme_bw()
ggplot(lake_daily, aes(DATE, lake)) + geom_line() + labs(y = "Lake Erie level (m, IGLD)") + theme_bw()

write_csv(gpp,        file.path(out_dir, "modis_gpp_8day.csv"))
write_csv(lake_daily, file.path(out_dir, "lake_erie_daily.csv"))
