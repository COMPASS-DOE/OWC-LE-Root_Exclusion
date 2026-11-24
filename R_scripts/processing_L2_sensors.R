# Roberta Peixoto 2026-11-24
# processing_L2_sensors.R  -  COMPASS L2 v2-1 sensors at OWC: weather, groundwater, soil (TEROS)
#
# Saves (hourly, Etc/GMT+5, whole record 2022-2025; filter to the study window in the analysis):
# Raw data: COMPASS L2 v2-1 is NOT in this repository (compass_dir).

require(pacman)
p_load(tidyverse, lubridate, arrow, here, plantecophys)

# Paths and choices
compass_dir <- "D:/COMPASS_L2-v2-1"
out_dir <- here::here("data", "L2_sensors", "processing_outputs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

site  <- "OWC"
years <- 2022:2025
tz_l  <- "Etc/GMT+5"

l2_files <- list.files(compass_dir, recursive = TRUE, full.names = TRUE, pattern = "\\.parquet$")

# Reads all L2 files of site/plots/variables, hourly mean per Plot, Sensor and variable
read_l2 <- function(vars, plots = "[A-Z]+") {
  f <- l2_files[grepl(paste0("^", site, "_(", paste(plots, collapse = "|"), ")_(",
                             paste(years, collapse = "|"), ")_(", paste(vars, collapse = "|"), ")_L2"),
                      basename(l2_files))]
  message(length(f), " files: ", paste(vars, collapse = ", "))
  purrr::map_dfr(f, read_parquet) %>%
    transmute(Site, Plot, Sensor_ID = as.character(Sensor_ID), research_name,
              TIMESTAMP = force_tz(TIMESTAMP, tz_l), Value) %>%
    filter(year(TIMESTAMP) %in% years)
}

# 1. Weather station (15 min -> hourly) 
# rain: sum per hour (mm); others: hourly mean. VPD from RH and air temperature (kPa).
# wx-slr-fd15 = solar radiation flux density (W m-2); wx-par-den15 = PAR (umol m-2 s-1).
wx <- read_l2(c("wx-rain15", "wx-rh15", "wx-tempavg15", "wx-slr-fd15", "wx-par-den15",
                "wx-windspeed15", "wx-barpress15")) %>%
  group_by(TIMESTAMP, research_name) %>%                       # one station per site
  summarise(Value = mean(Value, na.rm = TRUE), .groups = "drop") %>%
  mutate(datetime = floor_date(TIMESTAMP, "hour")) %>%
  group_by(datetime, research_name) %>%
  summarise(Value = if (first(research_name) == "wx-rain15") sum(Value, na.rm = TRUE) else mean(Value, na.rm = TRUE),
            n = sum(!is.na(Value)), .groups = "drop") %>%
  mutate(Value = if_else(n == 0 | !is.finite(Value), NA_real_, Value)) %>%
  select(-n) %>%
  pivot_wider(names_from = research_name, values_from = Value) %>%
  rename(any_of(c(rain_sum = "wx-rain15", rh = "wx-rh15", tair = "wx-tempavg15",
                  solar_wm2 = "wx-slr-fd15", par = "wx-par-den15",
                  wind_ms = "wx-windspeed15", pressure = "wx-barpress15"))) %>%
  mutate(VPD_hour = RHtoVPD(RH = rh, TdegC = tair))
summary(wx)

# 2. Groundwater (wells). gw-wl-below-surface: water level relative to the ground (m, negative = below).
# Also dissolved oxygen (gw-rdo-conc) and pH/ORP (gw-ph-orp), useful for CH4 (redox).
gw <- read_l2(c("gw-wl-below-surface", "gw-rdo-conc", "gw-ph-orp")) %>%
  mutate(datetime = floor_date(TIMESTAMP, "hour")) %>%
  group_by(datetime, Site, Plot, research_name) %>%
  summarise(Value = median(Value, na.rm = TRUE), .groups = "drop") %>%
  filter(is.finite(Value)) %>%
  pivot_wider(names_from = research_name, values_from = Value) %>%
  rename(any_of(c(gwl_median = "gw-wl-below-surface", do_mgL = "gw-rdo-conc",
                  ph = "gw-ph", orp_mV = "gw-orp")))
gw %>% group_by(Plot) %>% summarise(first = min(datetime), last = max(datetime), n = n()) %>% print()

# 3. Soil TEROS (all plots: UP surface vs 30 cm drying; TR and W linked to the water table)
soil <- read_l2(c("soil-vwc-10cm", "soil-vwc-30cm", "soil-temp-10cm", "soil-temp-30cm",
                  "soil-EC-10cm", "soil-EC-30cm")) %>%
  mutate(datetime = floor_date(TIMESTAMP, "hour"),
         variable = str_extract(research_name, "vwc|temp|EC"),
         depth_cm = as.integer(str_extract(research_name, "\\d+(?=cm)"))) %>%
  group_by(datetime, Site, Plot, Sensor_ID, variable, depth_cm) %>%
  summarise(value = mean(Value, na.rm = TRUE), .groups = "drop") %>%
  filter(is.finite(value))
soil %>% group_by(Plot, variable, depth_cm) %>%
  summarise(n_sensors = n_distinct(Sensor_ID), first = min(datetime), last = max(datetime), .groups = "drop") %>%
  print(n = Inf)

# 4. Save 
write_parquet(wx,   file.path(out_dir, "wx_hourly.parquet"),   compression = "zstd")
write_parquet(gw,   file.path(out_dir, "gw_hourly.parquet"),   compression = "zstd")
write_parquet(soil, file.path(out_dir, "soil_hourly.parquet"), compression = "zstd")
