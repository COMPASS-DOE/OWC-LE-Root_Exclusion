#Roberta Peixoto 2026 - 11 - 23

# processing_sapflow.R  -  COMPASS synoptic sap flux (Granier thermal dissipation probes)
# Adapted from sapflow_processing.qmd 

# Reads the COMPASS L2 v2-1 sap flow (2.5 cm) and weather (RH, air temperature) files,
# computes dTmax with the Baseliner-style night criteria, converts to sap flux density F
# (Granier), applies the QA/QC cascade and saves:
# data/sapflow/processing_outputs/sapflow_hourly_sensor.parquet  (F per sensor and hour)

# Raw data: COMPASS L2 v2-1 is NOT in this repository; only whoever regenerates the tables
# needs it, at compass_dir.
# Design tables: Design/sapflow_daily_issues.csv (technician notes), Design/sapflow_manual_checks.csv
# (days removed after visual check), Design/sapflow_sensors.csv (tree and species of each sensor).

require(pacman)
p_load(tidyverse, lubridate, arrow, here, slider, plantecophys)

#Paths and choices
# COMPASS L2 v2-1 folder (not in the repository)
compass_dir <- "D:/COMPASS_L2-v2-1"

out_dir <- here::here("data", "sap_flow", "processing_outputs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

sites <- c("OWC")          
plots <- c("UP", "TR")
years <- 2022:2025
data_start <- as.POSIXct("2022-01-01 00:00:00", tz = "Etc/GMT+5")
data_end   <- as.POSIXct("2026-01-01 00:00:00", tz = "Etc/GMT+5")

# 1. Read sap flow and weather from L2
l2_files <- list.files(compass_dir, recursive = TRUE, full.names = TRUE, pattern = "\\.parquet$")
pick <- function(vars, plot_pattern) {
  l2_files[grepl(paste0("^(", paste(sites, collapse = "|"), ")_", plot_pattern, "_(",
                        paste(years, collapse = "|"), ")_(", paste(vars, collapse = "|"), ")"),
                 basename(l2_files))]
}
sap_files <- pick("sapflow-3\\.5cm", paste0("(", paste(plots, collapse = "|"), ")"))
wx_files  <- pick(c("wx-rh15", "wx-tempavg15"), "[A-Z]+")
message(length(sap_files), " sap flow files | ", length(wx_files), " weather files")

# Sap flow (mV), 15 min, one row per sensor; deep sensors (xxD) are not used
sap <- purrr::map_dfr(sap_files, read_parquet) %>%
  transmute(Site, Plot, Sensor_ID = str_pad(as.character(Sensor_ID), 2, pad = "0"),
            TIMESTAMP = force_tz(TIMESTAMP, "Etc/GMT+5"), sapflow = Value) %>%
  filter(!str_detect(Sensor_ID, "D$"), TIMESTAMP >= data_start, TIMESTAMP < data_end) %>%
  distinct(Site, Plot, Sensor_ID, TIMESTAMP, .keep_all = TRUE)

# Weather station, one value per site and timestamp
wx <- purrr::map_dfr(wx_files, read_parquet) %>%
  mutate(TIMESTAMP = force_tz(TIMESTAMP, "Etc/GMT+5")) %>%
  group_by(Site, TIMESTAMP, research_name) %>%
  summarise(Value = mean(Value, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = research_name, values_from = Value)

# Tree (Tag) and species of each sap flow sensor, from the tree inventory.
# Trees with a sensor have the sensor number in the Sapflow column ("N" = no sensor).
species <- read_csv(here::here("data", "sap_flow", "LE_TreePlotInventory_OWC.csv"),
                    col_types = cols(.default = "c")) %>%
  filter(!is.na(Sapflow), Sapflow != "N") %>%
  transmute(Site, Plot, Sensor_ID = str_pad(Sapflow, 2, pad = "0"), Tree_ID = Tag, Species) %>%
  distinct(Site, Plot, Sensor_ID, .keep_all = TRUE)
message(nrow(species), " sap flow trees in the inventory (sensors without a match keep Species = NA)")

sf <- sap %>%
  left_join(wx, by = c("Site", "TIMESTAMP")) %>%
  left_join(species, by = c("Site", "Plot", "Sensor_ID")) %>%
  mutate(vpd = RHtoVPD(RH = `wx-rh15`, TdegC = `wx-tempavg15`))   # kPa
message(nrow(sf), " sap flow rows; VPD missing in ", round(100 * mean(is.na(sf$vpd)), 1), "%")

# 2. Technician notes: sensor/power problems set to NA 
# PORT = sensor number(s) separated by commas, or ALL (every sensor of that Site/Plot)
bad <- read_csv(here::here("Design", "sapflow_daily_issues.csv"), col_types = cols(.default = "c")) %>%
  mutate(start_dt = as.POSIXct(paste(S_DATE, trimws(gsub(" PM| AM", "", S_TIME))),
                               format = "%m/%d/%Y %H:%M:%S", tz = "Etc/GMT+5"),
         end_dt   = as.POSIXct(paste(END, trimws(gsub(" PM| AM", "", E_TIME))),
                               format = "%m/%d/%Y %H:%M:%S", tz = "Etc/GMT+5"),
         PORT = strsplit(PORT, ",\\s*")) %>%
  unnest(PORT) %>%
  mutate(PORT = trimws(PORT), PORT = if_else(PORT == "ALL", PORT, str_pad(PORT, 2, pad = "0")))

bad_rows <- sf %>% mutate(row_id = row_number()) %>%
  inner_join(bad, by = c("Site" = "SITE", "Plot" = "PLOT"), relationship = "many-to-many") %>%
  filter(PORT == "ALL" | PORT == Sensor_ID, TIMESTAMP >= start_dt, TIMESTAMP <= end_dt) %>%
  pull(row_id) %>% unique()
sf$sapflow[bad_rows] <- NA_real_
message(length(bad_rows), " readings set to NA from technician notes")

# 3. dTmax (Baseliner-style night criteria)
# Night = 21:00-06:00, assigned to the following day (day_of_night). A night reading is used
# for dTmax only if: sapflow not NA, VPD < 0.2 kPa, and not a spike (|Z| <= 2.5 within the night).
# dTmax = max accepted reading of the night (mV).
ZSCORE_THRESHOLD   <- 2.5
VPD_THRESHOLD      <- 0.2    # kPa
WINDOW_WIDTH_DTMAX <- 14     # days

sf <- sf %>%
  mutate(day = as_date(TIMESTAMP), hour = hour(TIMESTAMP),
         day_of_night = if_else(hour >= 21, day + 1, day), #COMPASS standard is 23 - 1 am but our data shows different pattern
         is_night = hour >= 21 | hour <= 6) %>% # COMPASS standard is 23 - 4 am but our data shows different pattern 
  group_by(Site, Plot, Sensor_ID, day_of_night) %>%
  mutate(sap_zscore = if_else(is_night,
                              (sapflow - mean(sapflow[is_night], na.rm = TRUE)) /
                                sd(sapflow[is_night], na.rm = TRUE), NA_real_),
         flag_spike = is_night & !is.na(sap_zscore) & abs(sap_zscore) > ZSCORE_THRESHOLD) %>%
  ungroup()

# Justification for the dTmax night window: % of readings with dT > dTmax_final by hour,
# computed with the original 23:00-04:00 window (grey) vs the adopted 21:00-06:00 window (blue)

# sf %>% filter(!is.na(sapflow), !is.na(dTmax_final), sapflow >= 0.1) %>%
#   group_by(Plot, hour) %>%
#   summarise(pct_invalid = 100 * mean(sapflow > dTmax_final), .groups = "drop") %>%
#   ggplot(aes(hour, pct_invalid, colour = Plot)) +
#   annotate("rect", xmin = 20.5, xmax = 23.5, ymin = -Inf, ymax = Inf, fill = "steelblue", alpha = 0.15) +
#   annotate("rect", xmin = -0.5, xmax = 6.5,  ymin = -Inf, ymax = Inf, fill = "steelblue", alpha = 0.15) +
#   annotate("rect", xmin = 22.5, xmax = 23.5, ymin = -Inf, ymax = Inf, fill = "grey40", alpha = 0.25) +
#   annotate("rect", xmin = -0.5, xmax = 4.5,  ymin = -Inf, ymax = Inf, fill = "grey40", alpha = 0.25) +
#   geom_line(linewidth = 0.8) + geom_point() +
#   scale_x_continuous(breaks = 0:23) +
#   labs(x = "Hour of day", y = "Readings with dT > dTmax (%)",
#        title = "dT exceeds dTmax just outside the 23:00-04:00 window",
#        subtitle = "Grey: original night window | Blue: adopted window (21:00-06:00)") +
#   theme_bw()

dtmax_night <- sf %>%
  filter(is_night) %>%
  mutate(accept = !is.na(sapflow) & !is.na(vpd) & vpd < VPD_THRESHOLD & !flag_spike) %>%
  group_by(Site, Plot, Sensor_ID, day_of_night) %>%
  summarise(dTmax = if (any(accept)) max(sapflow[accept]) else NA_real_, .groups = "drop")

# dTmax stability: nights with dTmax more than 2 SD from the 14-day rolling mean are removed,
# and gaps are filled with the nearest accepted dTmax (dTmax_final).
dtmax_night <- dtmax_night %>%
  arrange(Site, Plot, Sensor_ID, day_of_night) %>%
  group_by(Site, Plot, Sensor_ID) %>%
  mutate(roll_mean = slide_dbl(dTmax, ~ mean(.x, na.rm = TRUE), .before = 7, .after = 7),
         roll_sd   = slide_dbl(dTmax, ~ sd(.x, na.rm = TRUE),   .before = 7, .after = 7),
         dTmax_final = if_else(!is.na(roll_sd) & abs(dTmax - roll_mean) <= 2 * roll_sd, dTmax, NA_real_)) %>%
  fill(dTmax_final, .direction = "downup") %>%
  ungroup() %>%
  select(Site, Plot, Sensor_ID, day_of_night, dTmax, dTmax_final)

# 4. Granier sap flux density 
# F = 0.011899 * (dTmax / dT - 1)^1.231   (cm3 cm-2 s-1)
# NOTE: F uses the nightly dTmax_final (14-day stable, gap-filled)
# TODO - CHECK WITH NATE: readings with dT > dTmax_final give non-finite F (~3.6% of rows,
# mostly 20:00-07:00; UP sensors 06-08 also ~8-16% at midday). They are currently REMOVED
# by is.finite(F). Alternative: set F = 0 (dT >= dTmax = no flow), which keeps the zeros (2.3% TR and 6.4% UP)
# in hourly/daily means:
#   mutate(F = if_else(sapflow > dTmax_final, 0, 0.011899 * ((dTmax_final / sapflow) - 1)^1.231))
# Question: remove these readings, or treat them as zero flow?
sf <- sf %>%
  left_join(dtmax_night, by = c("Site", "Plot", "Sensor_ID", "day_of_night")) %>%
  mutate(F = 0.011899 * ((dTmax_final / sapflow) - 1)^1.231)
n0 <- nrow(sf)

# 5. QA/QC cascade on mV 
# Step 1 - bad mV days: >= 10% of the day's readings near zero (below the tree's 1st percentile
#          of positive mV) AND any F > 0.005 that day.
near_zero <- sf %>% filter(sapflow > 0) %>% group_by(Site, Plot, Sensor_ID) %>%
  summarise(near_zero_threshold = quantile(sapflow, 0.01, na.rm = TRUE), .groups = "drop")
bad_mV_days <- sf %>%
  left_join(near_zero, by = c("Site", "Plot", "Sensor_ID")) %>%
  filter(is.finite(F)) %>%
  group_by(Site, Plot, Sensor_ID, day) %>%
  summarise(flag_bad_mV = mean(sapflow <= near_zero_threshold, na.rm = TRUE) >= 0.1 &
              any(F > 0.005), .groups = "drop") %>%
  filter(flag_bad_mV)
sf <- sf %>% anti_join(bad_mV_days, by = c("Site", "Plot", "Sensor_ID", "day"))
n1 <- nrow(sf)

# Step 2 - sharp drop: reading below 10% of the median of the previous 2 h (8 readings),
#          when that median is >= 0.01 mV.
sf <- sf %>%
  arrange(Site, Plot, Sensor_ID, TIMESTAMP) %>%
  group_by(Site, Plot, Sensor_ID) %>%
  mutate(mV_prev = slide_dbl(sapflow, ~ median(.x, na.rm = TRUE), .before = 8, .after = -1)) %>%
  ungroup() %>%
  filter(!coalesce(mV_prev >= 0.01 & sapflow / mV_prev < 0.10, FALSE)) %>%
  select(-mV_prev)
n2 <- nrow(sf)

# Step 3 - sustained fault: sensor daily median below 10% of the plot daily median (>= 3 sensors)
#          for 2 or more consecutive days.
fault_days <- sf %>%
  filter(sapflow > 0) %>%
  group_by(Site, Plot, day) %>%
  mutate(plot_median = median(sapflow, na.rm = TRUE), n_sensors = n_distinct(Sensor_ID)) %>%
  group_by(Site, Plot, Sensor_ID, day) %>%
  summarise(ratio = median(sapflow, na.rm = TRUE) / first(plot_median),
            n_sensors = first(n_sensors), .groups = "drop") %>%
  filter(n_sensors >= 3, !is.na(ratio)) %>%
  arrange(Site, Plot, Sensor_ID, day) %>%
  group_by(Site, Plot, Sensor_ID) %>%
  mutate(is_low = ratio < 0.10, streak_id = cumsum(is_low != lag(is_low, default = FALSE))) %>%
  filter(is_low) %>%
  group_by(Site, Plot, Sensor_ID, streak_id) %>%
  filter(n() >= 2) %>%
  ungroup() %>%
  select(Site, Plot, Sensor_ID, day)
# (streaks are runs of consecutive low days present in the data, as in the original qmd)
sf <- sf %>% anti_join(fault_days, by = c("Site", "Plot", "Sensor_ID", "day"))
n3 <- nrow(sf)

# Why rows are lost in step 4 (per Site/Plot): low mV, missing dTmax, or other non-finite F
sf %>% group_by(Site, Plot) %>%
  summarise(n = n(),
            mV_below_0.1 = sum(sapflow < 0.1, na.rm = TRUE),
            mV_NA        = sum(is.na(sapflow)),
            dTmax_NA     = sum(is.na(dTmax) & !is.na(sapflow)),
            dTmaxfinal_NA = sum(is.na(dTmax_final) & !is.na(sapflow)),
            vpd_NA       = sum(is.na(vpd)),
            F_bad_other  = sum(!is.finite(F) & !is.na(dTmax) & !is.na(sapflow)),
            .groups = "drop") %>% print()

# Step 4 - mV below 0.1 (real zero or bad sensor), and F not finite
sf <- sf %>% filter(!coalesce(sapflow < 0.1, FALSE), is.finite(F))
n4 <- nrow(sf)

# Step 5 - days removed after visual check
manual_days <- read_csv(here::here("Design", "sapflow_manual_checks.csv"), col_types = cols(.default = "c")) %>%
  mutate(start_day = as_date(parse_date_time(start_day, c("ymd", "mdy"))),
         end_day   = as_date(parse_date_time(end_day,   c("ymd", "mdy"))),
         Sensor_ID = str_pad(Sensor_ID, 2, pad = "0"))
if (anyNA(manual_days$start_day) || anyNA(manual_days$end_day))
  stop("Dates not read in sapflow_manual_checks.csv: check the start_day/end_day columns")
manual_days <- manual_days %>%
  mutate(day = map2(start_day, end_day, seq, by = "day")) %>%
  unnest(day) %>%
  distinct(Site, Plot, Sensor_ID, day)

sf <- sf %>% anti_join(manual_days, by = c("Site", "Plot", "Sensor_ID", "day"))
n5 <- nrow(sf)

message("QA/QC rows: start ", n0, " | bad mV days -", n0 - n1, " | sharp drop -", n1 - n2,
        " | sustained fault -", n2 - n3, " | mV < 0.1 or F not finite -", n3 - n4,
        " | manual -", n4 - n5, " | final ", n5)

# MAD outlier step NOT applied (it removed too much real signal).

# 6. Hourly mean per sensor and unreliable periods 
# Winter months (sensors not maintained): Nov, Dec, Jan; and bad data from 14 Oct 2024 until
# the end of that winter.
sap_bad_start <- as.POSIXct("2024-10-14 00:00:00", tz = "Etc/GMT+5")
sap_bad_end   <- as.POSIXct("2025-06-01 00:00:00", tz = "Etc/GMT+5")

sapflow_hourly_sensor <- sf %>%
  mutate(datetime = floor_date(TIMESTAMP, "hour")) %>%
  group_by(datetime, Site, Plot, Sensor_ID, Tree_ID, Species) %>%
  summarise(F = mean(F), sapflow_mV = mean(sapflow), n_readings = n(), .groups = "drop") %>%
  mutate(F = if_else(month(datetime) %in% c(11, 12, 1) |
                       (datetime >= sap_bad_start & datetime < sap_bad_end), NA_real_, F))

# 7. Daily proxy per Site and Plot
# sap_mid: midday median (11-14 h); sap_sum: daily integral (cm3 cm-2 d-1).
# Median across sensors, complete days only (>= 20 h). F >= 0.1 treated as sensor error.
sapflow_daily_plot <- sapflow_hourly_sensor %>%
  mutate(F = if_else(F >= 0.1, NA_real_, F), DATE = as_date(datetime), hr = hour(datetime)) %>%
  group_by(Site, Plot, DATE, Sensor_ID) %>%
  summarise(sap_mid = median(F[hr %in% 11:14], na.rm = TRUE),
            sap_sum = sum(F, na.rm = TRUE) * 3600,
            n_hours = sum(is.finite(F)), .groups = "drop") %>%
  filter(n_hours >= 20, is.finite(sap_mid)) %>%
  group_by(Site, Plot, DATE) %>%
  summarise(sap_mid = median(sap_mid), sap_sum = median(sap_sum), n_sensors = n(), .groups = "drop")

# 8. Check and save 
sapflow_hourly_sensor %>% group_by(Site, Plot) %>%
  summarise(n_sensors = n_distinct(Sensor_ID), first = min(datetime), last = max(datetime),
            pct_valid = round(100 * mean(!is.na(F)), 1), .groups = "drop") %>% print()

# Daily mean F per sensor over the whole record
sapflow_hourly_sensor %>% filter(is.finite(F), F < 0.1) %>%
  group_by(Plot, Sensor_ID, DATE = as_date(datetime)) %>%
  summarise(F = mean(F), .groups = "drop") %>%
  ggplot(aes(DATE, F, colour = Sensor_ID)) +
  geom_line(linewidth = 0.3) +
  facet_wrap(~ Plot, ncol = 1) +
  labs(x = NULL, y = expression(F~(cm^3~cm^-2~s^-1)), title = "Daily mean sap flux density per sensor") +
  theme_bw()

# Mean diel cycle per plot, growing season (May-Sep)
sapflow_hourly_sensor %>% filter(is.finite(F), F < 0.1, month(datetime) %in% 5:9) %>%
  group_by(Plot, hr = hour(datetime)) %>%
  summarise(med = median(F), q25 = quantile(F, 0.25), q75 = quantile(F, 0.75), .groups = "drop") %>%
  ggplot(aes(hr, med, colour = Plot, fill = Plot)) +
  geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.2, colour = NA) +
  geom_line(linewidth = 0.8) +
  scale_x_continuous(breaks = seq(0, 23, 3)) +
  labs(x = "Hour of day", y = expression(F~(cm^3~cm^-2~s^-1)),
       title = "Diel cycle of sap flux density (May-Sep, median and IQR)") +
  theme_bw()

write_parquet(sapflow_hourly_sensor, file.path(out_dir, "sapflow_hourly_sensor.parquet"), compression = "zstd")
write_parquet(sapflow_daily_plot,    file.path(out_dir, "sapflow_daily_plot.parquet"),    compression = "zstd")
# Read with arrow::read_parquet(here::here("data", "sap_flow", "processing_outputs", "sapflow_hourly_sensor.parquet"))
