
# processing_flux.R  -  OWC root exclusion: LI-COR multiplexer soil fluxes (CH4, CO2)

#Roberta Peixoto 2026 - 11 - 23

# Reads the raw multiplexer CSVs, combines them, sets the time zone, removes duplicates,
# applies the flux QA/QC and saves one final table:
#   soil_flux_licor/processing_outputs/flux_licor_qc_port.parquet
# (one row per chamber measurement; CH4 and CO2 flux, R2, SWC_1, TS_1, PORT, Treat, season)
#
# Raw data: the multiplexer CSVs are NOT in this repository (too large). They are backed up
# on the COMPASS drive; only whoever regenerates the final table needs them, at multiplexer_dir.
# Treatment of each port: Design/port_treatment.csv


require(pacman)
p_load(tidyverse, lubridate, data.table, arrow, here)

# Paths
# Raw multiplexer CSVs (backup on the COMPASS drive; not in the repository)
multiplexer_dir <- "C:/Users/rbitten/OneDrive - University of Toledo/Documents/Multiplexer/OWC_root_exclusion/multiplexer_data"

# Design table (PORT -> treatment) and output folder, inside the repository
port_file <- here::here("Design", "port_treatment.csv")
out_dir   <- here::here("data", "soil_flux_licor", "processing_outputs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# 1. Read and combine the raw CSVs 
# Each file has header lines on top; the column names are usually on the 2nd line. Each
# file gets the names from its OWN header line, and files are stacked by column name: some files have the
# columns in another order (with the names of the first file, DOY ended up in PORT).
# Only the raw LI-COR *_dense_summary.csv files. Excluded: the zipped backup folder and
# old combined data frames saved as df_*.csv. The same file can be in more than one
# download folder: the copies are removed later (duplicates, step 3).
csv_files <- list.files(multiplexer_dir, pattern = "_dense_summary\\.csv$", full.names = TRUE, recursive = TRUE)
csv_files <- csv_files[!grepl("multiplexer_zipped_backup", csv_files) &
                         !grepl("^df_", basename(csv_files))]
message(length(csv_files), " multiplexer files")

data_list <- list(); error_files <- c(); header_list <- list()
for (file in csv_files) {
  tryCatch({
    if (file.info(file)$size > 0) {
      # Header = first of the first 30 lines with DATE, PORT and FCH4 (line 2 in most files;
      # the first file of a deployment can have extra lines on top).
      top <- readLines(file, n = 30, warn = FALSE)
      h <- which(grepl("\\bDATE\\b", top) & grepl("\\bPORT\\b", top) & grepl("\\bFCH4\\b", top))[1]
      if (is.na(h)) {
        message("Skipped (no line with DATE, PORT and FCH4 names): ", file)
        error_files <<- c(error_files, file)
      } else {
        # all columns read as text (numbers are converted in step 2)
        data <- read_csv(file, skip = h - 1, col_names = TRUE, name_repair = "minimal",
                         col_types = cols(.default = col_character()), show_col_types = FALSE)
        headers <- trimws(names(data))
        empty <- is.na(headers) | headers == ""
        headers[empty] <- paste0("empty_", which(empty))          # unnamed (empty) columns
        colnames(data) <- headers
        data <- data[grepl("^[0-9]{8}$", trimws(data$DATE)), ]   # drops the units line and repeated headers
        if (nrow(data) > 0) {
          data_list[[length(data_list) + 1]] <- data
          header_list[[file]] <- paste(headers[!empty], collapse = ",")
        }
      }
    } else {
      message("Empty file: ", file)
    }
  }, error = function(e) {
    message("Error processing file: ", file, ": ", e$message)
    error_files <<- c(error_files, file)
  })
}
if (length(error_files) > 0) { message("Files with errors:"); print(error_files) }

# Files whose column names or order differ from the first file (stacked by name below)
hdr <- tibble(file = basename(names(header_list)), header = unlist(header_list))
hdr_diff <- hdr %>% filter(header != first(header))
message(nrow(hdr_diff), " of ", nrow(hdr), " files with columns in another order or with other names")
if (nrow(hdr_diff) > 0) print(hdr_diff$file)

complete_data <- data.table::rbindlist(data_list, use.names = TRUE, fill = TRUE)

# Keep the named columns, in the order of the first file (the empty columns are dropped)
col_keep <- strsplit(hdr$header[1], ",")[[1]]
data_combined <- complete_data[, col_keep, with = FALSE]

# 2. Date, time, time zone and treatment 
port_treatment <- read_csv(port_file, col_types = cols(PORT = col_character(), Treat = col_character(),
                                                       .default = col_character()))

df_fluxroot2024 <- data_combined %>%
  mutate(across(c(3, 5:51), ~ na_if(as.numeric(.), -9999)),
         TIMESTAMP = as.POSIXct(paste0(DATE, TIME), format = "%Y%m%d%H%M%S"),
         DATE = as.Date(DATE, format = "%Y%m%d"),
         TIME = format(strptime(TIME, format = "%H%M%S"), "%H:%M:%S"),
         PORT = as.character(PORT)) %>%
  left_join(port_treatment %>% select(PORT, Treat), by = "PORT") %>%
  mutate(Treat = as.factor(Treat))

if (anyNA(df_fluxroot2024$Treat))
  warning("Ports without a treatment in Design/port_treatment.csv: ",
          paste(unique(df_fluxroot2024$PORT[is.na(df_fluxroot2024$Treat)]), collapse = ", "))

# Logger clock = EST without daylight saving (UTC-5)
# Logger clock = EST without daylight saving (UTC-5)
df_fluxroot2024b <- df_fluxroot2024 %>%
  mutate(TIMESTAMP = force_tz(TIMESTAMP, tzone = "Etc/GMT+5"),
         datetime  = floor_date(TIMESTAMP, unit = "hour"))

# 3. Duplicates (same data saved in more than one file)
n_dup  <- df_fluxroot2024b %>% count(TIMESTAMP, PORT) %>% filter(n > 1) %>% nrow()
n_2prt <- df_fluxroot2024b %>% distinct(TIMESTAMP, PORT) %>% count(TIMESTAMP) %>% filter(n > 1) %>% nrow()
message(n_dup, " TIMESTAMP x PORT saved more than once | ",
        n_2prt, " timestamps with two different ports (should be 0)")

two_ports <- df_fluxroot2024b %>% distinct(TIMESTAMP, PORT, .keep_all = TRUE) %>%
  group_by(TIMESTAMP) %>% filter(n() > 1) %>% ungroup()
if (nrow(two_ports) > 0) print(two_ports %>% select(TIMESTAMP, PORT, FCH4, FCO2, FCH4_R2))

df_fluxroot2024b <- df_fluxroot2024b %>%
  arrange(TIMESTAMP) %>%
  distinct(TIMESTAMP, PORT, .keep_all = TRUE)

# 4. QA/QC
# Initial chamber CO2 > 550 ppm (flushing too short): kept only if FCO2 R2 >= 0.8.
num_obs_CO2_DRY <- sum(df_fluxroot2024b$CO2_DRY > 550, na.rm = TRUE)
num_obs_both    <- sum(df_fluxroot2024b$CO2_DRY > 550 & df_fluxroot2024b$FCO2_R2 > 0.8, na.rm = TRUE)
cat("Initial CO2 > 550:", round(100 * num_obs_CO2_DRY / sum(!is.na(df_fluxroot2024b$CO2_DRY)), 2), "% of obs;",
    round(100 * num_obs_both / num_obs_CO2_DRY, 2), "% of those with FCO2 R2 > 0.8\n")

# Outliers (median +- 2.5 MAD) are flagged only, not removed (CH4: could be bubbles or high fluxes).
mad_outlier <- function(x, ndev = 2.5) {
  xmed <- median(x, na.rm = TRUE); xmad <- mad(x, na.rm = TRUE)
  xout <- abs(x - xmed) / xmad > ndev
  message("Median = ", signif(xmed, 4), ", MAD = ", signif(xmad, 4), ", outliers = ",
          sum(xout, na.rm = TRUE), " of ", length(xout))
  xout
}

df_root_clean <- df_fluxroot2024b %>%
  mutate(month = floor_date(DATE, "month")) %>%
  group_by(month, Treat) %>%
  mutate(FLAG_CO2_Outliers = mad_outlier(FCO2),
         FLAG_CH4_Outliers = mad_outlier(FCH4)) %>%
  ungroup() %>%
  select(-month) %>%
  filter(FCH4 != -9999,
         FCO2 != -9999,
         FCO2 > 0,
         !(CO2_DRY > 550 & FCO2_R2 < 0.8),
         FCO2_R2 > 0,
         FCH4_R2 > 0)

df_root_clean$SWC_1[df_root_clean$SWC_1 == -9999] <- NA
df_root_clean$TS_1[df_root_clean$TS_1 == -9999]   <- NA
df_root_clean$SWC_1[df_root_clean$SWC_1 < 0.03]   <- 0.03   # below the detection limit

# Soil temperature sensor failure in July (TS_1 <= 1): set to NA within the failure window
# (2024-07-11 10:00 to 2024-07-31 15:00), flagged in TS_1_flag.
ts_fail <- df_root_clean %>%
  filter(month(datetime) == 7, TS_1 <= 1) %>%
  summarise(first_bad = min(datetime), last_bad = max(datetime), n = n())
print(ts_fail)

# 5. Final table 
flux_licor_qc_port <- as_tibble(df_root_clean) %>%
  select(DATE, TIME, DOY, PORT, Treat, FCH4, FCH4_R2, FCO2, FCO2_R2, SWC_1, TS_1, TIMESTAMP, datetime,FLAG_CH4_Outliers, FLAG_CO2_Outliers) %>%
  filter(!is.na(FCH4)) %>%
  mutate(Treat = as.factor(Treat),
         season = case_when(
           DATE >= as.Date("2024-03-20") & DATE < as.Date("2024-06-22") ~ "Spring",
           DATE >= as.Date("2024-06-22") & DATE < as.Date("2024-09-24") ~ "Summer",
           DATE >= as.Date("2024-09-24") & DATE < as.Date("2024-12-22") ~ "Fall",
           DATE >= as.Date("2024-12-22") & DATE < as.Date("2025-03-20") ~ "Winter",
           DATE >= as.Date("2025-03-20") & DATE < as.Date("2025-06-22") ~ "Spring",
           DATE >= as.Date("2025-06-22") & DATE < as.Date("2025-09-24") ~ "Summer"),
         TS_1_flag = if_else(datetime >= ts_fail$first_bad & datetime <= ts_fail$last_bad & TS_1 <= 1,
                             "sensor_failure_july", "ok", missing = "ok"),
         TS_1 = if_else(TS_1_flag == "sensor_failure_july", NA_real_, TS_1))

message(nrow(flux_licor_qc_port), " measurements | ", min(flux_licor_qc_port$DATE), " to ",
        max(flux_licor_qc_port$DATE), " | time zone ", tz(flux_licor_qc_port$datetime))
flux_licor_qc_port %>% count(PORT, Treat) %>% print()

# 6. Save 
# Parquet: ~4 MB instead of ~20 MB as CSV, and it keeps the time zone (Etc/GMT+5), dates
# and factors. Read it with arrow::read_parquet(here::here("data", "soil_flux_licor", "processing_outputs", "flux_licor_qc_port.parquet")).
out_file <- file.path(out_dir, "flux_licor_qc_port.parquet")
arrow::write_parquet(flux_licor_qc_port, out_file, compression = "zstd")
message("Saved ", out_file, " (", round(file.size(out_file) / 1e6, 1), " MB)")
