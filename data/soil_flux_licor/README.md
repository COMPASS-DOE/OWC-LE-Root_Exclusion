# soil_flux_licor

Soil CH4 and CO2 fluxes from the LI-COR multiplexer chambers, OWC root exclusion experiment
(TR plot; ports 1-4 Control, 5-8 Root-free; see `Design/port_treatment.csv`).

## Raw data

The raw multiplexer CSVs are not in this repository (too large). They are backed up on the
COMPASS drive. Only `R_scripts/processing_flux.R` reads them, to regenerate the final table.

## processing_outputs/

| File | Content |
|---|---|
| `flux_licor_qc_port.parquet` | One row per chamber measurement after QA/QC: DATE, TIME, DOY, PORT, Treat, FCH4 and FCO2 (umol m-2 s-1), FCH4_R2, FCO2_R2, SWC_1 (m3 m-3), TS_1 (°C), TIMESTAMP and datetime (Etc/GMT+5, EST without daylight saving), season, TS_1_flag |

Read it in R with `arrow::read_parquet("soil_flux_licor/processing_outputs/flux_licor_qc_port.parquet")`.

## QA/QC (processing_flux.R)

- Duplicated measurements (same data saved in more than one file) removed.
- Removed: -9999 values, FCO2 <= 0, R2 <= 0, and initial chamber CO2 > 550 ppm with FCO2 R2 < 0.8. 
- Outliers (median ± 2.5 MAD) were checked but not removed (CH4 may be bubbles or high fluxes).
- SWC_1 below 0.03 (detection limit) set to 0.03.
- TS_1 <= 1 °C during the July 2024 sensor failure set to NA (TS_1_flag = "sensor_failure_july").
