# Mortality inputs are not distributed

Obtain the source life tables from the Human Mortality Database and prepare the four files described in the [application README](../README.md):

- `Data_Female.RData`, containing `Data_q`.
- `Data_Male.RData`, containing `Data_q`.
- `Data_raw_female.Rdata`, containing `Data_female`.
- `Data_raw_male.Rdata`, containing `Data_male`.

The application README specifies dimensions, annual ordering, probability/age grids, and the matching country-code order. The historical life-table-to-quantile preprocessing script is unavailable; exact refitting requires the prepared inputs with that original processing.

Place prepared files here locally, or set `MORTALITY_DATA_DIR` to an external directory. No data, saved fits, or downloaded files are included in the public distribution.
