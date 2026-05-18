# Land Use Data for ESMs - Changelog

The same REMIND-MAgPIE runs are used for the downscaled land use dataset
described here and for the corresponding downscaled emissions dataset.
Find the changelog for the general setup including model versions and
settings, as well as for the emissions dataset at
<https://github.com/IAMconsortium/concordia/blob/main/CHANGELOG.md>.

## \## unreleased

## 2024-08-19

- changed harmonization interpolation function, mostly reduces
  difference at end of harmonization
- replace fulwd/rndwd NAs with country shares instead of 0.5
- all consistency check results are now written to a logfile
- calculate low to high resolution mapping (fill missing with nearest
  neighbor)
- data now matches land/sea mask of LUH2v2h/CMIP6
- fixed wood harvest is share of cell area
- fixed bug where total areas of a land type would not match before and
  after downscaling
- fixed: primf/primn expansion is no longer happening
- LUH data on wood harvest and fertilizer is now aggregated over n years
  to match input data semantics
- fixed a bug where wrong reference year was used for downscaling
- wood harvest area is now consistent with land
- nonland/management data is now downscaled using already downscaled
  land data as reference

### uploaded_2024-10-28

- use yearly LUH data from 1995-2015 as extrapolation base
- year 2015 is now included in the reported data

## 2024-04-25

- fixed flipped latitude
- removed bounds_time
- report data from 2015 to 2100 (as opposed to from 1995)
- added units to management and transitions nc files

## 2024-03-22

- fixed wood harvest area (unit) conversion bug
- removed NAs
- renamed longitude/latitude to lon/lat
- added bounds variables for lon/lat/time
- added more metadata, incl. harmonization period

## 2023-12-08

- dataset for 1995-2100 in 5-year timesteps for all REMIND-MAgPIE
  scenarios used for the emissions
- added metadata to nc files
- management.nc: added fulwd, rndwd
- transitions.nc: reporting all LUH variables incl. wood harvest
  variables
- all transitions are gross transitions
- transitions in 2100 are net-zero gross transitions (to be used for
  extended (2100-2300) scenarios)
- cover all LUH cells (e.g. Greenland was missing before)

## 2023-10-11

- dataset for 1995-2100 in 5-year timesteps for one scenario
- states.nc: reporting all LUH variables except secma, secmb
- management.nc: reporting all LUH variables except fulwd, rndwd, combf,
  flood, fharv_c3per, fharv_c4per
- additionally reporting 2nd gen biofuel (crpbf2_c3per, crpbf2_c4per)
  and managed forests (manaf) in management.nc
