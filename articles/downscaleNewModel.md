# How to downscale data from a new model

## Starting Point

You have regional (i.e., non-gridded) land use data produced by an
integrated assessment model (we use “CoolNewModel” as placeholder name
here) that should be harmonized and downscaled using LUH3 as a
historical high-resolution reference dataset.

This vignette describes how to modify mrdownscale to process
CoolNewModel land data. It focuses on land variables representing the
distribution of land area across the different land use types
(e.g. urban, primary forest) only and does not cover downscaling of
non-land variables (e.g. fertilizer), even though mrdownscale is capable
of handling those (as implemented for MAgPIE). If you have questions or
run into trouble, please contact the mrdownscale developers (mail
address in the DESCRIPTION file).

## Overview

To integrate a new model into mrdownscale, you will need to create and
modify five key components. The table below summarizes each step;
detailed instructions for each follow in the sections below.

| Step | Component               | File                                              | Purpose                                                            |
|------|-------------------------|---------------------------------------------------|--------------------------------------------------------------------|
| 1    | Reference mapping       | `inst/extdata/referenceMappings/coolnewmodel.csv` | Map CoolNewModel variables to fine-grained reference variables     |
| 2    | Read function           | `R/readCoolNewModel.R`                            | Load your model’s data and region mapping via the madrat framework |
| 3    | Land input calculation  | `R/calcLandInput.R`                               | Read and validate your data for consistency                        |
| 4    | Resolution mapping      | `R/calcResolutionMapping.R`                       | Map from regional level to grid level                              |
| 5    | Categories mapping tool | `R/toolLandCategoriesMapping.R`                   | Register your reference mapping in the mapping tool                |

## Requirements and Preparation

### Files Needed

Before starting, prepare the following files:

- **CoolNewModel dataset file** - Your regional land use data in a
  format R can read (e.g. CSV)
- **Region mapping file** - Defines which countries/grid cells belong to
  which region
- **Variable mapping file** - Maps CoolNewModel variables to
  fine-grained reference variables (see examples in
  `inst/extdata/referenceMappings`)

### Software Setup

#### Install R and configure your environment

1.  Install R from <https://cran.rstudio.com/>

2.  Open R and find your home directory:

``` r
normalizePath("~")
```

3.  Create or edit a file called `.Rprofile` in your home directory with
    the following content:

``` r
options(repos = c(runiverse = "https://pik-piam.r-universe.dev",
                  CRAN = "https://cran.rstudio.com/"))
options(MADRAT_MAINFOLDER = file.path(normalizePath("~"), "madrat_mainfolder"))
```

This configuration enables you to install mrdownscale from the PIK piam
r-universe repository (it is not on the central CRAN repository) and
tells madrat where to store its data.

#### Install mrdownscale and development tools

Start a fresh R session and install the required packages:

``` r
install.packages(c("mrdownscale", "pkgload"))
```

### Getting Started with the Integration

We recommend creating a fork of the [mrdownscale
repository](https://github.com/pik-piam/mrdownscale), then cloning your
fork to a local folder. Apply all changes described in the following
sections to that local folder. To test your changes, always start a
fresh R session in that folder and run
[`pkgload::load_all()`](https://pkgload.r-lib.org/reference/load_all.html)
to load your modified mrdownscale directly from the source files.

For a concrete example of integrating a new model, see [this
PR](https://github.com/pik-piam/mrdownscale/pull/38), which integrated
the COFFEE model into mrdownscale.

## Integration Steps

### Step 1: Create Reference Mapping

A reference mapping file defines how CoolNewModel variables map to
fine-grained reference variables. This mapping enables mrdownscale to
translate CoolNewModel variables to any supported target dataset format.

Create a CSV file named `coolnewmodel.csv` in
`inst/extdata/referenceMappings/`. The file should have two columns:
CoolNewModel variable names and the corresponding reference variable
names.

If your model does not cover all land categories (for example, if urban
areas are not included), map missing categories to a variable called
“rest”. This ensures the total land area remains consistent.

Here’s an example reference mapping from the COFFEE model:

``` r
cat(readLines(system.file("extdata/referenceMappings/coffee.csv",
                          package = "mrdownscale"), n = 30), sep = "\n")
#> reference;data
#> primf;Land_Cover_Forest_Primary
#> secdf;Land_Cover_Forest_Secondary
#> forestry_ndc;Land_Cover_Forest_Secondary
#> forestry_aff_secdf;Land_Cover_Forest_Secondary
#> forestry_aff_pltns;Land_Cover_Forest_Planted
#> forestry_plant;Land_Cover_Forest_Planted
#> primn;Land_Cover_Other_Natural
#> secdn;Land_Cover_Other_Natural
#> pastr;Land_Cover_Pasture
#> range;Land_Cover_Pasture
#> urban;rest
#> c3per_irrigated_biofuel_1st_gen;rest
#> c3per_rainfed_biofuel_1st_gen;rest
#> c4per_irrigated_biofuel_1st_gen;rest
#> c4per_rainfed_biofuel_1st_gen;rest
#> c3ann_irrigated_biofuel_1st_gen;rest
#> c3ann_rainfed_biofuel_1st_gen;rest
#> c4ann_irrigated_biofuel_1st_gen;rest
#> c4ann_rainfed_biofuel_1st_gen;rest
#> c3nfx_irrigated_biofuel_1st_gen;rest
#> c3nfx_rainfed_biofuel_1st_gen;rest
#> begr, irrigated;Land_Cover_Cropland_Grassy
#> begr, rainfed;Land_Cover_Cropland_Grassy
#> betr, irrigated;Land_Cover_Cropland_Woody
#> betr, rainfed;Land_Cover_Cropland_Woody
#> alfalfa, irrigated;rest
#> alfalfa, rainfed;rest
#> Almonds, with shell, irrigated;Land_Cover_Cropland_Nuts
#> Almonds, with shell, rainfed;Land_Cover_Cropland_Nuts
```

For detailed technical information about variable mapping, see the
documentation of the `calcLandInputRecategorized` function:

``` r
?mrdownscale:::calcLandInputRecategorized
```

### Step 2: Create Read Function

mrdownscale uses the madrat framework to handle data loading. The madrat
framework automatically manages file organization, caching, and
reproducibility. Before reading your data, you must set up the file
structure and create a read function.

#### Set up the data folder structure

First, check where madrat stores source data:

``` r
madrat::getConfig("sourcefolder", verbose = FALSE)
#> [1] "/tmp/RtmpeTjCMX/madrat/sources"
```

In this folder, create a new subfolder called `CoolNewModel` and place
your data files and region mapping file there.

#### Create the read function

Write a function named `readCoolNewModel` in a new file
`R/readCoolNewModel.R`. This function must read two different data
types, specified via the `subtype` parameter:

- `subtype = "data"` - reads your land use data and returns it as a data
  frame
- `subtype = "resolutionMapping"` - reads your region mapping file and
  returns it as a data frame

Here’s an example read function to use as a template:

``` r
mrdownscale:::readWITCH
#> function (subtype = "data") 
#> {
#>     if (subtype == "data") {
#>         x <- utils::read.csv("db_ssp2_M_luh3.csv")
#>         return(list(x = x, class = "data.frame", unit = paste(unique(x$units), 
#>             collapse = ", "), description = "WITCH data"))
#>     }
#>     else if (subtype == "resolutionMapping") {
#>         mapping <- utils::read.csv("resolution_mapping_witch.csv")
#>         regions <- unique(mapping$witch17)
#>         mapping <- mapping[, setdiff(colnames(mapping), "lowRes")]
#>         addId <- data.frame(witch17 = regions, lowRes = paste0(regions, 
#>             ".", seq_along(regions)))
#>         mapping <- merge(mapping, addId, "witch17")
#>         return(list(x = mapping, class = "data.frame", description = "WITCH resolution mapping"))
#>     }
#>     else {
#>         stop("Unexpected subtype, only data and resolutionMapping are accepted")
#>     }
#> }
#> <bytecode: 0x55a1226c5968>
#> <environment: namespace:mrdownscale>
```

**Important:** The function name must match your source folder name
exactly (with “read” prefix). So for “CoolNewModel” folder, name it
`readCoolNewModel`.

#### Test the read function

Start a fresh R session and verify your setup works:

``` r
pkgload::load_all("path/to/mrdownscale")
x <- readSource("CoolNewModel", subtype = "data")
print(x)
map <- readSource("CoolNewModel", subtype = "resolutionMapping")
print(map)
```

Note that `readCoolNewModel` should never be called directly. Instead,
use [`readSource()`](https://rdrr.io/pkg/madrat/man/readSource.html),
which handles changing to the correct folder, caching results, and
loading cached data when nothing has changed.

### Step 3: Add to calcLandInput

Add a conditional branch for your model in `R/calcLandInput.R`:

``` r
if (input == "magpie") {
  # existing code ...
} else if (input == "CoolNewModel") {
  x <- readSource("CoolNewModel", subtype = "data")
  # apply validations here (see below)
}
```

In this section, apply the following consistency checks and
transformations to ensure your data meets mrdownscale’s requirements:

1.  **No missing values** - Ensure there are no NA or NaN values in your
    data
2.  **Non-negative values** - All variables must be \>= 0
3.  **Non-overlapping categories** - Variables should not overlap
    (e.g. avoid having both “forest” and “primary forest” as separate
    categories that sum to more than total forest area)
4.  **Constant total area over time** - The sum of all variables should
    remain constant across all years. If needed, add a “rest” category
    to account for missing land use types
5.  **Convert to magpie format** - Use `as.magpie(x)` to convert your
    data to the magpie object format
6.  **Match variable names** - Ensure all variable names exactly match
    those in your reference mapping file

### Step 4: Add to calcResolutionMapping

Add a conditional branch for your model in `R/calcResolutionMapping.R`:

``` r
if (input == "magpie") {
  # existing code ...
} else if (input == "CoolNewModel") {
  x <- readSource("CoolNewModel", subtype = "resolutionMapping")
  # create mapping to target grid here
}
```

This function must create a mapping from your regional level to the
target grid level. If you have a country-to-region mapping, you can use
the existing magpie mapping, which includes country-to-grid mappings.
This was done for COFFEE, see the `input == "coffee"` if-branch in
`R/calcResolutionMapping.R`.

### Step 5: Register in toolLandCategoriesMapping

Add your reference mapping file to the tool function in
`R/toolLandCategoriesMapping.R`:

1.  Place your variable mapping CSV file in
    `inst/extdata/referenceMappings/` (you created this in Step 1)
2.  Add a new condition to the if-else statement in
    `toolLandCategoriesMapping.R` that returns your mapping for your
    model, following the pattern used for other models

## Running the Harmonization

Once all integration steps are complete, you can calculate harmonized
data on your model’s native regional resolution:

``` r
harmonized <- calcOutput("LandHarmonized",
                         input = "CoolNewModel",
                         target = "luh3",
                         harmonizationPeriod = c(2025, 2050),
                         harmonization = "fadeForest")
write.magpie(harmonized, "harmonized.csv")
```

This produces a CSV file containing your data harmonized to LUH3 at your
model’s regional resolution. The harmonization period specifies the time
range to harmonize, and the harmonization method determines how your
model data is modified to smoothly continue from historical
observations.

## Running the Full Pipeline - Harmonization and Downscaling

To run the complete pipeline including both harmonization and
downscaling to grid resolution:

``` r
pathToTgz <- retrieveData("SCENARIOMIP", input = "CoolNewModel")
untar(pathToTgz)
```

After this completes successfully, your working directory should contain
NetCDF files (`.nc`) with the harmonized and downscaled land use data on
grid resolution.

## Checklist

Before running harmonization or downscaling, verify that you have
completed the following:

- Created reference mapping file:
  `inst/extdata/referenceMappings/coolnewmodel.csv`
- Created read function: `R/readCoolNewModel.R` (reads both “data” and
  “resolutionMapping” subtypes)
- Tested read function with
  `readSource("CoolNewModel", subtype = "data")` and
  `readSource("CoolNewModel", subtype = "resolutionMapping")`
- Added branch to `R/calcLandInput.R` with transformations, so all
  consistency checks pass
- Data has no missing values, all values \>= 0, non-overlapping
  categories, and matches reference mapping variable names
- Added branch to `R/calcResolutionMapping.R` to map regions to grid
  cells
- Added condition to `R/toolLandCategoriesMapping.R` to register the
  reference mapping
