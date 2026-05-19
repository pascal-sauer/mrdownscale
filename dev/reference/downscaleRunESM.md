# downscaleRunESM

Downscale MAgPIE results from the given folder to 0.25 degree resolution
in LUH2 format for ESMs. The resulting tgz file will be generated in the
given folder.

## Usage

``` r
downscaleRunESM(outputdir, revision = NULL, scenario = NULL, ...)
```

## Arguments

- outputdir:

  path to a folder containing fulldata.gdx and clustermap\_\*.rds,
  resulting tgz will be written here

- revision:

  passed on to retrieveData, default: current date

- scenario:

  passed on to retrieveData, default: slightly modified folder name

- ...:

  additional arguments passed on to retrieveData

## Value

Invisibly, the path to the newly created tgz archive.

## Author

Pascal Sauer
