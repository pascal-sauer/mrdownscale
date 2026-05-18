# downscaleRun

Downscale MAgPIE results from the given folder. The resulting tgz file
will be generated in the given folder.

## Usage

``` r
downscaleRun(outputfolder, ...)
```

## Arguments

- outputfolder:

  path to a folder containing fulldata.gdx and clustermap\_\*.rds,
  resulting tgz will be written here

- ...:

  arguments passed on to retrieveData, e.g. model ("DOWNSCALEDMAGPIE" or
  "ESM"), rev, harmonizationPeriod

## Value

Invisibly, the path to the newly created tgz archive.

## Author

Pascal Sauer
