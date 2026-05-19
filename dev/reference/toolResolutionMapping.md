# toolResolutionMapping

See description of [`calcResolutionMapping`](calcResolutionMapping.md).
Here we are assuming target resolution is finer than what mapping
already provides.

## Usage

``` r
toolResolutionMapping(mapping, targetGrid)
```

## Arguments

- mapping:

  a data.frame with columns x, y, lowRes

- targetGrid:

  a terra SpatRaster with the target resolution

## Value

a data.frame with columns x, y, lowRes, country

## Author

Pascal Sauer
