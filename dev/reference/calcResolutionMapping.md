# calcResolutionMapping

Calculate a complete mapping from low (input dataset,
clusters/countries/regions) to high resolution (target dataset, grid).
As a basis the mapping from the low resolution
clusters/countries/regions to grid cells is used. Cells which are
present in that mapping, but not in the target dataset are discarded.
Cells which are present in the target dataset, but not in the mapping
are added using a nearest-neighbor approach: These cells are mapped to
the same low resolution cluster/country/region as the closest cell which
is already present in the mapping.

## Usage

``` r
calcResolutionMapping(input, target)
```

## Arguments

- input:

  name of an input dataset, options: "magpie", "witch"

- target:

  character, the target dataset

## Value

a list including a data.frame with columns x, y, lowRes, country

## Author

Pascal Sauer
