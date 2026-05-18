# calcLandInput

Prepare the land input data for the category mapping, checking data for
consistency before returning. All "Land" functions deal with area data,
as opposed to "Nonland" functions which deal with non-area data such as
the amount of applied fertilizer. These are treated differently, because
for area data other constraints apply, e.g. the total area must be
constant over time.

## Usage

``` r
calcLandInput(input)
```

## Arguments

- input:

  name of an input dataset, options: "magpie", "witch"

## Value

land input data

## Details

input = "magpie": includes the land use categories past (pasture,
including rangeland), forestry (managed forest plantations), primforest,
secdforest, urban, other (other land) and many specific crop types.
Furthermore, 1st gen biofuel is added and filled with zeros. 1st gen
biofuel is only modeled implicitly in magpie via demand, and because of
trade it is unclear on what area 1st gen biofuel is grown, also 1st gen
biofuel is quickly phased out in magpie, so we fill biofuel_1st_gen with
zeros and rely on the harmonization to produce a plausible 1st gen
biofuel time series. input = "witch": includes a subset of LUH land use
categories: primf, primn, secdn, pastr, c4per, pltns, secdf,
c3ann_irrigated, c3ann_rainfed These are given as shares. A "rest"
category is added so shares sum up to 1.

## Author

Jan Philipp Dietrich, Pascal Sauer
