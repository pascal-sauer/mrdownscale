# calcLandInputRecategorized

Computes the land input data in target land categories. Splitting of
land categories is performed under use of internal land weights
reflecting the prevalence of a certain land category in the given area.

## Usage

``` r
calcLandInputRecategorized(input, target)
```

## Arguments

- input:

  name of an input dataset, options: "magpie", "witch"

- target:

  name of the land target source to be used

## Details

Mapping from input to target categories is achieved via a merge of a
land input mapping to reference categories and a mapping between land
target categories and the same reference categories. Thereby a new
source or new target can be supported by supplying a map of that new
input and/or target to the reference categories.

input = "witch": The "rest" category added in calcLandInput is
disaggregated into all missing land variables using the reference
dataset included in mrdownscale. Thus the following variables are added,
but were not at all part of the input scenario data: range, urban,
c4ann\_\*, c3per\_\*, c3nfx\_\* This allows the harmonization and
downscaling pipeline to continue with a complete dataset, but these
variables should not be used or reported (at the very least they have to
be checked extensively).

## Author

Jan Philipp Dietrich, Pascal Sauer
