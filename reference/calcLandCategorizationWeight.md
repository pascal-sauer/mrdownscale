# calcLandCategorizationWeight

Computes weights for a weighted category mapping and aggregates the
weights to the spatial resolution described by the arguments `geometry`
and `crs`.

## Usage

``` r
calcLandCategorizationWeight(map, geometry, crs)
```

## Arguments

- map:

  a map in form of a data.frame containing a mapping between reference
  categories (column name "reference"), input categories (column name
  "dataInput") and merged categories (column name "merge") for the data
  dimension of a magpie object

- geometry:

  the geometry of the magpie object for which the categories should be
  mapped as given in the geometry attribute of a magpie object with
  geometry information `attr(x, "geometry")`.

- crs:

  the coordinate reference system as returned by `attr(x, "crs")` from a
  magpie object with coordinates information.

## Note

This calc-function has a rather unusual shape in that sense that the
arguments to be provided are not simple configuration settings but
rather relative complex. For this kind of implementation it is typically
advised to use tool functions. The reason that a calc function is used
in this particular case instead is, that the inputs of this function
usually do not change over a long period of time so that caching of
results becomes key for the overall performance of the data processing,
which is available for calc- but not tool-functions. If tool-functions
might support caching in the future as well as a conversion to a tool
function might become a sensible option.

## Author

Jan Philipp Dietrich
