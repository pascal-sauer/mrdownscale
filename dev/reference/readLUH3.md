# readLUH3

Read LUH3 data. For the states subtype, the secma and secmb categories
are removed. For the management subtype, only the categories cpbf1,
cpbf2, rndwd, fulwd, fertl and irrig are read. For the transitions
subtype, only the wood harvest categories bioh and harv are read. To
match magpie semantics years are shifted by 1 when reading transitions.
The LUH3 nc files have day-based time, which is converted to years.

## Usage

``` r
readLUH3(subtype, subset)
```

## Arguments

- subtype:

  one of states, management, transitions, cellArea

- subset:

  which years to read

## Value

data read from LUH3 historic nc files as SpatRaster

## Author

Pascal Sauer
