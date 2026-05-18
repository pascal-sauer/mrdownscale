# readLUH2v2h

Read LUH2v2h data. For the states subtype, the secma and secmb
categories are removed. For the management subtype, only the categories
crpbf, rndwd, fulwd, fertl and irrig are read. For the transitions
subtype, only the wood harvest categories bioh and harv are read. To
match magpie semantics years are shifted by 1 when reading transitions.

## Usage

``` r
readLUH2v2h(subtype)
```

## Arguments

- subtype:

  one of states, management, transitions, cellArea
