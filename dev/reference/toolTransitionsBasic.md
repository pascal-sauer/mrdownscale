# toolTransitionsBasic

tool function to extract transitions between categories from a land data
set with at least 2 time steps. The approach is rather simplistic by
assuming that expansion happens proportionally in all affected classes
(equal preference of transitions across all categories).

## Usage

``` r
toolTransitionsBasic(x, gross = FALSE)
```

## Arguments

- x:

  magpie data set containing land data with at least two time steps to
  extract net transitions from

- gross:

  either boolean or a magpie object containing bidirectional transition
  shares relative to the area of the involved land pools (transition
  divided by the area of the land pool in the "from" sub dimension). If
  set to FALSE only net transitions will be returned. If set to TRUE an
  internal gross transition estimate based on average gross transitions
  in LUH2 in the period from 1995 to 2015 will be used.

## Details

In addition to the net effect it can also estimate gross transition. For
that purpose a reference data set containing bidirectional transition
shares must be provided.

If the time step length is longer than 1 year the returned object
contains reference years for each period which can be repeated to
retrieve the full transition between two time periods, e.g. if you
provide two time steps 2000 and 2005 the return value will be the
transition for year 2001. Repeating the same transition also in 2002,
2003, 2004 and 2005 will give the full transition from 2000 to 2005.

## Author

Jan Philipp Dietrich
