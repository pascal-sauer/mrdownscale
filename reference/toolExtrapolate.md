# toolExtrapolate

Extrapolate a dataset into the future. A linear model is fitted for each
combination of spatial entity and category, and then used to predict the
value of the last requested extrapolation year. If the linear model is
not significant (p \> 0.05) the historical mean is used instead. A
spline-based interpolation is used to create a smooth transition from
historical values to the predicted value.

## Usage

``` r
toolExtrapolate(x, years)
```

## Arguments

- x:

  A magpie object with "year" as the temporal dimension and without any
  NAs

- years:

  A vector of years to extrapolate to

## Value

A magpie object like x but with the extrapolated years only

## Author

Pascal Sauer
