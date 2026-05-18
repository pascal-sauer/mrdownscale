# toolGetSmallerArea

The function creates a transition data set based on a state dataset
which contains for every possible connection (land type A \<-\> land
type B) the area of the smaller of these land types. This information is
relevant if transitional shares should be calculated based on the
smaller area of the two (assuming that the smaller area is determining
how much transitions there will be).

## Usage

``` r
toolGetSmallerArea(states)
```

## Arguments

- states:

  magpie dataset containing states information

## Value

A land type x land type data set containing for each possible
combination the smaller area of the two land types

## Author

Jan Philipp Dietrich
