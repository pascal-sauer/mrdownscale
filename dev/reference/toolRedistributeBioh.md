# toolRedistributeBioh

Move the bioh of clusters without wood harvest area to the clusters that
still have some. Their recalculated bioh (harvest area times harvest
rate) is zero, so rates times land cannot equal target bioh and would be
stranded. Per subcategory and year, the stranded bioh is added to
clusters with recalculated bioh \> 0 in the same region proportional to
their recalculated bioh, falling back to the global pool if a region has
no such clusters. Bioh that cannot be moved anywhere is reported and
dropped, so the returned bioh is zero everywhere where the recalculated
bioh is zero.

## Usage

``` r
toolRedistributeBioh(bioh, biohCalculated)
```

## Arguments

- bioh:

  bioh in kg C yr-1

- biohCalculated:

  bioh in in kg C yr-1 recalculated from harvest rate times harvest area

## Value

bioh in kg C yr-1 with the stranded bioh moved to the receivers, zero on
clusters without harvest area, totals preserved unless bioh had to be
dropped

## Author

Pascal Sauer
