# toolGetHarmonizer

Get a harmonizer function by name.

## Usage

``` r
toolGetHarmonizer(harmonizerName, constantTotal)
```

## Arguments

- harmonizerName:

  name of a harmonizer function, currently offset, fade, fadeForest,
  absoluteChanges

- constantTotal:

  TRUE for area conservative land data, FALSE otherwise, only used by
  the absoluteChanges harmonizer, ignored by the others

## Value

harmonizer function

## See also

[`toolHarmonizeOffset`](toolHarmonizeOffset.md),
[`toolHarmonizeFade`](toolHarmonizeFade.md),
[`toolHarmonizeFadeForest`](toolHarmonizeFadeForest.md),
[`toolHarmonizeAbsoluteChanges`](toolHarmonizeAbsoluteChanges.md)

## Author

Pascal Sauer
