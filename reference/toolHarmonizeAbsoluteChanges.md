# toolHarmonizeAbsoluteChanges

Tool function for creating a harmonized data set by applying the
absolute changes of the input data to the target data: up to the
harmonization year the target data is used, afterwards the difference of
the input data to the input data of the harmonization year is added to
the target data of the harmonization year.

## Usage

``` r
toolHarmonizeAbsoluteChanges(xInput, xTarget, harmonizationPeriod)
```

## Arguments

- xInput:

  input data as magpie object

- xTarget:

  target data as magpie object

- harmonizationPeriod:

  Two identical integer values, the year the absolute changes of the
  input data are applied to, must be present in both input and target
  data

## Value

harmonized data set as magpie object with data from target for years up
to and including the harmonization year and absolute changes from input
relative to the harmonization year afterwards

## Details

Negative values are handled within groups of related categories (forest
& other land; all cropland types; pasture & rangeland), so that they are
compensated by the other categories of the group. For crops, negative
values are first compensated by the corresponding rainfed/irrigated twin
of the same crop, which keeps that pair's total area, including biofuel
types like c3ann_rainfed_biofuel_1st_gen. Crop categories without their
twin are compensated by the whole group. Groups with a negative total
area are set to zero. Afterwards, if any negatives remain, all
categories except urban are scaled down to keep the total area constant.

## Author

Pascal Sauer
