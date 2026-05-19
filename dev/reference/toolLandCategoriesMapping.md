# toolLandCategoriesMapping

Computes a land category mapping between input and target data set via
retrieving mappings from these sources to an internally defined
reference categorization and merging them to a direct input-to-target
mapping. Mappings to reference categories are stored internally in the
package and have to be added to it if new input and/or target data
should be supported.

## Usage

``` r
toolLandCategoriesMapping(input, target)
```

## Arguments

- input:

  name of an input dataset, options: "magpie", "witch"

- target:

  name of the land target source to be used

## Author

Jan Philipp Dietrich
