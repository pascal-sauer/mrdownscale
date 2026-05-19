# toolReplaceExpansion

Expansion from one timestep to the next of one land type is replaced
with another land type.

## Usage

``` r
toolReplaceExpansion(
  x,
  from,
  to,
  ...,
  noteThreshold = 10^-10,
  warnThreshold = 10^-5,
  level = 1
)
```

## Arguments

- x:

  a magpie object

- from:

  name of a land category, e.g. "primf" Expansion of this category
  happening from one timestep to the next will be replaced. If missing
  in x return x unchanged

- to:

  name of another land category, e.g. "secdf" Expansion of 'from' will
  be replaced with expansion of 'to'

- ...:

  not used, will throw an error if supplied

- noteThreshold:

  expansion greater than this will trigger a note that expansion was
  replaced

- warnThreshold:

  expansion greater than this will trigger a warning that expanding
  considerably

- level:

  passed to toolStatusMessage

## Value

a magpie object with expansion of 'from' replaced by 'to'

## Author

Pascal Sauer
