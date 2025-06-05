# Downscale and harmonize land use data using high resolution
    reference data

R package **mrdownscale**, version **0.37.1**

[![CRAN status](https://www.r-pkg.org/badges/version/mrdownscale)](https://cran.r-project.org/package=mrdownscale) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.11244475.svg)](https://doi.org/10.5281/zenodo.11244475) [![R build status](https://github.com/pik-piam/mrdownscale/workflows/check/badge.svg)](https://github.com/pik-piam/mrdownscale/actions) [![codecov](https://codecov.io/gh/pik-piam/mrdownscale/branch/master/graph/badge.svg)](https://app.codecov.io/gh/pik-piam/mrdownscale) [![r-universe](https://pik-piam.r-universe.dev/badges/mrdownscale)](https://pik-piam.r-universe.dev/builds)

## Purpose and Functionality

Downscale and harmonize land use data (e.g. MAgPIE) using
    high resolution reference data (e.g. LUH2v2h).


## Installation

For installation of the most recent package version an additional repository has to be added in R:

```r
options(repos = c(CRAN = "@CRAN@", pik = "https://rse.pik-potsdam.de/r/packages"))
```
The additional repository can be made available permanently by adding the line above to a file called `.Rprofile` stored in the home folder of your system (`Sys.glob("~")` in R returns the home directory).

After that the most recent version of the package can be installed using `install.packages`:

```r 
install.packages("mrdownscale")
```

Package updates can be installed using `update.packages` (make sure that the additional repository has been added before running that command):

```r 
update.packages()
```

## Tutorial

The package comes with a vignette describing the basic functionality of the package and how to use it. You can load it with the following command (the package needs to be installed):

```r
vignette("basicUsage") # Basic Usage of mrdownscale
```

## Questions / Problems

In case of questions / problems please contact Pascal Sauer <pascal.sauer@pik-potsdam.de>.

## Citation

To cite package **mrdownscale** in publications use:

Sauer P, Dietrich J (2025). "mrdownscale: Downscale and harmonize land use data using high resolution reference data." doi:10.5281/zenodo.11244475 <https://doi.org/10.5281/zenodo.11244475>, Version: 0.37.1, <https://github.com/pik-piam/mrdownscale>.

A BibTeX entry for LaTeX users is

 ```latex
@Misc{,
  title = {mrdownscale: Downscale and harmonize land use data using high resolution
    reference data},
  author = {Pascal Sauer and Jan Philipp Dietrich},
  doi = {10.5281/zenodo.11244475},
  date = {2025-06-05},
  year = {2025},
  url = {https://github.com/pik-piam/mrdownscale},
  note = {Version: 0.37.1},
}
```
