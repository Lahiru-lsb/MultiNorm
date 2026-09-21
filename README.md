# MultiNorm

**MultiNorm** is an R package for normalization of **multi-slice spatial transcriptomics (ST) data**.

The method decomposes observed gene-expression variation into biological and technical components. MultiNorm is designed to reduce unwanted technical variation while preserving biological spatial structure across slices.

It uses pre-identified biological regions, referred to as **niches**, to encourage consistency of average biological signals across matched regions in different slices.

## Overview

MultiNorm fits a negative binomial model that includes:

-   spatial biological effects;
-   library-size-associated technical effects;
-   slice-specific technical intercept deviations; and
-   cross-slice constraints on matched biological niches.

Biological and technical spatial effects are represented using two-dimensional spline bases. The cross-slice constraint is applied to the biological component so that shared biological structure is preserved during normalization.

## Installation

You can install the development version of MultiNorm directly from GitHub:

``` r
# Install remotes if needed
install.packages("remotes")

# Install MultiNorm
remotes::install_github("Lahiru-lsb/MultiNorm")
```

Then load the package:

``` r
library(MultiNorm)
```

## Input

The main function, `MultiNorm()`, expects a list of `SpatialExperiment` objects:

``` r
spe_list
```

Each element of the list should represent one spatial transcriptomics slice.

Each `SpatialExperiment` object must contain a spot-level niche annotation in `colData`. Niches should represent biologically comparable regions that can be matched across slices.

For example, if the niche labels are stored in:

``` r
colData(spe)$niche
```

then use:

``` r
niche_col = "niche"
```

## Quick start

A basic MultiNorm analysis can be run as follows:

``` r
library(MultiNorm)

result <- MultiNorm(
  spe_list = spe_list,
  niche_col = "niche"
)
```

By default, the returned object is a list of normalized `SpatialExperiment` objects:

``` r
norm_spe <- result
```

## Using selected niches

By default, MultiNorm uses all shared niches available across slices.

To use only selected niches:

``` r
result <- MultiNorm(
  spe_list = spe_list,
  niche_col = "niche",
  preferred_niches = c("A", "B", "C")
)
```

## Subsampling spots

For large datasets, a proportion of spots can be used for model fitting using `spot_frac`.

For example, to fit the model using 50% of spots:

``` r
result <- MultiNorm(
  spe_list = spe_list,
  niche_col = "niche",
  spot_frac = 0.5,
  seed = 123,
  verbose = TRUE
)
```

`spot_frac` must be between 0 and 1.

## Controlling spatial flexibility

The flexibility of the biological and technical spatial effects can be controlled using `df_bio` and `df_tech`.

``` r
result <- MultiNorm(
  spe_list = spe_list,
  niche_col = "niche",
  df_bio = c(8, 8),
  df_tech = c(4, 4)
)
```

The default values are:

``` r
df_bio  = c(6, 6)
df_tech = c(3, 3)
```

## Regularization

MultiNorm includes regularization parameters for the biological and technical spline components.

``` r
result <- MultiNorm(
  spe_list = spe_list,
  niche_col = "niche",
  lambdaC = 10000,
  lambda1 = 0.01,
  lambda2 = 2
)
```

The main parameters are:

-   `lambdaC`: strength of the cross-slice biological consistency constraint;
-   `lambda1`: ridge regularization for biological spline coefficients;
-   `lambda2`: ridge regularization for technical spline coefficients.

For most analyses, the default values can be used as a starting point.

## Returning fitted model parameters

To return fitted model parameters together with the normalized data:

``` r
result <- MultiNorm(
  spe_list = spe_list,
  niche_col = "niche",
  return_fit = TRUE
)

norm_spe <- result$norm_spe
fit <- result$fit
```

When `return_fit = TRUE`, the returned object also contains:

-   `norm_spe`: normalized `SpatialExperiment` objects;
-   `fit`: fitted model parameters;
-   `call`: matched function call;
-   `Eta_full_list`: fitted full linear predictors for each slice;
-   `percentiles`: fitted negative-binomial percentiles used for the PAC transformation.

## Main function

The main normalization function is:

``` r
MultiNorm()
```

To view its complete documentation in R:

``` r
?MultiNorm
```

or:

``` r
help("MultiNorm", package = "MultiNorm")
```

The full function includes additional options for convergence control, damping, initialization, contrast specification, and model fitting. These advanced arguments are documented in the R help page.

## Example workflow

A typical workflow is:

``` r
library(MultiNorm)

# 1. Prepare a list of SpatialExperiment objects
spe_list <- list(slice1, slice2, slice3)

# 2. Run MultiNorm
norm_spe <- MultiNorm(
  spe_list = spe_list,
  niche_col = "niche"
)

# 3. Access a normalized slice
norm_slice1 <- norm_spe[[1]]
```

The normalized `SpatialExperiment` objects can then be used in downstream analyses such as spatial visualization, spatially variable gene detection, clustering, or other multi-slice analyses.

## Method summary

For each gene, MultiNorm models expression using a negative binomial framework that separates biological and technical sources of variation.

Conceptually, the method decomposes expression into:

``` text
Observed expression
        |
        +-- Biological spatial smooth variation
        |
        +-- Library-size-associated smooth technical variation
        |
        +-- Slice-specific technical variation
```

Cross-slice constraints are applied to matched niches to preserve shared biological structure while normalization targets unwanted technical effects.

## Repository

Source code and development versions of MultiNorm are available at:

<https://github.com/Lahiru-lsb/MultiNorm>

## Issues

If you encounter a bug or have a feature request, please open an issue in the GitHub repository.

## 
