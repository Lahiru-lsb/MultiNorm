# Normalize multislice spatial transcriptomics (ST) data

Performs normalization of multislice ST data using a negative binomial
model. The method aims to reduce unwanted technical variation while
preserving biologically meaningful spatial expression patterns and
variations within and across tissue slices.

## Usage

``` r
MultiNorm(
  spe_list,
  niche_col,
  preferred_niches = NULL,
  spot_frac = 1,
  verbose = FALSE,
  df_bio = c(6, 6),
  df_tech = c(3, 3),
  lambdaC = 10000,
  lambda1 = 0.01,
  lambda2 = 2,
  max_outer = 20,
  max_inner = 100,
  tol_disp = 1e-04,
  tol_param = 1e-05,
  psi_damp = 1,
  param_damp = 0.2,
  seed = 1,
  alpha_init = NULL,
  Beta_init = NULL,
  Gamma_init = NULL,
  Psi_init = NULL,
  zeta_init = NULL,
  rho_init = NULL,
  A = NULL,
  return_fit = FALSE
)
```

## Arguments

- spe_list:

  List of `SpatialExperiment` objects, one per slice. Each object must
  contain a column with spot-level niche annotations.

- niche_col:

  Character or integer. Column name or index in `colData` specifying
  niche labels. Niches represent biologically similar clusters shared
  across slices.

- preferred_niches:

  Optional vector specifying a subset of niche labels. Default is
  `NULL`, in which case all shared niches are used.

- spot_frac:

  Numeric in `(0, 1]`. Proportion of spots to sample for model fitting.
  Default is `1`.

- verbose:

  Logical. If `TRUE`, prints sampling details. Default is `FALSE`.

- df_bio:

  Numeric vector of length 2 specifying degrees of freedom for the
  two-dimensional spline basis used to model biological spatial effects.
  Default is `c(6, 6)`.

- df_tech:

  Numeric vector of length 2 specifying degrees of freedom for the
  two-dimensional spline basis used to model library-size-associated
  technical effects. Default is `c(3, 3)`.

- lambdaC:

  Numeric. Penalty coefficient enforcing similarity of biological
  effects across matched niches and slices. Default is `10000`.

- lambda1:

  Numeric. Ridge regularization parameter for biological spline
  coefficients. Default is `0.01`.

- lambda2:

  Numeric. Ridge regularization parameter for technical spline
  coefficients. Default is `2`.

- max_outer:

  Integer. Maximum number of outer iterations for dispersion updates.
  Default is `20`.

- max_inner:

  Integer. Maximum number of inner IRWLS iterations for model parameter
  updates. Default is `100`.

- tol_disp:

  Numeric. Convergence tolerance for dispersion updates. Default is
  `1e-4`.

- tol_param:

  Numeric. Convergence tolerance for inner-loop parameter updates.
  Default is `1e-5`.

- psi_damp:

  Numeric in `(0, 1]`. Damping factor for dispersion updates. Default is
  `1`.

- param_damp:

  Numeric in `(0, 1]`. Damping factor for IRWLS parameter updates.
  Default is `0.2`.

- seed:

  Integer. Random seed used for spot subsampling. Default is `1`.

- alpha_init:

  Optional numeric vector of initial values for slice-specific global
  library-size effects. Should have length equal to the number of
  slices. Default is `NULL`.

- Beta_init:

  Optional list of initial biological spline coefficient matrices. Each
  element should have dimension `Kb x G`, where `Kb` is the number of
  biological spline basis functions and `G` is the number of genes.
  Default is `NULL`.

- Gamma_init:

  Optional list of initial technical spline coefficient matrices. Each
  element should have dimension `Kt x G`, where `Kt` is the number of
  technical spline basis functions. Default is `NULL`.

- Psi_init:

  Optional list of initial dispersion vectors. Each element should have
  length `G`. Default is `NULL`.

- zeta_init:

  Optional numeric vector of initial values for global gene-specific
  intercepts. Should have length `G`. Default is `NULL`.

- rho_init:

  Optional matrix of initial values for slice-deviation contrast
  parameters. Should have dimension `(S - 1) x G`, where `S` is the
  number of slices and `G` is the number of genes. Default is `NULL`.

- A:

  Optional sum-to-zero contrast matrix for slice-specific technical
  deviations. Should have dimension `S x (S - 1)` and satisfy
  `colSums(A) == 0`. If `NULL`, a default contrast matrix is constructed
  internally.

- return_fit:

  Logical. If `TRUE`, returns fitted model parameters along with
  normalized data. Default is `FALSE`.

## Value

If `return_fit = FALSE`, returns a list of `SpatialExperiment` objects
with normalized values stored by `norm_spe`. If `return_fit = TRUE`,
returns a list containing:

- `norm_spe`: normalized `SpatialExperiment` objects.

- `fit`: fitted model parameters.

- `call`: matched function call.

- `Eta_full_list`: fitted full linear predictors by slice.

- `percentiles`: fitted NB percentiles used for PAC transformation.

## Details

The method fits a negative binomial model with spatial biological
effects, library-size-associated technical effects, and slice-specific
technical intercept deviations. Biological and technical effects are
modeled using two-dimensional spline bases. Cross-slice biological
consistency is encouraged through matched niche-level constraints
applied only to the biological spline component.

The model includes a global gene-specific intercept and a slice-specific
technical intercept deviation represented using sum-to-zero contrasts.
This allows slice-level technical shifts not explained by library size
to be modeled explicitly, while preserving the biological niche-matching
constraint on the biological component.

Parameter estimation is performed using iteratively reweighted least
squares (IRWLS), with outer-loop updates for negative binomial
dispersion parameters.

## Examples

``` r
if (FALSE) { # \dontrun{
# Basic usage
result <- MultiNorm(
  spe_list,
  niche_col = "niche"
)

# Access normalized output
norm_spe <- result

# Subsample spots for faster fitting
result <- MultiNorm(
  spe_list,
  niche_col = "niche",
  spot_frac = 0.5,
  seed = 123,
  verbose = TRUE
)

# Tune spline flexibility and regularization
result <- MultiNorm(
  spe_list,
  niche_col = "niche",
  df_bio = c(8, 8),
  df_tech = c(4, 4),
  lambda1 = 0.05,
  lambda2 = 2.0
)

# Use a subset of niches
result <- MultiNorm(
  spe_list,
  niche_col = "niche",
  preferred_niches = c("A", "B", "C")
)

# Return fitted model parameters
result <- MultiNorm(
  spe_list,
  niche_col = "niche",
  return_fit = TRUE
)

norm_spe <- result$norm_spe
fit <- result$fit
} # }
```
