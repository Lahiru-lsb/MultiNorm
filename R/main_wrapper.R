#' Normalize multislice spatial transcriptomics (ST) data
#'
#' Performs normalization of multislice ST data using a negative binomial model.
#' The method aims to reduce unwanted technical variation while preserving biologically
#' meaningful spatial expression patterns and variations within and across tissue slices.
#'
#' @param spe_list List of \code{SpatialExperiment} objects, one per slice.
#'   Each object must contain a column with spot-level niche annotations.
#'
#' @param niche_col Character or integer. Column name or index in \code{colData}
#'   specifying niche labels. Niches represent biologically similar clusters
#'   shared across slices.
#'
#' @param preferred_niches Optional vector specifying a subset of niche labels.
#'   Default is \code{NULL}, in which case all shared niches are used.
#'
#' @param spot_frac Numeric in \code{(0, 1]}. Proportion of spots to sample for
#'   model fitting. Default is \code{1}.
#'
#' @param verbose Logical. If \code{TRUE}, prints sampling details.
#'   Default is \code{FALSE}.
#'
#' @param df_bio Numeric vector of length 2 specifying degrees of freedom for
#'   the two-dimensional spline basis used to model biological spatial effects.
#'   Default is \code{c(6, 6)}.
#'
#' @param df_tech Numeric vector of length 2 specifying degrees of freedom for
#'   the two-dimensional spline basis used to model library-size-associated
#'   technical effects. Default is \code{c(3, 3)}.
#'
#' @param lambdaC Numeric. Penalty coefficient enforcing similarity of biological
#'   effects across matched niches and slices. Default is \code{10000}.
#'
#' @param lambda1 Numeric. Ridge regularization parameter for biological spline
#'   coefficients. Default is \code{0.01}.
#'
#' @param lambda2 Numeric. Ridge regularization parameter for technical spline
#'   coefficients. Default is \code{2}.
#'
#' @param max_outer Integer. Maximum number of outer iterations for dispersion
#'   updates. Default is \code{20}.
#'
#' @param max_inner Integer. Maximum number of inner IRWLS iterations for model
#'   parameter updates. Default is \code{100}.
#'
#' @param tol_disp Numeric. Convergence tolerance for dispersion updates.
#'   Default is \code{1e-4}.
#'
#' @param tol_param Numeric. Convergence tolerance for inner-loop parameter
#'   updates. Default is \code{1e-5}.
#'
#' @param psi_damp Numeric in \code{(0, 1]}. Damping factor for dispersion
#'   updates. Default is \code{1}.
#'
#' @param param_damp Numeric in \code{(0, 1]}. Damping factor for IRWLS
#'   parameter updates. Default is \code{0.2}.
#'
#' @param seed Integer. Random seed used for spot subsampling.
#'   Default is \code{1}.
#'
#' @param alpha_init Optional numeric vector of initial values for slice-specific
#'   global library-size effects. Should have length equal to the number of
#'   slices. Default is \code{NULL}.
#'
#' @param Beta_init Optional list of initial biological spline coefficient
#'   matrices. Each element should have dimension \code{Kb x G}, where
#'   \code{Kb} is the number of biological spline basis functions and \code{G}
#'   is the number of genes. Default is \code{NULL}.
#'
#' @param Gamma_init Optional list of initial technical spline coefficient
#'   matrices. Each element should have dimension \code{Kt x G}, where
#'   \code{Kt} is the number of technical spline basis functions.
#'   Default is \code{NULL}.
#'
#' @param Psi_init Optional list of initial dispersion vectors. Each element
#'   should have length \code{G}. Default is \code{NULL}.
#'
#' @param zeta_init Optional numeric vector of initial values for global
#'   gene-specific intercepts. Should have length \code{G}. Default is
#'   \code{NULL}.
#'
#' @param rho_init Optional matrix of initial values for slice-deviation contrast
#'   parameters. Should have dimension \code{(S - 1) x G}, where \code{S} is the
#'   number of slices and \code{G} is the number of genes. Default is
#'   \code{NULL}.
#'
#' @param A Optional sum-to-zero contrast matrix for slice-specific technical
#'   deviations. Should have dimension \code{S x (S - 1)} and satisfy
#'   \code{colSums(A) == 0}. If \code{NULL}, a default contrast matrix is
#'   constructed internally.
#'
#' @param return_fit Logical. If \code{TRUE}, returns fitted model parameters
#'   along with normalized data. Default is \code{FALSE}.
#'
#' @return If \code{return_fit = FALSE}, returns a list of
#'   \code{SpatialExperiment} objects with normalized values stored by
#'   \code{norm_spe}. If \code{return_fit = TRUE}, returns a list containing:
#'   \itemize{
#'     \item \code{norm_spe}: normalized \code{SpatialExperiment} objects.
#'     \item \code{fit}: fitted model parameters.
#'     \item \code{call}: matched function call.
#'     \item \code{Eta_full_list}: fitted full linear predictors by slice.
#'     \item \code{percentiles}: fitted NB percentiles used for PAC transformation.
#'   }
#'
#' @details
#' The method fits a negative binomial model with spatial varying biological effects,
#' library size-associated technical effects, and slice-specific technical
#' effect. Biological and technical effects are modeled using
#' two-dimensional spline bases. Cross-slice biological consistency is encouraged
#' through matched niche-level constraints applied only to the biological spline
#' component.
#'
#'
#' Parameter estimation is performed using iteratively reweighted least squares
#' (IRWLS), with outer-loop updates for negative binomial dispersion parameters.
#'
#' @examples
#' \dontrun{
#' # Basic usage
#' result <- MultiNorm(
#'   spe_list,
#'   niche_col = "niche"
#' )
#'
#' # Access normalized output
#' norm_spe <- result
#'
#' # Subsample spots for faster fitting
#' result <- MultiNorm(
#'   spe_list,
#'   niche_col = "niche",
#'   spot_frac = 0.5,
#'   seed = 123,
#'   verbose = TRUE
#' )
#'
#' # Tune spline flexibility and regularization
#' result <- MultiNorm(
#'   spe_list,
#'   niche_col = "niche",
#'   df_bio = c(8, 8),
#'   df_tech = c(4, 4),
#'   lambda1 = 0.05,
#'   lambda2 = 2.0
#' )
#'
#' # Use a subset of niches
#' result <- MultiNorm(
#'   spe_list,
#'   niche_col = "niche",
#'   preferred_niches = c("A", "B", "C")
#' )
#'
#' # Return fitted model parameters
#' result <- MultiNorm(
#'   spe_list,
#'   niche_col = "niche",
#'   return_fit = TRUE
#' )
#'
#' norm_spe <- result$norm_spe
#' fit <- result$fit
#' }
#'
#' @export

MultiNorm <- function(
    spe_list,
    niche_col,
    preferred_niches = NULL,
    spot_frac      = 1,
    verbose        = FALSE,
    df_bio         = c(6, 6),
    df_tech        = c(3, 3),
    lambdaC        = 10000,
    lambda1        = 0.01,
    lambda2        = 2,
    max_outer      = 20,
    max_inner      = 100,
    tol_disp       = 1e-4,
    tol_param      = 1e-5,
    psi_damp       = 1,
    param_damp     = 0.2,
    seed           = 1,
    alpha_init     = NULL,
    Beta_init      = NULL,
    Gamma_init     = NULL,
    Psi_init       = NULL,
    zeta_init      = NULL,
    rho_init       = NULL,
    A              = NULL,
    return_fit     = FALSE
) {

  # Step 1: Fit
  message("Loading v2")

  model_fit <- fit_spatial_nb_multislice(
    spe_list   = spe_list,
    niche_col  = niche_col,
    preferred_niches = preferred_niches,
    spot_frac  = spot_frac,
    verbose   = verbose,
    df_bio     = df_bio,
    df_tech    = df_tech,
    lambdaC    = lambdaC,
    lambda1    = lambda1,
    lambda2    = lambda2,
    max_outer  = max_outer,
    max_inner  = max_inner,
    tol_disp   = tol_disp,
    tol_param  = tol_param,
    psi_damp   = psi_damp,
    param_damp = param_damp,
    seed       = seed,
    alpha_init = alpha_init,
    Beta_init  = Beta_init,
    Gamma_init = Gamma_init,
    Psi_init   = Psi_init,
    zeta_init  = zeta_init,
    rho_init   = rho_init,
    A          = A
  )

  # Step 2: Predict
  out <- predict_logPAC_multislice(
    fit      = model_fit$fit,
    df_bio   = df_bio,
    df_tech  = df_tech,
    spe_list = spe_list
  )

  run_config <- list(
    call = match.call(),
    params = list(
      niche_col = niche_col,
      preferred_niches = preferred_niches,
      spot_frac = spot_frac,
      df_bio = df_bio,
      df_tech = df_tech,
      lambdaC = lambdaC,
      lambda1 = lambda1,
      lambda2 = lambda2,
      max_outer = max_outer,
      max_inner = max_inner,
      tol_disp = tol_disp,
      tol_param = tol_param,
      psi_damp = psi_damp,
      param_damp = param_damp,
      seed = seed
    ),
    timestamp = Sys.time()
  )

  run_config$session <- utils::sessionInfo()

  # Step 3: restore original coords
  spe_list <- lapply(spe_list, restore_coords)

  # Step 4: attach normalized data back to SPE
  spe_list_norm <- norm_spe(spe_list, out$out_list)

  for (i in seq_along(spe_list_norm)) {
    metadata(spe_list_norm[[i]])$MultiNorm <- run_config
  }

  if (return_fit) {
    fit <- model_fit$fit

    fit$U_all <- out$U_all
    fit$V_all <- out$V_all
    fit$bio   <- out$bio
    fit$tech  <- out$tech

    return(list(
      norm_spe = spe_list_norm,
      fit = fit,
      call = model_fit$call,
      Eta_full_list = out$Eta_full_list,
      percentiles = out$percentiles
    ))
  }

  spe_list_norm
}
