
#' Normalize multislice spatial transcriptomics (ST) data
#'
#' Performs normalization of multislice ST data using
#' spatially dependent, spot- and gene-specific size factors.
#'
#' @param spe_list List of \code{SpatialExperiment} objects, one per slice.
#'   Each object must contain a column with spot-level niche annotations.
#'
#' @param niche_col Character or integer. Column name or index in \code{colData}
#'   specifying niche labels. Niches represent biologically similar clusters
#'   shared across slices.
#'
#' @param preferred_niches Optional vector specifying a subset of niche labels.
#'   Default is \code{NULL} (all niches are used).
#'
#' @param spot_frac Numeric in (0,1]. Proportion of spots to sample for
#'   model fitting. Default is \code{1}.
#'
#' @param verbose Logical. If \code{TRUE}, prints sampling details.
#'   Default is \code{FALSE}.
#'
#' @param df_bio Numeric vector of length 2 specifying degrees of freedom
#'   for the thin-plate spline modeling biological spatial effects.
#'   Default is \code{c(6, 6)}.
#'
#' @param df_tech Numeric vector of length 2 specifying degrees of freedom
#'   for the thin-plate spline modeling technical (library size) effects.
#'   Default is \code{c(3, 3)}.
#'
#' @param lambdaC Numeric. Penalty coefficient enforcing similarity of
#'   biological effects across slices at the niche level.
#'   Default is \code{1e4}.
#'
#' @param lambda1 Numeric. Regularization parameter for gene-specific biological coefficients.
#'   Default is \code{0.001}.
#'
#' @param lambda2 Numeric. Regularization parameter for gene-specific technical
#'   (library size) coefficients.
#'   Default is \code{1.9}.
#'
#' @param max_outer Integer. Maximum number of outer iterations.
#'   Default is \code{30}.
#'
#' @param max_inner Integer. Maximum number of inner IRWLS iterations.
#'   Default is \code{200}.
#'
#' @param tol_disp Numeric. Convergence tolerance for dispersion updates.
#'   Default is \code{1e-4}.
#'
#' @param tol_param Numeric. Convergence tolerance for inner-loop parameter updates.
#'   Default is \code{1e-5}.
#'
#' @param psi_damp Numeric in (0,1]. Damping factor (specifying the multiplicative factor to decrease IRLS step by when log-likelihood diverges) for dispersion updates to
#'   stabilize optimization. Default is \code{0.5}.
#'
#' @param param_damp Numeric in (0,1]. Damping factor (specifying the multiplicative factor to decrease IRWLS step by when log-likelihood diverges) for IRWLS parameter updates
#'   (applies to all parameters except dispersion). Default is \code{0.2}.
#'
#' @param seed Integer. Random seed used for subsampling.
#'   Default is \code{1}.
#'
#' @param alpha_init Optional numeric vector of initial values for slice specific global library size effect.
#'   Default is \code{NULL}.
#'
#' @param Beta_init Optional list of initial values for biological coefficients.
#'   Default is \code{NULL}.
#'
#' @param Gamma_init Optional list of initial values for technical coefficients.
#'   Default is \code{NULL}.
#'
#' @param Psi_init Optional list of initial values for dispersion parameters.
#'   Default is \code{NULL}.
#' @param return_fit Logical. If TRUE, returns fitted model parameters along with normalized data.
#' Defaults is \code{FALSE}
#'
#' @return A list  of SpatialExperiment  containing normalized outputs stored in 'multinorm' and optionally fitted model parameters.
#'
#' @details
#' The method fits a negative binomial model with spatially varying effects.
#' Biological and technical effects are modeled using thin-plate splines, and
#' cross-slice consistency is enforced via niche-level constraints. Parameter
#' estimation is performed using iteratively reweighted least squares (IRWLS).
#'
#' @examples
#' \dontrun{
#' # ------------------------------------------------------------
#' # Expected input structure
#' # spe_list is a list of SpatialExperiment objects (one per slice)
#' # Each must contain a 'niche' column in colData
#' # ------------------------------------------------------------
#' # spe_list <- list(spe_slice1, spe_slice2)
#' # colData(spe_list[[1]])$niche <- factor(...)
#'
#' # ------------------------------------------------------------
#' # Basic usage: minimal inputs
#' # ------------------------------------------------------------
#' result <- MultiNorm(
#'   spe_list,
#'   niche_col = "niche"
#' )
#'
#' # Access normalized output
#' norm_spe <- result
#'
#' # ------------------------------------------------------------
#' # Adjust sampling and ensure reproducibility
#' # Subsample spots for faster fitting
#' # ------------------------------------------------------------
#' result <- MultiNorm(
#'   spe_list,
#'   niche_col = "niche",
#'   spot_frac = 0.5,
#'   seed = 123,
#'   verbose = TRUE
#' )
#'
#' # ------------------------------------------------------------
#' # Model tuning: spline flexibility and regularization
#' # ------------------------------------------------------------
#' result <- MultiNorm(
#'   spe_list,
#'   niche_col = "niche",
#'   df_bio = c(8, 8),
#'   df_tech = c(4, 4),
#'   lambda1 = 0.05,
#'   lambda2 = 2.0
#' )
#'
#' # ------------------------------------------------------------
#' # Use a subset of niches
#' # ------------------------------------------------------------
#' result <- MultiNorm(
#'   spe_list,
#'   niche_col = "niche",
#'   preferred_niches = c("A", "B", "C")
#' )
#'
#' # ------------------------------------------------------------
#' # Return fitted model parameters
#' # ------------------------------------------------------------
#' result <- MultiNorm(
#'   spe_list,
#'   niche_col = "niche",
#'   return_fit = TRUE
#' )
#'
#' # Access outputs when return_fit = TRUE
#' norm_spe <- result$norm_spe
#' fit <- result$fit
#' }
#' @export

MultiNorm <- function(
    spe_list,
    niche_col,
    preferred_niches = NULL,
    spot_frac      = 1,
    verbose        =FALSE,
    df_bio         = c(6, 6),
    df_tech        = c(3, 3),
    lambdaC        = 10000,
    lambda1        = 0.001,
    lambda2        = 1.9,
    max_outer      = 30,
    max_inner      = 200,
    tol_disp       = 1e-4,
    tol_param      = 1e-5,
    psi_damp       = 0.5,
    param_damp     = 0.2,
    seed           = 1,
    alpha_init     = NULL,
    Beta_init      = NULL,
    Gamma_init     = NULL,
    Psi_init       = NULL,
    return_fit     = FALSE
) {

  # Step 1: Fit
  message('Loading v1')
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
    Psi_init   = Psi_init
  )

  # Step 2: Predict
  out <- predict_logPAC_multislice(
    fit      = model_fit$fit,
    df_bio   = df_bio,
    df_tech  = df_tech,
    spe_list = spe_list #original data
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

  # optional
  run_config$session <- utils::sessionInfo()

  # Step 3: restore original coords
  spe_list <- lapply(spe_list, restore_coords)

  # Step 4: attach normalized data back to SPE
  spe_list_norm <- norm_spe(spe_list, out$out_list)

  for (i in seq_along(spe_list_norm)) {
    metadata(spe_list_norm[[i]])$MultiNorm <- run_config
  }

  # Step 4: return
  if (return_fit) {
    fit <- model_fit$fit

    fit$U_all   <- out$U_all
    fit$V_all   <- out$V_all
    fit$bio  <- out$bio
    fit$tech <- out$tech

    return(list(
      norm_spe = spe_list_norm,
      fit = fit,
      call = model_fit$call,
      Eta_full_list = out$Eta_full_list,
      percentiles = out$percentiles
    ))
  }

  return(spe_list_norm)
}
