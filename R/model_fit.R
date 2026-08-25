fit_spatial_nb_multislice <- function(
    spe_list,
    niche_col,
    preferred_niches = NULL,

    spot_frac = 0.1,
    verbose = FALSE,

    df_bio = c(6, 6),
    df_tech = c(3, 3),

    lambdaC = 1e4,
    lambda1 = 1e-2,
    lambda2 = 1.9,

    max_outer = 30,
    max_inner = 200,
    tol_disp = 1e-4,
    tol_param = 1e-4,
    psi_damp = 0.3,
    param_damp = 0.2,
    seed = 1,

    alpha_init = NULL,
    Beta_init = NULL,
    Gamma_init = NULL,
    Psi_init = NULL,
    zeta_init = NULL,
    rho_init = NULL,
    A = NULL
) {

  start_time <- Sys.time()
  on.exit({
    end_time <- Sys.time()
    runtime <- as.numeric(difftime(end_time, start_time, units = "mins"))
    message("Total runtime (mins): ", round(runtime, 2))
  })

  stopifnot(length(df_bio) == 2, length(df_tech) == 2)
  set.seed(seed)

  S <- length(spe_list)
  if (S < 2) stop("Need at least two slices.")

  gene_names <- lapply(spe_list, rownames)
  if (!all(vapply(gene_names[-1], identical, logical(1), gene_names[[1]]))) {
    stop("All slices must have identical genes in the same order.")
  }

  if (is.numeric(niche_col)) {
    niche_col <- colnames(SummarizedExperiment::colData(spe_list[[1]]))[niche_col]
  }

  spe_list <- lapply(spe_list, function(spe) {

    # Store original coordinates
    coords_orig <- SpatialExperiment::spatialCoords(spe)
    # select only x y coords
    coords_orig <- coords_orig[,1:2]
    S4Vectors::metadata(spe)$orig_coords <- coords_orig

    # Scale full-data coordinates before sampling
    coords_scaled <- apply(coords_orig, 2, function(x) {
      rng <- max(x) - min(x)
      if (rng == 0) stop("Coordinate column has zero range.")
      (x - min(x)) / rng - 0.5
    })

    SpatialExperiment::spatialCoords(spe) <- coords_scaled

    # Store niche column consistently
    spe$niche <- spe[[niche_col]]

    # Precompute full-data LS before sampling
    if (!"LS" %in% colnames(SummarizedExperiment::colData(spe))) {
      cts_full <- SummarizedExperiment::assay(spe, "counts")
      LS_full <- as.numeric(MatrixGenerics::colSums(cts_full))
      LS_full <- LS_full / mean(LS_full)

      SummarizedExperiment::colData(spe)$LS <- LS_full
    }

    spe
  })

  # Build full spline bases before sampling
  full_basis_list <- lapply(spe_list, function(spe) {

    coords <- SpatialExperiment::spatialCoords(spe)

    U_full <- sp2d(
      coords[, 1],
      coords[, 2],
      dfx = df_bio[1],
      dfy = df_bio[2]
    )

    V_full <- sp2d(
      coords[, 1],
      coords[, 2],
      dfx = df_tech[1],
      dfy = df_tech[2]
    )

    list(
      U_full = U_full,
      V_full = V_full
    )
  })

  if (!is.null(spot_frac) && spot_frac < 1) {
    spe_list <- lapply(spe_list, function(spe) {
      sample_spe_proportional(
        spe,
        niche_col = niche_col,
        sample.p  = spot_frac,
        verbose   = verbose,
        seed      = seed
      )
    })
  } else {
    message("Skipping sampling (spot_frac >= 1 or NULL)")
  }

  # Extract sampled fitting data
  data_list <- lapply(spe_list, get_slice_data)

  # Subset full spline bases to sampled spots for fitting
  U_list <- vector("list", S)
  V_list <- vector("list", S)

  for (s in seq_len(S)) {

    sampled_idx <- S4Vectors::metadata(spe_list[[s]])$sampled_idx

    if (is.null(sampled_idx)) {
      sampled_idx <- seq_len(ncol(spe_list[[s]]))
    }

    U_list[[s]] <- full_basis_list[[s]]$U_full[sampled_idx, , drop = FALSE]
    V_list[[s]] <- full_basis_list[[s]]$V_full[sampled_idx, , drop = FALSE]
  }

  niche_sets <- lapply(data_list, function(d) unique(d$niche))
  all_niches <- sort(unique(unlist(niche_sets)))
  all_niches <- all_niches[!is.na(all_niches)]

  if (!is.null(preferred_niches)) {
    message("User-defined niche set will be used")
    all_niches <- intersect(all_niches, preferred_niches)

    if (length(all_niches) == 0) {
      stop("None of the preferred niches are present in the data.")
    }
  }

  # Biological niche means, no intercept
  ubar_list <- lapply(seq_len(S), function(s)
    niche_means(
      U_list[[s]],
      data_list[[s]]$niche,
      all_niches
    )
  )

  Y_list <- lapply(spe_list, function(spe)
    as.matrix(SummarizedExperiment::assay(spe, "counts"))
  )

  LS_list <- lapply(data_list, function(d) d$LS)

  fit <- fit_all_genes_nb_multislice(
    Y_list = Y_list,
    U_list = U_list,
    V_list = V_list,
    LS_list = LS_list,
    ubar_list = ubar_list,
    lambdaC = lambdaC,
    lambda1 = lambda1,
    lambda2 = lambda2,
    tol_disp = tol_disp,
    max_outer = max_outer,
    max_inner = max_inner,
    psi_damp = psi_damp,
    param_damp = param_damp,
    tol_param = tol_param,
    alpha_init = alpha_init,
    Beta_init = Beta_init,
    Gamma_init = Gamma_init,
    Psi_init = Psi_init,
    zeta_init = zeta_init,
    rho_init = rho_init,
    A = A
  )

  # Store full bases for prediction
  fit$U_full_list <- lapply(full_basis_list, `[[`, "U_full")
  fit$V_full_list <- lapply(full_basis_list, `[[`, "V_full")

  list(
    fit = fit,
    genes = rownames(spe_list[[1]]),
    niches = all_niches,
    slices = paste0("slice", seq_len(S)),
    call = match.call()
  )
}
