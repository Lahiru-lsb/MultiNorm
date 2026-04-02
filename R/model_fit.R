
fit_spatial_nb_multislice <- function(
    spe_list,               # list of SpatialExperiment objects
    niche_col,
    preferred_niches = NULL,

    spot_frac = 0.1,
    verbose=FALSE,

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
    param_damp=0.2,
    seed = 1,

    alpha_init = NULL,
    Beta_init = NULL,
    Gamma_init = NULL,
    Psi_init = NULL
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

  # --- gene consistency check ---
  gene_names <- lapply(spe_list, rownames)
  if (!all(vapply(gene_names[-1], identical, logical(1), gene_names[[1]]))) {
    stop("All slices must have identical genes in the same order.")
  }

  # ---------- normalize coords + attach niches ----------
  spe_list <- lapply(spe_list, function(spe) {

    coords <- spatialCoords(spe)
    # store originals in metadata
    metadata(spe)$orig_coords <- coords

    coords <- apply(coords, 2, function(x)
      (x - min(x)) / (max(x) - min(x)) - 0.5
    )

    #cat("coords range:", range(coords[,1]), "\n")
    # --- update in place ---
    spatialCoords(spe) <- coords

    # attach niche
    spe$niche <- spe[[niche_col]]
    return(spe)

  })

  # # ---------- drop zero-library spots ----------
  # spe_list <- lapply(spe_list, function(spe)
  #   spe[, MatrixGenerics::colSums(SummarizedExperiment::assay(spe, "counts")) > 0]
  # )
  #
  # # ---------- subsample spots ----------
  if (!is.null(spot_frac) && spot_frac < 1) {
    spe_list <- lapply(spe_list, function(spe) {
      sample_spe_proportional(
        spe,
        niche_col = niche_col,
        sample.p  = spot_frac,
        verbose   = verbose,
        seed = seed
      )
    })
  } else {
    message("Skipping sampling (spot_frac >= 1 or NULL)")
  }

  # ---------- extract slice data ----------
  data_list <- lapply(spe_list, get_slice_data)

  # ---------- spline bases ----------
  U_list <- lapply(data_list, function(d)
    sp2d(d$x, d$y, df_bio[1], df_bio[2])
  )

  V_list <- lapply(data_list, function(d)
    sp2d(d$x, d$y, df_tech[1], df_tech[2])
  )

  niche_sets <- lapply(data_list, function(d) unique(d$niche))
  all_niches <- sort(unique(unlist(niche_sets)))
  #remove NA niche if present
  all_niches <- all_niches[!is.na(all_niches)]

  if (!is.null(preferred_niches)) {
    message("User-defined niche set will be used")
    all_niches <- intersect(all_niches, preferred_niches)

    if (length(all_niches) == 0) {
      stop("None of the preferred niches are present in the data.")
    }
  }

  # ---------- niche means ----------
  Ubeta_list <- lapply(U_list, function(U) cbind(1, U))

  ubar_list <- lapply(seq_len(S), function(s)
    niche_means(Ubeta_list[[s]],
                data_list[[s]]$niche,
                all_niches)
  )

  # ---------- counts + LS ----------
  Y_list  <- lapply(spe_list, function(spe)
    as.matrix(SummarizedExperiment::assay(spe, "counts"))
  )

  LS_list <- lapply(data_list, function(d) d$LS)

  # ---------- model fit ----------
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
    param_damp=param_damp,
    tol_param = tol_param,
    alpha_init = alpha_init,
    Beta_init = Beta_init,
    Gamma_init = Gamma_init,
    Psi_init = Psi_init
  )

  # ---------- return ----------
  list(
    fit = fit,
    genes = rownames(spe_list[[1]]),
    niches = all_niches,
    slices = paste0("slice", seq_len(S)),
    call = match.call()
  )
}
