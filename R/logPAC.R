#----optimized logpac and post estimation----
compute_logPAC <- function(
    Y,              # cells x genes
    Mu_full,
    Mu_bio,
    Psi,            # cells x genes (psi_sg expanded to cells)
    slice_index     # length Ncells
) {

  Y_mat   <- as.matrix(Y)
  Mu_mat  <- as.matrix(Mu_full)
  Mu_bio  <- as.matrix(Mu_bio)
  Psi_mat <- as.matrix(Psi)

  Psi_w <- Psi_mat

  slices <- unique(slice_index)

  # -------- slice-level winsorization --------
  for (s in slices) {

    rows <- which(slice_index == s)

    # take one representative row (since psi identical within slice)
    psi_vec <- Psi_mat[rows[1], ]

    psi_max <- exp(median(log(psi_vec)) + 3 * mad(log(psi_vec)))

    # Vectorized: broadcast psi_max across all rows
    Psi_w[rows, ] <- pmin(Psi_mat[rows, ], psi_max)
  }
  #--------  transform  --------
  inv_Psi_w <- 1 / Psi_w

  lb <- pnbinom(Y_mat - 1, mu = Mu_mat, size = inv_Psi_w)
  ub <- dnbinom(Y_mat, mu = Mu_mat, size = inv_Psi_w) + lb
  pU <- (lb + ub) / 2
  pU <- pmax(pmin(pU, 0.999), 0.001)

  # -------- biological back-transform --------
  norm=log2(qnbinom(pU, mu = Mu_bio, size = 1 / Psi_w) + 1)

  return(
    list(percentiles=pU,
         norm=norm)
  )
}


#post estimations
predict_logPAC_multislice <- function(
    fit,
    spe_list,
    df_bio  = c(6, 6),
    df_tech = c(3, 3)
) {

  S <- length(spe_list)
  nm <- names(spe_list)

  Eta_full_list <- vector("list", S)
  Eta_bio_list  <- vector("list", S)
  Eta_tech_list <- vector("list", S)
  Psi_list_exp  <- vector("list", S)
  Y_list        <- vector("list", S)
  U_full_list   <- vector("list", S)
  V_full_list   <- vector("list", S)
  n_list        <- integer(S)

  names(Eta_full_list) <- nm
  names(Eta_bio_list)  <- nm
  names(Eta_tech_list) <- nm
  names(Psi_list_exp)  <- nm
  names(Y_list)        <- nm
  names(U_full_list)   <- nm
  names(V_full_list)   <- nm
  names(n_list)        <- nm

  for (s in seq_len(S)) {

    spe <- spe_list[[s]]

    d <- get_slice_data(spe)

    U <- fit$U_full_list[[s]]
    V <- fit$V_full_list[[s]]

    if (is.null(U) || is.null(V)) {
      stop("Full spline bases are missing from fit.")
    }

    if (nrow(U) != ncol(SummarizedExperiment::assay(spe, "counts"))) {
      stop("nrow(U_full) does not match number of spots in prediction SPE.")
    }

    if (nrow(V) != ncol(SummarizedExperiment::assay(spe, "counts"))) {
      stop("nrow(V_full) does not match number of spots in prediction SPE.")
    }

    U_full_list[[s]] <- U
    V_full_list[[s]] <- V

    blk <- build_blocks(d$LS, U, V)

    X_beta  <- blk$X_beta      # n_s × Kb, no intercept
    X_gamma <- blk$X_gamma     # n_s × Kt

    logLS <- log(d$LS)
    n_s   <- length(d$LS)
    n_list[s] <- n_s

    alpha_s <- fit$alpha[s]
    zeta    <- fit$zeta
    delta_s <- fit$delta[s, ]

    Beta_s  <- fit$Beta[[s]]
    Gamma_s <- fit$Gamma[[s]]

    Eta_bio <-
      matrix(zeta, nrow = n_s, ncol = length(zeta), byrow = TRUE) +
      X_beta %*% Beta_s

    Eta_tech <-
      matrix(delta_s, nrow = n_s, ncol = length(delta_s), byrow = TRUE) +
      alpha_s * logLS +
      X_gamma %*% Gamma_s

    Eta_full <- Eta_bio + Eta_tech

    Eta_full <- pmin(pmax(Eta_full, -30), 30)
    Eta_bio  <- pmin(pmax(Eta_bio,  -30), 30)

    Eta_full_list[[s]] <- Eta_full
    Eta_bio_list[[s]]  <- Eta_bio
    Eta_tech_list[[s]] <- Eta_tech

    Psi_list_exp[[s]] <-
      matrix(rep(fit$Psi[[s]], each = n_s), nrow = n_s)

    Y_list[[s]] <- t(as.matrix(assay(spe, "counts")))
  }

  N <- sum(n_list)
  G <- ncol(Eta_full_list[[1]])

  Eta_full <- matrix(0, N, G)
  Eta_bio  <- matrix(0, N, G)
  Eta_tech <- matrix(0, N, G)
  Psi      <- matrix(0, N, G)
  Y        <- matrix(0, N, G)

  row_idx <- c(0, cumsum(n_list))

  for (s in seq_len(S)) {
    rows <- (row_idx[s] + 1):row_idx[s + 1]

    Eta_full[rows, ] <- Eta_full_list[[s]]
    Eta_bio[rows, ]  <- Eta_bio_list[[s]]
    Eta_tech[rows, ] <- Eta_tech_list[[s]]
    Psi[rows, ]      <- Psi_list_exp[[s]]
    Y[rows, ]        <- Y_list[[s]]
  }

  Mu_full <- exp(Eta_full)
  Mu_bio  <- exp(Eta_bio)

  slice_index <- rep(seq_len(S), times = n_list)

  norm_logPAC <- compute_logPAC(
    Y           = Y,
    Mu_full     = Mu_full,
    Mu_bio      = Mu_bio,
    Psi         = Psi,
    slice_index = slice_index
  )

  dimnames(norm_logPAC$norm) <- dimnames(Y)

  out_list <- vector("list", S)

  for (s in seq_len(S)) {
    rows <- (row_idx[s] + 1):row_idx[s + 1]
    out_list[[s]] <- norm_logPAC$norm[rows, , drop = FALSE]
  }

  names(out_list) <- names(spe_list)

  list(
    out_list      = out_list,
    Eta_full_list = Eta_full_list,
    percentiles   = norm_logPAC$percentiles,
    V_all         = V_full_list,
    U_all         = U_full_list,
    bio           = Eta_bio,
    tech          = Eta_tech,
    Mu_full       = Mu_full,
    Mu_bio        = Mu_bio
  )
}
