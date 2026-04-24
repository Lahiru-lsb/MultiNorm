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
  nm<- names(spe_list)
  # ---- Pre-allocate lists ----
  Eta_full_list <- vector("list", S)
  Eta_bio_list  <- vector("list", S)
  Psi_list_exp  <- vector("list", S)
  Y_list        <- vector("list", S)
  U_full_list   <- vector("list", S)
  V_full_list   <- vector("list", S)
  n_list        <- integer(S)  # Use integer for indices

  # ---- Assign names ----
  names(Eta_full_list) <- nm
  names(Eta_bio_list)  <- nm
  names(Psi_list_exp)  <- nm
  names(Y_list)        <- nm
  names(U_full_list)   <- nm
  names(V_full_list)   <- nm
  names(n_list)        <- nm

  # ================= LOOP OVER SLICES =================
  for (s in seq_len(S)) {

    spe <- spe_list[[s]]

    # ---- extract slice data ----
    d <- get_slice_data(spe)

    # ---- spline bases ----
    U <- sp2d(d$x, d$y, dfx = df_bio[1], dfy = df_bio[2])
    V <- sp2d(d$x, d$y, dfx = df_tech[1], dfy = df_tech[2])

    #store to output U and V
    V_full_list[[s]] <- V
    U_full_list[[s]] <- U

    blk <- build_blocks(d$LS, U, V)

    X_beta  <- blk$X_beta
    X_gamma <- blk$X_gamma

    logLS <- log(d$LS)
    n_s   <- length(d$LS)
    n_list[s] <- n_s

    # ---- parameters from fit ----
    alpha_s <- fit$alpha[s]
    Beta_s  <- fit$Beta[[s]]
    Gamma_s <- fit$Gamma[[s]]

    # ---- full mean ----
    Eta_full <-
      alpha_s * logLS +
      X_beta  %*% Beta_s +
      X_gamma %*% Gamma_s

    Eta_full <- pmin(pmax(Eta_full, -30), 30)
    # Eta_full_max <- matrixStats::colMedians(Eta_full) +
    #   4 * matrixStats::colMads(Eta_full)
    #
    # Eta_full <- pmin(
    #   Eta_full,
    #   matrix(Eta_full_max, nrow(Eta_full), ncol(Eta_full), byrow = TRUE)
    # )

    Mu_full <- exp(Eta_full)

    # ---- biological mean ----
    Eta_bio <- X_beta %*% Beta_s
    Eta_bio <- pmin(pmax(Eta_bio, -30), 30)
    # Eta_bio_max <- matrixStats::colMedians(Eta_bio) +
    #   4 * matrixStats::colMads(Eta_bio)
    #
    # Eta_bio <- pmin(
    #   Eta_bio,
    #   matrix(Eta_bio_max, nrow(Eta_bio), ncol(Eta_bio), byrow = TRUE)
    # )

    Mu_bio <- exp(Eta_bio)


    # ---- store ----
    Eta_full_list[[s]] <- Eta_full
    Eta_bio_list[[s]]  <- Eta_bio

    # ---- expand dispersion ----
    Psi_list_exp[[s]] <-
      matrix(rep(fit$Psi[[s]], each = n_s), nrow = n_s)

    # ---- counts ----
    Y_list[[s]] <- t(as.matrix(assay(spe, "counts")))
  }

  # ---------------- stack (optimized) ----------------
  N <- sum(n_list)
  G <- ncol(Eta_full_list[[1]])

  # Pre-allocate stacked matrices
  Eta_full <- matrix(0, N, G)
  Eta_bio  <- matrix(0, N, G)
  Psi      <- matrix(0, N, G)
  Y        <- matrix(0, N, G)

  # Pre-compute row indices
  row_idx <- c(0, cumsum(n_list))

  # Fill stacked matrices
  for (s in seq_len(S)) {
    rows <- (row_idx[s] + 1):row_idx[s + 1]
    Eta_full[rows, ] <- Eta_full_list[[s]]
    Eta_bio[rows, ]  <- Eta_bio_list[[s]]
    Psi[rows, ]      <- Psi_list_exp[[s]]
    Y[rows, ]        <- Y_list[[s]]
  }

  # Apply exp in-place
  Mu_full <- exp(Eta_full)
  Mu_bio  <- exp(Eta_bio)

  # Create slice index
  slice_index <- rep(seq_len(S), times = n_list)

  # ---------------- PAC normalization ----------------
  norm_logPAC <- compute_logPAC(
    Y       = Y,
    Mu_full = Mu_full,
    Mu_bio  = Mu_bio,
    Psi     = Psi,
    slice_index = slice_index
  )

  dimnames(norm_logPAC$norm) <- dimnames(Y)

  # ---------------- split back to slices ----------------
  out_list <- vector("list", S)

  for (s in seq_len(S)) {
    rows <- (row_idx[s] + 1):row_idx[s + 1]
    out_list[[s]] <- norm_logPAC$norm[rows, , drop = FALSE]
  }

  names(out_list) <-names(spe_list)
  list(
    out_list=out_list,
    Eta_full_list=Eta_full_list,
    percentiles=norm_logPAC$percentiles,
    V_all=V_full_list,
    U_all=U_full_list

  )
  #out_list
}
