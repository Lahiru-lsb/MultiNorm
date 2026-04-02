#inner loop optimized
irwls_multislice_backfit <- function(
    Y_list,        # list of (G × n_s) count matrices
    U_list,        # list of biological spline matrices (n_s × Kb)
    V_list,        # list of technical spline matrices (n_s × Kt)
    LS_list,       # list of library size vectors
    Psi_list,      # list of dispersion vectors (length G)
    Dmat,
    DTD,           # D^T D constraint matrix
    Pbeta,
    lambdaC,
    lambda1,
    lambda2,
    tol_param = 1e-4,
    max_iter = 50,
    param_damp = 0.5,
    alpha_init=NULL,
    Beta_init = NULL,
    Gamma_init = NULL

) {

  S <- length(Y_list)
  G <- nrow(Y_list[[1]])
  n_list <- sapply(Y_list, ncol)
  N <- sum(n_list)

  # ---- Pre-compute static quantities ----
  logLS_list <- lapply(LS_list, log)
  Ubeta_list <- lapply(U_list, function(U) cbind(1, U))

  # Pre-compute LV matrices once (they don't change)
  LV_list <- vector("list", S)
  for (s in seq_len(S)) {
    LV_list[[s]] <- logLS_list[[s]] * V_list[[s]]  # Vectorized scaling
  }

  # ---- Pre-allocate stacked matrices ----
  row_idx <- c(0, cumsum(n_list))

  Y_stack <- matrix(0, N, G)
  Psi_stack <- matrix(0, N, G)

  for (s in seq_len(S)) {
    rows <- (row_idx[s] + 1):row_idx[s + 1]
    Y_stack[rows, ] <- t(Y_list[[s]])
    Psi_stack[rows, ] <- matrix(rep(Psi_list[[s]], each = n_list[s]),
                                nrow = n_list[s])
  }

  # ---- Initialize ----

  if (is.null(alpha_init)) {
    alpha <- rep(0,times=S)
  } else {
    alpha <- alpha_init
  }

  Kb <- ncol(Ubeta_list[[1]])

  if (is.null(Beta_init)) {
    Beta_list <- lapply(seq_len(S), function(s) matrix(0.01, Kb, G))
  } else {
    #message("Using user-provided beta")
    Beta_list <- Beta_init
  }

  if (is.null(Gamma_init)) {
    Gamma_list <- lapply(V_list, function(V) matrix(0, ncol(V), G))
  } else {
    #message("Using user-provided gamma")
    Gamma_list <- Gamma_init
  }

  # Pre-allocate working arrays
  Eta <- matrix(0, N, G)
  Mu <- matrix(0, N, G)
  Z_stack <- matrix(0, N, G)
  W_stack <- matrix(0, N, G)

  # ---- Initial predictor ----
  for (s in seq_len(S)) {
    rows <- (row_idx[s] + 1):row_idx[s + 1]
    Eta[rows, ] <- alpha[s] * logLS_list[[s]] +
      Ubeta_list[[s]] %*% Beta_list[[s]] +
      LV_list[[s]] %*% Gamma_list[[s]]
  }

  Mu[] <- exp(Eta)

  ll_old <- -Inf

  # ================= IRWLS LOOP =================
  for (iter in seq_len(max_iter)) {

    # ---- Working quantities (in-place) ----
    Z_stack[] <- Eta + (Y_stack - Mu) / Mu
    W_stack[] <- Mu / (1 + Psi_stack * Mu)

    # ---- Split by slice ----
    Z_list <- vector("list", S)
    W_list <- vector("list", S)

    for (s in seq_len(S)) {
      rows <- (row_idx[s] + 1):row_idx[s + 1]
      Z_list[[s]] <- Z_stack[rows, , drop = FALSE]
      W_list[[s]] <- rowMeans(W_stack[rows, , drop = FALSE])
    }



    # ---- Alpha update ----
    alpha_new <- update_alpha_multislice(
      Z_list, W_list,
      logLS_list,
      Ubeta_list, Beta_list,
      V_list, Gamma_list
    )

    alpha <- (1 - param_damp) * alpha + param_damp * alpha_new
    #alpha <- alpha_new

    # ---- Gamma update ----
    Gamma_list_new <- update_gamma_multislice(
      Z_list, W_list,
      alpha, logLS_list,
      Ubeta_list, Beta_list,
      LV_list,
      lambda_gamma = lambda2
    )

    # Apply cap in-place
    for (s in seq_len(S)) {
      Gamma_list_new[[s]] <- cap_coeff(Gamma_list_new[[s]])
    }

    # ---- Damping (in-place) ----
    for (s in seq_len(S)) {
      Gamma_list[[s]] <- (1 - param_damp) * Gamma_list[[s]] + param_damp * Gamma_list_new[[s]]
    }


    # ---- Beta update ----
    Beta_list_new <- update_beta_multislice(
      Z_list, W_list,
      alpha, logLS_list,
      Ubeta_list, V_list,
      Gamma_list,
      DTD, Pbeta,
      lambdaC, lambda1
    )

    for (s in seq_len(S)) {
      Beta_list_new[[s]] <- cap_coeff(Beta_list_new[[s]])
    }

    # ---- Damping (in-place) ----
    for (s in seq_len(S)) {
      Beta_list[[s]] <- (1 - param_damp) * Beta_list[[s]] + param_damp * Beta_list_new[[s]]
    }


    # ---- Recompute Eta and Mu ----
    for (s in seq_len(S)) {
      rows <- (row_idx[s] + 1):row_idx[s + 1]
      Eta[rows, ] <- alpha[s] * logLS_list[[s]] +
        Ubeta_list[[s]] %*% Beta_list[[s]] +
        LV_list[[s]] %*% Gamma_list[[s]]
    }
    Mu[] <- exp(Eta)



    # ---- Penalized likelihood ----
    ll_new <- sum(penalized_ll(
      Y = Y_stack,
      Mu = Mu,
      Psi = Psi_stack,
      Beta = do.call(rbind, Beta_list),
      Gamma = do.call(rbind, Gamma_list),
      Cmat = Dmat,
      Pbeta = Pbeta,
      lambdaC = lambdaC,
      lambda1 = lambda1,
      lambda2 = lambda2,
      n_cells = N
    ))

    message(sprintf(
      "Iter %d | log-likelihood: %.4f | Mu range: [%.2f, %.2f]",
      iter, ll_new , min(Mu), max(Mu)
    ))

    # ---- Convergence check ----
    if (iter > 1 && abs(ll_new - ll_old) / (abs(ll_old)) < tol_param) {
      message(sprintf("Converged at iteration %d", iter))
      break
    }

    ll_old <- ll_new
  }

  list(
    alpha = alpha,
    Beta  = Beta_list,
    Gamma = Gamma_list,
    Mu    = Mu,
    Eta   = Eta,
    ll    = ll_old,
    iterations = iter
  )
}



