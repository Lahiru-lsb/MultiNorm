
irwls_multislice_backfit <- function(
    Y_list,
    U_list,
    V_list,
    LS_list,
    Psi_list,
    Dmat,
    DTD,
    lambdaC,
    lambda1,
    lambda2,
    tol_param = 1e-4,
    max_iter = 50,
    param_damp = 0.5,
    alpha,
    Beta_list,
    Gamma_list,
    zeta = NULL,
    rho = NULL,
    A = NULL
) {

  S <- length(Y_list)
  G <- nrow(Y_list[[1]])
  n_list <- sapply(Y_list, ncol)
  N <- sum(n_list)

  logLS_list <- lapply(LS_list, log)

  delta <- A %*% rho

  LV_list <- vector("list", S)
  for (s in seq_len(S)) {
    LV_list[[s]] <- logLS_list[[s]] * V_list[[s]]
  }

  row_idx <- c(0, cumsum(n_list))

  Y_stack <- matrix(0, N, G)
  Psi_stack <- matrix(0, N, G)

  for (s in seq_len(S)) {
    rows <- (row_idx[s] + 1):row_idx[s + 1]
    Y_stack[rows, ] <- t(Y_list[[s]])
    Psi_stack[rows, ] <- matrix(
      rep(Psi_list[[s]], each = n_list[s]),
      nrow = n_list[s]
    )
  }

  Eta   <- matrix(0, N, G)
  Mu    <- matrix(0, N, G)
  Z_stack <- matrix(0, N, G)
  W_stack <- matrix(0, N, G)
  Eta_try <- matrix(0, N, G)   # declared once, reused every attempt

  for (s in seq_len(S)) {
    rows <- (row_idx[s] + 1):row_idx[s + 1]

    Eta[rows, ] <-
      matrix(zeta, nrow = n_list[s], ncol = G, byrow = TRUE) +
      matrix(delta[s, ], nrow = n_list[s], ncol = G, byrow = TRUE) +
      alpha[s] * logLS_list[[s]] +
      U_list[[s]] %*% Beta_list[[s]] +
      LV_list[[s]] %*% Gamma_list[[s]]
  }

  Mu[] <- exp(Eta)

  ll_old <- sum(penalized_ll(
    Y = Y_stack,
    Mu = Mu,
    Psi = Psi_stack,
    Beta = do.call(rbind, Beta_list),
    Gamma = do.call(rbind, Gamma_list),
    Cmat = Dmat,
    lambdaC = lambdaC,
    lambda1 = lambda1,
    lambda2 = lambda2,
    n_cells = N
  ))

  for (iter in seq_len(max_iter)) {

    Z_stack[] <- Eta + (Y_stack - Mu) / Mu
    W_stack[] <- Mu / (1 + Psi_stack * Mu)

    Z_list <- vector("list", S)
    W_list <- vector("list", S)

    for (s in seq_len(S)) {
      rows <- (row_idx[s] + 1):row_idx[s + 1]
      Z_list[[s]] <- Z_stack[rows, , drop = FALSE]
      W_list[[s]] <- rowMeans(W_stack[rows, , drop = FALSE])
    }

    alpha_best <- alpha
    Beta_best  <- Beta_list
    Gamma_best <- Gamma_list
    zeta_best  <- zeta
    rho_best   <- rho
    delta_best <- delta
    Eta_best   <- Eta
    Mu_best    <- Mu
    ll_best    <- ll_old

    # ---- All _new values computed once, outside the attempt loop ----
    zeta_new <- update_zeta(
      Z_list, W_list,
      logLS_list,
      U_list, Beta_list,
      LV_list, Gamma_list,
      alpha,
      delta
    )

    rho_new <- update_rho(
      Z_list, W_list,
      logLS_list,
      U_list, Beta_list,
      LV_list, Gamma_list,
      alpha,
      zeta_new,
      A
    )

    alpha_new <- update_alpha_multislice(
      Z_list, W_list,
      logLS_list,
      U_list, Beta_list,
      V_list, Gamma_list,
      zeta_new,
      A %*% rho_new
    )

    Gamma_list_new <- update_gamma_multislice(
      Z_list, W_list,
      alpha_new,
      logLS_list,
      U_list,
      Beta_list,
      LV_list,
      zeta_new,
      A %*% rho_new,
      lambda_gamma = lambda2
    )

    for (s in seq_len(S)) {
      Gamma_list_new[[s]] <- cap_coeff(Gamma_list_new[[s]])
    }

    Beta_list_new <- update_beta_multislice(
      Z_list, W_list,
      alpha_new,
      logLS_list,
      U_list,
      V_list,
      Gamma_list_new,
      zeta_new,
      A %*% rho_new,
      DTD,
      lambdaC,
      lambda1
    )

    for (s in seq_len(S)) {
      Beta_list_new[[s]] <- cap_coeff(Beta_list_new[[s]])
    }

    # ---- Attempt loop: ----
    step_factor <- 1

    for (attempt in 1:4) {

      damp <- param_damp * step_factor

      zeta_try  <- (1 - damp) * zeta  + damp * zeta_new
      rho_try   <- (1 - damp) * rho   + damp * rho_new
      delta_try <- A %*% rho_try
      alpha_try <- (1 - damp) * alpha + damp * alpha_new

      # zeta_try  <- zeta_new    #(1 - damp) * zeta  + damp * zeta_new
      # rho_try   <- rho_new     #(1 - damp) * rho   + damp * rho_new
      # delta_try <- A %*% rho_try
      # alpha_try <- alpha_new   #(1 - damp) * alpha + damp * alpha_new

      Gamma_list_try <- vector("list", S)
      for (s in seq_len(S)) {
        Gamma_list_try[[s]] <-
          (1 - damp) * Gamma_list[[s]] + damp * Gamma_list_new[[s]]
      }

      Beta_list_try <- vector("list", S)
      for (s in seq_len(S)) {
        Beta_list_try[[s]] <-
          (1 - damp) * Beta_list[[s]] + damp * Beta_list_new[[s]]
      }

      for (s in seq_len(S)) {
        rows <- (row_idx[s] + 1):row_idx[s + 1]

        Eta_try[rows, ] <-
          matrix(zeta_try, nrow = n_list[s], ncol = G, byrow = TRUE) +
          matrix(delta_try[s, ], nrow = n_list[s], ncol = G, byrow = TRUE) +
          alpha_try[s] * logLS_list[[s]] +
          U_list[[s]] %*% Beta_list_try[[s]] +
          LV_list[[s]] %*% Gamma_list_try[[s]]
      }

      Mu_try <- exp(Eta_try)

      ll_try <- sum(penalized_ll(
        Y = Y_stack,
        Mu = Mu_try,
        Psi = Psi_stack,
        Beta = do.call(rbind, Beta_list_try),
        Gamma = do.call(rbind, Gamma_list_try),
        Cmat = Dmat,
        lambdaC = lambdaC,
        lambda1 = lambda1,
        lambda2 = lambda2,
        n_cells = N
      ))

      if (ll_try >= ll_best) {

        alpha_best <- alpha_try
        Beta_best  <- Beta_list_try
        Gamma_best <- Gamma_list_try
        zeta_best  <- zeta_try
        rho_best   <- rho_try
        delta_best <- delta_try
        Eta_best   <- Eta_try
        Mu_best    <- Mu_try
        ll_best    <- ll_try

        if (attempt > 1) message("Step halved ", attempt - 1, " times")
        break

      } else {
        step_factor <- step_factor / 2
      }
    }

    alpha <- alpha_best
    Beta_list <- Beta_best
    Gamma_list <- Gamma_best
    zeta <- zeta_best
    rho <- rho_best
    delta <- delta_best
    Eta <- Eta_best
    Mu <- Mu_best
    ll_new <- ll_best

    message(sprintf(
      "Iter %d | log-likelihood: %.2f | Mu range: [%.2f, %.2f]",
      iter, ll_new, min(Mu), max(Mu)
    ))

    if (iter > 1 && abs(ll_new - ll_old) / abs(ll_old) < tol_param) {
      message(sprintf("Converged at iteration %d", iter))
      break
    }

    ll_old <- ll_new
  }

  Eta_nointercept <- matrix(0, N, G)

  for (s in seq_len(S)) {
    rows <- (row_idx[s] + 1):row_idx[s + 1]

    Eta_nointercept[rows, ] <-
      alpha[s] * logLS_list[[s]] +
      U_list[[s]] %*% Beta_list[[s]] +
      LV_list[[s]] %*% Gamma_list[[s]]
  }

  Mu_nointercept <- exp(Eta_nointercept)

  Eta_bio <- matrix(0, N, G)

  for (s in seq_len(S)) {
    rows <- (row_idx[s] + 1):row_idx[s + 1]

    Eta_bio[rows, ] <-
      matrix(zeta, nrow = n_list[s], ncol = G, byrow = TRUE) +
      U_list[[s]] %*% Beta_list[[s]]
  }

  Mu_bio <- exp(Eta_bio)

  list(
    alpha = alpha,
    zeta  = zeta,
    rho   = rho,
    delta = delta,
    Beta  = Beta_list,
    Gamma = Gamma_list,
    Mu    = Mu,
    Mu_bio = Mu_bio,
    Mu_nointercept = Mu_nointercept,
    Eta   = Eta,
    ll    = ll_old,
    iterations = iter
  )
}



# irwls_multislice_backfit1 <- function(
#     Y_list,
#     U_list,
#     V_list,
#     LS_list,
#     Psi_list,
#     Dmat,
#     DTD,
#     lambdaC,
#     lambda1,
#     lambda2,
#     tol_param = 1e-4,
#     max_iter = 50,
#     param_damp = 0.5,
#     alpha,
#     Beta_list,
#     Gamma_list,
#     zeta = NULL,
#     rho = NULL,
#     A = NULL,
#     ll_start = NULL
# ) {
#
#   S <- length(Y_list)
#   G <- nrow(Y_list[[1]])
#   n_list <- sapply(Y_list, ncol)
#   N <- sum(n_list)
#
#   logLS_list <- lapply(LS_list, log)
#
#   delta <- A %*% rho
#
#   LV_list <- vector("list", S)
#   for (s in seq_len(S)) {
#     LV_list[[s]] <- logLS_list[[s]] * V_list[[s]]
#   }
#
#   row_idx <- c(0, cumsum(n_list))
#
#   Y_stack <- matrix(0, N, G)
#   Psi_stack <- matrix(0, N, G)
#
#   for (s in seq_len(S)) {
#     rows <- (row_idx[s] + 1):row_idx[s + 1]
#
#     Y_stack[rows, ] <- t(Y_list[[s]])
#
#     Psi_stack[rows, ] <- matrix(
#       rep(Psi_list[[s]], each = n_list[s]),
#       nrow = n_list[s]
#     )
#   }
#
#   Eta     <- matrix(0, N, G)
#   Mu      <- matrix(0, N, G)
#   Z_stack <- matrix(0, N, G)
#   W_stack <- matrix(0, N, G)
#   Eta_try <- matrix(0, N, G)
#
#   for (s in seq_len(S)) {
#     rows <- (row_idx[s] + 1):row_idx[s + 1]
#
#     Eta[rows, ] <-
#       matrix(zeta, nrow = n_list[s], ncol = G, byrow = TRUE) +
#       matrix(delta[s, ], nrow = n_list[s], ncol = G, byrow = TRUE) +
#       alpha[s] * logLS_list[[s]] +
#       U_list[[s]] %*% Beta_list[[s]] +
#       LV_list[[s]] %*% Gamma_list[[s]]
#   }
#
#   Mu[] <- exp(Eta)
#
#   # ------------------------------------------------------------
#   # Baseline likelihood under the fixed Psi for this inner loop.
#   # If outer loop supplies ll_start, use it.
#   # Otherwise compute it here.
#   # ------------------------------------------------------------
#
#   if (is.null(ll_start)) {
#
#     ll_old <- sum(penalized_ll(
#       Y = Y_stack,
#       Mu = Mu,
#       Psi = Psi_stack,
#       Beta = do.call(rbind, Beta_list),
#       Gamma = do.call(rbind, Gamma_list),
#       Cmat = Dmat,
#       lambdaC = lambdaC,
#       lambda1 = lambda1,
#       lambda2 = lambda2,
#       n_cells = N
#     ))
#
#   } else {
#
#     ll_old <- ll_start
#   }
#
#   converged <- FALSE
#
#   for (iter in seq_len(max_iter)) {
#
#     Z_stack[] <- Eta + (Y_stack - Mu) / Mu
#     W_stack[] <- Mu / (1 + Psi_stack * Mu)
#
#     Z_list <- vector("list", S)
#     W_list <- vector("list", S)
#
#     for (s in seq_len(S)) {
#       rows <- (row_idx[s] + 1):row_idx[s + 1]
#
#       Z_list[[s]] <- Z_stack[rows, , drop = FALSE]
#
#       # Common slice/cell weights across genes
#       W_list[[s]] <- rowMeans(W_stack[rows, , drop = FALSE])
#     }
#
#     alpha_best <- alpha
#     Beta_best  <- Beta_list
#     Gamma_best <- Gamma_list
#     zeta_best  <- zeta
#     rho_best   <- rho
#     delta_best <- delta
#     Eta_best   <- Eta
#     Mu_best    <- Mu
#     ll_best    <- ll_old
#
#     # ------------------------------------------------------------
#     # Compute full proposed update once.
#     # Damping/step-halving happens only in the attempt loop.
#     # ------------------------------------------------------------
#
#     zeta_new <- update_zeta(
#       Z_list, W_list,
#       logLS_list,
#       U_list, Beta_list,
#       LV_list, Gamma_list,
#       alpha,
#       delta
#     )
#
#     rho_new <- update_rho(
#       Z_list, W_list,
#       logLS_list,
#       U_list, Beta_list,
#       LV_list, Gamma_list,
#       alpha,
#       zeta_new,
#       A
#     )
#
#     delta_new <- A %*% rho_new
#
#     alpha_new <- update_alpha_multislice(
#       Z_list, W_list,
#       logLS_list,
#       U_list, Beta_list,
#       V_list, Gamma_list,
#       zeta_new,
#       delta_new
#     )
#
#     Gamma_list_new <- update_gamma_multislice(
#       Z_list, W_list,
#       alpha_new,
#       logLS_list,
#       U_list,
#       Beta_list,
#       LV_list,
#       zeta_new,
#       delta_new,
#       lambda_gamma = lambda2
#     )
#
#     for (s in seq_len(S)) {
#       Gamma_list_new[[s]] <- cap_coeff(Gamma_list_new[[s]])
#     }
#
#     Beta_list_new <- update_beta_multislice(
#       Z_list, W_list,
#       alpha_new,
#       logLS_list,
#       U_list,
#       V_list,
#       Gamma_list_new,
#       zeta_new,
#       delta_new,
#       DTD,
#       lambdaC,
#       lambda1
#     )
#
#     for (s in seq_len(S)) {
#       Beta_list_new[[s]] <- cap_coeff(Beta_list_new[[s]])
#     }
#
#     # ------------------------------------------------------------
#     # Line search / damping.
#     # If param_damp = 1, first attempt is full update.
#     # If param_damp = 0.5, first attempt is half update.
#     # ------------------------------------------------------------
#
#     step_factor <- 1
#     accepted <- FALSE
#
#     for (attempt in 1:4) {
#
#       damp <- param_damp * step_factor
#
#       zeta_try  <- (1 - damp) * zeta  + damp * zeta_new
#       rho_try   <- (1 - damp) * rho   + damp * rho_new
#       delta_try <- A %*% rho_try
#       alpha_try <- (1 - damp) * alpha + damp * alpha_new
#
#       Gamma_list_try <- vector("list", S)
#       Beta_list_try  <- vector("list", S)
#
#       for (s in seq_len(S)) {
#         Gamma_list_try[[s]] <-
#           (1 - damp) * Gamma_list[[s]] + damp * Gamma_list_new[[s]]
#
#         Beta_list_try[[s]] <-
#           (1 - damp) * Beta_list[[s]] + damp * Beta_list_new[[s]]
#       }
#
#       for (s in seq_len(S)) {
#         rows <- (row_idx[s] + 1):row_idx[s + 1]
#
#         Eta_try[rows, ] <-
#           matrix(zeta_try, nrow = n_list[s], ncol = G, byrow = TRUE) +
#           matrix(delta_try[s, ], nrow = n_list[s], ncol = G, byrow = TRUE) +
#           alpha_try[s] * logLS_list[[s]] +
#           U_list[[s]] %*% Beta_list_try[[s]] +
#           LV_list[[s]] %*% Gamma_list_try[[s]]
#       }
#
#       Mu_try <- exp(Eta_try)
#
#       ll_try <- sum(penalized_ll(
#         Y = Y_stack,
#         Mu = Mu_try,
#         Psi = Psi_stack,
#         Beta = do.call(rbind, Beta_list_try),
#         Gamma = do.call(rbind, Gamma_list_try),
#         Cmat = Dmat,
#         lambdaC = lambdaC,
#         lambda1 = lambda1,
#         lambda2 = lambda2,
#         n_cells = N
#       ))
#
#       if (ll_try >= ll_best) {
#
#         alpha_best <- alpha_try
#         Beta_best  <- Beta_list_try
#         Gamma_best <- Gamma_list_try
#         zeta_best  <- zeta_try
#         rho_best   <- rho_try
#         delta_best <- delta_try
#
#         # Use explicit copy to avoid accidental reuse issues
#         Eta_best <- Eta_try + 0
#         Mu_best  <- Mu_try + 0
#
#         ll_best <- ll_try
#         accepted <- TRUE
#
#         if (attempt > 1) {
#           message("Step halved ", attempt - 1, " times")
#         }
#
#         break
#
#       } else {
#
#         step_factor <- step_factor / 2
#       }
#     }
#
#     if (!accepted) {
#       message("No improving inner step accepted; keeping previous parameters")
#     }
#
#     alpha <- alpha_best
#     Beta_list <- Beta_best
#     Gamma_list <- Gamma_best
#     zeta <- zeta_best
#     rho <- rho_best
#     delta <- delta_best
#     Eta <- Eta_best
#     Mu <- Mu_best
#     ll_new <- ll_best
#
#     message(sprintf(
#       "Iter %d | log-likelihood: %.2f | Mu range: [%.2f, %.2f]",
#       iter, ll_new, min(Mu), max(Mu)
#     ))
#
#     rel_change <- abs(ll_new - ll_old) / max(1, abs(ll_old))
#
#     if (iter > 1 && rel_change < tol_param) {
#       message(sprintf("Converged at inner iteration %d", iter))
#       converged <- TRUE
#       break
#     }
#
#     ll_old <- ll_new
#   }
#
#   Eta_nointercept <- matrix(0, N, G)
#
#   for (s in seq_len(S)) {
#     rows <- (row_idx[s] + 1):row_idx[s + 1]
#
#     Eta_nointercept[rows, ] <-
#       alpha[s] * logLS_list[[s]] +
#       U_list[[s]] %*% Beta_list[[s]] +
#       LV_list[[s]] %*% Gamma_list[[s]]
#   }
#
#   Mu_nointercept <- exp(Eta_nointercept)
#
#   Eta_bio <- matrix(0, N, G)
#
#   for (s in seq_len(S)) {
#     rows <- (row_idx[s] + 1):row_idx[s + 1]
#
#     Eta_bio[rows, ] <-
#       matrix(zeta, nrow = n_list[s], ncol = G, byrow = TRUE) +
#       U_list[[s]] %*% Beta_list[[s]]
#   }
#
#   Mu_bio <- exp(Eta_bio)
#
#   list(
#     alpha = alpha,
#     zeta  = zeta,
#     rho   = rho,
#     delta = delta,
#     Beta  = Beta_list,
#     Gamma = Gamma_list,
#     Mu    = Mu,
#     Mu_bio = Mu_bio,
#     Mu_nointercept = Mu_nointercept,
#     Eta   = Eta,
#     ll    = ll_old,
#     iterations = iter,
#     converged = converged
#   )
# }
