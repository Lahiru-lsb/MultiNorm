fit_all_genes_nb_multislice <- function(
    Y_list,            # list of (G × n_s) count matrices
    U_list,            # list of biological spline matrices, no intercept
    V_list,            # list of technical spline matrices
    LS_list,           # list of library size vectors
    ubar_list,         # list of niche means, Kb × M, no intercept
    lambdaC = 1e4,
    lambda1 = 0.02,
    lambda2 = 1.9,
    tol_disp = 1e-4,
    tol_param = 1e-4,
    max_outer = 30,
    max_inner = 200,
    psi_damp = 0.3,
    param_damp = 0.2,
    alpha_init = NULL,
    Beta_init  = NULL,
    Gamma_init = NULL,
    Psi_init   = NULL,
    zeta_init  = NULL,
    rho_init   = NULL,
    A = NULL
) {

  S <- length(Y_list)
  G <- nrow(Y_list[[1]])
  n_list <- sapply(Y_list, ncol)
  N <- sum(n_list)

  Kb_list <- sapply(ubar_list, nrow)
  Kb <- Kb_list[1]

  row_idx <- c(0, cumsum(n_list))

  if (is.null(A)) {
    A <- make_sumzero_contrast(S)
  }

  if (is.null(alpha_init)) {
    alpha <- rep(0, S)
  } else {
    alpha <- alpha_init
  }

  if (is.null(Beta_init)) {
    Beta_list <- lapply(seq_len(S), function(s) matrix(0.0, Kb, G))
  } else {
    Beta_list <- Beta_init
  }

  if (is.null(Gamma_init)) {
    Gamma_list <- lapply(V_list, function(V) matrix(0, ncol(V), G))
  } else {
    Gamma_list <- Gamma_init
  }

  if (is.null(Psi_init)) {
    Psi_list <- lapply(seq_len(S), function(s) rep(0, G))
  } else {
    Psi_list <- Psi_init
  }

  if (is.null(zeta_init)) {
    zeta <- rep(0.1, G)
  } else {
    zeta <- zeta_init
  }

  if (is.null(rho_init)) {
    rho <- matrix(0.1, nrow = S - 1, ncol = G)
  } else {
    rho <- rho_init
  }

  delta <- A %*% rho

  ll_old <- -Inf
  ll_trace <- numeric(max_outer)

  # ---- Build biological constraint matrix D ----

  M <- ncol(ubar_list[[1]])
  niche_names <- paste0("m", 1:M)

  niche_present <- lapply(ubar_list, function(U) {
    colSums(abs(U)) > 0
  })

  D_rows <- list()
  row_names <- character()

  for (m in seq_len(M)) {

    for (s in seq_len(S - 1)) {

      if (niche_present[[s]][m] &&
          niche_present[[s + 1]][m]) {

        blocks <- vector("list", S)

        for (j in seq_len(S)) {

          if (j == s) {
            blocks[[j]] <- t(ubar_list[[j]][, m, drop = FALSE])
          } else if (j == s + 1) {
            blocks[[j]] <- -t(ubar_list[[j]][, m, drop = FALSE])
          } else {
            blocks[[j]] <- matrix(0, 1, Kb_list[j])
          }
        }

        D_rows[[length(D_rows) + 1]] <- do.call(cbind, blocks)

        row_names <- c(
          row_names,
          sprintf("%s_s%d_vs_s%d", niche_names[m], s, s + 1)
        )
      }
    }
  }

  Dmat <- do.call(rbind, D_rows)
  rownames(Dmat) <- row_names

  DTD <- crossprod(Dmat)

  Y_stack <- do.call(rbind, lapply(Y_list, t))
  Psi_stack <- matrix(0, N, G)

  # ================= OUTER LOOP =================

  for (k in seq_len(max_outer)) {

    cat("\n=== Outer iteration", k, "===\n")

    Psi_list_old <- Psi_list

    ir <- irwls_multislice_backfit(
      Y_list = Y_list,
      U_list = U_list,
      V_list = V_list,
      LS_list = LS_list,
      Psi_list = Psi_list,
      Dmat = Dmat,
      DTD = DTD,
      lambdaC = lambdaC,
      lambda1 = lambda1,
      lambda2 = lambda2,
      max_iter = max_inner,
      param_damp = param_damp,
      tol_param = tol_param,
      alpha = alpha,
      Beta_list = Beta_list,
      Gamma_list = Gamma_list,
      zeta = zeta,
      rho = rho,
      A = A
    )

    alpha <- ir$alpha
    zeta <- ir$zeta
    rho <- ir$rho
    delta <- ir$delta
    Beta_list <- ir$Beta
    Gamma_list <- ir$Gamma

    Mu <- ir$Mu
    Mu_nointercept <- ir$Mu_nointercept

    # ---- Dispersion update ----

    Psi_new <- vector("list", S)

    for (s in seq_len(S)) {
      rows <- (row_idx[s] + 1):row_idx[s + 1]

      Psi_new[[s]] <- update_psi_edgeR(
        Ys  = Y_list[[s]],
        Mus = Mu[rows, , drop = FALSE]
      )

      # protect against invalid values
      Psi_new[[s]][!is.finite(Psi_new[[s]])] <- Psi_list_old[[s]][!is.finite(Psi_new[[s]])]

      # lower bound
      Psi_new[[s]] <- pmax(Psi_new[[s]], 1e-8)

      # winsorize on log scale
      log_psi <- log(Psi_new[[s]])

      log_med <- median(log_psi, na.rm = TRUE)
      log_mad <- mad(log_psi, na.rm = TRUE)

      log_upper <- log_med + 3 * log_mad
      log_lower <- log_med - 3 * log_mad

      Psi_new[[s]] <- exp(pmin(pmax(log_psi, log_lower), log_upper))
    }

    for (s in seq_len(S)) {
      Psi_list[[s]] <-
        (1 - psi_damp) * Psi_list_old[[s]] +
        psi_damp * Psi_new[[s]]
    }

    for (s in seq_len(S)) {
      rows <- (row_idx[s] + 1):row_idx[s + 1]

      Psi_stack[rows, ] <- matrix(
        rep(Psi_list[[s]], each = n_list[s]),
        nrow = n_list[s]
      )
    }

    ll_new <- sum(penalized_ll(
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

    if (k > 1 && ll_new < ll_old) {

      max_halving <- 3
      halving_successful <- FALSE

      for (h in seq_len(max_halving)) {

        step_size <- 0.5^h
        cat("  Step halving attempt \n")

        Psi_list_halved <- lapply(seq_len(S), function(s) {
          Psi_list_old[[s]] + step_size * (Psi_new[[s]] - Psi_list_old[[s]])
        })

        for (s in seq_len(S)) {
          rows <- (row_idx[s] + 1):row_idx[s + 1]

          Psi_stack[rows, ] <- matrix(
            rep(Psi_list_halved[[s]], each = n_list[s]),
            nrow = n_list[s]
          )
        }

        ll_halved <- sum(penalized_ll(
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

        cat("Halved likelihood =", round(ll_halved, 4), "\n")

        if (ll_halved > ll_old) {
          ll_new <- ll_halved
          Psi_list <- Psi_list_halved
          halving_successful <- TRUE
          break
        }
      }

      if (!halving_successful) {
        cat("  All step halving attempts failed \n")

        Psi_list <- Psi_list_old
        ll_new <- ll_old

        for (s in seq_len(S)) {
          rows <- (row_idx[s] + 1):row_idx[s + 1]

          Psi_stack[rows, ] <- matrix(
            rep(Psi_list[[s]], each = n_list[s]),
            nrow = n_list[s]
          )
        }
      }
    }

    cat("log-likelihood =", round(ll_new, 4), "\n")
    ll_trace[k] <- ll_new

    if (k > 1 && abs(ll_new - ll_old) / abs(ll_old) < tol_disp) {
      cat("Converged at outer iteration", k, "\n")
      break
    }

    ll_old <- ll_new
  }

  list(
    alpha = alpha,
    zeta = zeta,
    rho = rho,
    delta = delta,
    Beta = Beta_list,
    Gamma = Gamma_list,
    LS_list = LS_list,
    Psi = Psi_list,
    Eta = ir$Eta,
    Mu = Mu,
    Mu_bio = ir$Mu_bio,
    Mu_nointercept = ir$Mu_nointercept,
    U_list = U_list,
    V_list = V_list,
    likelihood = ll_trace[seq_len(k)],
    Dmat = Dmat,
    constraints = Dmat %*% do.call(rbind, Beta_list),
    A = A,
    iterations = k
  )
}


# fit_all_genes_nb_multislice_1 <- function(
#     Y_list,            # list of (G × n_s) count matrices
#     U_list,            # list of biological spline matrices, no intercept
#     V_list,            # list of technical spline matrices
#     LS_list,           # list of library size vectors
#     ubar_list,         # list of niche means, Kb × M, no intercept
#     lambdaC = 1e4,
#     lambda1 = 0.02,
#     lambda2 = 1.9,
#     tol_disp = 1e-4,
#     tol_param = 1e-4,
#     max_outer = 30,
#     max_inner = 200,
#     psi_damp = 0.3,
#     param_damp = 0.2,
#     alpha_init = NULL,
#     Beta_init  = NULL,
#     Gamma_init = NULL,
#     Psi_init   = NULL,
#     zeta_init  = NULL,
#     rho_init   = NULL,
#     A = NULL
# ) {
#
#   S <- length(Y_list)
#   G <- nrow(Y_list[[1]])
#   n_list <- sapply(Y_list, ncol)
#   N <- sum(n_list)
#
#   Kb_list <- sapply(ubar_list, nrow)
#   Kb <- Kb_list[1]
#
#   row_idx <- c(0, cumsum(n_list))
#
#   if (is.null(A)) {
#     A <- make_sumzero_contrast(S)
#   }
#
#   if (is.null(alpha_init)) {
#     alpha <- rep(0, S)
#   } else {
#     alpha <- alpha_init
#   }
#
#   if (is.null(Beta_init)) {
#     Beta_list <- lapply(seq_len(S), function(s) matrix(0.0, Kb, G))
#   } else {
#     Beta_list <- Beta_init
#   }
#
#   if (is.null(Gamma_init)) {
#     Gamma_list <- lapply(V_list, function(V) matrix(0, ncol(V), G))
#   } else {
#     Gamma_list <- Gamma_init
#   }
#
#   if (is.null(Psi_init)) {
#     Psi_list <- lapply(seq_len(S), function(s) rep(0, G))
#   } else {
#     Psi_list <- Psi_init
#   }
#
#   if (is.null(zeta_init)) {
#     zeta <- rep(0.1, G)
#   } else {
#     zeta <- zeta_init
#   }
#
#   if (is.null(rho_init)) {
#     rho <- matrix(0.1, nrow = S - 1, ncol = G)
#   } else {
#     rho <- rho_init
#   }
#
#   delta <- A %*% rho
#
#   logLS_list <- lapply(LS_list, log)
#
#   LV_list <- vector("list", S)
#   for (s in seq_len(S)) {
#     LV_list[[s]] <- logLS_list[[s]] * V_list[[s]]
#   }
#
#   # ---- helper: compute Eta and Mu from current regression parameters ----
#   compute_eta_mu <- function(alpha, Beta_list, Gamma_list, zeta, delta) {
#
#     Eta <- matrix(0, N, G)
#
#     for (s in seq_len(S)) {
#       rows <- (row_idx[s] + 1):row_idx[s + 1]
#
#       Eta[rows, ] <-
#         matrix(zeta, nrow = n_list[s], ncol = G, byrow = TRUE) +
#         matrix(delta[s, ], nrow = n_list[s], ncol = G, byrow = TRUE) +
#         alpha[s] * logLS_list[[s]] +
#         U_list[[s]] %*% Beta_list[[s]] +
#         LV_list[[s]] %*% Gamma_list[[s]]
#     }
#
#     list(
#       Eta = Eta,
#       Mu = exp(Eta)
#     )
#   }
#
#   # ---- initial Mu needed for first dispersion update ----
#   em0 <- compute_eta_mu(
#     alpha = alpha,
#     Beta_list = Beta_list,
#     Gamma_list = Gamma_list,
#     zeta = zeta,
#     delta = delta
#   )
#
#   Eta <- em0$Eta
#   Mu  <- em0$Mu
#
#   ll_old <- -Inf
#   ll_trace <- numeric(max_outer)
#
#   # ---- Build biological constraint matrix D ----
#
#   M <- ncol(ubar_list[[1]])
#   niche_names <- paste0("m", 1:M)
#
#   niche_present <- lapply(ubar_list, function(U) {
#     colSums(abs(U)) > 0
#   })
#
#   D_rows <- list()
#   row_names <- character()
#
#   for (m in seq_len(M)) {
#
#     for (s in seq_len(S - 1)) {
#
#       if (niche_present[[s]][m] &&
#           niche_present[[s + 1]][m]) {
#
#         blocks <- vector("list", S)
#
#         for (j in seq_len(S)) {
#
#           if (j == s) {
#             blocks[[j]] <- t(ubar_list[[j]][, m, drop = FALSE])
#           } else if (j == s + 1) {
#             blocks[[j]] <- -t(ubar_list[[j]][, m, drop = FALSE])
#           } else {
#             blocks[[j]] <- matrix(0, 1, Kb_list[j])
#           }
#         }
#
#         D_rows[[length(D_rows) + 1]] <- do.call(cbind, blocks)
#
#         row_names <- c(
#           row_names,
#           sprintf("%s_s%d_vs_s%d", niche_names[m], s, s + 1)
#         )
#       }
#     }
#   }
#
#   Dmat <- do.call(rbind, D_rows)
#   rownames(Dmat) <- row_names
#
#   DTD <- crossprod(Dmat)
#
#   Y_stack <- do.call(rbind, lapply(Y_list, t))
#   Psi_stack <- matrix(0, N, G)
#
#   # ================= OUTER LOOP =================
#
#   for (k in seq_len(max_outer)) {
#
#     cat("\n=== Outer iteration", k, "===\n")
#
#
#     # 1. UPDATE DISPERSION FIRST, USING CURRENT ACCEPTED Mu
#
#
#     Psi_list_old <- Psi_list
#     Psi_new <- vector("list", S)
#
#     for (s in seq_len(S)) {
#       rows <- (row_idx[s] + 1):row_idx[s + 1]
#
#       Psi_new[[s]] <- update_psi_edgeR(
#         Ys  = Y_list[[s]],
#         Mus = Mu[rows, , drop = FALSE]
#       )
#     }
#
#     for (s in seq_len(S)) {
#       Psi_list[[s]] <-
#         (1 - psi_damp) * Psi_list_old[[s]] +
#         psi_damp * Psi_new[[s]]
#     }
#
#     # Optional: build Psi_stack only for baseline reporting.
#     # irwls_multislice_backfit() will rebuild internally.
#     for (s in seq_len(S)) {
#       rows <- (row_idx[s] + 1):row_idx[s + 1]
#
#       Psi_stack[rows, ] <- matrix(
#         rep(Psi_list[[s]], each = n_list[s]),
#         nrow = n_list[s]
#       )
#     }
#
#     ll_baseline <- sum(penalized_ll(
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
#     cat("baseline after Psi update =", round(ll_baseline, 4), "\n")
#
#
#     # 2. RUN INNER IRWLS/BACKFITTING WITH UPDATED Psi FIXED
#
#
#     ir <- irwls_multislice_backfit(
#       Y_list = Y_list,
#       U_list = U_list,
#       V_list = V_list,
#       LS_list = LS_list,
#       Psi_list = Psi_list,
#       Dmat = Dmat,
#       DTD = DTD,
#       lambdaC = lambdaC,
#       lambda1 = lambda1,
#       lambda2 = lambda2,
#       max_iter = max_inner,
#       param_damp = param_damp,
#       tol_param = tol_param,
#       alpha = alpha,
#       Beta_list = Beta_list,
#       Gamma_list = Gamma_list,
#       zeta = zeta,
#       rho = rho,
#       A = A,
#       ll_start = ll_baseline
#     )
#
#
#     # 3. ACCEPT INNER-LOOP RESULTS
#
#
#     alpha <- ir$alpha
#     zeta <- ir$zeta
#     rho <- ir$rho
#     delta <- ir$delta
#     Beta_list <- ir$Beta
#     Gamma_list <- ir$Gamma
#
#     Eta <- ir$Eta
#     Mu <- ir$Mu
#     Mu_nointercept <- ir$Mu_nointercept
#
#     ll_new <- ir$ll
#
#     cat("log-likelihood =", round(ll_new, 4), "\n")
#     ll_trace[k] <- ll_new
#
#
#     # 4. OUTER CONVERGENCE CHECK
#
#
#     if (k > 1 && abs(ll_new - ll_old) / max(1, abs(ll_old)) < tol_disp) {
#       cat("Converged at outer iteration", k, "\n")
#       break
#     }
#
#     ll_old <- ll_new
#   }
#
#   list(
#     alpha = alpha,
#     zeta = zeta,
#     rho = rho,
#     delta = delta,
#     Beta = Beta_list,
#     Gamma = Gamma_list,
#     LS_list = LS_list,
#     Psi = Psi_list,
#     Eta = Eta,
#     Mu = Mu,
#     Mu_bio = ir$Mu_bio,
#     Mu_nointercept = Mu_nointercept,
#     U_list = U_list,
#     V_list = V_list,
#     likelihood = ll_trace[seq_len(k)],
#     Dmat = Dmat,
#     constraints = Dmat %*% do.call(rbind, Beta_list),
#     A = A,
#     iterations = k
#   )
# }
