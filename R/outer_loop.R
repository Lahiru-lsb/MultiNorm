#main code inner +outer optimized
fit_all_genes_nb_multislice <- function(
    Y_list,            # list of (G × n_s) count matrices
    U_list,            # list of biological spline matrices
    V_list,            # list of technical spline matrices
    LS_list,           # list of library size vectors
    ubar_list,         # list of niche means
    lambdaC = 1e4,
    lambda1 = 0.02,
    lambda2 = 1.9,
    tol_disp = 1e-4,
    tol_param = 1e-4,
    max_outer = 30,
    max_inner = 200,
    psi_damp = 0.3,
    param_damp=0.2,
    alpha_init = NULL,
    Beta_init  = NULL,
    Gamma_init = NULL,
    Psi_init   = NULL
) {

  S <- length(Y_list)
  G <- nrow(Y_list[[1]])
  n_list <- sapply(Y_list, ncol)
  N <- sum(n_list)
  Kb_list <- sapply(ubar_list, nrow)
  Kb <- Kb_list[1]  # Assuming all same size

  # ---- Pre-compute row indices ----
  row_idx <- c(0, cumsum(n_list))


  # ---- Initialize  ----
  if (is.null(alpha_init)) {
    alpha <- rep(0, S)
  } else {
    #message("Using user-provided alpha")
    alpha <- alpha_init
  }

  if (is.null(Beta_init)) {
    Beta_list <- lapply(seq_len(S), function(s) matrix(0.0, Kb, G))
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


  # ---- Initialize Psi ----
  if (is.null(Psi_init)) {
    Psi_list <- lapply(seq_len(S), function(s) rep(0, G))
  } else {
    #message("Using user-provided psi")
    Psi_list <- Psi_init
  }

  ll_old <- -Inf
  ll_trace <- numeric(max_outer)

  # ---- Build full constraint matrix D ----

  M <- ncol(ubar_list[[1]])
  niche_names <- paste0("m", 1:M)

  niche_present <- lapply(ubar_list, function(U)
    colSums(abs(U)) > 0
  )

  D_rows <- list()
  row_names <- character()

  for (m in seq_len(M)) {           # niche first

    for (s in seq_len(S-1)) {       # then slice pairs

      if (niche_present[[s]][m] &&
          niche_present[[s+1]][m]) {

        blocks <- vector("list", S)

        for (j in seq_len(S)) {

          if (j == s) {

            blocks[[j]] <- t(ubar_list[[j]][, m, drop = FALSE])

          } else if (j == s+1) {

            blocks[[j]] <- -t(ubar_list[[j]][, m, drop = FALSE])

          } else {

            blocks[[j]] <- matrix(0,1,Kb_list[j])

          }

        }

        D_rows[[length(D_rows)+1]] <- do.call(cbind, blocks)
        row_names <- c(
          row_names,
          sprintf("%s_s%d_vs_s%d", niche_names[m], s, s+1)
        )

      }

    }

  }

  Dmat <- do.call(rbind, D_rows)
  rownames(Dmat) <- row_names

  DTD <- crossprod(Dmat)

  # ---- Biological spline penalty ----

  Pblock <- diag(c(0, rep(1, Kb - 1)))
  Pbeta  <- as.matrix(Matrix::bdiag(replicate(S, Pblock, simplify = FALSE)))

  # ---- Pre-stack Y once (doesn't change) ----
  Y_stack <- do.call(rbind, lapply(Y_list, t))

  # ---- Pre-allocate Psi_stack ----
  Psi_stack <- matrix(0, N, G)

  # ================= OUTER LOOP =================
  for (k in seq_len(max_outer)) {

    cat("\n=== Outer iteration", k, "===\n")

    # ---- ADDED: Store previous Psi for step halving (dispersion only) ----
    Psi_list_old <- Psi_list

    # ---- IRWLS ----
    ir <- irwls_multislice_backfit(
      Y_list = Y_list,
      U_list = U_list,
      V_list = V_list,
      LS_list = LS_list,
      Psi_list = Psi_list,
      Dmat = Dmat,
      DTD = DTD,
      Pbeta = Pbeta,
      lambdaC = lambdaC,
      lambda1 = lambda1,
      lambda2 = lambda2,
      max_iter = max_inner,
      param_damp=param_damp,
      tol_param = tol_param,
      alpha = alpha,
      Beta_list = Beta_list,
      Gamma_list = Gamma_list
    )

    alpha <- ir$alpha
    Beta_list  <- ir$Beta
    Gamma_list <- ir$Gamma


    Mu <- ir$Mu
    Mu_nointercept<-ir$Mu_nointercept

    # ---- Dispersion update ----
    Psi_new <- vector("list", S)

    for (s in seq_len(S)) {
      rows <- (row_idx[s] + 1):row_idx[s + 1]
      Psi_new[[s]] <- update_psi_edgeR(
        Ys  = Y_list[[s]],
        Mus = Mu_nointercept[rows, , drop = FALSE] # set mu without intercept
      )
    }

    # compute_one <- function(s, Y_list, Mu, row_idx) {
    #   rows <- (row_idx[s] + 1):row_idx[s + 1]
    #   update_psi_edgeR(
    #     Ys  = Y_list[[s]],
    #     Mus = Mu[rows, , drop = FALSE]
    #   )
    # }
    #
    # Psi_new <- future.apply::future_lapply(
    #   seq_len(S),
    #   compute_one,
    #   Y_list = Y_list,
    #   Mu = Mu,
    #   row_idx = row_idx,
    #   future.seed = TRUE
    # )
    # ---- Damped update (in-place) ----
    for (s in seq_len(S)) {
      Psi_list[[s]] <- (1 - psi_damp) * Psi_list[[s]] +  psi_damp* Psi_new[[s]]
    }

    # ---- Update Psi_stack in-place ----
    for (s in seq_len(S)) {
      rows <- (row_idx[s] + 1):row_idx[s + 1]
      Psi_stack[rows, ] <- matrix(rep(Psi_list[[s]], each = n_list[s]),
                                  nrow = n_list[s])
    }

    # ---- Likelihood ----
    # ll_new <- sum(ll_gene(
    #   Y = Y_stack,
    #   Mu = Mu,
    #   Psi = Psi_stack,
    #   n_cells = N
    # ))

    # ---- Penalized likelihood (aligned with inner loop) ----
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

    # ---- ADDED: Step halving for dispersion only if likelihood decreased ----
    if (k > 1 && ll_new < ll_old) {
      #cat("Likelihood decreased! Attempting step halving for dispersion...\n")

      max_halving <- 3
      halving_successful <- FALSE

      for (h in seq_len(max_halving)) {
        step_size <- 0.25^h
        cat("  Step halving attempt \n")

        # Interpolate only Psi between old and new values
        Psi_list_halved <- lapply(seq_len(S), function(s) {
          Psi_list_old[[s]] + step_size * (Psi_new[[s]] - Psi_list_old[[s]])
        })

        # Update Psi_stack with halved values
        for (s in seq_len(S)) {
          rows <- (row_idx[s] + 1):row_idx[s + 1]
          Psi_stack[rows, ] <- matrix(rep(Psi_list_halved[[s]], each = n_list[s]),
                                      nrow = n_list[s])
        }

        # Recompute likelihood (Mu unchanged, only Psi halved)
        ll_halved <- sum(penalized_ll(
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

        cat("Halved likelihood =", round(ll_halved, 4), "\n")

        if (ll_halved > ll_old) {
          #cat("  Step halving successful!\n")
          ll_new <- ll_halved
          Psi_list <- Psi_list_halved
          halving_successful <- TRUE
          break
        }
      }

      # If all step halving attempts failed, revert Psi to previous values
      if (!halving_successful) {
        cat("  All step halving attempts failed \n")
        Psi_list <- Psi_list_old
        ll_new <- ll_old

        # Restore Psi_stack with old values
        for (s in seq_len(S)) {
          rows <- (row_idx[s] + 1):row_idx[s + 1]
          Psi_stack[rows, ] <- matrix(rep(Psi_list[[s]], each = n_list[s]),
                                      nrow = n_list[s])
        }
      }
    }

    cat("log-likelihood =", round(ll_new,4), "\n")
    ll_trace[k] <- ll_new

    # ---- Convergence check ----
    if (k > 1 && abs(ll_new - ll_old) / (abs(ll_old) ) < tol_disp) {
      cat("Converged at outer iteration", k, "\n")
      break
    }

    ll_old <- ll_new
  }

  list(
    alpha = alpha,
    Beta  = Beta_list,
    Gamma = Gamma_list,
    LS_list =LS_list,
    Psi   = Psi_list,
    Eta = ir$Eta,
    Mu    = Mu,
    U_list = U_list,
    V_list = V_list,
    likelihood = ll_trace[seq_len(k)],
    Dmat=Dmat,
    constraints = Dmat %*% do.call(rbind, Beta_list),
    iterations = k
  )
}





# fit_all_genes_nb_multislice1 <- function(
#     Y_list,            # list of (G × n_s) count matrices
#     U_list,            # list of biological spline matrices
#     V_list,            # list of technical spline matrices
#     LS_list,           # list of library size vectors
#     ubar_list,         # list of niche means
#     lambdaC = 1e4,
#     lambda1 = 0.02,
#     lambda2 = 1.9,
#     tol_disp = 1e-4,
#     tol_param = 1e-4,
#     max_outer = 30,
#     max_inner = 200,
#     psi_damp = 0.3,
#     param_damp=0.2,
#     alpha_init = NULL,
#     Beta_init  = NULL,
#     Gamma_init = NULL,
#     Psi_init   = NULL
# ) {
#
#   S <- length(Y_list)
#   G <- nrow(Y_list[[1]])
#   n_list <- sapply(Y_list, ncol)
#   N <- sum(n_list)
#   Kb_list <- sapply(ubar_list, nrow)
#   Kb <- Kb_list[1]  # Assuming all same size
#
#   # ---- Pre-compute row indices ----
#   row_idx <- c(0, cumsum(n_list))
#
#
#   # ---- Initialize  ----
#   if (is.null(alpha_init)) {
#     alpha <- rep(0, S)
#   } else {
#     #message("Using user-provided alpha")
#     alpha <- alpha_init
#   }
#
#   if (is.null(Beta_init)) {
#     Beta_list <- lapply(seq_len(S), function(s) matrix(0.01, Kb, G))
#   } else {
#     #message("Using user-provided beta")
#     Beta_list <- Beta_init
#   }
#
#   if (is.null(Gamma_init)) {
#     Gamma_list <- lapply(V_list, function(V) matrix(0, ncol(V), G))
#   } else {
#     #message("Using user-provided gamma")
#     Gamma_list <- Gamma_init
#   }
#
#
#   # ---- Initialize Psi ----
#   if (is.null(Psi_init)) {
#     Psi_list <- lapply(seq_len(S), function(s) rep(0.1, G))
#   } else {
#     #message("Using user-provided psi")
#     Psi_list <- Psi_init
#   }
#
#   ll_old <- -Inf
#   ll_trace <- numeric(max_outer)
#
#   # ---- Build full constraint matrix D ----
#
#   M <- ncol(ubar_list[[1]])
#   #niche names
#   niche_names <- paste0("m", 1:M)
#
#   niche_present <- lapply(ubar_list, function(U)
#     colSums(abs(U)) > 0
#   )
#
#   D_rows <- list()
#
#   for (m in seq_len(M)) {           # niche first
#
#     for (s in seq_len(S-1)) {       # then slice pairs
#
#       if (niche_present[[s]][m] &&
#           niche_present[[s+1]][m]) {
#
#         blocks <- vector("list", S)
#
#         for (j in seq_len(S)) {
#
#           if (j == s) {
#
#             blocks[[j]] <- t(ubar_list[[j]][, m, drop = FALSE])
#
#           } else if (j == s+1) {
#
#             blocks[[j]] <- -t(ubar_list[[j]][, m, drop = FALSE])
#
#           } else {
#
#             blocks[[j]] <- matrix(0,1,Kb_list[j])
#
#           }
#
#         }
#
#         D_rows[[length(D_rows)+1]] <- do.call(cbind, blocks)
#
#       }
#
#     }
#
#   }
#
#
#   Dmat <- do.call(rbind, D_rows)
#   DTD <- crossprod(Dmat)
#
#   # ---- Biological spline penalty ----
#
#   Pblock <- diag(c(0, rep(1, Kb - 1)))
#   Pbeta  <- as.matrix(Matrix::bdiag(replicate(S, Pblock, simplify = FALSE)))
#
#   # ---- Pre-stack Y once (doesn't change) ----
#   Y_stack <- do.call(rbind, lapply(Y_list, t))
#
#   # ---- Pre-allocate Psi_stack ----
#   Psi_stack <- matrix(0, N, G)
#
#   # ================= OUTER LOOP =================
#   for (k in seq_len(max_outer)) {
#
#     cat("\n=== Outer iteration", k, "===\n")
#
#     # ---- IRWLS ----
#     ir <- irwls_multislice_backfit(
#       Y_list = Y_list,
#       U_list = U_list,
#       V_list = V_list,
#       LS_list = LS_list,
#       Psi_list = Psi_list,
#       Dmat = Dmat,
#       DTD = DTD,
#       Pbeta = Pbeta,
#       lambdaC = lambdaC,
#       lambda1 = lambda1,
#       lambda2 = lambda2,
#       max_iter = max_inner,
#       param_damp=param_damp,
#       tol_param = tol_param,
#       alpha = alpha,
#       Beta_list = Beta_list,
#       Gamma_list = Gamma_list
#     )
#
#     alpha <- ir$alpha
#     Beta_list  <- ir$Beta
#     Gamma_list <- ir$Gamma
#
#
#     Mu <- ir$Mu
#
#     # ---- Dispersion update ----
#     Psi_new <- vector("list", S)
#
#     for (s in seq_len(S)) {
#       rows <- (row_idx[s] + 1):row_idx[s + 1]
#       Psi_new[[s]] <- update_psi_edgeR(
#         Ys  = Y_list[[s]],
#         Mus = Mu[rows, , drop = FALSE]
#       )
#     }
#
#     # ---- Damped update (in-place) ----
#     for (s in seq_len(S)) {
#       Psi_list[[s]] <- (1 - psi_damp) * Psi_list[[s]] +  psi_damp* Psi_new[[s]]
#     }
#
#     # ---- Update Psi_stack in-place ----
#     for (s in seq_len(S)) {
#       rows <- (row_idx[s] + 1):row_idx[s + 1]
#       Psi_stack[rows, ] <- matrix(rep(Psi_list[[s]], each = n_list[s]),
#                                   nrow = n_list[s])
#     }
#
#     # ---- Likelihood ----
#     # ll_new <- sum(ll_gene(
#     #   Y = Y_stack,
#     #   Mu = Mu,
#     #   Psi = Psi_stack,
#     #   n_cells = N
#     # ))
#
#     # ---- Penalized likelihood (aligned with inner loop) ----
#     ll_new <- sum(penalized_ll(
#       Y = Y_stack,
#       Mu = Mu,
#       Psi = Psi_stack,
#       Beta = do.call(rbind, Beta_list),
#       Gamma = do.call(rbind, Gamma_list),
#       Cmat = Dmat,
#       Pbeta = Pbeta,
#       lambdaC = lambdaC,
#       lambda1 = lambda1,
#       lambda2 = lambda2,
#       n_cells = N
#     ))
#
#     cat("log-likelihood =", round(ll_new,4), "\n")
#     ll_trace[k] <- ll_new
#
#     # ---- Convergence check ----
#     if (k > 1 && abs(ll_new - ll_old) / (abs(ll_old) ) < tol_disp) {
#       cat("Converged at outer iteration", k, "\n")
#       break
#     }
#
#     ll_old <- ll_new
#   }
#
#   list(
#     alpha = alpha,
#     Beta  = Beta_list,
#     Gamma = Gamma_list,
#     LS_list =LS_list,
#     Psi   = Psi_list,
#     Eta = ir$Eta,
#     Mu    = Mu,
#     U_list = U_list,
#     V_list = V_list,
#     likelihood = ll_trace[seq_len(k)],
#     Dmat=Dmat,
#     constraints = Dmat %*% do.call(rbind, Beta_list),
#     iterations = k
#   )
# }
