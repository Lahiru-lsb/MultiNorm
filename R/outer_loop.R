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

  # ---- Pre-compute row indices ----
  row_idx <- c(0, cumsum(n_list))


  # ---- Initialize alpha ----
  if (is.null(alpha_init)) {
    alpha <- rep(0, S)
  } else {
    #message("Using user-provided alpha")
    alpha <- alpha_init
  }

  # ---- Initialize Psi ----
  if (is.null(Psi_init)) {
    Psi_list <- lapply(seq_len(S), function(s) rep(0.1, G))
  } else {
    #message("Using user-provided psi")
    Psi_list <- Psi_init
  }

  ll_old <- -Inf
  ll_trace <- numeric(max_outer)

  # ---- Build full constraint matrix D ----
  Kb_list <- sapply(ubar_list, nrow)
  M <- ncol(ubar_list[[1]])

  niche_present <- lapply(ubar_list, function(U)
    colSums(abs(U)) > 0
  )

  D_rows <- list()

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

      }

    }

  }


  Dmat <- do.call(rbind, D_rows)
  DTD <- crossprod(Dmat)

  # ---- Biological spline penalty ----
  Kb <- Kb_list[1]  # Assuming all same size
  Pblock <- diag(c(0, rep(1, Kb - 1)))
  Pbeta  <- as.matrix(Matrix::bdiag(replicate(S, Pblock, simplify = FALSE)))

  # ---- Pre-stack Y once (doesn't change) ----
  Y_stack <- do.call(rbind, lapply(Y_list, t))

  # ---- Pre-allocate Psi_stack ----
  Psi_stack <- matrix(0, N, G)

  # ================= OUTER LOOP =================
  for (k in seq_len(max_outer)) {

    cat("\n=== Outer iteration", k, "===\n")

    # ---- IRWLS ----
    ir <- irwls_multislice_backfit(
      Y_list = Y_list,
      U_list = U_list,
      V_list = V_list,
      LS_list = LS_list,
      alpha_init = alpha,
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
      Beta_init = Beta_init,
      Gamma_init = Gamma_init
    )

    alpha <- ir$alpha
    Beta_list  <- ir$Beta
    Gamma_list <- ir$Gamma
    Beta_init  <- Beta_list
    Gamma_init <- Gamma_list

    Mu <- ir$Mu

    # ---- Dispersion update ----
    Psi_new <- vector("list", S)

    for (s in seq_len(S)) {
      rows <- (row_idx[s] + 1):row_idx[s + 1]
      Psi_new[[s]] <- update_psi_edgeR(
        Ys  = Y_list[[s]],
        Mus = Mu[rows, , drop = FALSE]
      )
    }

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
    ll_new <- sum(ll_gene(
      Y = Y_stack,
      Mu = Mu,
      Psi = Psi_stack,
      n_cells = N
    ))

    cat("log-likelihood =", round(ll_new,4), "\n")
    ll_trace[k] <- ll_new

    # ---- Convergence check ----
    if (k > 1 && abs(ll_new - ll_old) / (abs(ll_old) + 1) < tol_disp) {
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
