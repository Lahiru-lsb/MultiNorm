update_gamma_multislice <- function(
    Z_list,
    w_list,
    alpha,
    logLS_list,
    Ubeta_list,
    Beta_list,
    LV_list,
    lambda_gamma
) {

  S  <- length(Z_list)
  G  <- ncol(Z_list[[1]])
  Gamma_list <- vector("list", S)

  for (s in seq_len(S)) {

    # Residual for slice s
    R_s <-
      Z_list[[s]] -
      alpha[s] * logLS_list[[s]] -
      Ubeta_list[[s]] %*% Beta_list[[s]]

    V_s <- LV_list[[s]]
    w_s <- w_list[[s]]

    # Apply weights
    W_s <- scale_rows(V_s, w_s)

    # Ridge system
    A_s <- crossprod(W_s, V_s)
    diag(A_s) <- diag(A_s) + lambda_gamma

    rhs_s <- crossprod(W_s, R_s)

    # Solve via Cholesky
    Rchol <- chol(A_s)
    Gamma_list[[s]] <-
      backsolve(Rchol, forwardsolve(t(Rchol), rhs_s))
  }

  Gamma_list
}

update_alpha_multislice <- function(
    Z_list,        # list of (n_s × G) matrices
    W_list,        # list of (n_s × G) matrices OR (n_s) vectors
    logLS_list,    # list of length-n_s vectors
    Ubeta_list,    # list of (n_s × Kb+1) matrices
    Beta_list,     # list of (Kb × G) matrices
    V_list,        # list of (n_s × Kt) matrices
    Gamma_list     # list of (Kt × G) matrices
) {
  S <- length(Z_list)
  G <- ncol(Z_list[[1]])
  alpha_new <- numeric(S)

  for (s in seq_len(S)) {

    Zs    <- Z_list[[s]]
    Ws    <- W_list[[s]]        # length n
    logLS <- logLS_list[[s]]    # length n

    UB  <- Ubeta_list[[s]] %*% Beta_list[[s]]
    LVG <- (V_list[[s]] *logLS)     %*% Gamma_list[[s]]

    # collapse rows FIRST
    row_resid_sum <- rowSums(Zs - UB - LVG)

    # now weight by Ws * logLS
    num <- sum(Ws * logLS * row_resid_sum)

    den <- G * sum(Ws * (logLS^2))

    alpha_new[s] <- num / den
  }

  alpha_new
}



update_beta_multislice <- function(
    Z_list,        # list of (n_s × G)
    w_list,        # list of length-n_s vectors
    alpha,         # length S
    logLS_list,    # list of vectors
    Ubeta_list,    # list of (n_s × Kb+1)
    V_list,        # list of (n_s × Kt)
    Gamma_list,    # list of (Kt × G)
    DTD,           # constraint penalty matrix
    Pbeta,         # biological spline penalty
    lambdaC,
    lambda1
) {

  S  <- length(Z_list)
  G  <- ncol(Z_list[[1]])
  Kb <- ncol(Ubeta_list[[1]])

  # ---- build stacked residual ----
  R_stack <- vector("list", S)
  WU_stack <- vector("list", S)
  U_stack <- vector("list", S)

  for (s in seq_len(S)) {

    LV <- scale_rows(V_list[[s]], logLS_list[[s]])

    R_stack[[s]] <-
      Z_list[[s]] -
      alpha[s] * logLS_list[[s]] -
      LV %*% Gamma_list[[s]]

    WU_stack[[s]] <- scale_rows(Ubeta_list[[s]], w_list[[s]])
    U_stack[[s]]  <- Ubeta_list[[s]]
  }

  # ---- block matrices ----
  U_big  <- Matrix::bdiag(U_stack)
  WU_big <- Matrix::bdiag(WU_stack)

  # ---- normal equations ----
  A <-
    crossprod(WU_big, U_big) +
    lambdaC * DTD +
    lambda1 * Pbeta

  rhs <- crossprod(WU_big, do.call(rbind, R_stack))

  # ---- solve ----
  #Beta_stack <- as.matrix(solve(A, rhs))
  #Beta_stack <- as.matrix(Matrix::solve(A, rhs, sparse = TRUE))
  Rchol <- chol(A)
  Beta_stack <- backsolve(Rchol, forwardsolve(t(Rchol), rhs))

  # ---- split back into slices ----
  Beta_list <- vector("list", S)
  idx <- 0
  for (s in seq_len(S)) {
    Beta_list[[s]] <- Beta_stack[(idx + 1):(idx + Kb), , drop = FALSE]
    idx <- idx + Kb
  }

  Beta_list
}




