update_gamma_multislice <- function(
    Z_list,
    w_list,
    alpha,
    logLS_list,
    U_list,
    Beta_list,
    LV_list,
    zeta,
    delta,
    lambda_gamma
) {

  S  <- length(Z_list)
  G  <- ncol(Z_list[[1]])
  Gamma_list <- vector("list", S)

  for (s in seq_len(S)) {

    Zs    <- Z_list[[s]]
    ws    <- w_list[[s]]
    logLS <- logLS_list[[s]]
    n_s   <- nrow(Zs)

    # Residual excluding gamma
    R_s <-
      Zs -
      matrix(zeta, nrow = n_s, ncol = G, byrow = TRUE) -
      matrix(delta[s, ], nrow = n_s, ncol = G, byrow = TRUE) -
      alpha[s] * logLS -
      U_list[[s]] %*% Beta_list[[s]]

    V_s <- LV_list[[s]]

    # Apply common slice weights
    W_s <- scale_rows(V_s, ws)

    A_s <- crossprod(W_s, V_s)
    diag(A_s) <- diag(A_s) + lambda_gamma

    rhs_s <- crossprod(W_s, R_s)

    Rchol <- chol(A_s)
    Gamma_list[[s]] <-
      backsolve(Rchol, forwardsolve(t(Rchol), rhs_s))
  }

  Gamma_list
}

update_alpha_multislice <- function(
    Z_list,        # list of (n_s × G) matrices
    W_list,        # list of (n_s × G) matrices OR length-n_s vectors
    logLS_list,    # list of length-n_s vectors
    U_list,        # list of (n_s × Kb) biological spline matrices, no intercept
    Beta_list,     # list of (Kb × G) matrices
    V_list,        # list of (n_s × Kt) matrices
    Gamma_list,    # list of (Kt × G) matrices
    zeta,          # length-G vector
    delta          # S × G matrix
) {
  S <- length(Z_list)
  G <- ncol(Z_list[[1]])
  alpha_new <- numeric(S)

  for (s in seq_len(S)) {

    Zs    <- Z_list[[s]]        # n_s × G
    Ws    <- W_list[[s]]        # length n_s OR n_s × G
    logLS <- logLS_list[[s]]    # length n_s
    n_s   <- nrow(Zs)

    UB  <- U_list[[s]] %*% Beta_list[[s]]          # n_s × G
    LVG <- (V_list[[s]] * logLS) %*% Gamma_list[[s]] # n_s × G

    R <- Zs -
      matrix(zeta, nrow = n_s, ncol = G, byrow = TRUE) -
      matrix(delta[s, ], nrow = n_s, ncol = G, byrow = TRUE) -
      UB -
      LVG

    row_resid_sum <- rowSums(R)

    num <- sum(Ws * logLS * row_resid_sum)
    den <- G * sum(Ws * logLS^2)

    alpha_new[s] <- num / pmax(den, 1e-8)
  }

  alpha_new
}



update_beta_multislice <- function(
    Z_list,        # list of (n_s × G) matrices
    w_list,        # list of length-n_s vectors OR (n_s × G) matrices
    alpha,         # length-S vector
    logLS_list,    # list of length-n_s vectors
    U_list,        # list of (n_s × Kb) biological spline matrices, no intercept
    V_list,        # list of (n_s × Kt) matrices
    Gamma_list,    # list of (Kt × G) matrices
    zeta,          # length-G vector
    delta,         # S × G matrix
    DTD,           # SKb × SKb constraint penalty matrix
    lambdaC,
    lambda1
) {

  S  <- length(Z_list)
  G  <- ncol(Z_list[[1]])
  Kb <- ncol(U_list[[1]])

  R_stack  <- vector("list", S)
  WU_stack <- vector("list", S)
  U_stack  <- vector("list", S)

  for (s in seq_len(S)) {

    Zs    <- Z_list[[s]]        # n_s × G
    ws    <- w_list[[s]]        # length n_s OR n_s × G
    logLS <- logLS_list[[s]]    # length n_s
    n_s   <- nrow(Zs)

    LV <- scale_rows(V_list[[s]], logLS)  # n_s × Kt

    R_stack[[s]] <-
      Zs -
      matrix(zeta, nrow = n_s, ncol = G, byrow = TRUE) -
      matrix(delta[s, ], nrow = n_s, ncol = G, byrow = TRUE) -
      alpha[s] * logLS -
      LV %*% Gamma_list[[s]]

    if (is.matrix(ws)) {
      stop("Current beta update assumes common length-n_s weights. Use common weights or implement gene-specific block update.")
    }

    WU_stack[[s]] <- scale_rows(U_list[[s]], ws)
    U_stack[[s]]  <- U_list[[s]]
  }

  U_big  <- Matrix::bdiag(U_stack)
  WU_big <- Matrix::bdiag(WU_stack)

  A <-
    crossprod(WU_big, U_big) +
    lambdaC * DTD +
    lambda1 * Matrix::Diagonal(n = S * Kb)

  rhs <- crossprod(WU_big, do.call(rbind, R_stack))

  Rchol <- chol(A)
  Beta_stack <- backsolve(Rchol, forwardsolve(t(Rchol), rhs))

  Beta_list <- vector("list", S)
  idx <- 0

  for (s in seq_len(S)) {
    Beta_list[[s]] <- Beta_stack[(idx + 1):(idx + Kb), , drop = FALSE]
    idx <- idx + Kb
  }

  Beta_list
}



update_zeta <- function(
    Z_list,        # list of (n_s × G) matrices
    W_list,        # list of length-n_s vectors OR (n_s × G) matrices
    logLS_list,    # list of length-n_s vectors
    U_list,        # list of (n_s × Kb) biological spline matrices, no intercept
    Beta_list,     # list of (Kb × G) matrices
    LV_list,       # list of (n_s × Kt) matrices, already V_s * logLS_s
    Gamma_list,    # list of (Kt × G) matrices
    alpha,         # length-S vector
    delta          # S × G matrix
) {
  S <- length(Z_list)
  G <- ncol(Z_list[[1]])

  num <- numeric(G)
  den <- numeric(G)

  for (s in seq_len(S)) {
    Zs <- Z_list[[s]]
    Ws <- W_list[[s]]
    n_s <- nrow(Zs)

    R <- Zs -
      matrix(delta[s, ], nrow = n_s, ncol = G, byrow = TRUE) -
      alpha[s] * logLS_list[[s]] -
      U_list[[s]] %*% Beta_list[[s]] -
      LV_list[[s]] %*% Gamma_list[[s]]

    if (is.matrix(Ws)) {
      num <- num + colSums(Ws * R)
      den <- den + colSums(Ws)
    } else {
      num <- num + colSums(R * Ws)
      den <- den + sum(Ws)
    }
  }

  num / pmax(den, 1e-8)
}


update_rho <- function(
    Z_list,        # list of n_s × G working-response matrices
    W_list,        # list of length-n_s common weight vectors
    logLS_list,    # list of length-n_s log-library-size vectors
    U_list,        # list of n_s × Kb biological spline matrices
    Beta_list,     # list of Kb × G biological coefficient matrices
    LV_list,       # list of n_s × Kt technical design matrices, L_s V_s
    Gamma_list,    # list of Kt × G technical coefficient matrices
    alpha,         # length-S vector
    zeta,          # length-G vector
    A              # S × (S-1) contrast matrix
) {
  S <- length(Z_list)
  G <- ncol(Z_list[[1]])
  q <- ncol(A)  # q = S - 1

  # Build stacked residual matrix R and stacked contrast design Z_delta
  R_list <- vector("list", S)
  Zdelta_list <- vector("list", S)
  W_vec <- numeric(0)

  for (s in seq_len(S)) {

    Zs <- Z_list[[s]]      # n_s × G
    Ws <- W_list[[s]]      # length n_s

    if (is.matrix(Ws)) {
      stop("update_rho assumes common slice-specific weights.")
    }

    n_s <- nrow(Zs)
    a_s <- A[s, ]          # length q

    # R_s = Z_s - zeta - alpha_s logLS_s - U_s Beta_s - LV_s Gamma_s
    # R_s: n_s × G
    R_list[[s]] <- Zs -
      matrix(zeta, nrow = n_s, ncol = G, byrow = TRUE) -
      alpha[s] * logLS_list[[s]] -
      U_list[[s]] %*% Beta_list[[s]] -
      LV_list[[s]] %*% Gamma_list[[s]]

    # Z_delta_s = 1_{n_s} a_s^T
    # Z_delta_s: n_s × q
    Zdelta_list[[s]] <- matrix(1, nrow = n_s, ncol = 1) %*%
      matrix(a_s, nrow = 1)

    # W_vec stacks the diagonal weights
    # W_vec: length N after loop
    W_vec <- c(W_vec, Ws)
  }

  # Stack across slices
  # R:       N × G
  # Z_delta: N × q
  R <- do.call(rbind, R_list)
  Z_delta <- do.call(rbind, Zdelta_list)

  # Weighted normal-equation pieces
  # XtWX = Z_delta^T W Z_delta, q × q
  # XtWR = Z_delta^T W R,       q × G
  WZ_delta <- Z_delta * W_vec
  WR <- R * W_vec

  XtWX <- crossprod(Z_delta, WZ_delta)
  XtWR <- crossprod(Z_delta, WR)

  # rho: q × G
  solve(XtWX + diag(1e-8, q), XtWR)
}
