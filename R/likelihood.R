# ---------------- NB log-likelihood  ----------------
nb_loglik <- function(Y, Mu, Psi) {
  eps <- 1e-12
  dnbinom(
    x    = as.vector(Y),
    mu   = as.vector(Mu),
    size = 1 / pmax(as.vector(Psi), eps),
    log  = TRUE
  )
}

# ---------------- gene-wise NB log-likelihood ----------------
ll_gene <- function(Y, Mu, Psi, n_cells) {
  ll <- nb_loglik(Y, Mu, Psi)
  colSums(matrix(ll, nrow = n_cells))   # length = G
}

# ---------------- gene-wise penalized log-likelihood ----------------
penalized_ll <- function(
    Y, Mu, Psi,
    Beta, Gamma,
    Cmat, Pbeta,
    lambdaC, lambda1, lambda2,
    n_cells
) {
  # ---------- NB likelihood per gene ----------
  llg <- ll_gene(Y, Mu, Psi, n_cells)   # length G

  # ---------- niche-matching penalty (per gene) ----------
  CB <- as.matrix(Cmat %*% Beta)                   # (constraints × G)
  penC_g <- 0.5 * lambdaC * colSums(CB^2)

  # ---------- biological spline penalty (per gene) ----------
  PB <- Pbeta %*% Beta                  # (basis × G)
  penB_g <- 0.5 * lambda1 * colSums(Beta * PB)

  # ---------- technical spline ridge (per gene) ----------
  penG_g <- 0.5 * lambda2 * colSums(Gamma^2)

  # ---------- gene-wise penalized log-likelihood ----------
  llg - penC_g - penB_g - penG_g

}





