update_psi_edgeR <- function(Ys, Mus, max_cells = 1000) {
  # Ys  : G × n_cells
  # Mus : n_cells × G

  n_cells <- ncol(Ys)

  set.seed(1)
  if (n_cells > max_cells) {
    idx <- sample.int(n_cells, max_cells)  # numeric indices only
    Ys  <- Ys[, idx, drop = FALSE]
    Mus <- Mus[idx, , drop = FALSE]
    n_cells <- length(idx)
  }

  offs <- log(t(Mus))
  design <- matrix(1, nrow = n_cells, ncol = 1)

  stopifnot(
    nrow(design) == ncol(Ys),
    all(dim(offs) == dim(Ys))
  )

  disp <- edgeR::estimateDisp(
    y      = Ys,
    design = design,
    offset = offs,
    tagwise = TRUE,
    robust  = TRUE
  )$tagwise.dispersion

  disp[!is.finite(disp)] <- median(disp[is.finite(disp)])
  disp
}
