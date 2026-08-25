#blocks
build_blocks <- function(LS, U, V) {
  if (any(LS == 0)) {
    stop("LS contains zero values, cannot take log.")
  }

  logLS <- log(LS)

  list(
    # Global slice-specific library size effect
    X_alpha = matrix(logLS, ncol = 1),

    # Biological spline design matrix only; no intercept here
    X_beta  = U,

    # Technical spline design matrix scaled by log library size
    X_gamma = scale_rows(V, logLS)
  )
}


# Compute average biological spline basis within each matched niche
# Returns a K_b × M matrix (columns = niches)
niche_means <- function(U, niche_vec, niches) {

  K <- ncol(U)
  M <- length(niches)

  out <- matrix(0, K, M)
  colnames(out) <- niches

  for (m in seq_len(M)) {

    idx <- which(niche_vec == niches[m])

    if (length(idx) > 0) {
      out[, m] <- colMeans(U[idx, , drop = FALSE])
    }
  }

  out
}


#basis construction
sp2d <- function(x, y, dfx = 6, dfy = 6) {
  # Construct 1D natural spline bases for x and y
  bs.x <- splines::ns(x, df = dfx)
  bs.y <- splines::ns(y, df = dfy)

  # Initialize tensor-product spline design matrix
  out <- matrix(0, nrow = length(x), ncol = dfx * dfy)

  # Form tensor products B_i(x) * B_j(y)
  for (i in seq_len(dfx)) {
    for (j in seq_len(dfy)) {
      out[, (i - 1) * dfy + j] <- bs.x[, i] * bs.y[, j]
    }
  }

  # Column-center the spline basis (no scaling)
  scale(out, scale = FALSE)
}

scale_rows <- function(A, v) {
  # Multiply each row of A by the corresponding element of v
  # Assumes length(v) == nrow(A)
  A * as.numeric(v)
}


#get slice data
get_slice_data <- function(spe) {

  coords <- SpatialExperiment::spatialCoords(spe)
  cts    <- SummarizedExperiment::assay(spe, "counts")
  cd     <- SummarizedExperiment::colData(spe)

  if ("LS" %in% colnames(cd)) {

    LS <- as.numeric(cd$LS)
    message("LS extracted from colData")

    if (length(LS) != ncol(cts)) {
      stop("colData(spe)$LS length does not match number of columns.")
    }

  } else {
    LS <- as.numeric(MatrixGenerics::colSums(cts))
    LS <- LS/mean(LS)
  }

  list(
    x     = coords[, 1],
    y     = coords[, 2],
    cts   = cts,
    LS    = LS,
    niche = cd$niche
  )
}



#winsorize coefficients to prevent numerical explosions
cap_coeff <- function(M, k = 20) {

  med <- matrixStats::rowMedians(M)
  mad <- matrixStats::rowMads(M)

  hi  <- med + k * mad
  lo  <- med - k * mad

  hi_mat <- matrix(hi, nrow = nrow(M), ncol = ncol(M))
  lo_mat <- matrix(lo, nrow = nrow(M), ncol = ncol(M))

  pmin(pmax(M, lo_mat), hi_mat)
}

#
# cap_coeff1 <- function(M, k = 20, name = "coeff") {
#
#   med <- matrixStats::rowMedians(M, na.rm = TRUE)
#   mad <- matrixStats::rowMads(M, na.rm = TRUE)
#
#   hi <- med + k * mad
#   lo <- med - k * mad
#
#   hi_mat <- matrix(hi, nrow = nrow(M), ncol = ncol(M))
#   lo_mat <- matrix(lo, nrow = nrow(M), ncol = ncol(M))
#
#   M_cap <- pmin(pmax(M, lo_mat), hi_mat)
#
#   n_cap <- sum(M_cap != M, na.rm = TRUE)
#   prop_cap <- n_cap / sum(is.finite(M))
#
#   message(sprintf(
#     "%s capped: %d / %d = %.4f%%",
#     name, n_cap, sum(is.finite(M)), 100 * prop_cap
#   ))
#
#   M_cap
# }
#
#
# cap_coeff_no_intercept <- function(M, k = 20) {
#
#   M2 <- M
#   rows <- 2:nrow(M)
#
#   med <- matrixStats::rowMedians(M2[rows, , drop = FALSE])
#   mad <- matrixStats::rowMads(M2[rows, , drop = FALSE])
#
#   hi <- med + k * mad
#   lo <- med - k * mad
#
#   hi_mat <- matrix(hi, nrow = length(rows), ncol = ncol(M2))
#   lo_mat <- matrix(lo, nrow = length(rows), ncol = ncol(M2))
#
#   M2[rows, ] <- pmin(pmax(M2[rows, , drop = FALSE], lo_mat), hi_mat)
#
#   M2
# }





#assign normalized data back to spe object
norm_spe <- function(spe_list, norm_data) {

  for (nm in names(spe_list)) {
    spe_obj <- spe_list[[nm]]

    # Get counts dimensions once
    cts <- SummarizedExperiment::assay(spe_obj, 'counts')

    # Transpose and assign dimnames in one step
    data <- t(norm_data[[nm]])
    dimnames(data) <- dimnames(cts)

    # Assign and update list
    SummarizedExperiment::assay(spe_obj, "logcounts") <- data
    spe_list[[nm]] <- spe_obj
  }

  return(spe_list)
}

#assign back original coords
restore_coords <- function(spe) {
  if (!is.null(metadata(spe)$orig_coords)) {
    spatialCoords(spe) <- metadata(spe)$orig_coords
    metadata(spe)$orig_coords <- NULL  # optional cleanup
  }
  spe
}





#sampling
sample_spe_proportional <- function(spe,
                                    niche_col  = "common_niche",
                                    sample.p   = 0.5,
                                    seed       = NULL,
                                    verbose    = TRUE,
                                    return_all = FALSE) {


  cd <- SummarizedExperiment::colData(spe)

  if (!niche_col %in% colnames(cd)) {
    stop("'", niche_col, "' not found in colData(spe).")
  }

  if (!is.numeric(sample.p) || length(sample.p) != 1 || sample.p <= 0 || sample.p > 1) {
    stop("'sample.p' must be a single numeric value in (0, 1].")
  }

  if (!is.null(seed)) {
    set.seed(seed)
  }

  # ---- extract ----
  niche_vec <- cd[[niche_col]]
  niche_vec <- factor(niche_vec)
  niche_vec <- droplevels(niche_vec) #drop unused cluster assignments
  cells     <- colnames(spe)

  # ---- split indices ----

  idx_valid <- which(!is.na(niche_vec))
  idx_na    <- which(is.na(niche_vec))

  # ---- stratified sampling (valid niches) ----
  sampled_idx_valid <- unlist(lapply(split(idx_valid, niche_vec[idx_valid]), function(idx) {
    n <- length(idx)
    k <- max(1, floor(sample.p * n))
    sample(idx, size = k)
  }), use.names = FALSE)

  # ---- sample NAs if present ----
  if (length(idx_na) > 0) {
    k_na <- max(1, floor(sample.p * length(idx_na)))
    sampled_idx_na <- sample(idx_na, size = k_na)
    sampled_idx <- c(sampled_idx_valid, sampled_idx_na)
  } else {
    sampled_idx <- sampled_idx_valid
  }

  # ---- subset SPE (correct way) ----
  spe_sampled <- spe[, sampled_idx, drop = FALSE]
  S4Vectors::metadata(spe_sampled)$sampled_idx <- sampled_idx
  # ---- construct tracking df (optional output) ----
  cell_niche <- data.frame(
    cell  = cells,
    niche = niche_vec,
    stringsAsFactors = FALSE
  )

  sampled_cells <- data.frame(
    cell  = cells[sampled_idx],
    niche = niche_vec[sampled_idx],
    stringsAsFactors = FALSE
  )

  # ---- verbose diagnostics ----
  if (verbose) {
    orig_tab <- table(niche_vec, useNA = "always")
    new_tab  <- table(sampled_cells$niche, useNA = "always")

    orig_prop <- orig_tab / sum(orig_tab)
    new_prop  <- new_tab / sum(new_tab)

    # align names safely
    all_names <- union(names(orig_prop), names(new_prop))
    orig_prop <- orig_prop[match(all_names, names(orig_prop))]
    new_prop  <- new_prop[match(all_names, names(new_prop))]

    orig_prop[is.na(orig_prop)] <- 0
    new_prop[is.na(new_prop)]   <- 0

    diff <- abs(orig_prop - new_prop)

    result <- data.frame(
      niche        = all_names,
      prop         = orig_prop,
      sampled_prop = new_prop,
      diff         = diff,
      row.names    = NULL
    )

    cat("=== Sample SPE Proportional ===\n")
    cat("NAs detected:                ", length(idx_na), "\n")
    cat("Original cells:              ", length(cells), "\n")
    cat("Sampled cells:               ", length(sampled_idx), "\n")
    #cat("Max proportion diff:         ", round(max(diff), 5), "\n")
    #cat("Proportions preserved (<1%): ", all(diff < 0.01), "\n\n")
    #print(result)
  }

  # ---- return ----
  if (return_all) {
    return(list(
      spe_sampled   = spe_sampled,
      spe_original  = spe,
      sampled_cells = sampled_cells,
      cell_niche    = cell_niche
    ))
  } else {
    return(spe_sampled)
  }
}




make_sumzero_contrast <- function(S) {
  if (S < 2) {
    stop("S must be at least 2.")
  }

  rbind(
    diag(S - 1),
    -rep(1, S - 1)
  )
}

build_Z_delta <- function(A, n_list) {
  S <- length(n_list)

  if (nrow(A) != S) {
    stop("nrow(A) must equal number of slices.")
  }

  do.call(rbind, lapply(seq_len(S), function(s) {
    matrix(1, nrow = n_list[s], ncol = 1) %*%
      matrix(A[s, ], nrow = 1)
  }))
}
