#' Sensitivity biplot
#'
#' @param mds_matrices The output from make_mds_matrices.
#' @param dist_fns The output from make_dist_fns.
#' @param k The number of embedding dimensions for multi-dimensional scaling. Defaults to 2.
#' @param samples Which of the points to compute sensitivities for.
#'
#' @return A list. The ith element of the list is a p x k matrix, and
#' the (jl)th element of which describes the sensitivity of variable j
#' on axis l around point i.
compute_sensitivity <- function(mds_matrices, dist_fns, k, samples) {
    dist_matrix = mds_matrices$delta^(.5)
    biplot_list = list()
    for(j in samples) {
        biplot_at_j = matrix(0, nrow = ncol(mds_matrices$Y), ncol = k)
        for(k_idx in 1:k) {
            ## TODO: check which number this should be
            lambda = mds_matrices$Lambda[k_idx]
            ## TODO: rename , check that we're doing the right thing
            A = apply(mds_matrices$Y, 1, function(y) {
                deriv = as.matrix(dist_fns$dist_deriv(y, mds_matrices$Y[j,]))
                return(deriv)
            })
            M = -.5 * lambda^(-1) * sweep(A, 2, mds_matrices$X[,k_idx] * dist_matrix[,j], FUN = '*')
            biplot_at_j[,k_idx] = rowSums(M)
        }
        biplot_list[[j]] = biplot_at_j
    }
    return(biplot_list)
}

#' Computes MDS representation and other associated values
#'
#' @param Y A samples x variables data matrix
#' @param dist_fn A function that computes the distances between the rows of Y.
#'
#' @return 
make_mds_matrices <- function(Y, dist_fn) {
    n = nrow(Y)
    dist_output = dist_fn(Y)
    ## as.matrix allows us to handle the output from the 'dist'
    ## function as well as matrix-valued outputs
    ## delta is the matrix that contains the squared distances
    delta = as.matrix(dist_output)^2
    A = -.5 * delta
    ## there is a faster way to do this
    P = diag(1, n) - n^(-1) * matrix(1, nrow = n, ncol = n)
    B = P %*% A %*% P
    Beig = eigen(B, symmetric = TRUE)
    Beig$vectors = Beig$vectors[,1:(n-1)]
    Beig$values = Beig$values[1:(n-1)]
    X = Beig$vectors %*% diag(sqrt(Beig$values))
    out = list()
    out$d2 = diag(B)
    out$B = B
    out$A = A
    out$delta = delta
    out$X = X
    out$Lambda = Beig$values
    out$Y = Y
    return(out)
}

#' Creates distance function and corresponding derivative function
#'
#' @param dist Either a string or a function.
#' @param dist_deriv Either NULL or a function.
#'
#' @return A list containing two functions, dist_fn and dist_deriv
make_dist_fns <- function(dist_fn, dist_deriv) {
    if(typeof(dist_fn) == "closure" & typeof(dist_deriv) == "closure") {
        return(list(dist_fn = dist_fn, dist_deriv = dist_deriv))
    }
    if(dist_fn == "euclidean") {
        dist_fn = function(x) dist(x, method = "euclidean")
        return(list(dist_fn = dist_fn, dist_deriv = euclidean_dist_deriv))

    } else if(dist_fn == "manhattan") {
        dist_fn = function(x) dist(x, method = "manhattan")
        return(list(dist_fn = dist_fn, dist_deriv = manhattan_dist_deriv))
    } else if(dist_fn == "maximum") {
        dist_fn = function(x) dist(x, method = "maximum")
        return(list(dist_fn = dist_fn, dist_deriv = maximum_dist_deriv))
    } else {
        stop("Unsupported distance")
    }
}


#' Partial gradient for Euclidean distance
#'
#' @param x
#' @param y
#'
#' @return If x and y each have length p, the function returns a
#' p-vector with jth element equal to \frac{\partial}{\partial y_{j}}
#' d(x,y)
euclidean_dist_deriv <- function(x, y) {
    if(sum((y - x)^2) == 0) {
        return(rep(0, length(y)))
    }
    return((y - x) * (sum((y - x)^2))^(-.5))
}

#' Partial gradient for manhattan distance
#'
#' @param x
#' @param y
#'
#' @return If x and y each have length p, the function returns a
#' p-vector with jth element equal to \frac{\partial}{\partial y_{j}}
#' d(x,y)
manhattan_dist_deriv <- function(x, y) {
    derivs = ifelse(y < x, -1, 1)
    derivs[x == y] = 0
    return(derivs)
}

#' Partial gradient for max distance
#'
#' @param x
#' @param y
#'
#' @return If x and y each have length p, the function returns a
#' p-vector with jth element equal to \frac{\partial}{\partial y_{j}}
#' d(x,y)
maximum_dist_deriv <- function(x, y) {
    max_abs = max(abs(y - x))
    deriv = ifelse(y < x, -1, 1)
    deriv[abs(y - x) != max_abs] = 0
    return(deriv)
}

minkowski_dist_deriv <- function(x, y, p) {
    
}

#' Gradient for weighted Unifrac
#'
#' If d_{x,T}(y) denotes the weighted Unifrac distance between x and y
#' with tree T, computes the gradient of d_{x,T} evaluated at y.
#' @param x A p-vector.
#' @param y A p-vector.
#' @param A A p x B matrix, with A_{ij} = 1 indicating that taxon i
#' descends from branch j.
#'
#' @return A p-vector giving the gradient.
#' @export
wuf_dist_deriv <- function(x, y, A) {
    sy = sum(y)
    sx = sum(x)
    py = y / sy
    px = x / sx
    by = py %*% A
    bx = px %*% A
    bp_deriv = as.vector((1 - by)) / sy
    signs = ifelse(by > bx, 1, -1)
    grad = A %*% (bp_deriv * signs)
    return(grad)
}

#' Weighted Unifrac distance
#'
#' @param X A n x p matrix, n samples, p taxa.
#' @param A A p x B marix, A_{ij} = 1 indicating that taxon i descends
#' from branch j.
#'
#' @return An n x n distance matrix, with element i,j giving the
#' weighted Unifrac distance between X[i,] and X[j,].
#' @export
wuf_dist <- function(X, A) {
    bp_matrix = X %*% A
    return(dist(bp_matrix, method = "manhattan"))
}

#' Plots a sensitivity biplot
#'
#' @param Y A matrix or data frame with the data the distances were computed from.
#' @param sensitivity_list A list giving the output from compute_sensitivity
#'
#' @return
#' @export
sensitivity_biplot <- function(Y, dist, dist_deriv = NULL, k = 2, plotting_axes = 1:2,
                               samples = 1:nrow(Y), only_df = FALSE, ...) {
    dist_fns = make_dist_fns(dist, dist_deriv)
    mds_matrices = make_mds_matrices(Y, dist_fns$dist_fn)
    sensitivity_list = compute_sensitivity(mds_matrices, dist_fns, k = k, samples = samples)
    biplot_df = make_biplot_data_frame(sensitivity_list, mds_matrices,
        axes = plotting_axes, samples = samples)
    if(only_df) {
        return(biplot_df)
    }
    ggplot(biplot_df) +
        geom_segment(aes(x = x, y = y, xend = xend, yend = yend, color = variable)) +
        xlab("MDS1") + ylab("MDS2")

}

make_biplot_data_frame <- function(sensitivity_list, mds_matrices, axes = 1:2, samples) {
    n = nrow(mds_matrices$Y)
    p = ncol(mds_matrices$Y)
    varnames = get_variable_names(mds_matrices$Y)
    biplot_list = lapply(1:n, function(i) {
        if(i %in% samples) {
            return(make_one_sample_biplot_df(sensitivity = sensitivity_list[[i]][,axes],
                                             center = mds_matrices$X[i,axes],
                                             sample = i, varnames = varnames))

        } else {
            point = mds_matrices$X[i,axes]
            return(data.frame(x = point[1], y = point[2], xend = NA, yend = NA,
                              variable = NA, sample = i))
        }
    })
    biplot_df = Reduce(rbind, biplot_list)
    return(biplot_df)
}


get_variable_names <- function(m) {
    if(is.null(colnames(m))) {
        return(paste("Var", 1:ncol(m), sep = ""))
    }
    return(colnames(m))
}

make_one_sample_biplot_df <- function(sensitivity, center, sample, varnames) {
    endpoints = sweep(sensitivity, MARGIN = 2, STATS = center, FUN = '+')
    centers = matrix(center, nrow = nrow(sensitivity), ncol = ncol(sensitivity), byrow = TRUE)
    one_sample_df = data.frame(centers, endpoints, variable = varnames, sample = sample)
    colnames(one_sample_df)[1:4] = c("x", "y", "xend", "yend")
    return(one_sample_df)
}
