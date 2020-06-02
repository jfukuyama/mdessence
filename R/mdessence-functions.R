#' Sensitivity biplot at input data points
#'
#' @param mds_matrices The output from make_mds_matrices.
#' @param dist_fns The output from make_dist_fns.
#' @param k The number of embedding dimensions for multi-dimensional scaling. Defaults to 2.
#' @param samples Which of the points to compute sensitivities for.
#'
#' @return A list. The ith element of the list is a p x k matrix, and
#' the (jl)th element of which describes the sensitivity of variable j
#' on axis l around point i.
compute_sensitivity_samples <- function(mds_matrices, dist_fns, k, samples, ridge = 0) {
    dist_matrix = mds_matrices$delta^(.5)
    biplot_list = list()
    Ylambdainv = sweep(mds_matrices$Y[,1:k], MARGIN = 2,
                       STATS = mds_matrices$Lambda[1:k] + ridge, FUN = "/")
    for(j in samples) {
        dist_to_j = dist_matrix[,j]
        dist_jacobian = apply(mds_matrices$X, 1, function(x) {
            -dist_fns$dist_deriv(x, mds_matrices$X[j,])
        })
        Jd = sweep(dist_jacobian, MARGIN = 2, STATS = dist_to_j, FUN = "*")
        biplot_list[[j]] = Jd %*% Ylambdainv
    }
    return(biplot_list)
}

#' Sensitivity biplot ot new points
#'
#' @param mds_matrices The output from make_mds_matrices.
#' @param dist_fns The output from make_dist_fns.
#' @param k The number of embedding dimensions for multi-dimensional scaling. Defaults to 2.
#' @param new_points A list with new points to compute biplot axes for.
#' 
#' @return A list. The ith element of the list is a p x k matrix, and
#' the (jl)th element of which describes the sensitivity of variable j
#' on axis l around point i.
compute_sensitivity_new_points <- function(mds_matrices, dist_fns, k, new_points, n_random_points = NA, ridge = 0) {
    if(!is.na(n_random_points)) {
        new_points = replicate(n = n_random_points, {
            convex_comb = DirichletReg::rdirichlet(n = 1, alpha = rep(1, nrow(mds_matrices$X)))
            return(convex_comb %*% mds_matrices$X)
        },
        simplify = FALSE)
    }
    biplot_list = list()
    Ylambdainv = sweep(mds_matrices$Y[,1:k], MARGIN = 2,
                       STATS = mds_matrices$Lambda[1:k] + ridge, FUN = "/")

    for(i in 1:length(new_points)) {
        new_point = new_points[[i]]
        ## mds_matrices$X is the original data
        dist_to_new_point = apply(mds_matrices$X, 1, function(x) {
            as.matrix(dist_fns$dist_fn(rbind(x, new_point)))[1,2]
        })
        dist_jacobian = apply(mds_matrices$X, 1, function(x) {
            -dist_fns$dist_deriv(x, new_point)
        })
        embedding = .5 * dist_to_new_point %*% Ylambdainv
        axis_center = matrix(embedding, nrow = ncol(mds_matrices$X), ncol = k, byrow = TRUE)
        Jd = sweep(dist_jacobian, MARGIN = 2, STATS = dist_to_new_point, FUN = "*")
        biplot_axes = Jd %*% Ylambdainv
        biplot_df = data.frame(axis_center, axis_center + biplot_axes, biplot_axes)
        names(biplot_df) = c(paste0("Embedding", 1:k), paste0("AxisEnd", 1:k), paste0("Axis", 1:k))
        biplot_df$variable = colnames(mds_matrices$X)
        biplot_df$sample = i
        biplot_list[[i]] = biplot_df
    }
    return(Reduce(rbind, biplot_list))
}

#' Computes MDS representation and other associated values
#'
#' @param Y A samples x variables data matrix
#' @param dist_fn A function that computes the distances between the rows of Y.
#'
#' @return A list, containing
#' - delta: Matrix of squared distances.
#' - jdj: Row- and column-centered -.5 * delta
#' - d2: The diagonal elements of B
#' - Y: The embeddings of the samples in the MDS space.
#' - Lambda: The eigenvalues of B.
#' - X: The original data.
make_mds_matrices <- function(X, dist_fn) {
    n = nrow(X)
    dist_output = dist_fn(X)
    ## as.matrix allows us to handle the output from the 'dist'
    ## function as well as matrix-valued outputs
    ## delta is the matrix that contains the squared distances
    delta = as.matrix(dist_output)^2
    A = -.5 * delta
    ## there is a faster way to do this
    J = diag(1, n) - n^(-1) * matrix(1, nrow = n, ncol = n)
    jdj = J %*% A %*% J
    Beig = eigen(jdj, symmetric = TRUE)
    Beig$vectors = Beig$vectors[,1:(n-1)]
    Beig$values = Beig$values[1:(n-1)]
    Y = Beig$vectors %*% diag(sqrt(Beig$values))
    colnames(Y) = paste("Axis", 1:ncol(Y), sep = "")
    out = list()
    out$d2 = diag(jdj)
    out$jdj = jdj
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
    } else if(dist_fn == "jaccard") {
        dist_fn = function(x) {
            d = matrix(0, nrow = nrow(x), ncol = nrow(x))
            for(i in 1:(nrow(x) - 1)) {
                for(j in (i+1):nrow(x)) {
                    d[i,j] = jaccard_dis(x[i,], x[j,])
                }
            }
            d = d + t(d)
            return(d)
        }
        return(list(dist_fn = dist_fn, dist_deriv = jaccard_deriv))
    }  else {
        stop("Unsupported distance")
    }
}


#' Partial gradient for Euclidean distance
#'
#' @param x A p-vector.
#' @param y A p-vector.
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
#' @param x A p-vector.
#' @param y A p-vector.
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
#' @param x A p-vector.
#' @param y A p-vector.
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
wuf_dist_deriv <- function(x, y, A, L) {
    sy = sum(y)
    sx = sum(x)
    py = y / sy
    px = x / sx
    by = py %*% A
    bx = px %*% A
    bp_deriv = as.vector((1 - by)) / sy
    signs = ifelse(by > bx, 1, -1)
    grad = A %*% as.vector((L * (bp_deriv * signs)))
    return(as.vector(grad))
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
wuf_dist <- function(X, A, L) {
    bp_matrix = X %*% sweep(A, MARGIN = 2, STATS = L, FUN = "*")
    return(dist(bp_matrix, method = "manhattan"))
}

#' Unweighted Unifrac distance
#'
#' @param X A n x p matrix, n samples, p taxa.
#' @param A A p x B marix, A_{ij} = 1 indicating that taxon i descends
#' from branch j.
#' @param L A vector giving branch lengths.
#'
#' @return An n x n distance matrix, with element i,j giving the
#' weighted Unifrac distance between X[i,] and X[j,].
#' @export
uf_dist <- function(X, A, L) {
    bp_matrix = X %*% A
    br_ind_matrix = bp_matrix > 0
    br_ind_matrix_scaled = sweep(br_ind_matrix, MARGIN = 2, STATS = L, FUN = "*")
    return(dist(br_ind_matrix_scaled, method = "manhattan") / sum(L))
}

#' Generalized gradient for unweighted Unifrac
#'
#' If d_{x,T}(y) denotes the unweighted Unifrac distance between x and y
#' with tree T, computes the generalized gradient of d_{x,T} evaluated at y.
#' @param x A p-vector.
#' @param y A p-vector.
#' @param A A p x B matrix, with A_{ij} = 1 indicating that taxon i
#' descends from branch j.
#' @param L A vector of branch lengths.
#' @param delta_min Smallest perturbation.
#' @param positive Logical, positive perturbation?
#'
#' @return A p-vector giving the gradient.
#' @export
uf_dist_deriv <- function(x, y, A, L, delta_min = 1, positive = TRUE,
                          return_dist_difference = FALSE) {
    y = as.vector(y)
    x = as.vector(x)
    if(positive) {
        rho = ifelse(y == 0, 1, Inf)
    } else {
        rho = ifelse(y == 0, Inf, y)
    }
    dxy = as.numeric(uf_dist(rbind(x, y), A, L))
    deriv = sapply(seq_along(y), function(i) {
        if(rho[i] == Inf) return(0)
        if(positive) {
            delta = uf_delta(x, y, A, L, i, rho = rho[i])
        } else {
            delta = -uf_delta(x, y, A, L, i, rho = -rho[i])
        }
        if(return_dist_difference) {
            return((-1) * delta * (delta + 2 * dxy))
        } else {
            return(delta_min * rho[i]^(-1) * (-1) * delta * (delta + 2 * dxy))
        }
    })
    return(deriv)
}


#' Pre-computations for unifrac biplots
#'
#' Saves matrices that we need to access many times when computing the
#' unifrac biplot axes.
#'
#' @param X A samples by species matrix
#' @param A The branches by species ancestry matrix
#' @param L Branch-length vector
#'
#' @return A list, containing 
#' 
#' @export
uf_sensitivity_computations <- function(X, A, L) {
    A_lgc = A != 0
    ##X_lgc = as(X != 0, "lgCMatrix")
    X_lgc = X != 0
    sample_indicator_list = list()
    for(i in 1:nrow(X)) {
        sample_indicator_list[[i]] = Matrix::t(X_lgc[i,,drop=FALSE])
    }
    sample_species_combinations = expand.grid(1:nrow(X), 1:ncol(X))
    ## this is a list where sample_sensitivity_indicators[[j]][[k]]
    ## has A^{(jk)}, the matrix indicating whether the contribution of
    ## the bth branch to the unifrac distance between sample j and any
    ## other sample will change if taxon k is increased
    sample_sensitivity_indicators = list()
    for(i in 1:nrow(X)) {
        sample_sensitivity_indicators[[i]] = list()
    }
    for(i in 1:nrow(sample_species_combinations)) {
        sample = as.numeric(sample_species_combinations[i,1])
        species = as.numeric(sample_species_combinations[i,2])
        d = as.logical(A_lgc[,species]) & !as.logical(A_lgc %*% X_lgc[sample,])
        diagonal = Diagonal(n = nrow(A), x = d)
        sample_sensitivity_indicators[[sample]][[species]] = drop0(diagonal %*% A != 0)
    }
    uf_distances = dist(sweep((X_lgc %*% Matrix::t(A_lgc) > 0), MARGIN = 2, STATS = L, FUN = "*"), "manhattan") / sum(L)
    return(list(
        sample_indicator_list = sample_indicator_list,
        sample_sensitivity_indicators = sample_sensitivity_indicators,
        uf_distances = uf_distances,
        L = L,
        T = sum(L)
    ))

}

#' Computes biplot axes around a point
#'
#' @param X A samples by species data matrix.
#' @param sample_idx The point around which to compute biplot axes.
#' @param species_indices The species for which to compute biplot axes.
#' @param delta_min Minimum value for perturbation.
#' @param positive Logical, perturbation in the positive direction?
#' @param uf_cached The output from uf_sensitivity_computations
#' @param mds_matrices The output from make_mds_matrices
#'
#' @return A species by n_axes matrix with biplot axes.
#' @export
uf_sensitivity <- function(X, sample_idx, species_indices = 1:ncol(X),
                           delta_min = 1, positive = TRUE, uf_cached,
                           mds_matrices, n_axes = 2) {
    if(class(X) != "Matrix") {
        X = as(X, "Matrix")
    }
    generalized_gradients =
        make_uf_generalized_gradients(X, sample_idx, species_indices,
                                      delta_min, positive = positive,
                                      uf_cached = uf_cached)
    biplot_axes =
        t(.5 * diag(mds_matrices$Lambda[1:n_axes]^(-1)) %*%
            t(mds_matrices$X[,1:n_axes,drop=FALSE]) %*% generalized_gradients)
    return(biplot_axes)
}


explicit_generalized_gradients <- function(X, sample_idx, delta_min, positive, A, L) {
    if(positive) {
        rho = ifelse(X[sample_idx,] == 0, delta_min, Inf)
    } else {
        rho = -ifelse(X[sample_idx,] == 0, Inf, X[sample_idx,])
    }
    generalized_gradients = matrix(0, nrow(X), ncol(X))
    for(i in 1:nrow(X)) {
        for(k in 1:ncol(X)) {
            x = X[i,]
            y = X[sample_idx,]
            y_pert = y
            y_pert[k] = y_pert[k] + rho[k]
            dxy2 = uf_dist(rbind(x, y), A, L)^2
            dxypert2 = uf_dist(rbind(x, y_pert), A, L)^2
            generalized_gradients[i,k] = (delta_min * rho[k]^(-1) * (dxy2 - dxypert2))
        }
    }
    return(generalized_gradients)
}

#' Computes change in uf distance after perturbation
#'
#' Computes d(xi,xj+rho e_k) - d(xi,xj)
#'
#' @param uf_cache The output from uf_sensitivity_computations.
delta_from_A_sens <- function(uf_cache, i, j, k) {
    x_j = x_j_pert = uf_cache$sample_indicator_list[[j]]
    ## next five lines should be equivalent to x_j_pert[k,1] = TRUE
    if(all(x_j_pert@i != (k-1))) {
        x_j_pert@x = rep(TRUE, length(x_j_pert@x) + 1)
        x_j_pert@p[2] = x_j_pert@p[2] + 1L
        x_j_pert@i = sort(c(x_j_pert@i, as.integer(k-1)))
    }
    A_sens = uf_cache$sample_sensitivity_indicators[[j]][[k]]
    A_sens_x_i = A_sens %*% uf_cache$sample_indicator_list[[i]]
    A_sens_x_j = A_sens %*% x_j
    A_sens_x_j_pert = A_sens %*% x_j_pert
    #if(any(A_sens_x_i@x == 0))
    #    browser()
    #if(any(A_sens_x_j@x == 0))
    #    browser()
    #if(any(A_sens_x_j_pert@x == 0))
    #       browser
    ##d1 = sum(abs((A_sens_x_i > 0) - (A_sens_x_j > 0)) * uf_cache$L)
    ##d2 = sum(abs((A_sens_x_i> 0) - (A_sens_x_j_pert > 0)) * uf_cache$L)
    d1 = mdessence:::uf_num_from_indicators(A_sens_x_i, A_sens_x_j, uf_cache$L)
    d2 = mdessence:::uf_num_from_indicators(A_sens_x_i, A_sens_x_j_pert, uf_cache$L)
    (d2 - d1) / sum(uf_cache$T)
}

sym_diff <- function(a,b) unique(c(setdiff(a,b), setdiff(b,a)))

#' Unifrac numerator
#'
#' More efficient, less transparent version of numerator computation.
#'
#' @param sm1 A p x 1 sparse matrix of class "dgCMatrix". Assumed to
#' have zeros dropped.
#' @param sm2 A p x 1 sparse matrix of class "dgCMatrix" Assumed to
#' have zeros dropped.
#' @param L A numeric p-vector.
#' @return Sum of unique branch lengths.
uf_num_from_indicators <- function(sm1, sm2, L) {
    #sm1 = as(sm1, "dgTMatrix")
    #sm2 = as(sm2, "dgTMatrix")
    #sm1 = drop0(sm1)
    #sm2 = drop0(sm2)
    sum_indices = sym_diff(sm1@i, sm2@i)
    ## +1 because Matrix uses zero indexing
    sum(L[sum_indices + 1])
}


make_uf_generalized_gradients <- function(X, sample_idx, species_indices,
                                          delta_min, positive, uf_cached) {
    X_lgc = as(X != 0, "lgTMatrix")
    if(positive) {
        rho = ifelse(!X_lgc[sample_idx,], delta_min, Inf)
    } else {
        rho = -ifelse(!X_lgc[sample_idx,], Inf, X[sample_idx,])
    }
    sample_sensitivity_matrices = uf_cached$sample_sensitivity_matrices
    uf_distances = as.matrix(uf_cached$uf_distances)
    Xt_lgc = Matrix::t(X_lgc)
    generalized_gradients = matrix(0, nrow(X), length(species_indices))
    for(i in 1:nrow(X)) {
        for(k in 1:length(species_indices)) {
            species_idx = species_indices[k]
            if(rho[species_idx] == Inf) {
                generalized_gradients[i, k] = 0
            } else {
                dxy = uf_distances[i, sample_idx]
                delta = delta_from_A_sens(uf_cache, i, sample_idx, species_idx)
                generalized_gradients[i, species_idx] =
                    delta_min / rho[species_idx] * (-1) * delta * (delta + 2 * dxy)
            }
        }
    }
    return(generalized_gradients)
}
#' Computes delta for unifrac
#'
#' Computes d(x,y) - d(x,y + rho e_i) for d the unweighted Unifrac
#' distance.
#'
#' @param A A species by branches ancestry matrix.
#' @param L A branch length vector.
uf_delta <- function(x, y, A, L, i, rho) {
    i_s_ancestors = which(A[i,] > 0)
    A_sub = A[,i_s_ancestors,drop=FALSE]
    L_sub = L[i_s_ancestors]
    y_pert = y
    y_pert[i] = y_pert[i] + rho
    original = sum(abs(((t(x) %*% A_sub) > 0) - ((t(y) %*% A_sub) > 0)) * L_sub)
    ##new = sum(abs(((t(x) %*% A_sub) > 0) - (t(y_pert) %*% A_sub) > 0) * L_sub)
    new = sum(abs(((t(x) %*% A_sub) > 0) - ((t(y_pert) %*% A_sub) > 0)) * L_sub)
    delta = (original - new) / sum(L)
}

#' Derivative for Jaccard distance
#'
#' @param x A p-vector.
#' @param y A p-vector
#' @param delta_min Minimum value for the perturbation.
#' @param positive Positive perturbation?
#'
#' @return If x and y each have length p, the function returns a
#' p-vector with jth element equal to
#' \frac{delta_min}{rho} (d^2(x, y)- d^2(x, y + rho e_j))
jaccard_deriv <- function(x, y, delta_min = 1, positive = TRUE) {
    if(positive) {
        rho = ifelse(y == 0, 1, Inf)
    } else {
        rho = ifelse(y == 0, Inf, y)
    }
    dxy = jaccard_dis(x, y)^2
    deriv = sapply(seq_along(y), function(i) {
        if(rho[i] == Inf) return(0)
        y_pert = y
        if(positive) {
            y_pert[i] = y_pert[i] + rho[i]
        } else {
            y_pert[i] = y_pert[i] - rho[i]
        }
        return(delta_min * rho[i]^(-1) * (dxy - jaccard_dis(x,y_pert)^2))
    })
    return(deriv)
}

jaccard_dis <- function(x, y) {
    a = sum(x > 0 & y > 0)
    b = sum(x > 0 & y == 0)
    c = sum(x == 0 & y > 0)
    if(a + b + c == 0) {
        return(0)
    }
    return(1 - a / (a + b + c))
}

#' Plots a sensitivity biplot
#'
#' @param Y A matrix or data frame with the data the distances were computed from.
#' @param sensitivity_list A list giving the output from compute_sensitivity
#'
#' @return A plot.
#' @export
sensitivity_biplot <- function(X, dist, dist_deriv = NULL, k = 2, plotting_axes = 1:2,
                               n_random_points = NA, new_points = NA,
                               only_df = FALSE, scale = 1,
                               ridge = 0, ...) {
    dist_fns = make_dist_fns(dist, dist_deriv)
    mds_matrices = make_mds_matrices(X, dist_fns$dist_fn)
    sensitivity_list = compute_sensitivity_new_points(mds_matrices, dist_fns, k = k,
                                                      n_random_points = n_random_points,
                                                      new_points = new_points,
                                                      ridge = ridge)
    if(only_df) {
        return(sensitivity_list)
    }
    ggplot(biplot_df) +
        geom_segment(aes(x = x, y = y, xend = xend, yend = yend, color = variable)) +
        xlab("MDS1") + ylab("MDS2")

}

make_biplot_data_frame <- function(sensitivity_list, mds_matrices, axes = 1:2, samples, scale) {
    n = nrow(mds_matrices$X)
    p = ncol(mds_matrices$X)
    varnames = get_variable_names(mds_matrices$X)
    biplot_list = lapply(1:n, function(i) {
        if(i %in% samples) {
            return(make_one_sample_biplot_df(sensitivity = sensitivity_list[[i]][,axes],
                                             center = mds_matrices$Y[i,axes],
                                             sample = i, varnames = varnames,
                                             scale = scale))

        } else {
            point = mds_matrices$Y[i,axes]
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

make_one_sample_biplot_df <- function(sensitivity, center, sample, varnames, scale) {
    endpoints = sweep(sensitivity * scale, MARGIN = 2, STATS = center, FUN = '+')
    centers = matrix(center, nrow = nrow(sensitivity), ncol = ncol(sensitivity), byrow = TRUE)
    one_sample_df = data.frame(centers, endpoints, variable = varnames, sample = sample)
    colnames(one_sample_df)[1:4] = c("x", "y", "xend", "yend")
    return(one_sample_df)
}


correlation_biplot <- function(X, dist, k = 2, plotting_axes = 1:2) {
    mds_matrices = make_mds_matrices(X, dist)
    biplot_axes = cor(X, mds_matrices$Y)[,plotting_axes]
    colnames(biplot_axes) = sapply(plotting_axes, function(i) sprintf("Axis%i", i))
    return(biplot_axes)
}
