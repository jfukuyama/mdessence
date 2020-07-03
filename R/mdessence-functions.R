#' Local biplot at input data points
#'
#' @param mds_matrices The output from make_mds_matrices.
#' @param dist_fns The output from make_dist_fns.
#' @param k The number of embedding dimensions for multi-dimensional scaling. Defaults to 2.
#' @param samples Which of the points to compute sensitivities for.
#'
#' @return A list. The ith element of the list is a p x k matrix, and
#' the (jl)th element of which describes the LB axis of variable j
#' on axis l around point i.
compute_lb_samples <- function(mds_matrices, dist_fns, k, samples) {
    dist_matrix = mds_matrices$delta^(.5)
    biplot_list = list()
    Ylambdainv = sweep(mds_matrices$Y[,1:k], MARGIN = 2,
                       STATS = mds_matrices$Lambda[1:k], FUN = "/")
    for(j in samples) {
        dist_to_j = dist_matrix[,j]
        dist_jacobian = apply(mds_matrices$X, 1, function(x) {
            -dist_fns$dist_deriv(x, mds_matrices$X[j,])
        })
        Jd = sweep(dist_jacobian, MARGIN = 2, STATS = dist_to_j, FUN = "*")
        biplot_axes = Jd %*% Ylambdainv
        embedding = mds_matrices$Y[j,1:k]
        axis_center = matrix(embedding, nrow = ncol(mds_matrices$X), ncol = k, byrow = TRUE)
        biplot_df = data.frame(axis_center,  biplot_axes)
        names(biplot_df) = c(paste0("Embedding", 1:k), paste0("Axis", 1:k))
        biplot_df$variable = colnames(mds_matrices$X)
        biplot_df$sample = paste0("Original", j)
        biplot_list[[j]] = biplot_df
    }
    return(Reduce(rbind, biplot_list))
}

#' Local biplot at new points
#'
#' @param mds_matrices The output from make_mds_matrices.
#' @param dist_fns The output from make_dist_fns.
#' @param k The number of embedding dimensions for multi-dimensional scaling. Defaults to 2.
#' @param new_points A list with new points to compute biplot axes for.
#' 
#' @return A list. The ith element of the list is a p x k matrix, and
#' the (jl)th element of which describes the LB axis of variable j
#' on axis l around point i.
compute_lb_new_points <- function(mds_matrices, dist_fns, k, new_points, n_random_points, alpha = 1) {
    if(n_random_points > 0) {
        new_points = replicate(n = n_random_points, {
            convex_comb = DirichletReg::rdirichlet(n = 1, alpha = rep(alpha, nrow(mds_matrices$Y)))
            return(convex_comb %*% mds_matrices$X)
        },
        simplify = FALSE)
    }
    biplot_list = list()
    Ylambdainv = sweep(mds_matrices$Y[,1:k], MARGIN = 2,
                       STATS = mds_matrices$Lambda[1:k], FUN = "/")

    for(i in 1:length(new_points)) {
        new_point = new_points[[i]]
        dist_to_new_point = apply(mds_matrices$X, 1, function(x) {
            as.matrix(dist_fns$dist_fn(rbind(x, new_point)))[1,2]
        })
        dist_jacobian = apply(mds_matrices$X, 1, function(x) {
            -dist_fns$dist_deriv(x, new_point)
        })
        embedding = .5 * (diag(mds_matrices$jdj) - dist_to_new_point^2) %*% Ylambdainv
        axis_center = matrix(embedding, nrow = ncol(mds_matrices$X), ncol = k, byrow = TRUE)
        Jd = sweep(dist_jacobian, MARGIN = 2, STATS = dist_to_new_point, FUN = "*")
        biplot_axes = Jd %*% Ylambdainv
        biplot_df = data.frame(axis_center, biplot_axes)
        names(biplot_df) = c(paste0("Embedding", 1:k), paste0("Axis", 1:k))
        biplot_df$variable = colnames(mds_matrices$X)
        biplot_df$sample = paste0("New", i)
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

    } else if(dist_fn == "manhattan-pos") {
        dist_fn = function(x) dist(x, method = "manhattan")
        return(list(dist_fn = dist_fn, dist_deriv = manhattan_dist_deriv_pos))
    } else if(dist_fn == "manhattan-neg") {
        dist_fn = function(x) dist(x, method = "manhattan")
        return(list(dist_fn = dist_fn, dist_deriv = manhattan_dist_deriv_neg))
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
manhattan_dist_deriv_pos <- function(x, y) {
    derivs = ifelse(y < x, -1, 1)
    return(derivs)
}


#' Partial gradient for manhattan distance
#'
#' @param x A p-vector.
#' @param y A p-vector.
#'
#' @return If x and y each have length p, the function returns a
#' p-vector with jth element equal to \frac{\partial}{\partial y_{j}}
#' d(x,y)
manhattan_dist_deriv_neg <- function(x, y) {
    derivs = ifelse(y <= x, -1, 1)
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

#' Plots a local biplot
#'
#' @param Y A matrix or data frame with the data the distances were computed from.
#' @param lb_list A list giving the output from compute_lb
#'
#' @return A plot.
#' @export
local_biplot <- function(X, dist, dist_deriv = NULL, k = 2,
                         samples = c(),
                         n_random_points = 0,
                         new_points = list(),
                         alpha = 1) {
    dist_fns = make_dist_fns(dist, dist_deriv)
    mds_matrices = make_mds_matrices(X, dist_fns$dist_fn)
    lb_dfs = list()
    if(length(samples) > 0) {
        lb_dfs[["original"]] = compute_lb_samples(
            mds_matrices, dist_fns, k = k, samples = samples
        )
        
    }
    if(length(new_points) > 0 | n_random_points > 0) {
        lb_dfs[["new"]] = compute_lb_new_points(
            mds_matrices, dist_fns, k = k,
            n_random_points = n_random_points,
            new_points = new_points,
            alpha = alpha)
    }
    return(Reduce(rbind, lb_dfs))
}

#' Make a correlation biplot
#' @export
correlation_biplot <- function(X, dist, plotting_axes = 1:2) {
    mds_matrices = make_mds_matrices(X, dist)
    biplot_axes = cor(X, mds_matrices$Y)[,plotting_axes]
    colnames(biplot_axes) = sapply(plotting_axes, function(i) sprintf("Axis%i", i))
    return(biplot_axes)
}
