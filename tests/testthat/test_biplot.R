context("Derivative tests")
x = 1:5
y = 5:1
test_that("Derivative functions have right size output", {
    expect_equal(length(euclidean_dist_deriv(x, y)), length(x))
    expect_equal(length(manhattan_dist_deriv_pos(x, y)), length(x))
    expect_equal(length(manhattan_dist_deriv_neg(x, y)), length(x))
    expect_equal(length(maximum_dist_deriv(x, y)), length(x))
})

test_that("Max dist derivative gives correct results", {
    expect_equal(maximum_dist_deriv(x, y), c(1, 0, 0, 0, -1))
})

context("Distance function creation")
euc_dist_fns = make_dist_fns(dist_fn = "euclidean", dist_deriv = NULL)
manhattan_dist_fns = make_dist_fns(dist_fn = "manhattan-pos", dist_deriv = NULL)
max_dist_fns = make_dist_fns(dist_fn = "manhattan-pos", dist_deriv = NULL)
jaccard_dist_fns = make_dist_fns(dist_fn = "jaccard", dist_deriv = NULL)
test_that("make_dist_fns gives correct type of output", {
    expect_equal(typeof(euc_dist_fns), "list")
    expect_equal(length(euc_dist_fns), 2)
    expect_equal(typeof(euc_dist_fns$dist_fn), "closure")
    expect_equal(typeof(euc_dist_fns$dist_deriv), "closure")
    expect_equal(typeof(manhattan_dist_fns), "list")
    expect_equal(length(manhattan_dist_fns), 2)
    expect_equal(typeof(manhattan_dist_fns$dist_fn), "closure")
    expect_equal(typeof(manhattan_dist_fns$dist_deriv), "closure")
    expect_equal(typeof(max_dist_fns), "list")
    expect_equal(length(max_dist_fns), 2)
    expect_equal(typeof(max_dist_fns$dist_fn), "closure")
    expect_equal(typeof(max_dist_fns$dist_deriv), "closure")

})


context("MDS matrix creation functions")
X = matrix(rnorm(20), nrow = 5, ncol = 4)
mds_matrices = make_mds_matrices(X = X, euc_dist_fns$dist_fn)
test_that("make_mds_matrices works", {
    ## delta is the matrix of squared distances between the samples
    expect_equal(nrow(mds_matrices$delta), nrow(X))
    ## X is the input data
    expect_equal(X, mds_matrices$X)
    ## Lambda are the eigenvalues
    expect_equal(length(mds_matrices$Lambda), nrow(X) - 1)
    ## Y gives the embedding of the samples in the MDS space
    expect_equal(nrow(mds_matrices$Y), nrow(X))
    ## jdj is the row- and column-centered version of -.5 * delta
    expect_equal(rowSums(mds_matrices$jdj), rep(0, nrow(X)))
    expect_equal(colSums(mds_matrices$jdj), rep(0, nrow(X)))
})

