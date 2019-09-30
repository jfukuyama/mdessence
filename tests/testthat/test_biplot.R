context("Derivative tests")
x = 1:5
y = 5:1
test_that("Derivative functions have right size output", {
    expect_equal(length(euclidean_dist_deriv(x, y)), length(x))
    expect_equal(length(manhattan_dist_deriv(x, y)), length(x))
    expect_equal(length(maximum_dist_deriv(x, y)), length(x))
})

test_that("Max dist derivative gives correct results", {
    expect_equal(maximum_dist_deriv(x, y), c(1, 0, 0, 0, -1))
})

context("Distance function creation")
euc_dist_fns = make_dist_fns(dist_fn = "euclidean", dist_deriv = NULL)
manhattan_dist_fns = make_dist_fns(dist_fn = "manhattan", dist_deriv = NULL)
max_dist_fns = make_dist_fns(dist_fn = "manhattan", dist_deriv = NULL)
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
    expect_equal(typeof(jaccard_dist_fns), "list")
    expect_equal(length(jaccard_dist_fns), 2)
    expect_equal(typeof(jaccard_dist_fns$dist_fn), "closure")
    expect_equal(typeof(jaccard_dist_fns$dist_deriv), "closure")

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
    expect_equal(ncol(mds_matrices$Y), nrow(X) - 1)
    expect_equal(nrow(mds_matrices$Y), nrow(X))
    ## jdj is the row- and column-centered version of -.5 * delta
    expect_equal(rowSums(mds_matrices$jdj), rep(0, nrow(X)))
    expect_equal(colSums(mds_matrices$jdj), rep(0, nrow(X)))
})

context("Sensitivity computation")
sample_subset = c(1,5)
sensitivity_output = compute_sensitivity(mds_matrices, euc_dist_fns, k = 2, samples = 1:nrow(X))
test_that("Sensitivity computations are correct", {
    ## check that things are the right size
    expect_equal(length(sensitivity_output), nrow(X))
    expect_equal(nrow(sensitivity_output[[1]]), ncol(X))
    expect_equal(ncol(sensitivity_output[[1]]), 2)
    ## for the Euclidean distance, the sensitivities are the same at each sample
    expect_equal(sensitivity_output[[1]], sensitivity_output[[2]])
    expect_equal(sensitivity_output[[3]], sensitivity_output[[4]])
})
sensitivity_sample_output = compute_sensitivity(mds_matrices, euc_dist_fns, k = 2, samples = sample_subset)
test_that("'samples' argument to compute_sensitivity works", {
    expect_equal(sensitivity_sample_output[[1]], sensitivity_output[[1]])
    expect_equal(sensitivity_sample_output[[5]], sensitivity_output[[5]])
    expect_null(sensitivity_sample_output[[2]])
})

context("Making the biplot data frame")
biplot_df = make_biplot_data_frame(sensitivity_output, mds_matrices, samples = 1:nrow(X))
biplot_sub_df = make_biplot_data_frame(sensitivity_sample_output, mds_matrices, samples = sample_subset)
test_that("Biplot data frame made correctly", {
    expect_equal(ncol(biplot_df), 6)
    expect_equal(colnames(biplot_df), c("x", "y", "xend", "yend", "variable", "sample"))
    ## one row per variable per sample
    expect_equal(nrow(biplot_df), ncol(X) * nrow(X))
    ## one row per variable for the samples we compute biplot axes for
    ## plus one row for each of the samples we don't compute biplot
    ## axes for
    expect_equal(nrow(biplot_sub_df),
                 ncol(X) * length(sample_subset) + nrow(X) - length(sample_subset))
#    expect_equal(subset(biplot_sub_df, sample == 1), subset(biplot_df, sample == 1))
    expect_equal(matrix(subset(biplot_sub_df, sample == 5)),
                 matrix(subset(biplot_df, sample == 5)))
})

context("Discrete distance tests")
test_that("Jaccard computations are correct", {
    expect_equal(jaccard_dis(c(1,1), c(3,4)), 0)
    expect_equal(jaccard_dis(c(0,0), c(0,0)), 0)
    expect_equal(jaccard_dis(c(1,1), c(0,4)), .5)
})
test_that("Derivative function for categorical data works", {
    ## output is of the right size
    expect_equal(length(jaccard_deriv(x, y)), length(y))
    ## derivative in positive direction is zero if the coordinate is non-zero
    expect_equal(jaccard_deriv(c(0,1), c(3,4), delta_min = 1, positive = TRUE), c(0,0))
    ## derivative in the negative direction is non-zero if the coordinate is non-zero
    expect_true(all(jaccard_deriv(c(0,1), c(3,4), delta_min = 1, positive = FALSE) != 0))
    ## derivative in the negative direction is zero if the coordinate is zero
    expect_true(all(jaccard_deriv(c(0,1), c(0,0), delta_min = 1, positive = FALSE) == 0))
}
)

Y_count = matrix(rpois(30, lambda = 2), nrow = 6, ncol = 5)
mds_matrices_count = make_mds_matrices(Y_count, jaccard_dist_fns$dist_fn)
sensitivity_output_count = compute_sensitivity(mds_matrices_count,
    jaccard_dist_fns, k = 2, samples = 1:nrow(Y_count))
test_that("Sensitivity computations on categorical data work", {
    ## check that things are the right size
    expect_equal(length(sensitivity_output_count), nrow(Y_count))
    expect_equal(nrow(sensitivity_output_count[[1]]), ncol(Y_count))
    expect_equal(ncol(sensitivity_output_count[[1]]), 2)
})
