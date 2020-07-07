# Local biplots vignette



In this vignette, we demonstrate how to use the `mdessence` package to make local biplots.
The primary function in the package is `local_biplot`, and it can be used either with one of the built-in distances or with a distance function that the user provides.
Both methods are demonstrated.


We start off by loading some additional packages that we need.
`Matrix`, `phyloseq`, and `ape` are all related to the computations for the custom distance we will make.
`rmutil ` is used to create the simulated data, and the remaining packages are for plotting.
Some extra functions for simulating data are in `data-simulation-functions.R`, and some functions related to the custom distances are in `unifrac-distances-and-derivs.R`.

```r
library(mdessence)
library(Matrix)
library(phyloseq)
library(ape)
library(rmutil)
library(ggtree)
library(gridExtra)
library(viridis)
library(latex2exp)
source("data-simulation-functions.R")
source("unifrac-distances-and-derivs.R")
```

## Simulating the example data

We start off by making the example data.
The dataset is supposed to mimic microbiome data, so we have an abundance matrix with $2^6$ variables and 20 samples.
The samples fall into two groups of ten each.
The variables are supposed to represent bacterial taxa, and we create a phylogenetic tree describing the relationship between the taxa.


```r
set.seed(0)
example_data = make_example_data(n_per_group = 10, log_2_p = 6, deep_split_scaling = 3, overall_scaling = 1, s = .1)
example_data_phyloseq = make_phyloseq_from_example_data(example_data, mean_matrix = FALSE)
```

The idea behind the simulation is taxa in one half of the tree are over-represented in one group, and taxa in the other half of the tree are over-represented in the other group.
In addition, the specific taxa represented in each group are non-overlapping: if a taxon is represented in one group, its sister taxon (the closest taxon on the phylogenetic tree) is not represented in that group but it is represented in the other group.
This pattern can be seen more easily in a plot of the data matrix with the phylogenetic tree below.

```r
variable_palette = viridis(4)[c(1,2,4,3)]
variable_labels = c(TeX("(\\textbf{1}_{deep})_i = 1, (\\textbf{1}_{shallow})_i = 1"),
                       TeX("(\\textbf{1}_{deep})_i = 1, (\\textbf{1}_{shallow})_i = 0"),
                       TeX("(\\textbf{1}_{deep})_i = 0, (\\textbf{1}_{shallow})_i = 1"),
                       TeX("(\\textbf{1}_{deep})_i = 0, (\\textbf{1}_{shallow})_i = 0"))
names(variable_palette) = names(variable_labels) = c("TRUETRUE", "TRUEFALSE", "FALSETRUE", "FALSEFALSE")
p_tree = ggtree(example_data_phyloseq) +
    geom_tippoint(aes(color=var_description), size = 1.5) +
    scale_color_manual(values = variable_palette, labels = variable_labels, "Variable Type")

tip_positions = unique(subset(p_tree$data, isTip)[,c("label", "y")])
example_data_melted = melt(example_data$C)
abundance_and_position = merge(example_data_melted,
                               tip_positions,
                               by.x = "Var2", by.y = "label",
                               all.x = TRUE)
## a hack to make the colors show up better, set the zero values to NA
abundance_and_position_NAs = within(abundance_and_position, value[value == 0] <- NA)
p_tree + 
    geom_tile(aes(x = Var1 / 2 + 6, y = y, fill = value), color = "darkgray", size = .1,
              data = abundance_and_position_NAs) +
    scale_fill_distiller(palette = "Spectral", na.value = "white", "Abundance") +
    theme(legend.position = "right") +
    coord_flip() + scale_x_reverse() +
    guides(color = guide_legend(order = 1),
           fill_distiller = guide_legend(order = 2))
```

![plot of chunk plot-data-and-tree](figure/plot-data-and-tree-1.png)

## Creating local biplots with a built-in distance

Next we make local biplots with the Manhattan distance, which is built in to the package.
The first argument to the function is the data matrix, formatted as samples by variables.
The second argument is the type of distance to use: in this case, we are using the Manhattan/l1 distance, and the positive version of local biplot axes.
By default local biplot axes are computed at of the original samples, but this can be modified by the `samples`, `new`, or `random` arguments.

```r
sb_l1 = local_biplot(X = example_data$C, dist = "manhattan-pos")
```

The output from `local_biplot` is a data frame that gives information about the embedding of the samples with the specified distance and the local biplot axes at a collection of points.
The `Embedding1` and `Embedding2` columns record where the corresponding sample embeds in the MDS space (or, if using a different set of embedding dimensions, up to `Embeddingn`).
There is a good deal of replication in the `Embedding` columns: each sample shows up in multiple rows, and each row corresponding to the same sample will have the same values for `Embedding1` and `Embedding2`.

The local biplot axes are given in the `Axis1` and `Axis2` columns (or, if using a different set of embedding dimensions, up to `Axisn`).
Each local biplot axis corresponds to a variable (given in the `variable` column), and a vector in the input data space.
The `sample` column records the position in the input data space the local biplot axes correspond to: `Originaln` means that the local biplot axis for that row corresponds to the nth sample in the input dataset.
The function also has the capability to compute local biplot axes at other points in the data space, or at random points, in which case the `sample` column would be `Newn`, for the nth new sample.

```r
head(sb_l1)
```

```
##    Embedding1 Embedding2      Axis1      Axis2 variable    sample
## t1  -314.9535   -15.8234  1.5666791 -0.4101203       t1 Original1
## t2  -314.9535   -15.8234 -0.8744008 -0.8474948       t2 Original1
## t3  -314.9535   -15.8234  1.5666791 -0.4101203       t3 Original1
## t4  -314.9535   -15.8234 -0.8711629  0.8671601       t4 Original1
## t5  -314.9535   -15.8234  1.5666791 -0.4101203       t5 Original1
## t6  -314.9535   -15.8234 -0.9685534 -0.4386861       t6 Original1
```

We can first plot the MDS embeddings of the samples using the data contained in `sb_l1`:

```r
group_a = paste0("Original", 1:10)
group_labels = c("a","b")
names(group_labels) = c("TRUE", "FALSE")
ggplot(subset(sb_l1, variable == "t1")) +
    geom_point(aes(x = Embedding1, y = Embedding2, shape = sample %in% group_a)) +
    scale_shape(name = "Sample type", labels = group_labels) +
    ggtitle("Sample Embeddings\nMDS with Manhattan Distance") +
    xlab("Axis 1") + ylab("Axis 2") + coord_fixed() +
    theme(plot.title = element_text(size = 12))
```

![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3-1.png)


In this particular example, we have more information about the variables in the `example_data$variable_descriptions` data frame.
To plot this information along with the local biplot axes, we merge it with `sb_l1`.

```r
head(example_data$variable_descriptions)
```

```
##   variable shallow deep var_description
## 1       t1    TRUE TRUE        TRUETRUE
## 2       t2   FALSE TRUE       TRUEFALSE
## 3       t3    TRUE TRUE        TRUETRUE
## 4       t4   FALSE TRUE       TRUEFALSE
## 5       t5    TRUE TRUE        TRUETRUE
## 6       t6   FALSE TRUE       TRUEFALSE
```

```r
sb_and_var_descriptions_l1 = merge(sb_l1, example_data$variable_descriptions, by = "variable")
```

Then we can plot the local biplot axes.
In this case, we plot the set of local biplot axes corresponding to one sample around the point where that sample embeds in MDS space.
We also scale up the axes so that we can see them more clearly.

```r
ggplot(sb_and_var_descriptions_l1,
       aes(x = Embedding1, y = Embedding2,
           xend = Axis1 * 30 + Embedding1,
           yend = Axis2 * 30 + Embedding2,
           color = var_description)) +
    geom_segment(alpha = 1) +
    geom_point(aes(shape = sample %in% group_a),
               data = sb_and_var_descriptions_l1, color = "black") + 
    scale_color_manual(values = variable_palette, labels = variable_labels, "Variable Type") +
    scale_shape("Sample Type", labels = group_labels) +
    ggtitle("Local biplot axes for Manhattan distance") +
    xlab("Axis 1") + ylab("Axis 2") + coord_fixed() + xlim(c(-590, 590))
```

![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-5-1.png)

## Creating local biplots with a custom distance


We can also create local biplots with custom distances.
To do this, we need to provide a function that computes the distance and a function that computes the gradient of the distance.
The rules are:
- The distance function should take a matrix and compute the distance between the rows of the matrix. The output can be formatted either as a `dist` object or as a matrix.
- The gradient function should take two vectors, `x` and `y`, and compute $\frac{\partial}{\partial y_i}d(x, y)$. The output should be a vector of the same length as `x` or `y`.

In the example below, we create the distance function `wuf_dist_ours` and the gradient function `wuf_dist_deriv_ours`.
These are based on functions that are defined in `unifrac-distances-and-derivs.R`.
Those functions require a specific tree, and the definition below tailors the distance to the tree in our simulated dataset.

```r
A = makeDescendantMatrix(example_data$tree)
L = getBranchLengths(example_data$tree)
wuf_dist_ours = function(X) wuf_dist(X, A, L)
wuf_dist_deriv_ours = function(x, y) wuf_dist_deriv(x, y, A, L)
```

When using the `local_biplot` function with a custom distance, we again provide the data matrix, but instead of specifying the distance by name, we provide the functions we created as values for `dist` and `dist_deriv`.

```r
sb_wuf = local_biplot(X = example_data$C, dist = wuf_dist_ours, dist_deriv = wuf_dist_deriv_ours)
```

As before, we can use the output from the `local_biplot` function to plot the samples in the MDS embedding space.

```r
ggplot(subset(sb_wuf, variable == "t1")) +
    geom_point(aes(x = Embedding1, y = Embedding2, shape = sample %in% group_a)) +
    xlab("Axis 1") + ylab("Axis 2") + coord_fixed() + ylim(c(-1000, 800)) +
    scale_shape(name = "Sample type", labels = group_labels) +
    ggtitle("Sample Embeddings\nMDS with Weighted Unifrac") +
    theme(plot.title = element_text(size = 12))
```

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-8-1.png)

We again add the extra information in `example_data$variable_descriptions` to the local biplot axes in `sb_wuf`.

```r
sb_and_var_descriptions_wuf = merge(sb_wuf, example_data$variable_descriptions, by = "variable")
```

We plot the local biplot axes, again centering each set of local biplot axes at the position where the corresponding sample embeds in MDS space, and scaling up the axes so that they are visible.

```r
ggplot(sb_and_var_descriptions_wuf,
       aes(x = Embedding1, y = Embedding2,
           xend = Axis1 * 10000 + Embedding1,
           yend = Axis2 * 10000 + Embedding2,
           color = var_description)) +
    geom_segment() +
    geom_point(aes(shape = sample %in% group_a),
               data = sb_and_var_descriptions_wuf, color = "black") + 
    scale_color_manual(values = variable_palette, labels = variable_labels, "Variable Type") +
    scale_shape("Sample Type", labels = group_labels) +
    ggtitle("Local biplot axes for weighted UniFrac") +
    xlab("Axis 1") + ylab("Axis 2") + coord_fixed()
```

![plot of chunk unnamed-chunk-10](figure/unnamed-chunk-10-1.png)


## Interactive version

In some cases, a large number of samples and variables make it hard to visualize the biplot axes.
To help in this situation, we also provide a function that will create an interactive plot with the samples and the biplot axes.

The `interactive_biplot` function takes a data matrix, either a string describing the distance or a pair `dist` and `dist_deriv` (like the `local_biplot`) function.
In addition, it has arguments describing addition data about the samples and variables, and how to plot the samples and the biplot axes.

In the code below, we are creating an interactive biplot with the same simulated data as before and the Manhattan distance.
The `var_data` argument gives the extra information about the variables, and the `var_mapping` argument tells the function to plot points corresponding to the first and second local biplot axes, and to color the local biplot axes according to the value of the variable `shallow` in `example_data$variable_descriptions`.


```r
interactive_biplot(X = example_data$C, dist = "manhattan-pos",
                  var_mapping = aes_string(x = "Axis1", y = "Axis2", color = "shallow"),
                  var_data = example_data$variable_descriptions)
```

Analogously, in the code below, we are asking for an interactive plot of the local biplot axes with the weighted UniFrac distance, and we want the local biplot axes to be colored according to the value of the `deep` column in `example_data$variable_descriptions`.

```r
interactive_biplot(X = example_data$C, dist = wuf_dist_ours,
                  dist_deriv = wuf_dist_deriv_ours,
                  var_mapping = aes_string(x = "Axis1", y = "Axis2", color = "deep"),
                  var_data = example_data$variable_descriptions)
```

In both cases, running the code above will cause a browser window to open.
The left panel will have the sample embeddings, and clicking on one of the sample points on the left plot will lead to the biplot axes corresponding to that point to show up on the right panel.
