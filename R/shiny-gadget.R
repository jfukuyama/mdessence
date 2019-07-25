#' Interactive biplot
#'
#' @import shiny
#' @export
interactiveBiplot <- function(Y, dist, dist_deriv = NULL, k = 2, axes = 1:2, samples = NULL,
                              sample_data = NULL,
                              sample_mapping = aes_string(x = "Axis1", y = "Axis2"),
                              sample_facet = NULL,
                              var_data = NULL,
                              var_mapping = aes_string(x = "xend", y = "yend"),
                              layout = c(5,5), ...) {
    dist_fns = make_dist_fns(dist, dist_deriv)
    mds_matrices = make_mds_matrices(Y, dist_fns$dist_fn)
    embedding_and_sample_num = data.frame(mds_matrices$X, sample_num = 1:nrow(mds_matrices$X))
    sensitivity_list = compute_sensitivity(mds_matrices, dist_fns, k = k, samples = samples)
    varnames = get_variable_names(mds_matrices$Y)
    ui <- fluidPage(
        headerPanel("Interactive Sensitivity Biplot"),
        sidebarPanel(width = 1,
                     actionButton("done", "Done")),
        mainPanel(width = 11,
                  fluidPage(
                      fluidRow(
                          column(layout[1], plotOutput("plot_samples", click = "plot_click")),
                          column(layout[2], plotOutput("plot_variables")))
                  )
                  )
    )

    server <- function(input, output, session) {
        output$plot_samples = renderPlot({
            if(!is.null(sample_data))
                p = ggplot(data.frame(embedding_and_sample_num, sample_data), sample_mapping) +
                    geom_point()
                else
                    p = ggplot(embedding_and_sample_num, sample_mapping) +
                        geom_point()
                if(!is.null(sample_facet))
                    p = p + sample_facet
                p
        })
        output$plot_variables = renderPlot({
            biplot_center = nearPoints(embedding_and_sample_num, input$plot_click, maxpoints = 1)
            sample_num = biplot_center$sample_num
            if(!(sample_num %in% samples)) {
                samples <<- c(samples, sample_num)            
                sensitivity_list[[sample_num]] <<-
                    compute_sensitivity(mds_matrices, dist_fns, k, sample_num)[[sample_num]]
            }
            biplot_df = make_one_sample_biplot_df(
                sensitivity = sensitivity_list[[sample_num]][,axes],
                center = mds_matrices$X[sample_num, axes],
                sample = sample_num, varnames = varnames)
            if(!is.null(var_data))
                p = ggplot(data.frame(biplot_df, var_data), var_mapping) + geom_point()
            else
                p = ggplot(biplot_df, var_mapping) +
                    geom_point()
            p
        })

        observeEvent(input$done, {
            biplot_df = make_biplot_data_frame(sensitivity_list, mds_matrices, axes, samples)
            stopApp(biplot_df)
        })
    }
    
    runGadget(ui, server)
}
