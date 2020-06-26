#' Interactive biplot
#'
#' @import shiny
#' @export
interactiveBiplot <- function(X, dist, dist_deriv = NULL, k = 2, axes = 1:2, samples = NULL,
                              sample_data = NULL,
                              sample_mapping = aes_string(x = "Axis1", y = "Axis2"),
                              sample_facet = NULL,
                              var_data = NULL,
                              var_mapping = aes_string(x = "Axis1", y = "Axis2"),
                              layout = c(5,5), ...) {
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
        dist_fns = make_dist_fns(dist, dist_deriv)
        mds_matrices = make_mds_matrices(X, dist_fns$dist_fn)
        embedding_and_sample_name = data.frame(mds_matrices$Y, sample = paste0("Original", 1:nrow(mds_matrices$Y)))
        sensitivity_df = compute_sensitivity_samples(mds_matrices, dist_fns, k = k, samples = samples)
        varnames = get_variable_names(mds_matrices$X)
        output$plot_samples = renderPlot({
            if(!is.null(sample_data))
                p = ggplot(data.frame(embedding_and_sample_name, sample_data), sample_mapping) +
                    geom_point()
                else
                    p = ggplot(embedding_and_sample_name, sample_mapping) +
                        geom_point()
                if(!is.null(sample_facet))
                    p = p + sample_facet
                p
        })
        output$plot_variables = renderPlot({
            biplot_center = nearPoints(embedding_and_sample_name, input$plot_click, maxpoints = 1)
            sample_name = biplot_center$sample
            if(!(sample_name %in% sensitivity_df$sample)) {
                sample_num = as.numeric(gsub("Original", "", sample_name))
                sensitivity_df = rbind(sensitivity_df,
                                        compute_sensitivity_samples(mds_matrices, dist_fns, k, samples = sample_num))
            }
            if(!is.null(var_data)) {
                p = ggplot(data.frame(subset(sensitivity_df, sample == sample_name), var_data), var_mapping) +
                    geom_point()
            } else {
                p = ggplot(subset(sensitivity_df, sample == sample_name), var_mapping) +
                    geom_point()
            }
            p
        })

        observeEvent(input$done, {
            stopApp(sensitivity_df)
        })
    }
    
    runGadget(ui, server)
}
