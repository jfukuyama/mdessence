#' Interactive biplot
#'
#' @import shiny
#' @export
interactiveBiplot <- function(Y, dist, dist_deriv = NULL, k = 2, axes = 1:2, samples = NULL, ...) {
    dist_fns = make_dist_fns(dist, dist_deriv)
    mds_matrices = make_mds_matrices(Y, dist_fns$dist_fn)
    embedding_and_sample_num = data.frame(mds_matrices$X, sample_num = 1:nrow(mds_matrices$X))
    sensitivity_list = compute_sensitivity(mds_matrices, dist_fns, k = k, samples = samples)
    varnames = get_variable_names(mds_matrices$Y)
    ui <- fluidPage(
        headerPanel("Interactive Sensitivity Biplot"),
        sidebarPanel(width = 2,
                     actionButton("done", "Done")),
        mainPanel(width = 10,
                  fluidPage(
                      fluidRow(
                          column(5, plotOutput("plot_samples", click = "plot_click")),
                          column(5, plotOutput("plot_variables")))
                  )
                  )
    )

    server <- function(input, output, session) {
        output$plot_samples = renderPlot({
            ggplot(embedding_and_sample_num) + geom_point(aes(x = Axis1, y = Axis2))
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
            ggplot(biplot_df) + geom_point(aes(x = xend, y = yend))
        })

        observeEvent(input$done, {
            biplot_df = make_biplot_data_frame(sensitivity_list, mds_matrices, axes, samples)
            stopApp(biplot_df)
        })
    }
    
    runGadget(ui, server)
}
