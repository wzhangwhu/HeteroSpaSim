#' Launch the optional HeteroSpaSim Shiny interface
#'
#' @return A Shiny application object, launched interactively.
#' @export
run_heterospa_app <- function() {
    if (!requireNamespace("shiny", quietly = TRUE)) {
        stop("Install the shiny package to use the graphical interface.", call. = FALSE)
    }
    ui <- shiny::fluidPage(
        shiny::titlePanel("HeteroSpaSim"),
        shiny::sidebarLayout(
            shiny::sidebarPanel(
                shiny::selectInput("scenario", "Scenario", choices = c("easy", "moderate", "hard"), selected = "moderate"),
                shiny::checkboxGroupInput("modalities", "Modalities", choices = c("RNA", "ATAC", "Protein"), selected = c("RNA", "ATAC", "Protein")),
                shiny::numericInput("n_cells", "Latent cells", value = 2000, min = 100, step = 100),
                shiny::numericInput("seed", "Seed", value = 1, min = 1, step = 1),
                shiny::actionButton("run", "Simulate", class = "btn-primary"),
                shiny::hr(),
                shiny::downloadButton("download_rds", "Download RDS")
            ),
            shiny::mainPanel(
                shiny::tabsetPanel(
                    shiny::tabPanel("Tissue", shiny::plotOutput("tissue_plot", height = "650px")),
                    shiny::tabPanel("Modalities", shiny::plotOutput("modality_plot", height = "650px")),
                    shiny::tabPanel("Summary", shiny::tableOutput("summary_table"))
                )
            )
        )
    )
    server <- function(input, output, session) {
        sim <- shiny::eventReactive(input$run, {
            mods <- input$modalities
            shiny::validate(shiny::need(length(mods) >= 1L, "Select at least one modality."))
            simulate_spatial_multiomics(
                scenario = input$scenario,
                modalities = mods,
                n_cells = as.integer(input$n_cells),
                seed = as.integer(input$seed)
            )
        }, ignoreNULL = FALSE)
        output$tissue_plot <- shiny::renderPlot(plot(sim(), what = "domains"))
        output$modality_plot <- shiny::renderPlot(plot(sim(), what = "modalities"))
        output$summary_table <- shiny::renderTable(summary(sim()), digits = 3)
        output$download_rds <- shiny::downloadHandler(
            filename = function() paste0("HeteroSpaSim_", input$scenario, "_seed", input$seed, ".rds"),
            content = function(file) saveRDS(sim(), file)
        )
    }
    shiny::shinyApp(ui, server)
}
