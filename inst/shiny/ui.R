# BRA (BioRangeAnalyzer Shiny) - User Interface

# Create the UI function
shinyUI(
  navbarPage(
    title = "BRA (BioRangeAnalyzer Shiny)",
    id = "navbar",
    theme = "bootstrap",
    
    # CSS styling
    tags$head(
      tags$style(HTML("
        body {
          font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
          font-size: 14px;
        }
        .navbar {
          font-size: 16px !important;
        }
        .navbar-default .navbar-nav > li > a {
          font-size: 16px !important;
          font-weight: 500 !important;
          color: #ffffff !important;
          padding: 15px 20px !important;
        }
        .navbar-default .navbar-brand {
          font-size: 18px !important;
          font-weight: bold !important;
          color: #ffffff !important;
        }
        .navbar-default {
          background-color: #2c5aa0 !important;
          border-color: #2c5aa0 !important;
          box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }
        .navbar-default .navbar-nav > li > a:hover,
        .navbar-default .navbar-nav > li > a:focus {
          background-color: #1e3f5a !important;
          color: #ffffff !important;
        }
        .step-title {
          font-size: 28px;
          font-weight: bold;
          color: #2c5aa0;
          margin-bottom: 15px;
        }
        .guidance-text {
          background-color: #e8f4f8;
          padding: 15px;
          border-left: 5px solid #2c5aa0;
          margin-bottom: 20px;
          font-size: 14px;
          border-radius: 4px;
        }
        h4 {
          font-size: 18px;
          font-weight: 600;
          color: #2c5aa0;
        }
        .upload-status-success {
          background-color: #d4edda;
          border: 1px solid #c3e6cb;
          color: #155724;
          padding: 12px 15px;
          border-radius: 4px;
          margin-top: 10px;
          font-weight: 500;
        }
        .upload-status-error {
          background-color: #f8d7da;
          border: 1px solid #f5c6cb;
          color: #721c24;
          padding: 12px 15px;
          border-radius: 4px;
          margin-top: 10px;
          font-weight: 500;
        }
        .upload-status-info {
          background-color: #d1ecf1;
          border: 1px solid #bee5eb;
          color: #0c5460;
          padding: 12px 15px;
          border-radius: 4px;
          margin-top: 10px;
          font-weight: 500;
        }
      "))
    ),
    
    # Tab 1: Welcome
    tabPanel(
      "Welcome",
      icon = icon("home"),
      div(
        class = "container-fluid",
        style = "padding: 40px;",
        h1("BRA (BioRangeAnalyzer Shiny)"),
        h3("Ancestral Range Estimation and Biogeographic Analysis"),
        hr(),
        p("This application provides a user-friendly interface for:"),
        tags$ul(
          tags$li("Loading occurrence data and phylogenetic trees"),
          tags$li("Extrapolating species ranges using multiple methods"),
          tags$li("Visualizing distributions and ancestral ranges"),
          tags$li("Running BioGeoBEARS analyses"),
          tags$li("Exporting results in multiple formats (NEXUS, TNT, NDM)")
        ),
        hr(),
        h3("Getting Started"),
        p("Follow these steps:"),
        tags$ol(
          tags$li("Go to 'Data Input' and upload your occurrence data (CSV/TXT) and phylogenetic tree"),
          tags$li("Go to 'Tree Validation' to check your tree properties"),
          tags$li("Go to 'Range Extrapolation' to create species range maps"),
          tags$li("Go to 'Visualizations' to view your results"),
          tags$li("Go to 'Export Files' to download your results"),
          tags$li("Go to 'BioGeoBEARS Setup' to configure your analysis"),
          tags$li("Go to 'Analysis & Results' to run the analysis")
        )
      )
    ),
    
    # Tab 2: Data Input
    tabPanel(
      "Data Input",
      icon = icon("upload"),
      div(
        class = "container-fluid",
        style = "padding: 20px;",
        div(
          class = "step-title",
          "Step 1: Data Input"
        ),
        div(
          class = "guidance-text",
          p("Upload your occurrence data (CSV/TXT), phylogenetic tree (Newick or Nexus format), and study area shapefile (optional).")
        ),
        div(
          class = "row",
          div(
            class = "col-md-6",
            h4("Occurrence Data"),
            div(
              class = "guidance-text",
              p(strong("Required CSV/TXT Format:")),
              p("Your CSV/TXT file MUST have exactly 3 columns (in any order):")
            ),
            tags$ul(
              tags$li(strong("spp"), " - Species name (e.g., 'Homo_sapiens')"),
              tags$li(strong("long"), " - Longitude in decimal degrees (e.g., -50.5)"),
              tags$li(strong("lat"), " - Latitude in decimal degrees (e.g., -25.3)")
            ),
            p("Example CSV/TXT content:"),
            pre("spp,long,lat\nSpecies1,-50.5,-25.3\nSpecies1,-51.2,-26.1\nSpecies2,-48.3,-24.5"),
            fileInput("occurrence_file", "Choose CSV/TXT file:", accept = c("text/csv", "text/plain", ".csv", ".txt")),
            actionButton("load_occurrence", "Load Data", class = "btn btn-primary"),
            br(), br(),
            htmlOutput("data_status")
          ),
          div(
            class = "col-md-6",
            h4("Phylogenetic Tree"),
            p("Newick or Nexus format file"),
            fileInput("tree_file", "Choose tree file:", accept = c(".nwk", ".newick", ".nex", ".nexus", ".txt", ".tre")),
            actionButton("load_tree", "Load Tree", class = "btn btn-primary"),
            br(), br(),
            htmlOutput("tree_load_status")
          )
        ),
        div(
          class = "row",
          style = "margin-top: 30px;",
          div(
            class = "col-md-12",
            h4("Study Area Shapefile (Optional)"),
            div(
              class = "guidance-text",
              p(strong("Use this shapefile to:")),
              tags$ul(
                tags$li("Delimit where extrapolation runs"),
                tags$li("Filter occurrence points outside the study area"),
                tags$li("Create regular grid cells within the boundary")
              ),
              p("You can use a simple boundary polygon (e.g., Brazil outline only) to speed up analyses."),
              p(strong("IMPORTANT: Select ALL shapefile files at once!")),
              p("Your shapefile consists of multiple files. You MUST select all of them together:"),
              tags$ul(
                tags$li(strong(".shp"), " - Geometry (required)"),
                tags$li(strong(".shx"), " - Shape index (required)"),
                tags$li(strong(".dbf"), " - Attributes (required)"),
                tags$li(strong(".prj"), " - Projection (optional but recommended)")
              ),
              p("How to select multiple files: Hold Ctrl (or Cmd on Mac) and click each file, then click Open."),
              p("Example: Select America_Sul.shp, America_Sul.shx, America_Sul.dbf, America_Sul.prj all together.")
            ),
            fileInput("study_area_shapefile", "Choose ALL shapefile files together:", accept = c(".shp", ".dbf", ".shx", ".prj"), multiple = TRUE),
            actionButton("load_shapefile", "Load Shapefile", class = "btn btn-primary"),
            br(), br(),
            htmlOutput("shapefile_status")
          )
        )
      )
    ),
    
    # Tab 3: Tree Validation
    tabPanel(
      "Tree Validation",
      icon = icon("sitemap"),
      div(
        class = "container-fluid",
        style = "padding: 20px;",
        div(
          class = "step-title",
          "Step 2: Tree Validation"
        ),
        div(
          class = "guidance-text",
          p("Inspect your phylogenetic tree and check branch-length availability."),
          p("Important: for BioGeoBEARS analyses, use an ultrametric and time-calibrated (dated) tree.", style = "color: #8a6d3b; font-weight: 600;"),
          p("Tip: Use the checkbox to toggle node numbers. Zoom (mouse wheel), pan (drag), and hover to explore.", style = "font-size: 12px; color: #666;")
        ),
        div(
          class = "row",
          div(
            class = "col-md-12",
            h4("Tree Properties"),
            actionButton("validate_tree", "Inspect Tree", class = "btn btn-primary"),
            br(), br(),
            verbatimTextOutput("tree_validation_output")
          )
        ),
        br(),
        h4("Tree Visualization"),
        div(
          style = "margin-bottom: 10px;",
          checkboxInput("show_node_labels", "Show internal node numbers", value = TRUE),
          p("Tip: Uncheck to reduce visual clutter. Zoom with mouse wheel, pan with click+drag, hover to see node IDs.", style = "font-size: 11px; color: #666; margin-top: 5px;")
        ),
        ggiraph::girafeOutput("tree_plot", height = "650px"),
        br(),
        h4("Optional Clade Filter"),
        p("Optionally restrict downstream extrapolation to a single clade (all descendants of one internal node)."),
        checkboxInput("use_clade_filter", "Use a specific clade in next steps", value = FALSE),
        conditionalPanel(
          condition = "input.use_clade_filter == true",
          selectizeInput("clade_node_id", "Internal node ID:", choices = NULL, multiple = FALSE),
          verbatimTextOutput("clade_filter_status"),
          br(),
          h5("Clade Preview"),
          div(
            style = "margin-bottom: 10px;",
            checkboxInput("show_clade_node_labels", "Show internal node numbers", value = TRUE),
            p("Tip: Uncheck to reduce visual clutter. Zoom with mouse wheel, pan with click+drag, click node to see terminal taxa.", style = "font-size: 11px; color: #666; margin-top: 5px;")
          ),
          p("Interactive visualization of the selected clade:", style = "font-size: 11px; color: #666;"),
          ggiraph::girafeOutput("clade_tree_plot", height = "600px")
        ),

      )
    ),
    
    # Tab 4: Data Preprocessing
    tabPanel(
      "Data Preprocessing",
      icon = icon("filter"),
      div(
        class = "container-fluid",
        style = "padding: 20px;",
        div(
          class = "step-title",
          "Step 3: Data Preprocessing"
        ),
        div(
          class = "guidance-text",
          p("Remove problematic taxa (singletons and degenerate cases) before extrapolation.")
        ),
        div(
          class = "row",
          div(
            class = "col-md-12",
            h4("Filter Points Outside Study Area"),
            p("If you loaded a study area shapefile in Step 1, you can remove occurrence points that fall outside the shapefile boundary."),
            actionButton("filter_points_by_shapefile", "Remove Points Outside Shapefile", class = "btn btn-success"),
            br(), br(),
            verbatimTextOutput("shapefile_filter_output")
          )
        ),
        br(),
        div(
          class = "row",
          div(
            class = "col-md-12",
            h4("Problematic Taxa Detection"),
            p("The following taxa may cause problems in range extrapolation:"),
            tags$ul(
              tags$li("Singletons: Taxa with only 1 occurrence point"),
              tags$li("Degenerate: Taxa with 2 identical points (for convex hull)"),
              tags$li("Collinear: Taxa with all points on a line")
            ),
            actionButton("detect_problems", "Detect Problematic Taxa", class = "btn btn-info"),
            br(), br(),
            verbatimTextOutput("problem_detection_output")
          )
        ),
        br(),
        div(
          class = "row",
          div(
            class = "col-md-12",
            h4("Remove Problematic Taxa"),
            p("Choose what to remove, then apply the filter to the dataset (and tree tips, if provided):"),
            radioButtons(
              "problematic_removal_mode",
              "Removal strategy:",
              choices = list(
                "Remove singletons + doubletons" = "singletons_doubletons",
                "Remove only singletons" = "singletons_only",
                "Remove only doubletons" = "doubletons_only"
              ),
              selected = "singletons_doubletons"
            ),
            verbatimTextOutput("problematic_recommendation"),
            actionButton("remove_problems", "Remove Problematic Taxa", class = "btn btn-danger"),
            br(), br(),
            verbatimTextOutput("removal_output")
          )
        ),
        br(),
        div(
          class = "row",
          div(
            class = "col-md-12",
            h4("Remove Duplicates"),
            p("After detection, remove identical occurrence records within taxa (same spp, longitude and latitude):"),
            actionButton("remove_duplicates", "Remove Duplicates", class = "btn btn-warning"),
            br(), br(),
            verbatimTextOutput("duplicates_output")
          )
        ),
        br(),
        div(
          class = "row",
          div(
            class = "col-md-12",
            h4("Harmonize Tree and Data"),
            p("Align taxa between occurrence data and tree without applying singleton/doubleton removal."),
            actionButton("harmonize_tree_data", "Harmonize Tree <-> Data", class = "btn btn-default"),
            br(), br(),
            verbatimTextOutput("harmonization_output")
          )
        ),
        br(),
        div(
          class = "row",
          div(
            class = "col-md-6",
            h4("Cleaned Data Summary"),
            verbatimTextOutput("cleaned_data_summary")
          ),
          div(
            class = "col-md-6",
            h4("Updated Tree"),
            plotOutput("updated_tree_plot", height = "300px"),
            br(),
            downloadButton("download_pruned_tree", "Download Pruned Tree", class = "btn btn-success btn-sm")
          )
        )
      )
    ),
    
    # Tab 5: Range Extrapolation
    tabPanel(
      "Range Extrapolation",
      icon = icon("map"),
      div(
        class = "container-fluid",
        style = "padding: 20px;",
        div(
          class = "step-title",
          "Step 3: Range Extrapolation"
        ),
        div(
          class = "guidance-text",
          p("Choose how occurrence points will be converted into ranges/areas for matrix and map outputs."),
          p("Occurrence points mode has two workflows: (1) regular grid cells (requires grid resolution) or (2) direct counting in user irregular polygons (no regular grid).", style = "font-size: 12px; color: #666;")
        ),
        div(
          class = "row",
          div(
            class = "col-md-12",
            h4("Study Area (Shapefile)"),
            div(
              class = "guidance-text",
              p(strong("Shapefile 1 (required): study-area boundary")),
              p("Use this shapefile to delimit where extrapolation runs and where regular grid cells are created."),
              p("You can use a simple boundary polygon (e.g., Brazil outline only) to speed up analyses."),
              p(strong("IMPORTANT: Select ALL shapefile files at once!")),
              p("Your shapefile consists of multiple files. You MUST select all of them together:"),
              tags$ul(
                tags$li(strong(".shp"), " - Geometry (required)"),
                tags$li(strong(".shx"), " - Shape index (required)"),
                tags$li(strong(".dbf"), " - Attributes (required)"),
                tags$li(strong(".prj"), " - Projection (optional but recommended)")
              ),
              p(strong("How to select multiple files:"), " Hold Ctrl (or Cmd on Mac) and click each file, then click Open."),
              p("Example: Select America_Sul.shp, America_Sul.shx, America_Sul.dbf, America_Sul.prj all together.")
            ),
            fileInput("study_area_shapefile", "Choose ALL shapefile files together:", accept = c(".shp", ".dbf", ".shx", ".prj"), multiple = TRUE),
            htmlOutput("shapefile_status"),
            br()
          )
        ),
        div(
          class = "row",
          div(
            class = "col-md-8",
            h4("Extrapolation Method"),
            selectInput(
              "extrap_method",
              "Select method:",
              choices = list(
                "Occurrence points only" = "occurrence_only",
                "Buffer (circular buffers)" = "buffer",
                "Convex Hull (minimum convex polygon)" = "convex_hull",
                "Minimum Spanning Tree" = "mst"
              ),
              selected = "buffer"
            ),
            conditionalPanel(
              condition = "input.extrap_method == 'buffer'",
              numericInput("buffer_width", "Buffer width (km):", value = 100, min = 1),
              checkboxInput("use_mean_dist", "Use mean distance between points", value = FALSE)
            ),
            conditionalPanel(
              condition = "input.extrap_method == 'occurrence_only'",
              p("Build matrix from occurrence points. Choose one workflow below before running extrapolation:"),
              radioButtons(
                "points_occurrence_mode",
                "Occurrence workflow:",
                choices = list(
                  "Regular grid only (points -> grid cells)" = "regular_grid",
                  "Irregular polygons only (direct points -> polygons, no regular grid)" = "irregular_direct",
                  "Irregular polygons + regular grid (points -> grid -> polygons)" = "irregular_with_grid"
                ),
                selected = "regular_grid"
              ),
              conditionalPanel(
                condition = "input.points_occurrence_mode == 'irregular_direct' || input.points_occurrence_mode == 'irregular_with_grid'",
                p(strong("Shapefile 2 (optional): irregular polygons for counting/reporting units")),
                p("Examples: states, ecoregions, municipalities. These polygons define output columns/areas for this workflow."),
                fileInput("points_irregular_bins_shapefile", "Choose ALL irregular polygon files together:", accept = c(".shp", ".dbf", ".shx", ".prj"), multiple = TRUE),
                selectInput("points_irregular_bin_id_column", "Irregular polygon ID column:", choices = c("(Auto-detect)" = ""), selected = ""),
                p("Available columns are populated automatically after upload.", style = "font-size: 12px; color: #666;")
              ),
              conditionalPanel(
                condition = "input.points_occurrence_mode == 'regular_grid'",
                p("Regular grid workflow: occurrence points are assigned to grid cells based on the selected grid resolution.", style = "font-size: 12px; color: #666;")
              ),
              conditionalPanel(
                condition = "input.points_occurrence_mode == 'irregular_direct'",
                p("Direct irregular workflow: richness is computed by counting occurrence points directly inside each polygon.", style = "font-size: 12px; color: #666;"),
                p("In this mode, diversity by irregular polygons is computed automatically and the exported matrix basis becomes irregular polygons (direct point-in-polygon).", style = "font-size: 12px; color: #666;")
              ),
              conditionalPanel(
                condition = "input.points_occurrence_mode == 'irregular_with_grid'",
                p("Grid-overlap workflow: points are first assigned to regular grid cells; then species counts per polygon are computed from grid-polygon overlap.", style = "font-size: 12px; color: #666;"),
                p("In this mode, diversity by irregular polygons is computed automatically and the exported matrix basis becomes irregular polygons (derived from grid cells).", style = "font-size: 12px; color: #666;")
              )
            ),
            conditionalPanel(
              condition = "input.extrap_method == 'buffer' || input.extrap_method == 'convex_hull' || input.extrap_method == 'mst'",
              tags$hr(),
              h5("Optional: Diversity in Irregular Polygons"),
              p(strong("Shapefile 2 (optional): subdivisions for diversity and optional matrix aggregation outputs.")),
              p("Workflow: choose Buffer/Convex Hull/MST + enable this option + upload a SECOND shapefile with subdivisions + choose ID column + click Run Extrapolation once.", style = "font-size: 12px; color: #666;"),
              checkboxInput("enable_irregular_richness", "Compute diversity by irregular polygons", value = FALSE),
              conditionalPanel(
                condition = "input.enable_irregular_richness == true",
                p("Upload all subdivision shapefile files (.shp, .shx, .dbf, .prj):"),
                fileInput("irregular_richness_shapefile", "Choose ALL subdivision shapefile files together:", accept = c(".shp", ".dbf", ".shx", ".prj"), multiple = TRUE),
                selectInput("irregular_richness_id_column", "Subdivision ID column:", choices = c("(Auto-detect)" = ""), selected = ""),
                p("When this option is enabled, the output matrix is aggregated to these subdivisions (recommended for BioGeoBEARS).", style = "font-size: 12px; color: #666;"),
                p("If left empty, the app will auto-detect a suitable ID column.", style = "font-size: 12px; color: #666;")
              )
            ),
            conditionalPanel(
              condition = "input.extrap_method == 'mst'",
              h5("MST Analysis Level"),
               radioButtons("mst_level", "Select level of analysis:",
                            choices = list(
                             "All Taxa (Individual MSTs)" = "all",
                             "Single MST for All Taxa" = "single_all",
                             "Specific Node(s) (Ancestral Area)" = "node",
                             "Specific Node(s) + Extra Taxa" = "node_taxa",
                             "Specific Taxa" = "taxa"
                            ),
                            selected = "all"),
               
                conditionalPanel(
                 condition = "input.mst_level == 'node' || input.mst_level == 'node_taxa'",
                  p("View the tree below to identify node numbers."),
                  p("Node mode builds an ancestral MST focused on selected internal node(s), instead of separate MST tracks for all taxa.", style = "font-size: 12px; color: #666;"),
                  p("Node numbers are shown as high-contrast labels (white text on red circles).", style = "font-size: 12px; color: #666;"),
                  selectizeInput("mst_node_ids", "Select Internal Node ID(s):", choices = NULL, multiple = TRUE),
                  plotOutput("mst_tree_plot", height = "520px")
                ),
               
               conditionalPanel(
                 condition = "input.mst_level == 'taxa' || input.mst_level == 'node_taxa'",
                 selectizeInput("mst_taxa_select", "Select Taxa:", choices = NULL, multiple = TRUE),
                 p("In 'node + extra taxa' mode, selected taxa are added to the taxa descending from selected node(s).", style = "font-size: 12px; color: #666;")
               )
             ),
            conditionalPanel(
              condition = "input.extrap_method != 'occurrence_only' || input.points_occurrence_mode == 'regular_grid' || input.points_occurrence_mode == 'irregular_with_grid'",
              numericInput("grid_resolution", "Grid resolution (degrees):", value = 1, min = 0.1, step = 0.1)
            ),
            conditionalPanel(
              condition = "input.extrap_method == 'occurrence_only'",
              p("If you only want to visualize occurrence points in Step 4, use the button below (no extrapolation, no matrix rebuild).", style = "font-size: 12px; color: #666;"),
              actionButton("plot_points_only", "Plot Points Only (No Extrapolation)", class = "btn btn-info"),
              tags$span(" ")
            ),
            actionButton("run_extrap", "Run Extrapolation", class = "btn btn-primary")
          ),
          div(
            class = "col-md-4",
            h4("Status"),
            verbatimTextOutput("extrap_status"),
            h5("Analysis Log"),
            verbatimTextOutput("extrap_console_log"),
            downloadButton("download_extrap_log", "Download Log (.txt)")
          )
        )
      )
    ),
    
    # Tab 5: PAE-PCE Analysis
    tabPanel(
      "PAE-PCE (Optional)",
      icon = icon("project-diagram"),
      div(
        class = "container-fluid",
        style = "padding: 20px;",
        div(
          class = "step-title",
          "Optional Analysis: PAE-PCE"
        ),
        div(
          class = "guidance-text",
          p("Run Parsimony Analysis of Endemicity with Parsimony Constraint (PAE-PCE) only if you want a post-MST generalized track analysis.")
        ),
        
        div(
          class = "row",
          div(
            class = "col-md-4",
            h4("Analysis Parameters"),
            p("Optional module. Requires a presence-absence matrix generated by the Minimum Spanning Tree (MST) method.", style = "color: #d9534f; font-weight: bold;"),
            
            # Option to use custom data
            div(style = "background-color: #f5f5f5; padding: 10px; border-radius: 4px; margin-bottom: 15px;",
              checkboxInput("pae_use_custom_data", "Use your own data?", value = FALSE),
              conditionalPanel(
                condition = "input.pae_use_custom_data == true",
                fileInput("pae_custom_matrix", "Upload presence-absence matrix (CSV/TXT):", accept = c(".csv", ".txt")),
                fileInput("pae_custom_shapefile", "Upload shapefile (.shp + .dbf + .shx):", multiple = TRUE, accept = c(".shp", ".dbf", ".shx")),
                numericInput("pae_custom_resol_x", "Grid resolution X:", value = 1, min = 0.1),
                numericInput("pae_custom_resol_y", "Grid resolution Y:", value = 1, min = 0.1),
                div(style = "background-color: #e8f4f8; padding: 8px; border-left: 3px solid #2c5aa0; margin-top: 8px; font-size: 11px;",
                  strong("Matrix format:"), br(),
                  "First column: grid IDs (numeric or text)", br(),
                  "Other columns: species names (0/1 for absence/presence)", br(),
                  "Example: grid_1,sp1,sp2,sp3", br(),
                  "         1,1,0,1", br(),
                  "         2,1,1,0"
                )
              )
            ),
            
             numericInput("pae_n_iterations", "Number of iterations (N):", value = 10, min = 1, max = 100),
             numericInput("pae_seed", "Random seed (reproducibility):", value = 123, min = 1),
             
              checkboxInput("pae_grid_view", "View grid on map", value = TRUE),
              checkboxInput("pae_label_grid", "Show grid numbers", value = TRUE),
              checkboxInput("pae_sobrepo", "Superimpose on existing MSTs", value = FALSE),
              p("PAE-PCE uses the MST presence-absence matrix (pres_abs.txt), with the same study-area shapefile and grid resolution from Step 3.", style = "font-size: 12px; color: #666;"),
            
            h5("Zoom / Extent (Optional)"),
            p("Leave empty to use full extent", style = "font-size: 12px; color: #666;"),
            fluidRow(
              column(6, numericInput("pae_xmin", "X Min:", value = NA)),
              column(6, numericInput("pae_xmax", "X Max:", value = NA))
            ),
            fluidRow(
              column(6, numericInput("pae_ymin", "Y Min:", value = NA)),
              column(6, numericInput("pae_ymax", "Y Max:", value = NA))
            ),
            
            actionButton("run_pae_pce", "Run PAE-PCE Analysis", class = "btn btn-primary", style = "width: 100%; margin-top: 15px;")
          ),
          div(
            class = "col-md-8",
            h4("Results"),
            tabsetPanel(
              tabPanel(
                "Generalized Track Map",
                br(),
                plotOutput("pae_pce_plot", height = "600px")
              ),
              tabPanel(
                "Species List",
                br(),
                DT::dataTableOutput("pae_pce_table"),
                br(),
                downloadButton("download_pae_table", "Download Table (.csv)", class = "btn btn-success")
              ),
              tabPanel(
                "Analysis Log",
                br(),
                downloadButton("download_pae_input_matrix", "Download exact PAE input matrix", class = "btn btn-info"),
                br(), br(),
                verbatimTextOutput("pae_input_matrix_info"),
                br(),
                verbatimTextOutput("pae_pce_log")
              )
            )
          )
        )
      )
    ),
    
    # Tab 6: Visualizations
    tabPanel(
      "Visualizations",
      icon = icon("image"),
      div(
        class = "container-fluid",
        style = "padding: 20px;",
        div(
          class = "step-title",
          "Step 4: Visualizations"
        ),
        div(
          class = "guidance-text",
          p("View maps of species distributions and phylogenetic trees with ancestral ranges."),
          p("After running BioGeoBEARS (Step 7), choose a fitted model below to visualize ancestral states.", style = "font-size: 12px; color: #666;")
        ),
        selectInput("bgb_visual_model", "BioGeoBEARS model for visualization:", choices = c("(Run BioGeoBEARS first)" = ""), selected = ""),
        sliderInput(
          "ancestral_ranges_label_offset",
          "Ancestral ranges label offset:",
          min = 0,
          max = 1,
          value = 0.35,
          step = 0.05
        ),
        sliderInput(
          "ancestral_ranges_tip_cex_app",
          "Ancestral ranges (app) - terminal label size:",
          min = 0.3,
          max = 2,
          value = 0.95,
          step = 0.05
        ),
        sliderInput(
          "ancestral_ranges_state_cex_app",
          "Ancestral ranges (app) - node range-box size:",
          min = 0.3,
          max = 2,
          value = 1.0,
          step = 0.05
        ),
        sliderInput(
          "range_uncertainty_label_offset",
          "Range uncertainty label offset:",
          min = 0,
          max = 1,
          value = 0.2,
          step = 0.05
        ),
        sliderInput(
          "range_uncertainty_tip_cex_app",
          "Range uncertainty (app) - terminal label size:",
          min = 0.3,
          max = 2,
          value = 0.95,
          step = 0.05
        ),
        sliderInput(
          "range_uncertainty_state_cex_app",
          "Range uncertainty (app) - pie size:",
          min = 0.3,
          max = 2,
          value = 1.0,
          step = 0.05
        ),
        downloadButton("download_bgb_ancestral_pdf", "Download ancestral ranges + uncertainty (PDF)", class = "btn btn-success btn-sm"),
        br(), br(),
        
        tabsetPanel(
            tabPanel(
              "Distribution Map",
            br(),
            div(
              class = "row",
              style = "margin-bottom: 15px;",
              div(
                class = "col-md-12",
                h4("Map Controls"),
                sliderInput(
                  "polygon_opacity",
                  "Polygon Opacity:",
                  min = 0,
                  max = 1,
                  value = 0.5,
                  step = 0.1
                ),
                checkboxInput("viz_overlay_all_methods", "Overlay layers from all extrapolation methods", value = FALSE),
                actionButton("load_polygons", "Reload output layers", class = "btn btn-default btn-sm"),
                actionButton("clear_visual_outputs", "Clear Output Layers", class = "btn btn-warning btn-sm"),
                br(),
                p("Colors are automatically assigned by species/layer and applied immediately after extrapolation.", style = "font-size: 12px; color: #666;"),
                p("Grid/raster files are still generated for outputs; map controls focus on shapefile/vector layers.", style = "font-size: 12px; color: #666;")
              )
            ),
            uiOutput("distribution_map_panel")
            ),

            tabPanel(
              "Irregular Polygon Diversity",
              br(),
              p("How to use this panel: in Step 3, select Buffer/Convex Hull/MST, enable 'Compute diversity by irregular polygons', upload the SECOND subdivision shapefile, choose the subdivision ID column, then click Run Extrapolation. For points-only, select 'Occurrence points only', enable irregular polygons, upload the SECOND shapefile, choose ID column, and run once."),
              leaflet::leafletOutput("irregular_bins_map", height = "500px"),
              br(),
              DT::dataTableOutput("irregular_bins_table")
            ),
           
            tabPanel(
              "Phylogenetic Tree",
            br(),
            h4("Your Phylogenetic Tree"),
            plotOutput("phylo_tree_plot", height = "600px")
          ),
          
          tabPanel(
            "Ancestral Ranges",
            br(),
            h4("Ancestral Ranges"),
            p("Most likely ancestral ranges at each node (text labels)."),
            plotOutput("ancestral_ranges_plot", height = "650px"),
            br(),
            verbatimTextOutput("ancestral_ranges_status")
          ),
          
          tabPanel(
            "Range Uncertainty",
            br(),
            h4("Range Uncertainty"),
            p("Node-level uncertainty as pie charts of ancestral-range probabilities."),
            plotOutput("range_uncertainty_plot", height = "650px"),
            br(),
            verbatimTextOutput("range_uncertainty_status")
          )
        )
      )
    ),
    
    # Tab 6: Export Files
    tabPanel(
      "Export Files",
      icon = icon("download"),
      div(
        class = "container-fluid",
        style = "padding: 20px;",
        div(
          class = "step-title",
          "Step 5: Export Results"
        ),
        div(
          class = "guidance-text",
          p("Export your results in multiple formats."),
          p("For compatibility, area names are abbreviated (A, B, C, ...). Download the area-code mapping after export.", style = "font-size: 12px; color: #666;")
        ),

        downloadButton("download_area_code_mapping", "Download area-code mapping", class = "btn btn-info btn-sm"),
        br(), br(),
        verbatimTextOutput("area_code_mapping_status"),
        radioButtons(
          "matrix_export_basis",
          "Matrix basis for exports:",
          choices = list(
            "Current matrix from last extrapolation" = "current",
            "Regular-grid matrix (when available)" = "regular",
            "Irregular-polygon matrix (when available)" = "irregular"
          ),
          selected = "current"
        ),
        verbatimTextOutput("matrix_export_status"),
        
        tabsetPanel(
          tabPanel(
            "BioGeoBEARS",
            br(),
            p("Export presence-absence matrix in BioGeoBEARS format (.data)."),
            downloadButton("download_biogeobears", "Download BioGeoBEARS File", class = "btn btn-success"),
            br(), br(),
            verbatimTextOutput("biogeobears_export_status")
          ),
          
          tabPanel(
            "TNT Matrix",
            br(),
            p("Export presence-absence matrix in TNT format."),
            downloadButton("download_tnt_matrix", "Download TNT Matrix", class = "btn btn-success"),
            br(), br(),
            verbatimTextOutput("tnt_export_status")
          ),
          
          tabPanel(
            "NEXUS Export",
            br(),
            p("Export presence-absence matrix in NEXUS format."),
            selectInput(
              "nexus_method",
              "Select extrapolation method:",
              choices = list(
                "Occurrence Points (no extrapolation)" = "PTS",
                "Buffer" = "BUFF",
                "Convex Hull" = "MPC",
                "Minimum Spanning Tree" = "MST",
                "Irregular Bins" = "IRR_BINS"
              )
            ),
            downloadButton("download_nexus", "Download NEXUS File", class = "btn btn-success"),
            br(), br(),
            verbatimTextOutput("nexus_export_status")
          ),
          
          tabPanel(
            "TNT PAE-PCE Script",
            br(),
            p("Generate a script for TNT PAE-PCE analysis."),
            downloadButton("download_tnt_script", "Download TNT Script", class = "btn btn-success"),
            br(), br(),
            verbatimTextOutput("tnt_script_status")
          ),
          
          tabPanel(
            "NDM Export",
            br(),
            p("Export data in NDM/VNDM format."),
            downloadButton("download_ndm", "Download NDM File", class = "btn btn-success"),
            br(), br(),
            verbatimTextOutput("ndm_export_status")
          ),
          
          tabPanel(
            "Shapefiles",
            br(),
            p("Download generated shapefiles from extrapolation methods."),
            fluidRow(
              column(6,
                selectInput(
                  "shapefile_dir",
                  "Select Output Directory:",
                  choices = list(
                    "Buffer (out_buffers)" = "out_buffers",
                    "Convex Hull (out_MCP)" = "out_MCP",
                    "Minimum Spanning Tree (out_MST)" = "out_MST",
                    "Regular Grid (out_grid)" = "out_grid",
                    "Irregular Bins (out_irregular_bins)" = "out_irregular_bins"
                  )
                ),
                actionButton("refresh_shapefiles", "Refresh File List", icon = icon("sync"))
              ),
              column(6,
                uiOutput("shapefile_download_ui")
              )
            )
          )
        )
      )
    ),
    
    # Tab 7: BioGeoBEARS Setup
    tabPanel(
      "BioGeoBEARS Setup",
      icon = icon("cog"),
      div(
        class = "container-fluid",
        style = "padding: 20px;",
        div(
          class = "step-title",
          "Step 6: BioGeoBEARS Configuration"
        ),
        div(
          class = "guidance-text",
          p("Configure your BioGeoBEARS analysis."),
          p("You can use the matrix generated in this app (Step 3/5) or upload an external BioGeoBEARS .data file.", style = "font-size: 12px; color: #666;")
        ),

        div(
          class = "row",
          div(
            class = "col-md-6",
            h4("Geographic Matrix Source"),
            radioButtons(
              "bgb_geog_source",
              "Geographic data source:",
              choices = list(
                "Use matrix generated in this app" = "current",
                "Upload external BioGeoBEARS .data file" = "upload"
              ),
              selected = "current"
            ),
            conditionalPanel(
              condition = "input.bgb_geog_source == 'upload'",
              fileInput("bgb_geog_file", "Upload .data file:", accept = c(".data", ".txt"))
            )
          ),
          div(
            class = "col-md-6",
            h4("Tree Source"),
            radioButtons(
              "bgb_tree_source",
              "Tree source:",
              choices = list(
                "Use tree loaded in Step 1" = "step1",
                "Upload external Newick tree" = "upload"
              ),
              selected = "step1"
            ),
            conditionalPanel(
              condition = "input.bgb_tree_source == 'upload'",
              fileInput("bgb_tree_file", "Upload Newick tree:", accept = c(".tre", ".tree", ".newick", ".nwk", ".txt"))
            )
          )
        ),

        div(
          class = "row",
          div(
            class = "col-md-3",
            numericInput("bgb_max_range_size", "Max range size:", value = 3, min = 1, step = 1)
          ),
          div(
            class = "col-md-3",
            numericInput("bgb_min_branchlength", "Min branch length:", value = 0.000001, min = 0.000001, step = 0.000001)
          ),
          div(
            class = "col-md-3",
            numericInput("bgb_num_cores", "Number of cores:", value = 1, min = 1, step = 1)
          ),
          div(
            class = "col-md-3",
            selectInput("bgb_optimizer", "Optimizer:", choices = c("GenSA", "optimx"), selected = "GenSA")
          )
        ),
        checkboxInput("bgb_include_null_range", "Include null range", value = FALSE),
        verbatimTextOutput("bgb_input_status"),
        
        div(
          class = "row",
          div(
            class = "col-md-6",
            h4("Model Selection"),
            checkboxGroupInput(
              "bgb_models",
              "Select models to run:",
              choices = list(
                "DEC" = "DEC",
                "DEC+J" = "DECJ",
                "DIVALIKE" = "DIVALIKE",
                "DIVALIKE+J" = "DIVALIKEJ",
                "BAYAREALIKE" = "BAYAREALIKE",
                "BAYAREALIKE+J" = "BAYAREALIKEJ"
              ),
              selected = c("DEC", "DECJ")
            )
          ),
          div(
            class = "col-md-6",
            h4("Advanced Options"),
            checkboxInput("bgb_use_nested_starts", "Warm-start +J from nested model (d/e)", value = TRUE),
            checkboxInput("bgb_use_nested_starts_xwn", "Also warm-start x/w/n from nested model", value = FALSE),
            checkboxInput("use_distance_matrix", "Use distance matrix (x parameter)", value = FALSE),
            checkboxInput("use_dispersal_multiplier", "Use dispersal multiplier", value = FALSE),
            checkboxInput("use_time_slices", "Use time slices", value = FALSE)
          )
        ),
        
        hr(),
        
        h4("Upload Custom Matrices"),
        div(
          class = "guidance-text",
          p("Upload custom matrices for BioGeoBEARS analysis. Format: CSV or TXT with row and column names.")
        ),
        
        # Row 1: Distance and Dispersal Multipliers
        div(
          class = "row",
          div(
            class = "col-md-6",
            h5("Distance Matrix"),
            p("Pairwise distances between areas (for x parameter). Can include multiple time periods.", 
              style = "font-size: 12px; color: #666;"),
            fileInput(
              "distance_matrix_file",
              "Upload Distance Matrix (CSV/TXT):",
              accept = c(".csv", ".txt")
            ),
            actionButton("load_distance_matrix", "Load Matrix", class = "btn btn-sm btn-info"),
            br(),
            verbatimTextOutput("distance_matrix_status")
          ),
          div(
            class = "col-md-6",
            h5("Dispersal Multipliers"),
            p("Matrix of 0s and 1s restricting dispersal routes between areas (1=allowed, 0=blocked).", 
              style = "font-size: 12px; color: #666;"),
            fileInput(
              "dispersal_multipliers_file",
              "Upload Dispersal Multipliers (CSV/TXT):",
              accept = c(".csv", ".txt")
            ),
            actionButton("load_dispersal_multipliers", "Load Matrix", class = "btn btn-sm btn-info"),
            br(),
            verbatimTextOutput("dispersal_multipliers_status")
          )
        ),
        
        # Row 2: Environmental Distance and Areas-Allowed
        div(
          class = "row",
          div(
            class = "col-md-6",
            h5("Environmental Distance Matrix"),
            p("Environmental distances between areas (for n parameter). Same format as distance matrix.", 
              style = "font-size: 12px; color: #666;"),
            fileInput(
              "env_distance_matrix_file",
              "Upload Environmental Distance Matrix (CSV/TXT):",
              accept = c(".csv", ".txt")
            ),
            actionButton("load_env_distance_matrix", "Load Matrix", class = "btn btn-sm btn-info"),
            br(),
            verbatimTextOutput("env_distance_matrix_status")
          ),
          div(
            class = "col-md-6",
            h5("Areas-Allowed Matrix"),
            p("Defines which areas are available at each time period (1=available, 0=not available).", 
              style = "font-size: 12px; color: #666;"),
            fileInput(
              "areas_allowed_file",
              "Upload Areas-Allowed Matrix (CSV/TXT):",
              accept = c(".csv", ".txt")
            ),
            actionButton("load_areas_allowed", "Load Matrix", class = "btn btn-sm btn-info"),
            br(),
            verbatimTextOutput("areas_allowed_status")
          )
        ),
        
        # Row 3: Areas-Adjacency and Time Periods
        div(
          class = "row",
          div(
            class = "col-md-6",
            h5("Areas-Adjacency Matrix"),
            p("Defines which areas can coexist in species ranges (1=adjacent, 0=not adjacent).", 
              style = "font-size: 12px; color: #666;"),
            fileInput(
              "areas_adjacency_file",
              "Upload Areas-Adjacency Matrix (CSV/TXT):",
              accept = c(".csv", ".txt")
            ),
            actionButton("load_areas_adjacency", "Load Matrix", class = "btn btn-sm btn-info"),
            br(),
            verbatimTextOutput("areas_adjacency_status")
          ),
          div(
            class = "col-md-6",
            h5("Time Periods"),
            p("Vector of time periods in millions of years (Ma). One value per line.", 
              style = "font-size: 12px; color: #666;"),
            fileInput(
              "time_periods_file",
              "Upload Time Periods (TXT):",
              accept = c(".txt")
            ),
            actionButton("load_time_periods", "Load Periods", class = "btn btn-sm btn-info"),
            br(),
            verbatimTextOutput("time_periods_status")
          )
        ),
        
        # Summary
        div(
          class = "row",
          div(
            class = "col-md-12",
            h5("Loaded Matrices Summary"),
            verbatimTextOutput("matrices_summary"),
            br(),
            h5("Time-Stratified Validator"),
            verbatimTextOutput("bgb_strat_validation")
          )
        )
      )
    ),
    
    # Tab 8: Analysis & Results
    tabPanel(
      "Analysis & Results",
      icon = icon("bar-chart"),
      div(
        class = "container-fluid",
        style = "padding: 20px;",
        div(
          class = "step-title",
          "Step 7: Analysis & Results"
        ),
        div(
          class = "guidance-text",
          p("Run BioGeoBEARS analysis and view results.")
        ),
        
        actionButton("run_analysis", "Run BioGeoBEARS Analysis", class = "btn btn-primary btn-lg"),
        br(), br(),
        verbatimTextOutput("analysis_status"),
        br(),
        
        tabsetPanel(
          tabPanel(
            "Model Comparison",
            br(),
            h4("Model Comparison"),
            DT::dataTableOutput("model_comparison_table"),
            br(),
            downloadButton("download_model_comparison", "Download Model Comparison (.csv)", class = "btn btn-success"),
            br(), br(),
            downloadButton("download_bgb_models_rdata", "Download Fitted Models (.RData)", class = "btn btn-primary"),
            br(), br(),
            downloadButton("download_bgb_run_report", "Download Run Config Report (.txt)", class = "btn btn-info"),
            br(), br(),
            p("Tip: use this file as the model-selection report for reproducibility.", style = "font-size: 12px; color: #666;")
          ),
          
          tabPanel(
            "Likelihood Ratio Test",
            br(),
            h4("Likelihood Ratio Test (LRT)"),
            div(
              class = "row",
              div(
                class = "col-md-6",
                selectInput("lrt_model1", "Model 1:", choices = c("DEC", "DECJ")),
                selectInput("lrt_model2", "Model 2:", choices = c("DEC", "DECJ"))
              ),
              div(
                class = "col-md-6",
                actionButton("run_lrt", "Run LRT", class = "btn btn-primary")
              )
            ),
            br(),
            verbatimTextOutput("lrt_results")
          )
        )
      )
    ),
    
    # Tab 9: About
    tabPanel(
      "About",
      icon = icon("info-circle"),
      div(
        class = "container-fluid",
        style = "padding: 40px;",
        h2("About BRA (BioRangeAnalyzer Shiny)"),
        hr(),
        h3("Project"),
        p("BRA project page:"),
        p(tags$a(href = "https://jozecaricardo.github.io/biogeografiaAULAS/", "https://jozecaricardo.github.io/biogeografiaAULAS/", target = "_blank")),
        hr(),
        h3("Authors"),
        tags$ul(
          tags$li("José Ricardo Inacio Ribeiro (Universidade Federal do Pampa (UNIPAMPA), campus São Gabriel, Rio Grande do Sul State, Brazil)"),
          tags$li("Augusto Ferrari (Universidade Federal do Rio Grande (FURG), Rio Grande State, Brazil)")
        ),
        hr(),
        h3("Features"),
        tags$ul(
          tags$li("Range extrapolation (buffers, convex hulls, MST)"),
          tags$li("Phylogenetic tree plotting"),
          tags$li("Comprehensive model comparison"),
          tags$li("Statistical tests (AIC, AICc, LRT)"),
          tags$li("Export to multiple formats")
        ),
        hr(),
        h3("Citation"),
        p("If you use this application, please cite:"),
        code("BRA (BioRangeAnalyzer Shiny) (2025)"),
        hr(),
        h3("References"),
        tags$ul(
          tags$li("Matzke, N.J. (2013). BioGeoBEARS: BioGeography with Bayesian and Likelihood Evolutionary Analysis in R Scripts.")
        )
      )
    )
  )
)
