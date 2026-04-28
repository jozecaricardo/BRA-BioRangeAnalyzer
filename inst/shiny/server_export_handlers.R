# ===== EXPORT FILES TAB HANDLERS =====

# TNT Matrix Export
observeEvent(input$export_tnt_matrix, {
  if (is.null(data_store$pres_abs)) {
    output$tnt_matrix_status <- renderPrint({
      cat("Error: Please run range extrapolation first\n")
    })
  } else {
    tryCatch({
      # Create output directory
      if (!dir.exists("out")) {
        dir.create("out", recursive = TRUE)
      }
      
      # Call the tnt_matrix function
      tnt_matrix(data_store$pres_abs)
      
      # Store the file path for download
      data_store$tnt_matrix_file <- "out/pres_abs.tnt"
      
      output$tnt_matrix_status <- renderPrint({
        cat("✓ TNT matrix exported successfully!\n")
        cat("File: out/pres_abs.tnt\n")
        cat("Dimensions:", nrow(data_store$pres_abs), "taxa x", ncol(data_store$pres_abs), "species\n")
      })
    }, error = function(e) {
      output$tnt_matrix_status <- renderPrint({
        cat("Error exporting TNT matrix:\n")
        cat(e$message, "\n")
      })
    })
  }
})

# TNT PAE-PCE Script Generation
observeEvent(input$export_tnt_script, {
  if (is.null(data_store$tnt_matrix_file)) {
    output$tnt_script_status <- renderPrint({
      cat("Error: Please export TNT matrix first\n")
    })
  } else {
    tryCatch({
      # Create output directory
      if (!dir.exists("out_TNT")) {
        dir.create("out_TNT", recursive = TRUE)
      }
      
      # Call the generate_tnt_pae_pce function
      script_file <- generate_tnt_pae_pce(
        matrix_file = data_store$tnt_matrix_file,
        output_file = "PAE_PCE_tntAUTO.run",
        max_iterations = input$tnt_max_iterations,
        search_replicates = input$tnt_search_replicates,
        search_method = input$tnt_search_method,
        output_dir = "out_TNT/"
      )
      
      # Store the file path for download
      data_store$tnt_script_file <- script_file
      
      output$tnt_script_status <- renderPrint({
        cat("✓ TNT PAE-PCE script generated successfully!\n")
        cat("File:", script_file, "\n")
        cat("Parameters:\n")
        cat("  - Max iterations:", input$tnt_max_iterations, "\n")
        cat("  - Search replicates:", input$tnt_search_replicates, "\n")
        cat("  - Search method:", input$tnt_search_method, "\n")
      })
    }, error = function(e) {
      output$tnt_script_status <- renderPrint({
        cat("Error generating TNT script:\n")
        cat(e$message, "\n")
      })
    })
  }
})

# NDM Export
observeEvent(input$export_ndm, {
  if (is.null(data_store$occurrence)) {
    output$ndm_status <- renderPrint({
      cat("Error: Please load occurrence data first\n")
    })
  } else {
    tryCatch({
      # Create output directory
      if (!dir.exists("out_NDM")) {
        dir.create("out_NDM", recursive = TRUE)
      }
      
      # Prepare input based on selected mode
      if (input$ndm_input_type == "occurrence") {
        # Use occurrence data directly
        input_data <- data_store$occurrence
      } else {
        # Use shapefile directory
        input_data <- input$ndm_shapefile_dir
      }
      
      # Call the toNDM function
      if (input$ndm_input_type == "shapefile") {
        toNDM(
          input_data = input_data,
          shapefile_dir = input$ndm_shapefile_dir,
          resolution = c(input$ndm_resolution_x, input$ndm_resolution_y),
          output_file = input$ndm_output_file,
          output_dir = "out_NDM/"
        )
      } else {
        toNDM(
          input_data = input_data,
          output_file = input$ndm_output_file,
          output_dir = "out_NDM/"
        )
      }
      
      # Store the file path for download
      data_store$ndm_file <- file.path("out_NDM", input$ndm_output_file)
      
      output$ndm_status <- renderPrint({
        cat("✓ NDM file exported successfully!\n")
        cat("File:", data_store$ndm_file, "\n")
        if (input$ndm_input_type == "shapefile") {
          cat("Grid resolution:", input$ndm_resolution_x, "x", input$ndm_resolution_y, "degrees\n")
        }
      })
    }, error = function(e) {
      output$ndm_status <- renderPrint({
        cat("Error exporting NDM file:\n")
        cat(e$message, "\n")
      })
    })
  }
})

# Display available files for download
output$available_files <- renderUI({
  files <- c()
  
  if (!is.null(data_store$tnt_matrix_file) && file.exists(data_store$tnt_matrix_file)) {
    files <- c(files, data_store$tnt_matrix_file)
  }
  if (!is.null(data_store$tnt_script_file) && file.exists(data_store$tnt_script_file)) {
    files <- c(files, data_store$tnt_script_file)
  }
  if (!is.null(data_store$ndm_file) && file.exists(data_store$ndm_file)) {
    files <- c(files, data_store$ndm_file)
  }
  
  if (length(files) == 0) {
    return(p("No files exported yet. Export files above to see them here."))
  }
  
  # Create download buttons for each file
  lapply(files, function(file) {
    div(
      downloadButton(
        paste0("download_", gsub("[^[:alnum:]]", "_", file)),
        label = basename(file),
        class = "btn btn-sm btn-info"
      ),
      span(paste0(" (", file, ")"), style = "margin-left: 10px;")
    )
  })
})

# Download individual files
observe({
  files <- c()
  
  if (!is.null(data_store$tnt_matrix_file) && file.exists(data_store$tnt_matrix_file)) {
    files <- c(files, data_store$tnt_matrix_file)
  }
  if (!is.null(data_store$tnt_script_file) && file.exists(data_store$tnt_script_file)) {
    files <- c(files, data_store$tnt_script_file)
  }
  if (!is.null(data_store$ndm_file) && file.exists(data_store$ndm_file)) {
    files <- c(files, data_store$ndm_file)
  }
  
  for (file in files) {
    local({
      file_path <- file
      button_id <- paste0("download_", gsub("[^[:alnum:]]", "_", file_path))
      
      output[[button_id]] <- downloadHandler(
        filename = function() {
          basename(file_path)
        },
        content = function(con) {
          file.copy(file_path, con)
        }
      )
    })
  }
})

# Download all exports as ZIP
output$download_all_exports <- downloadHandler(
  filename = function() {
    paste0("biogeoshiny_exports_", Sys.Date(), ".zip")
  },
  content = function(con) {
    # Collect all exported files
    files_to_zip <- c()
    
    if (!is.null(data_store$tnt_matrix_file) && file.exists(data_store$tnt_matrix_file)) {
      files_to_zip <- c(files_to_zip, data_store$tnt_matrix_file)
    }
    if (!is.null(data_store$tnt_script_file) && file.exists(data_store$tnt_script_file)) {
      files_to_zip <- c(files_to_zip, data_store$tnt_script_file)
    }
    if (!is.null(data_store$ndm_file) && file.exists(data_store$ndm_file)) {
      files_to_zip <- c(files_to_zip, data_store$ndm_file)
    }
    
    if (length(files_to_zip) > 0) {
      # Create a temporary directory
      tmpdir <- tempdir()
      
      # Copy files to temporary directory
      for (file in files_to_zip) {
        file.copy(file, file.path(tmpdir, basename(file)), overwrite = TRUE)
      }
      
      # Create ZIP file
      zip(con, files = file.path(tmpdir, basename(files_to_zip)))
    }
  }
)
