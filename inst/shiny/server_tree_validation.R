# ===== TREE VALIDATION TAB HANDLERS =====

# Validate tree
observeEvent(input$validate_tree, {
  if (is.null(data_store$tree)) {
    output$tree_validation_info <- renderPrint({
      cat("ERROR: No tree loaded!\n")
      cat("Please load a phylogenetic tree in the Data Input tab first.\n")
    })
  } else {
    tryCatch({
      tree <- data_store$tree
      
      # Perform validation checks
      n_taxa <- length(tree$tip.label)
      has_branch_lengths <- !is.null(tree$edge.length)
      is_ultrametric <- ape::is.ultrametric(tree)
      
      # Calculate tree statistics
      if (has_branch_lengths) {
        max_branch_length <- max(tree$edge.length, na.rm = TRUE)
        min_branch_length <- min(tree$edge.length, na.rm = TRUE)
        mean_branch_length <- mean(tree$edge.length, na.rm = TRUE)
      }
      
      # Calculate root to tip distances
      root_to_tip <- ape::node.depth.edgelength(tree)[1:n_taxa]
      max_root_to_tip <- max(root_to_tip, na.rm = TRUE)
      min_root_to_tip <- min(root_to_tip, na.rm = TRUE)
      
      output$tree_validation_info <- renderPrint({
        cat("=== TREE VALIDATION REPORT ===\n\n")
        cat("✓ Tree loaded successfully!\n\n")
        
        cat("BASIC INFORMATION:\n")
        cat("  Number of taxa:", n_taxa, "\n")
        cat("  Number of internal nodes:", tree$Nnode, "\n")
        cat("  Total branches:", nrow(tree$edge), "\n\n")
        
        cat("BRANCH LENGTHS:\n")
        if (has_branch_lengths) {
          cat("  ✓ Branch lengths present\n")
          cat("  Max branch length:", round(max_branch_length, 4), "\n")
          cat("  Min branch length:", round(min_branch_length, 4), "\n")
          cat("  Mean branch length:", round(mean_branch_length, 4), "\n\n")
        } else {
          cat("  ✗ No branch lengths\n\n")
        }
        
        cat("ULTRAMETRIC STATUS:\n")
        if (is_ultrametric) {
          cat("  ✓ Tree is ULTRAMETRIC\n")
          cat("  All tips are equidistant from the root\n")
          cat("  Root to tip distance:", round(max_root_to_tip, 4), "\n\n")
        } else {
          if (has_branch_lengths) {
            cat("  ✗ Tree is NOT ultrametric\n")
            cat("  Max root to tip distance:", round(max_root_to_tip, 4), "\n")
            cat("  Min root to tip distance:", round(min_root_to_tip, 4), "\n")
            cat("  Difference:", round(max_root_to_tip - min_root_to_tip, 4), "\n")
            cat("  → You can convert it using the controls on the right\n\n")
          } else {
            cat("  ✗ Cannot determine ultrametric status (no branch lengths)\n\n")
          }
        }
        
        cat("TAXA NAMES:\n")
        if (n_taxa <= 20) {
          cat("  ", paste(tree$tip.label, collapse = ", "), "\n\n")
        } else {
          cat("  First 10 taxa:", paste(tree$tip.label[1:10], collapse = ", "), "\n")
          cat("  ... and", n_taxa - 10, "more taxa\n\n")
        }
        
        cat("RECOMMENDATIONS:\n")
        if (!has_branch_lengths) {
          cat("  ⚠ Add branch lengths to your tree for phylogenetic analyses\n")
        }
        if (!is_ultrametric && has_branch_lengths) {
          cat("  ⚠ For BioGeoBEARS: Convert to ultrametric using the controls on the right\n")
        }
        if (is_ultrametric) {
          cat("  ✓ Tree is ready for BioGeoBEARS analysis\n")
        }
      })
    }, error = function(e) {
      output$tree_validation_info <- renderPrint({
        cat("ERROR validating tree:\n")
        cat(e$message, "\n")
      })
    })
  }
})

# Convert tree to ultrametric
observeEvent(input$convert_to_ultrametric, {
  if (is.null(data_store$tree)) {
    output$ultrametric_status <- renderPrint({
      cat("ERROR: No tree loaded!\n")
    })
  } else {
    tryCatch({
      tree <- data_store$tree
      
      # Check if tree has branch lengths
      if (is.null(tree$edge.length)) {
        output$ultrametric_status <- renderPrint({
          cat("ERROR: Tree has no branch lengths!\n")
          cat("Cannot ultrametricize a tree without branch lengths.\n")
        })
      } else {
        # Convert to ultrametric
        method <- input$ultrametric_method_validation
        
        if (method == "equal") {
          # Equal branch lengths
          tree_ultrametric <- tree
          tree_ultrametric$edge.length <- rep(1, length(tree$edge.length))
        } else if (method == "scale") {
          # Scale proportionally
          tree_ultrametric <- ape::chronos(tree, model = "relaxed", quiet = TRUE)
        } else if (method == "penalized") {
          # Penalized least squares
          tree_ultrametric <- ape::chronos(tree, model = "strict", quiet = TRUE)
        }
        
        # Store the converted tree
        data_store$tree_ultrametric <- tree_ultrametric
        
        # Check if conversion was successful
        is_now_ultrametric <- ape::is.ultrametric(tree_ultrametric)
        
        output$ultrametric_status <- renderPrint({
          if (is_now_ultrametric) {
            cat("✓ Tree successfully converted to ultrametric!\n\n")
            cat("Method:", method, "\n")
            cat("Original root-to-tip distances:\n")
            
            root_to_tip_orig <- ape::node.depth.edgelength(tree)[1:length(tree$tip.label)]
            cat("  Max:", round(max(root_to_tip_orig), 4), "\n")
            cat("  Min:", round(min(root_to_tip_orig), 4), "\n")
            cat("  Difference:", round(max(root_to_tip_orig) - min(root_to_tip_orig), 4), "\n\n")
            
            cat("New root-to-tip distances:\n")
            root_to_tip_new <- ape::node.depth.edgelength(tree_ultrametric)[1:length(tree_ultrametric$tip.label)]
            cat("  Max:", round(max(root_to_tip_new), 4), "\n")
            cat("  Min:", round(min(root_to_tip_new), 4), "\n")
            cat("  Difference:", round(max(root_to_tip_new) - min(root_to_tip_new), 4), "\n\n")
            
            cat("✓ This tree is now ready for BioGeoBEARS analysis!\n")
          } else {
            cat("✗ Conversion failed or tree is still not ultrametric.\n")
            cat("Try a different method.\n")
          }
        })
      }
    }, error = function(e) {
      output$ultrametric_status <- renderPrint({
        cat("ERROR converting tree:\n")
        cat(e$message, "\n")
      })
    })
  }
})

# Tree validation plot
output$tree_validation_plot <- renderPlot({
  if (is.null(data_store$tree)) {
    plot(1, type = "n", main = "Phylogenetic Tree", xlab = "", ylab = "", axes = FALSE)
    text(1, 1, "Load a phylogenetic tree first", cex = 1.5, col = "red")
  } else {
    tryCatch({
      tree <- data_store$tree
      
      # Create a nice plot
      plot(tree, 
           main = "Phylogenetic Tree",
           cex = 0.8,
           label.offset = 0.01)
      
      # Add branch length information if available
      if (!is.null(tree$edge.length)) {
        edgelabels(round(tree$edge.length, 3), 
                   cex = 0.6, 
                   bg = "white",
                   adj = c(0.5, -0.2))
      }
      
      # Add title with tree information
      n_taxa <- length(tree$tip.label)
      has_branch_lengths <- !is.null(tree$edge.length)
      is_ultrametric <- ape::is.ultrametric(tree)
      
      status_text <- paste("Taxa:", n_taxa, 
                          "| Branch lengths:", ifelse(has_branch_lengths, "Yes", "No"),
                          "| Ultrametric:", ifelse(is_ultrametric, "Yes", "No"))
      
      mtext(status_text, side = 1, line = 3, cex = 0.8)
      
    }, error = function(e) {
      plot(1, type = "n", main = "Error plotting tree", xlab = "", ylab = "", axes = FALSE)
      text(1, 1, paste("Error:", e$message), cex = 1, col = "red")
    })
  }
})
