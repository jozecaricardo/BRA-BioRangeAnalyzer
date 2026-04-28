#' generate_tnt_pae_pce
#'
#' Creates a TNT script file for automated PAE-PCE (Parsimony Analysis of Endemicity
#' with Parsimony Constraint) analysis with multiple rounds of character removal.
#'
#' @param matrix_file Character. Path to the TNT matrix file (default: "out_TNT/pres_abs.tnt")
#' @param output_file Character. Name of the output TNT script file (default: "PAE_PCE_tntAUTO.run")
#' @param max_iterations Integer. Maximum number of PAE-PCE rounds to perform (default: 10)
#' @param search_replicates Integer. Number of search replicates to perform (default: 4)
#' @param search_method Character. Search method to use: "traditional" or "new_technology" (default: "new_technology")
#' @param output_dir Character. Directory for output files (default: "out_TNT/")
#'
#' @details
#' This function generates a TNT script that performs iterative PAE-PCE analysis:
#'
#' **Search Methods:**
#' - `"traditional"`: Uses traditional search with TBR branch swapping
#'   - Command: `mult 1000 = hold 100 tbr; bb=tbr;`
#' - `"new_technology"`: Uses new technology search with sectorial searches
#'   - Command: `xmult= hits 7 rep 5 rss level 10;`
#'
#' **Iterative Process:**
#' 1. Performs parsimony analysis
#' 2. Calculates CI (Consistency Index) and RI (Retention Index) for all characters
#' 3. Identifies non-homoplastic synapomorphies (CI=1, RI=1)
#' 4. Deactivates these characters cumulatively
#' 5. Repeats until no more synapomorphies are found or max_iterations is reached
#'
#' **Output Files (per round N):**
#' - `Round_OutputN.txt`: Complete log of the round
#' - `Round_CI_RI_TableN.txt`: CI/RI table for all characters
#' - `Round_SynapomorphiesN.txt`: List of non-homoplastic synapomorphies
#' - `Round_MPTN.tre`: Most parsimonious trees
#' - `Round_ConsensusN.tmp`: Consensus tree (TNT format)
#' - `Round_TreeN.tre`: Consensus tree (Newick format)
#' - `Round_MapN.svg`: Tree visualization with mapped synapomorphies
#'
#' **Summary Files:**
#' - `PAE_PCE_Summary.txt`: Summary of all rounds
#' - `resul_CI_RI.log`: CI/RI statistics for all rounds
#'
#' @return Invisibly returns the path to the generated script file
#'
#' @examples
#' \dontrun{
#' # Generate script with default settings (4 replicates, new technology)
#' generate_tnt_pae_pce()
#'
#' # Generate script with 20 replicates and traditional search
#' generate_tnt_pae_pce(
#'   search_replicates = 20,
#'   search_method = "traditional"
#' )
#'
#' # Generate script with custom settings
#' generate_tnt_pae_pce(
#'   matrix_file = "data/my_matrix.tnt",
#'   output_file = "my_pae_pce.run",
#'   max_iterations = 15,
#'   search_replicates = 10,
#'   search_method = "new_technology"
#' )
#' }
#'
#' @export
generate_tnt_pae_pce <- function(matrix_file = "out_TNT/pres_abs.tnt",
                                 output_file = "PAE_PCE_tntAUTO.run",
                                 max_iterations = 10,
                                 search_replicates = 4,
                                 search_method = "new_technology",
                                 output_dir = "out_TNT/") {

  ######################
  ## Input validation ##
  ######################

  if (!file.exists(matrix_file)) {
    stop(paste0("Error: Matrix file not found: ", matrix_file))
  }

  valid_methods <- c("traditional", "new_technology")
  if (!(search_method %in% valid_methods)) {
    stop(paste0("Error: search_method must be one of: ", paste(valid_methods, collapse = ", ")))
  }

  if (!is.numeric(max_iterations) || max_iterations < 1) {
    stop("Error: max_iterations must be a positive integer")
  }

  if (!is.numeric(search_replicates) || search_replicates < 1) {
    stop("Error: search_replicates must be a positive integer")
  }

  ######################
  ## Generate script ###
  ######################

  message("Generating TNT PAE-PCE script based on template...")
  message(paste0("  Matrix: ", matrix_file))
  message(paste0("  Output: ", output_file))
  message(paste0("  Max iterations: ", max_iterations))
  message(paste0("  Search replicates: ", search_replicates))
  message(paste0("  Search method: ", search_method))

  # Open connection to write script
  con <- file(output_file, "w", encoding = "UTF-8")

  # Write header (lines 1-56 from template)
  cat("macro= ;\n", file = con)
  cat("clb;\n", file = con)
  cat("mxram 1000;\n", file = con)
  cat("macro * 15 1000;\n", file = con)
  cat("macfloat 3;\n", file = con)
  cat("piwe-;\n", file = con)
  cat("\n", file = con)
  cat("report-;\n", file = con)
  cat("\n", file = con)
  cat(paste0("proc ", matrix_file, " ;\n"), file = con)
  cat("\n", file = con)
  cat("hold 99999;\n", file = con)
  cat("\n", file = con)
  cat(paste0("log ", output_dir, "PAE_PCE_Summary.txt ;\n"), file = con)
  cat("quote\n", file = con)
  cat("================================================================================\n", file = con)
  cat("           partially AUTOMATED PAE-PCE ANALYSIS - SUMMARY\n", file = con)
  cat("================================================================================\n", file = con)
  cat(paste0("Matrix: ", matrix_file, "\n"), file = con)
  cat(paste0("Max iterations: ", max_iterations, " (by author)\n"), file = con)
  cat("================================================================================\n", file = con)
  cat(";\n", file = con)
  cat("log/ ;\n", file = con)
  cat("\n", file = con)
  cat("\n", file = con)
  cat("\n", file = con)
  cat(paste0("log+ ", output_dir, "resul_CI_RI.log;\n"), file = con)
  cat("quote RUN, CI, RI, tree_length;\n", file = con)
  cat("log/;\n", file = con)
  cat("\n", file = con)
  cat("var: removed_count;\n", file = con)
  cat("\n", file = con)
  cat("/* ========================================================================== */\n", file = con)
  cat("/* MAIN ITERATION LOOP                                                        */\n", file = con)
  cat("/* ========================================================================== */\n", file = con)
  cat("\n", file = con)
  cat("/* Define variables for automation */\n", file = con)
  cat("var: total_removed conta contagem;\n", file = con)
  cat("var: theminchar themaxchar ci_val[(nchar + 1 )] ri_val[(nchar + 1)] rc_val[(nchar + 1)] char_i num_chars;\n", file = con)
  cat("var: THEMINtree THEMAXtree; \n", file = con)
  cat("var: THIStree CItree RItree RCtree nome list_idx;\n", file = con)
  cat("\n", file = con)
  cat("/* [OK] NEW: Global array to store ALL deactivated characters */\n", file = con)
  cat("var: all_deactivated[(nchar + 1)] total_deactivated;\n", file = con)
  cat("\n", file = con)
  cat("set num_chars nchar;\n", file = con)
  cat("set total_removed 0;\n", file = con)
  cat("set total_deactivated 0;\n", file = con)
  cat("\n", file = con)
  cat("outgroup ROOT;\n", file = con)
  cat("\n", file = con)
  cat("set theminchar 0;\n", file = con)
  cat("set themaxchar 1;\n", file = con)
  cat("\n", file = con)
  cat("cnames=;\n", file = con)
  cat("\n", file = con)

  # Generate rounds (replicate the structure search_replicates times)
  for (round_num in 1:search_replicates) {

    # Round header
    cat("/* --------------ROUNDS------------------- */\n", file = con)
    cat("\n", file = con)
    cat("quote\n", file = con)
    cat("============================================================================\n", file = con)
    cat(paste0("                ROUND ", round_num, " - Starting iteration\n"), file = con)
    cat("============================================================================\n", file = con)
    cat(";\n", file = con)
    cat("\n", file = con)
    cat(paste0("log ", output_dir, "Round_Output", round_num, ".txt;\n"), file = con)
    cat("\n", file = con)
    cat("quote\n", file = con)
    cat("============================================================================\n", file = con)
    cat(paste0("                ROUND ", round_num, " - PAE-PCE Iteration\n"), file = con)
    cat("============================================================================\n", file = con)
    cat(";\n", file = con)
    cat("\n", file = con)

    # Clear memory
    cat("/* Clear memory */\n", file = con)
    cat("keep 0 ;\n", file = con)
    cat("ttags - ;\n", file = con)
    cat("rseed 0 ; \n", file = con)
    cat("\n", file = con)

    # Parsimony search - DIFFERENT BASED ON METHOD
    if (search_method == "traditional") {
      cat("/* Traditional search */  \n", file = con)
      cat("quote Running parsimony analysis... ;\n", file = con)
      cat("mult 1000 = hold 100 tbr;\n", file = con)
      cat("bb=tbr;\n", file = con)
      cat("reroot [;\n", file = con)
    } else {
      cat("/* New tech search */  \n", file = con)
      cat("quote Running parsimony analysis... ;\n", file = con)
      cat("xmult= hits 7 rep 5 rss level 10;\n", file = con)
      cat("reroot [;\n", file = con)
    }
    cat("\n", file = con)

    # Save most parsimonious trees
    cat("/* Save most parsimonious trees */\n", file = con)
    cat(paste0("tsave *", output_dir, "Round_MPT", round_num, ".tre ;\n"), file = con)
    cat("save ;\n", file = con)
    cat("tsave / ;\n", file = con)
    cat("quote Most parsimonious trees saved ;\n", file = con)
    cat("\n", file = con)

    # Calculate strict consensus
    cat("/* Calculate strict consensus */\n", file = con)
    cat("if ((ntrees+1) > 1)\n", file = con)
    cat("    quote Calculating strict consensus... ;\n", file = con)
    cat("    nelsen*;      /*  calculate strict consensus   */\n", file = con)
    cat("    tchoose {strict};\n", file = con)
    cat(paste0("    tsave ", output_dir, "Round_Consensus", round_num, ".tmp ;\n"), file = con)
    cat("    save;\n", file = con)
    cat("    tsave/;      /*  close file with saved consensus   */\n", file = con)
    cat("    tt-;\n", file = con)
    cat("    ttags = ;\n", file = con)
    cat("    apo ] ;\n", file = con)
    cat(paste0("    ttag & ", output_dir, "Round_Map", round_num, ".svg ;\n"), file = con)
    cat("    taxname=;\n", file = con)
    cat(paste0("    export - ", output_dir, "Round_Tree", round_num, ".tre ;\n"), file = con)
    cat("    tt-;    \n", file = con)
    cat("    quote Consensus tree saved ;\n", file = con)
    cat("    keep 0;\n", file = con)
    cat("else\n", file = con)
    cat("    quote Last tree... ;\n", file = con)
    cat(paste0("    tsave ", output_dir, "Round_Consensus", round_num, ".tmp ;\n"), file = con)
    cat("    save 0;\n", file = con)
    cat("    tsave/;      /*  close file with saved consensus   */\n", file = con)
    cat("    tt-;\n", file = con)
    cat("    ttags = ;\n", file = con)
    cat("    apo ] ;\n", file = con)
    cat(paste0("    ttag & ", output_dir, "Round_Map", round_num, ".svg ;\n"), file = con)
    cat("    taxname=;\n", file = con)
    cat(paste0("    export - ", output_dir, "Round_Tree", round_num, ".tre ;\n"), file = con)
    cat("    tt-;    \n", file = con)
    cat("    quote Last tree saved ;\n", file = con)
    cat("    keep 0;\n", file = con)
    cat("end;    \n", file = con)
    cat("\n", file = con)

    # Calculate and save complete CI/RI statistics
    cat("/* ====================================================================== */\n", file = con)
    cat("/* CALCULATE AND SAVE COMPLETE CI/RI STATISTICS                          */\n", file = con)
    cat("/* ====================================================================== */\n", file = con)
    cat("\n", file = con)
    cat("quote\n", file = con)
    cat("----------------------------------------------------------------------------\n", file = con)
    cat("Tree Statistics:\n", file = con)
    cat("----------------------------------------------------------------------------\n", file = con)
    cat(";\n", file = con)
    cat("\n", file = con)
    cat("/* get the tree */\n", file = con)
    cat(paste0("short ", output_dir, "Round_Consensus", round_num, ".tmp;\n"), file = con)
    cat("\n", file = con)
    cat("set THEMINtree minsteps;\n", file = con)
    cat("set THEMAXtree maxsteps;\n", file = con)
    cat("set THIStree length[0];\n", file = con)
    cat("\n", file = con)
    cat(paste0("log+ ", output_dir, "Round_CI_RI_Table", round_num, ".txt ;\n"), file = con)
    cat("quote\n", file = con)
    cat("============================================================================\n", file = con)
    cat(paste0("        ROUND ", round_num, " - COMPLETE CI AND RI TABLE FOR ALL CHARACTERS\n"), file = con)
    cat("============================================================================\n", file = con)
    cat("\n", file = con)
    cat("Character | CI    | RI    | RC    | Status\n", file = con)
    cat("----------|-------|-------|-------|----------------------------------\n", file = con)
    cat(";\n", file = con)
    cat("\n", file = con)

    # Loop through all characters and print CI/RI
    cat("/* Loop through all characters and print CI/RI */\n", file = con)
    cat("loop 0 ('num_chars' - 1) ;    \n", file = con)
    cat("    set char_i length[0 #1];\n", file = con)
    cat("\n", file = con)
    cat("    if ('char_i' > 0)\n", file = con)
    cat("        set themaxchar maxsteps[#1] ;\n", file = con)
    cat("        set theminchar minsteps[#1] ;\n", file = con)
    cat("\n", file = con)
    cat("        /* Calculate CI */\n", file = con)
    cat("        if ('char_i' == 0)\n", file = con)
    cat("            set ci_val[#1] 1 ;\n", file = con)
    cat("        else\n", file = con)
    cat("            set ci_val[#1] 'theminchar'/'char_i' ;\n", file = con)
    cat("        end;\n", file = con)
    cat("\n", file = con)
    cat("        /* Calculate RI */\n", file = con)
    cat("        if (('themaxchar'-'theminchar') == 0)\n", file = con)
    cat("            set ri_val[#1] 1 ;\n", file = con)
    cat("        else\n", file = con)
    cat("            set ri_val[#1] ('themaxchar'-'char_i')/('themaxchar'-'theminchar') ;\n", file = con)
    cat("        end;\n", file = con)
    cat("\n", file = con)
    cat("        /* Calculate RC */\n", file = con)
    cat("        set rc_val[#1]  'ri_val[#1]'*'ci_val[#1]';\n", file = con)
    cat("\n", file = con)
    cat("        /* Print character info */\n", file = con)
    cat("        if ('ci_val[#1]' == 1)\n", file = con)
    cat("            if ('ri_val[#1]' == 1)\n", file = con)
    cat("                quote #1 | 'ci_val[#1]' | 'ri_val[#1]' | 'rc_val[#1]' | NON-HOMOPLASTIC SYNAPOMORPHY ;\n", file = con)
    cat("            else\n", file = con)
    cat("                quote #1 | 'ci_val[#1]' | 'ri_val[#1]' | 'rc_val[#1]' | Homoplastic ;\n", file = con)
    cat("            end;\n", file = con)
    cat("        else\n", file = con)
    cat("            quote #1 | 'ci_val[#1]' | 'ri_val[#1]' | 'rc_val[#1]' | Homoplastic ;\n", file = con)
    cat("        end;\n", file = con)
    cat("    else\n", file = con)
    cat("        /* Character is INACTIVE - skip */\n", file = con)
    cat("        quote #1 | INACTIVE | INACTIVE | INACTIVE | Character deactivated in previous round;\n", file = con)
    cat("    end;\n", file = con)
    cat("stop;\n", file = con)
    cat("\n", file = con)

    # Tree statistics
    cat("/* tree statistics (usando valores capturados anteriormente) */\n", file = con)
    cat("set CItree 'THEMINtree'/'THIStree';\n", file = con)
    cat("set RItree ('THEMAXtree'-'THIStree')/('THEMAXtree'-'THEMINtree');\n", file = con)
    cat("set RCtree 'RItree'*'CItree';\n", file = con)
    cat("\n", file = con)
    cat(paste0("log+ ", output_dir, "resul_CI_RI.log;\n"), file = con)
    cat("\n", file = con)
    cat("quote\n", file = con)
    cat("============================================================================\n", file = con)
    cat(paste0("RUN_", round_num, ", 'CItree', 'RItree', 'RCtree', 'THIStree';\n"), file = con)
    cat("\n", file = con)
    cat("log/ ;\n", file = con)
    cat("\n", file = con)

    # Identify and list non-homoplastic synapomorphies
    cat("/* ====================================================================== */\n", file = con)
    cat("/* IDENTIFY AND LIST NON-HOMOPLASTIC SYNAPOMORPHIES                      */\n", file = con)
    cat("/* ====================================================================== */\n", file = con)
    cat("\n", file = con)
    cat(paste0("log+ ", output_dir, "Round_Synapomorphies", round_num, ".txt ;\n"), file = con)
    cat("quote\n", file = con)
    cat("============================================================================\n", file = con)
    cat(paste0("ROUND ", round_num, " - NON-HOMOPLASTIC SYNAPOMORPHIES (CI=1, RI=1)\n"), file = con)
    cat("============================================================================\n", file = con)
    cat("\n", file = con)
    cat("These characters (taxa/species) will be removed for the next round:\n", file = con)
    cat("\n", file = con)
    cat("Character(Taxon) | Run\n", file = con)
    cat("------------------------------------------------------------    \n", file = con)
    cat(";\n", file = con)
    cat("\n", file = con)
    cat("set removed_count 0;\n", file = con)
    cat("set conta 0;\n", file = con)
    cat("set contagem 0;\n", file = con)
    cat("\n", file = con)

    # Loop to identify and list synapomorphies
    cat("/* Loop to identify and list synapomorphies */\n", file = con)
    cat("loop 0 ('num_chars' - 1) ;\n", file = con)
    cat("    set char_i length[0 #1];\n", file = con)
    cat("\n", file = con)
    cat("    if ('char_i' > 0)\n", file = con)
    cat("        set themaxchar maxsteps[#1] ;\n", file = con)
    cat("        set theminchar minsteps[#1] ;\n", file = con)
    cat("\n", file = con)
    cat("        /* Calculate CI */\n", file = con)
    cat("        if ('char_i' == 0)\n", file = con)
    cat("            set ci_val[#1] 1 ;\n", file = con)
    cat("        else\n", file = con)
    cat("            set ci_val[#1] 'theminchar'/'char_i' ;\n", file = con)
    cat("        end ;\n", file = con)
    cat("\n", file = con)
    cat("        /* Calculate RI */\n", file = con)
    cat("        if ( ('themaxchar'-'theminchar') == 0 )\n", file = con)
    cat("            set ri_val[#1] 1 ;\n", file = con)
    cat("        else\n", file = con)
    cat("            set ri_val[#1] ('themaxchar'-'char_i')/('themaxchar'-'theminchar') ;\n", file = con)
    cat("        end ;\n", file = con)
    cat("\n", file = con)
    cat("        set nome #1;\n", file = con)
    cat("\n", file = con)
    cat("        /* If CI=1 AND RI=1, it's a non-homoplastic synapomorphy */\n", file = con)
    cat("        if ('ci_val[#1]' == 1)\n", file = con)
    cat("            if ('ri_val[#1]' == 1)\n", file = con)
    cat(paste0("                quote Character 'nome' | Run ", round_num, ";\n"), file = con)
    cat("                set conta ++;\n", file = con)
    cat("                set removed_count ++ ;\n", file = con)
    cat("            else\n", file = con)
    cat("                quote It is a homoplastic synapomorphy;\n", file = con)
    cat("            end;\n", file = con)
    cat("        else\n", file = con)
    cat("            quote It is a homoplastic synapomorphy;\n", file = con)
    cat("        end;\n", file = con)
    cat("    else\n", file = con)
    cat("        /* Character is INACTIVE - skip */\n", file = con)
    cat("        quote Character #1 is INACTIVE (deactivated in previous round);\n", file = con)
    cat("    end;            \n", file = con)
    cat("stop;\n", file = con)
    cat("\n", file = con)

    # If synapomorphies found, deactivate them
    cat("if ('removed_count' > 0)\n", file = con)
    cat("    var: caraExtract['conta'];\n", file = con)
    cat("    loop 0 ('num_chars' - 1);\n", file = con)
    cat("        set char_i length[0 #1];\n", file = con)
    cat("        if ('char_i' > 0)\n", file = con)
    cat("            if ('ci_val[#1]' == 1)\n", file = con)
    cat("                if ('ri_val[#1]' == 1)\n", file = con)
    cat("                    set caraExtract['contagem'] #1;\n", file = con)
    cat("                    set all_deactivated['total_deactivated'] #1;\n", file = con)
    cat("                    set total_deactivated ++;\n", file = con)
    cat("                    set contagem ++;\n", file = con)
    cat("                    set total_removed ++;\n", file = con)
    cat("                else\n", file = con)
    cat("                    quote Not counted;\n", file = con)
    cat("                end;\n", file = con)
    cat("            else\n", file = con)
    cat("                quote Not counted;\n", file = con)
    cat("            end;\n", file = con)
    cat("        else\n", file = con)
    cat("            /* Character is INACTIVE - skip */\n", file = con)
    cat("            quote Character #1 is INACTIVE - not counted;\n", file = con)
    cat("        end;        \n", file = con)
    cat("    stop;\n", file = con)
    cat("\n", file = con)
    cat("    quote\n", file = con)
    cat("    ----------------------------------------------------------------------------\n", file = con)
    cat("    Total synapomorphies found: 'removed_count'\n", file = con)
    cat("    Characters removed:;\n", file = con)
    cat("    loop 0 ('conta' - 1) ;\n", file = con)
    cat("        quote character 'caraExtract[#1]';\n", file = con)
    cat("    stop;\n", file = con)
    cat("    \n", file = con)
    cat("    quote\n", file = con)
    cat("    ============================================================================\n", file = con)
    cat("    ;\n", file = con)
    cat("\n", file = con)
    cat("\n", file = con)
    cat("    set list_idx 0;\n", file = con)
    cat("    loop 0 ('total_deactivated' - 1) ;\n", file = con)
    cat("        ccode] 'all_deactivated['list_idx']';\n", file = con)
    cat("        quote Deactivating character 'all_deactivated['list_idx']';\n", file = con)
    cat("        set list_idx ++;\n", file = con)
    cat("    stop;\n", file = con)
    cat("\n", file = con)
    cat("    /* [OK] PROTECTION 2: Update num_chars after deactivation */\n", file = con)
    cat("    set num_chars nchar;\n", file = con)
    cat("\n", file = con)
    cat("    quote\n", file = con)
    cat("    ----------------------------------------------------------------------------\n", file = con)
    cat(paste0("    Round ", round_num, " Summary:\n"), file = con)
    cat("    Characters removed this round: 'all_deactivated[0-('total_deactivated' - 1) &45]'\n", file = con)
    cat("    Total characters removed so far: 'total_removed'\n", file = con)

    # Add additional info for rounds after the first
    if (round_num > 1) {
      cat("    Total cumulative deactivated: 'total_deactivated'\n", file = con)
      cat("    Active characters remaining: 'num_chars'\n", file = con)
    }

    cat("    ----------------------------------------------------------------------------\n", file = con)
    cat("    ;\n", file = con)
    cat("    log/;\n", file = con)
    cat("\n", file = con)
    cat(paste0("    log+ ", output_dir, "PAE_PCE_Summary.txt;\n"), file = con)
    cat("\n", file = con)
    cat("    quote\n", file = con)
    cat("    ================================================================================\n", file = con)
    cat(paste0("                        REACHED ITERATION NUMBER ", round_num, "\n"), file = con)
    cat("    ================================================================================\n", file = con)
    cat(paste0("    Maximum iterations (", max_iterations, ")\n"), file = con)
    cat("    Total characters removed: 'total_removed'\n", file = con)
    cat("    Total cumulative deactivated: 'total_deactivated'\n", file = con)
    cat("\n", file = con)

    # Different message for last round vs. intermediate rounds
    if (round_num == search_replicates) {
      cat("    Note: Analysis has stopped due to iteration limit.\n", file = con)
    } else {
      cat("    Note: Analysis MIGHT BE stopped due to iteration limit.\n", file = con)
    }

    cat("    Consider increasing max_iterations if more rounds are needed.\n", file = con)
    cat("    ================================================================================\n", file = con)
    cat("    ;\n", file = con)
    cat("    log/;\n", file = con)

    # If no synapomorphies found, stop
    cat("else\n", file = con)
    cat("    log/;\n", file = con)
    cat("    quote\n", file = con)
    cat("    ============================================================================\n", file = con)
    cat(paste0("    STOPPING: No more non-homoplastic synapomorphies found in Round ", round_num, "\n"), file = con)
    cat("    ============================================================================\n", file = con)
    cat("    ;\n", file = con)
    cat("    \n", file = con)
    cat(paste0("    log+ ", output_dir, "PAE_PCE_Summary.txt;\n"), file = con)
    cat("    quote\n", file = con)
    cat("    ================================================================================\n", file = con)
    cat("                        ANALYSIS COMPLETED SUCCESSFULLY\n", file = con)
    cat("    ================================================================================\n", file = con)
    cat(paste0("    Stopped at iteration: ", round_num, "\n"), file = con)
    cat("    Reason: No more non-homoplastic synapomorphies (CI=1, RI=1)\n", file = con)
    cat("    Total characters removed: 'total_removed'\n", file = con)
    cat("    Total cumulative deactivated: 'total_deactivated'\n", file = con)
    cat("    ================================================================================\n", file = con)
    cat("    ;\n", file = con)
    cat("    log/;\n", file = con)
    cat("    silent-;\n", file = con)
    cat("    report=;\n", file = con)
    cat("    proc/;\n", file = con)
    cat("end;\n", file = con)

    # Add separator between rounds (except after last round)
    if (round_num < search_replicates) {
      cat("\n", file = con)
      cat(paste0("/* ---------------- RUN ", round_num + 1, " ----------------- */\n"), file = con)
      cat("set theminchar 0;\n", file = con)
      cat("set themaxchar 1;\n", file = con)
      cat("\n", file = con)
      cat("cnames=;\n", file = con)
      cat("\n", file = con)
    }
  }

  # Final closing
  cat("silent-;\n", file = con)
  cat("report=;\n", file = con)
  cat("proc/;", file = con)

  # Close connection
  close(con)

  message(paste0("[OK] TNT script generated successfully: ", output_file))
  message(paste0("   Total rounds: ", search_replicates))
  message(paste0("   Script will stop automatically when no more synapomorphies are found"))

  invisible(output_file)
}
