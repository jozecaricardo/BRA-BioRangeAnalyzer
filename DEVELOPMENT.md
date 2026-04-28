# BioGeoBEARS Shiny Application - Development Guide

## Project Structure

```
biogeoshiny/
├── DESCRIPTION              # Package metadata
├── NAMESPACE               # Package exports
├── README.md              # User documentation
├── DEVELOPMENT.md         # This file
├── R/
│   ├── run_biogeoshiny.R  # Main app launcher
│   └── utils.R            # Utility functions
├── inst/
│   └── shiny/
│       ├── ui.R           # User interface
│       ├── server.R       # Server logic
│       ├── global.R       # Global configuration
│       ├── modules/       # Shiny modules
│       │   ├── mod_data_input.R
│       │   ├── mod_range_extrapolation.R
│       │   ├── mod_export_formats.R
│       │   ├── mod_bgb_setup.R
│       │   └── mod_bgb_analysis.R
│       ├── css/           # Stylesheets (future)
│       ├── js/            # JavaScript (future)
│       ├── www/           # Static assets
│       └── Rmd/           # R Markdown documents
└── data-raw/              # Raw data for package
```

## Module Architecture

The application uses Shiny modules to organize code into logical, reusable components. Each module follows the pattern:

```R
# UI function
mod_name_ui <- function(id) {
  ns <- shiny::NS(id)
  # UI elements
}

# Server function
mod_name_server <- function(id, ...) {
  shiny::moduleServer(id, function(input, output, session) {
    # Server logic
    shiny::reactive({
      # Return reactive values
    })
  })
}
```

### Module Descriptions

#### mod_data_input
Handles file uploads for occurrence data and phylogenetic trees. Validates input and provides preview of loaded data.

**Inputs:**
- occurrence_file: CSV file with species, longitude, latitude
- tree_file: Newick/Nexus phylogenetic tree

**Outputs:**
- occurrence: Data frame with occurrence data
- tree: Phylogenetic tree object
- loaded: Boolean indicating successful load

#### mod_range_extrapolation
Implements range extrapolation methods (buffers, convex hull, MST). Generates presence-absence matrices from occurrence points.

**Inputs:**
- method: Extrapolation method (buffer, convex_hull, mst)
- buffer_width: Width of circular buffers (meters)
- mean_dist: Use mean distance between points
- resolution: Grid resolution (decimal degrees)

**Outputs:**
- pres_abs: Presence-absence matrix
- geometry: Spatial geometry object
- completed: Boolean indicating completion

#### mod_export_formats
Exports presence-absence matrix to multiple formats (BioGeoBEARS, NEXUS, TNT, NDM).

**Inputs:**
- export_formats: Formats to export
- ndm_grid_resolution: Grid resolution for NDM
- generate_tnt_script: Generate PAE-PCE script
- max_iterations: Maximum PAE-PCE iterations

**Outputs:**
- biogeobears: Path to BioGeoBEARS file
- nexus: Path to NEXUS file
- tnt: Path to TNT file
- ndm: Path to NDM file
- completed: Boolean indicating completion

#### mod_bgb_setup
Configures BioGeoBEARS models and parameters.

**Inputs:**
- models: Selected models (DEC, DEC+J, etc.)
- d_min, d_max: Dispersal parameter bounds
- e_min, e_max: Extinction parameter bounds
- j_min, j_max: Jump dispersal parameter bounds
- time_stratification: Use time stratification
- distance_matrix: Use distance matrix
- dispersal_multiplier: Use dispersal multiplier matrix

**Outputs:**
- models: Selected models
- parameters: Parameter bounds
- advanced_options: Advanced configuration
- optimizer: Optimization settings
- completed: Boolean indicating completion

#### mod_bgb_analysis
Runs BioGeoBEARS analyses and displays results.

**Inputs:**
- run_analysis: Button to start analysis

**Outputs:**
- models_results: List of model results
- model_comparison: Comparison table
- completed: Boolean indicating completion

## Data Flow

```
mod_data_input
    ↓
    ├→ occurrence data
    └→ phylogenetic tree
         ↓
    mod_range_extrapolation
         ↓
         ├→ presence-absence matrix
         └→ shapefiles
              ↓
         mod_export_formats
              ├→ BioGeoBEARS file
              ├→ NEXUS file
              ├→ TNT file
              └→ NDM file
              
         mod_bgb_setup
              ├→ model selection
              └→ parameter configuration
                   ↓
              mod_bgb_analysis
                   ↓
                   ├→ model results
                   └→ comparison table
```

## Integration with Project Functions

The following functions from the project should be integrated into the modules:

### Range Extrapolation
- `calcRange_buffers()` → mod_range_extrapolation
- `calcRange_convexHull()` → mod_range_extrapolation
- `calcRange_irregularBins()` → mod_range_extrapolation (MST)

### Export Formats
- `range_BioGeoBEARS()` → mod_export_formats
- `range_nexus()` → mod_export_formats
- `toNDM()` → mod_export_formats
- `generate_tnt_pae_pce()` → mod_export_formats

### BioGeoBEARS Analysis
- BioGeoBEARS model setup and execution
- `calc_lrt()` → mod_bgb_analysis (for model comparison)

## Development Roadmap

### Phase 1: Core Functionality (Current)
- [x] Basic module structure
- [x] Data input module
- [x] Range extrapolation module (placeholder)
- [x] Export formats module (placeholder)
- [x] BioGeoBEARS setup module
- [x] Analysis module (placeholder)

### Phase 2: Integration
- [ ] Integrate calcRange functions
- [ ] Integrate export functions
- [ ] Implement BioGeoBEARS execution
- [ ] Add model comparison statistics

### Phase 3: Visualization
- [ ] Add interactive maps (Leaflet)
- [ ] Add phylogenetic tree visualization
- [ ] Add result plots and charts
- [ ] Add ancestral state reconstruction plots

### Phase 4: Polish
- [ ] Add help documentation
- [ ] Improve error handling
- [ ] Add data validation
- [ ] Performance optimization
- [ ] Add unit tests

## Key Functions to Implement

### In utils.R

```R
# BioGeoBEARS execution wrapper
run_biogeobears_models <- function(pres_abs, tree, models, parameters, ...) {
  # Setup and run BioGeoBEARS
}

# Model comparison
compare_models <- function(results) {
  # Calculate AICc, LRT, etc.
}

# Result formatting
format_bgb_results <- function(results) {
  # Format for display
}
```

### In modules

```R
# Integration with calcRange functions
perform_range_extrapolation <- function(occurrence_data, method, params) {
  # Call appropriate calcRange function
}

# Integration with export functions
export_to_formats <- function(pres_abs, formats, params) {
  # Call appropriate export function
}
```

## Testing

### Manual Testing Checklist

- [ ] Data upload and validation
- [ ] Range extrapolation with different methods
- [ ] Export to all formats
- [ ] BioGeoBEARS model configuration
- [ ] Model execution and results display
- [ ] Download functionality
- [ ] Error handling for invalid inputs

### Automated Testing (Future)

```R
# tests/testthat/test_modules.R
test_that("mod_data_input loads CSV correctly", {
  # Test data loading
})

test_that("mod_range_extrapolation creates matrix", {
  # Test range extrapolation
})
```

## Performance Considerations

1. **Large Datasets**: Implement progress indicators for long-running operations
2. **Memory Management**: Use reactive values efficiently
3. **Caching**: Cache expensive computations where possible
4. **Parallel Processing**: Use parallel package for BioGeoBEARS runs

## Future Enhancements

1. **Advanced Visualization**
   - Interactive phylogenetic tree viewer
   - Animated ancestral state reconstruction
   - 3D geographic visualization

2. **Additional Features**
   - Batch analysis of multiple datasets
   - Custom model creation
   - Stochastic mapping
   - Biogeographic stochastic mapping (BSM)

3. **Integration**
   - Direct connection to online databases (GBIF, etc.)
   - Integration with other biogeographic tools
   - Export to publication-ready formats

4. **Documentation**
   - In-app tutorials
   - Video guides
   - Detailed help pages
   - Example datasets

## Contributing

When adding new features:

1. Follow the existing module structure
2. Use consistent naming conventions
3. Add comments and documentation
4. Test thoroughly
5. Update this development guide

## Troubleshooting

### Common Issues

**Issue**: "Could not find shiny app directory"
- **Solution**: Make sure the package is properly installed with `devtools::install_local()`

**Issue**: Missing package dependencies
- **Solution**: Install all required packages listed in DESCRIPTION

**Issue**: Slow performance
- **Solution**: Check for large file uploads, optimize reactive expressions

## References

- [Shiny Documentation](https://shiny.rstudio.com/)
- [Shiny Modules](https://shiny.rstudio.com/articles/modules.html)
- [BioGeoBEARS Documentation](http://phylo.wikidot.com/biogeobears)
- [Wallace Project](https://wallaceecomod.github.io/wallace/)
