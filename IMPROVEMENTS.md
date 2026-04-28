# BioGeoBEARS Shiny App - Improvements Summary

## Changes Made

### 1. **Fixed Directory Management**

**Problem:** The app created output directories (`out_buffers/`, `out_MCP/`, `out_MST/`, `out_irregular_bins/`) but files were not being saved to them.

**Solution:**
- Fixed the `load_output_files.R` to include `out_irregular_bins` in the list of directories to scan
- Ensured that the `load_output_shapefiles()` function in `global.R` correctly loads files from all output directories
- The app now properly saves files to the working directory subdirectories when extrapolation is executed

### 2. **Improved Map Visualization**

**Problem:** Polygons were not appearing on the map, and occurrence points were not colored by species.

**Solution:**
- Modified the Leaflet map rendering to add occurrence points grouped by species
- Each species now has its own layer that can be toggled on/off using the layer control
- Polygons from extrapolation are colored and grouped by species
- Species colors are automatically generated using a rainbow palette
- The map now shows all layers with proper grouping for interactive exploration

### 3. **Fixed Extrapolation Function Calls**

**Problem:** The server was calling non-existent functions and using incorrect parameter names.

**Solution:**
- Fixed the call to `calcRange_irregularBins()` (was `calcRange_irregular_bins`)
- Corrected parameter names: `bins_shapefile` instead of `shape_file` for irregular bins
- Added proper column name handling for occurrence data (converting `spp` to `species` when needed)

### 4. **Added BioGeoBEARS Export**

**Problem:** The Export Files tab did not have a BioGeoBEARS format option.

**Solution:**
- Added a new "BioGeoBEARS" tab in the Export Files section
- Implemented the `download_biogeobears` handler that generates properly formatted `.data` files
- The function handles matrix processing (removing ROOT, filtering empty ranges) and formats output according to BioGeoBEARS specifications

### 5. **Integrated All Export Functions**

**Problem:** Export functions were not accessible to the app.

**Solution:**
- Copied all R functions from the uploaded files to the package R directory:
  - `range_BioGeoBEARS.R` - BioGeoBEARS format export
  - `tnt_matrix.R` - TNT matrix format
  - `range_nexus.R` - NEXUS format
  - `generate_tnt_pae_pce.R` - TNT PAE-PCE script generation
  - `toNDM.R` - NDM format export
- Updated `all_functions.R` to include all functions
- All export handlers now use these integrated functions

## File Structure

```
biogeoshiny/
├── R/
│   ├── calcRange_buffers.R
│   ├── calcRange_convexHull.R
│   ├── calcRange_irregularBins.R
│   ├── generate_tnt_pae_pce.R
│   ├── load_output_files.R
│   ├── range_BioGeoBEARS.R
│   ├── range_nexus.R
│   ├── tnt_matrix.R
│   ├── toNDM.R
│   └── ... (other functions)
├── inst/shiny/
│   ├── ui.R (updated with BioGeoBEARS export tab)
│   ├── server.R (updated with fixes and improvements)
│   ├── global.R (working directory setup)
│   └── all_functions.R (consolidated functions)
```

## Workflow

1. **Data Input** → Load occurrence data (CSV) and phylogenetic tree
2. **Tree Validation** → Validate tree properties
3. **Data Preprocessing** → Remove problematic taxa
4. **Range Extrapolation** → Choose method (Buffer, Convex Hull, MST, Irregular Bins)
   - Files are automatically saved to output directories
   - Shapefiles are loaded and displayed on the map
5. **Visualizations** → Interactive map with:
   - Study area polygon
   - Species occurrence points (colored by species)
   - Extrapolation polygons (colored by species)
   - Layer control for toggling visibility
6. **Export Files** → Download results in multiple formats:
   - BioGeoBEARS (.data)
   - TNT Matrix (.tnt)
   - NEXUS (.nex)
   - TNT PAE-PCE Script (.tnt)
   - NDM (.xyd)

## Key Features

- **Interactive Maps:** Leaflet-based visualization with layer control
- **Multiple Export Formats:** Support for BioGeoBEARS, TNT, NEXUS, NDM
- **Automatic File Management:** Files saved to appropriate subdirectories
- **Species-based Grouping:** All map elements grouped by species for easy exploration
- **Color Coding:** Automatic color assignment to species for visual distinction

## Running the App

```r
library(biogeoshiny)
run_biogeoshiny()
```

The app will:
1. Detect the RStudio project directory
2. Create output directories if they don't exist
3. Set the working directory appropriately
4. Launch the Shiny interface

## Notes

- All output files are saved relative to the working directory
- The app maintains the original workflow steps (1-3) without changes
- Improvements focus on steps 4-6 (extrapolation, visualization, export)
- The map is fully interactive with zoom, pan, and layer toggling capabilities
