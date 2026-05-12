# Debug Notes - Session 2 (May 12, 2026)

## Problem: Irregular Polygon Diversity - Mismatch between n_occurrences and n_species

### Current Status
- ✅ Filtered shapefiles correctly: 14 species maintained (only those in occ_data after preprocessing)
- ✅ n_occurrences is being calculated for some polygons (> 0)
- ❌ BUT: Some polygons show `n_occurrences > 0` with `n_species = 0` and empty `species_list`

### Root Cause Analysis
The issue is that:
1. Occurrence points are being counted correctly in some polygons
2. BUT the species from those occurrence points are NOT being added to the `species_list`
3. This suggests a problem in how occurrence points are being joined to species names

### Code Location
- File: `/mnt/desktop/BRA/inst/shiny/server.R`
- Function: `compute_irregular_richness_from_layers()` (lines 1960-2095)
- Specific issue: Lines 2043-2062 (occurrence counting logic)

### Debug Output from Last Run
```
[DEBUG] occ rows: 91 
[DEBUG] occ taxa: 14 
[DEBUG] valid_species: Belostoma_estevezae, Belostoma_lariversi_var1, ... (14 species)
[DEBUG] filtered shapefiles count: 14 
[DEBUG] filtered shapefiles names: MST_mst_Belostoma_estevezae, ... (14 species)
```

### Next Steps for Tomorrow
1. Add debug logs inside `compute_irregular_richness_from_layers()` to see:
   - How many occurrence points are in each polygon
   - What species_id values are being assigned to those points
   - Why species_list is empty when n_occurrences > 0

2. Check if the occurrence points are being properly converted to sf objects with correct CRS

3. Verify that the `st_join()` between bins and occurrence points is working correctly

4. Ensure that occurrence point species names are being properly added to the final `species_list`

### Related Commits
- `a617808`: Filter shapefiles to only include species in preprocessed occurrence data
- `138eb43`: Extract species name from shapefile names for proper filtering

### Test Data
- 91 occurrence records
- 14 taxa
- Arctic Subregion: Shows n_occurrences > 0 but n_species = 0 (example of the bug)
