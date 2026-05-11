# Debug Notes - Remove Duplicates Issue

## Current Status
- **Error Location**: Remove Duplicates step (line ~1604 in server.R)
- **Error Message**: "Occurrence data must include columns spp, long and lat."
- **Issue**: The `occ_data` arriving at Remove Duplicates is missing `long` and/or `lat` columns

## Investigation Done
1. ✅ Fixed `remove_duplicates` to use `c("spp", "long", "lat")` directly instead of `intersect()`
2. ✅ Removed duplicate `occ_data <- data_store$occurrence` assignment
3. ✅ Added guarantee checks in extrapolation to catch missing columns

## Root Cause (Not Yet Found)
One of these preprocessing steps is dropping the `long`/`lat` columns:
1. Filter Points Outside Study Area
2. Detect Singletons/Doubletons
3. Remove Problematic Taxa
4. Harmonize Tree ↔ Data

## Next Steps (Tomorrow)
1. Add debug log at line ~1604 in Remove Duplicates to print:
   - `names(occ_data)` - what columns exist
   - `nrow(occ_data)` - how many rows
   - `head(occ_data)` - first few rows

2. Then trace backwards through each preprocessing step to find which one drops columns

3. Fix the offending step

## Code Locations
- Remove Duplicates: line 1595-1655
- Detect Singletons: line ~1350-1455
- Remove Problematic: line 1457-1593
- Harmonize Tree: line 1657-1760
- Filter Points: line ~1200-1350

## Test Data
User is testing with:
- 91 occurrence rows
- 14 taxa
- Tree with 14 tips
