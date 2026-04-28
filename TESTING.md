# BioGeoBEARS Shiny Application - Testing and Validation Guide

## Testing Overview

This document provides comprehensive testing procedures for the BioGeoBEARS Shiny application, including manual testing checklists, automated testing setup, and validation procedures.

## Manual Testing Checklist

### Phase 1: Installation and Startup

- [ ] Install all dependencies without errors
- [ ] Install biogeoshiny package successfully
- [ ] Load library without errors: `library(biogeoshiny)`
- [ ] Run app without errors: `run_biogeoshiny()`
- [ ] App opens in web browser
- [ ] All tabs are visible and accessible
- [ ] Help panel displays correctly

### Phase 2: Data Input Module (Step 1)

#### CSV Upload
- [ ] Upload valid CSV with correct columns (species, longitude, latitude)
- [ ] Preview shows first 5 rows correctly
- [ ] Data loads successfully
- [ ] Success message appears
- [ ] Reject CSV with missing columns (show error)
- [ ] Reject empty CSV (show error)
- [ ] Reject malformed CSV (show error)

#### Tree File Upload
- [ ] Upload valid Newick tree file
- [ ] Tree loads successfully
- [ ] Number of taxa displayed correctly
- [ ] Success message appears
- [ ] Reject invalid tree format (show error)
- [ ] Reject empty tree file (show error)

#### Data Validation
- [ ] Coordinates are within valid range (-180 to 180 longitude, -90 to 90 latitude)
- [ ] Species names are preserved correctly
- [ ] No data loss during upload
- [ ] Large files (>10MB) handled appropriately

### Phase 3: Range Extrapolation Module (Step 2)

#### Buffer Method
- [ ] Select buffer method
- [ ] Set buffer width parameter
- [ ] Toggle "use mean distance" option
- [ ] Set grid resolution
- [ ] Run extrapolation without errors
- [ ] Presence-absence matrix generated
- [ ] Matrix dimensions correct (areas × species)
- [ ] Matrix values are binary (0 or 1)
- [ ] Progress indicator shows during computation

#### Convex Hull Method
- [ ] Select convex hull method
- [ ] Set grid resolution
- [ ] Run extrapolation without errors
- [ ] Presence-absence matrix generated
- [ ] Matrix quality reasonable

#### MST Method
- [ ] Select MST method
- [ ] Set grid resolution
- [ ] Run extrapolation without errors
- [ ] Presence-absence matrix generated
- [ ] Matrix quality reasonable

#### Error Handling
- [ ] Show error if no data loaded
- [ ] Show error for invalid parameters
- [ ] Show error for invalid grid resolution
- [ ] Graceful handling of edge cases

### Phase 4: Export Formats Module (Step 3)

#### BioGeoBEARS Export
- [ ] Select BioGeoBEARS format
- [ ] Run export without errors
- [ ] File generated correctly
- [ ] File format matches specification
- [ ] Download button appears
- [ ] Downloaded file opens correctly

#### NEXUS Export
- [ ] Select NEXUS format
- [ ] Run export without errors
- [ ] File generated correctly
- [ ] File format matches NEXUS specification
- [ ] Download button appears
- [ ] Downloaded file opens correctly

#### TNT Export
- [ ] Select TNT format
- [ ] Run export without errors
- [ ] File generated correctly
- [ ] File format matches TNT specification
- [ ] Download button appears
- [ ] Downloaded file opens correctly

#### NDM Export
- [ ] Select NDM format
- [ ] Set grid resolution
- [ ] Run export without errors
- [ ] File generated correctly
- [ ] File format matches XYD specification
- [ ] Download button appears
- [ ] Downloaded file opens correctly

#### Multiple Format Export
- [ ] Select multiple formats
- [ ] All exports complete successfully
- [ ] All download buttons appear
- [ ] All files download correctly

### Phase 5: BioGeoBEARS Setup Module (Step 4)

#### Model Selection
- [ ] Select DEC model
- [ ] Select DEC+J model
- [ ] Select DIVALIKE model
- [ ] Select DIVALIKE+J model
- [ ] Select BAYAREALIKE model
- [ ] Select BAYAREALIKE+J model
- [ ] Select multiple models
- [ ] Deselect models
- [ ] Error if no models selected

#### Parameter Configuration
- [ ] Set d (dispersal) bounds
- [ ] Set e (extinction) bounds
- [ ] Set j (jump dispersal) bounds
- [ ] Parameters within valid ranges
- [ ] Min < Max for all parameters
- [ ] Error if Min >= Max

#### Advanced Options
- [ ] Toggle time stratification
- [ ] Upload time stratification file
- [ ] Toggle distance matrix
- [ ] Upload distance matrix file
- [ ] Toggle dispersal multiplier
- [ ] Upload dispersal multiplier file
- [ ] All options work independently

#### Optimization Settings
- [ ] Select optimizer (GenSA, optim)
- [ ] Set maximum iterations
- [ ] Valid iteration range

### Phase 6: Analysis Module (Step 5)

#### Model Execution
- [ ] Run single model without errors
- [ ] Run multiple models without errors
- [ ] Progress indicator shows
- [ ] Results display after completion
- [ ] Success message appears

#### Results Display
- [ ] Model comparison table displays
- [ ] All columns present (Model, LnL, d, e, j, AICc)
- [ ] Values formatted correctly
- [ ] Table is sortable
- [ ] Table is scrollable for many models

#### Model Details
- [ ] Individual model details tab accessible
- [ ] All parameters displayed
- [ ] Values formatted correctly
- [ ] Details for all models present

#### Download Results
- [ ] Download button present
- [ ] Results CSV downloads
- [ ] CSV opens correctly
- [ ] Data matches displayed table

### Phase 7: User Interface

#### Navigation
- [ ] All tabs accessible
- [ ] Tab switching smooth
- [ ] No data loss when switching tabs
- [ ] Back button works (if applicable)

#### Help System
- [ ] Help panel visible on each tab
- [ ] Help text relevant to current step
- [ ] Help text clear and informative
- [ ] Help panel can be collapsed/expanded

#### Responsive Design
- [ ] App works on desktop (1920×1080)
- [ ] App works on laptop (1366×768)
- [ ] App works on tablet (768×1024)
- [ ] Text readable on all screen sizes
- [ ] Buttons clickable on all screen sizes

#### Error Messages
- [ ] Error messages clear and helpful
- [ ] Error messages appear in appropriate location
- [ ] Error messages don't prevent further use
- [ ] Warning messages distinct from errors

### Phase 8: Performance

#### Load Time
- [ ] App loads in < 5 seconds
- [ ] Tabs switch in < 1 second
- [ ] No lag when typing in inputs

#### File Upload
- [ ] Small files (< 1MB) upload instantly
- [ ] Medium files (1-10MB) upload in < 5 seconds
- [ ] Large files (10-100MB) upload with progress indicator

#### Computation
- [ ] Range extrapolation completes in reasonable time
- [ ] Export formats complete in < 10 seconds
- [ ] BioGeoBEARS setup completes instantly
- [ ] Analysis completes in reasonable time

#### Memory Usage
- [ ] App doesn't consume excessive RAM
- [ ] No memory leaks during extended use
- [ ] Multiple analyses don't cause crashes

## Automated Testing

### Unit Tests

```R
# tests/testthat/test_utils.R
test_that("calc_lrt calculates correctly", {
  # Test LRT calculation
})

test_that("calc_aicc calculates correctly", {
  # Test AICc calculation
})

test_that("validate_occurrence_data rejects invalid data", {
  # Test validation
})
```

### Integration Tests

```R
# tests/testthat/test_modules.R
test_that("mod_data_input loads CSV correctly", {
  # Test data loading
})

test_that("mod_range_extrapolation creates matrix", {
  # Test range extrapolation
})

test_that("mod_export_formats exports correctly", {
  # Test export functions
})
```

### Running Tests

```R
# Install testing packages
install.packages("testthat")

# Run all tests
devtools::test()

# Run specific test file
devtools::test_file("tests/testthat/test_utils.R")

# Run with coverage
devtools::test_coverage()
```

## Data Validation

### Occurrence Data Validation

```R
# Valid data
valid_data <- data.frame(
  species = c("Sp1", "Sp1", "Sp2"),
  longitude = c(-60.5, -58.7, -60.5),
  latitude = c(-3.2, -5.3, -3.2)
)

# Invalid: missing column
invalid_data1 <- data.frame(
  species = c("Sp1", "Sp1"),
  longitude = c(-60.5, -58.7)
)

# Invalid: NA values
invalid_data2 <- data.frame(
  species = c("Sp1", "Sp1", NA),
  longitude = c(-60.5, -58.7, -60.5),
  latitude = c(-3.2, -5.3, -3.2)
)

# Invalid: out of range coordinates
invalid_data3 <- data.frame(
  species = c("Sp1", "Sp1"),
  longitude = c(-60.5, 200),
  latitude = c(-3.2, -5.3)
)
```

### Tree File Validation

```R
# Valid Newick tree
valid_tree <- "(Sp1:1.0,Sp2:1.0):0.0;"

# Invalid: malformed
invalid_tree1 <- "(Sp1:1.0,Sp2:1.0"

# Invalid: empty
invalid_tree2 <- ""

# Invalid: wrong format
invalid_tree3 <- "Sp1,Sp2,Sp3"
```

### Presence-Absence Matrix Validation

```R
# Valid matrix
valid_matrix <- matrix(
  c(1,0,1,0, 0,1,1,1, 1,1,0,0),
  nrow = 3,
  ncol = 4,
  dimnames = list(
    c("Area1", "Area2", "Area3"),
    c("Sp1", "Sp2", "Sp3", "Sp4")
  )
)

# Check properties
all(valid_matrix %in% c(0, 1))  # Binary values
nrow(valid_matrix) > 0          # Has areas
ncol(valid_matrix) > 0          # Has species
```

## Regression Testing

### Test Cases to Run After Each Update

1. **Data Input**
   - Upload sample CSV
   - Upload sample tree
   - Verify data loads correctly

2. **Range Extrapolation**
   - Run buffer method
   - Run convex hull method
   - Run MST method
   - Verify matrix generation

3. **Export Formats**
   - Export to all formats
   - Verify file integrity
   - Verify downloads

4. **BioGeoBEARS Setup**
   - Configure all models
   - Set parameters
   - Verify configuration saved

5. **Analysis**
   - Run single model
   - Run multiple models
   - Verify results display

## Performance Benchmarks

### Target Performance Metrics

| Operation | Target Time | Acceptable Range |
|-----------|------------|-----------------|
| App startup | 3-5 seconds | < 10 seconds |
| Tab switch | < 1 second | < 2 seconds |
| CSV upload (1MB) | < 1 second | < 5 seconds |
| CSV upload (10MB) | 2-5 seconds | < 30 seconds |
| Range extrapolation | 5-30 seconds | < 60 seconds |
| Export formats | < 5 seconds | < 10 seconds |
| BioGeoBEARS analysis | 30-120 seconds | < 300 seconds |

## Browser Compatibility

### Tested Browsers

- [ ] Chrome (latest)
- [ ] Firefox (latest)
- [ ] Safari (latest)
- [ ] Edge (latest)

### Features to Test in Each Browser

- [ ] File upload works
- [ ] File download works
- [ ] Tables render correctly
- [ ] Buttons are clickable
- [ ] Text is readable
- [ ] No console errors

## Accessibility Testing

- [ ] Tab navigation works
- [ ] Keyboard shortcuts work
- [ ] Screen reader compatible
- [ ] Color contrast adequate
- [ ] Font sizes readable
- [ ] Error messages clear

## Security Testing

- [ ] File upload validation
- [ ] Input sanitization
- [ ] SQL injection prevention (if applicable)
- [ ] XSS prevention
- [ ] CSRF protection (if applicable)

## Deployment Testing

### Pre-Deployment Checklist

- [ ] All tests pass
- [ ] No console errors
- [ ] Performance acceptable
- [ ] Documentation complete
- [ ] Example data provided
- [ ] Installation guide tested
- [ ] Troubleshooting guide accurate

### Post-Deployment Verification

- [ ] App accessible from target URL
- [ ] File uploads work
- [ ] Downloads work
- [ ] All features functional
- [ ] Performance acceptable
- [ ] No error logs

## Test Report Template

```markdown
# Test Report - BioGeoBEARS Shiny v0.1.0

**Date**: [Date]
**Tester**: [Name]
**Environment**: [OS, R version, Browser]

## Summary
- Total Tests: X
- Passed: X
- Failed: X
- Skipped: X

## Issues Found
1. [Issue description]
   - Severity: [Critical/High/Medium/Low]
   - Steps to reproduce: [Steps]
   - Expected: [Expected behavior]
   - Actual: [Actual behavior]

## Recommendations
- [Recommendation 1]
- [Recommendation 2]

## Sign-off
- [Tester name]
- [Date]
```

## Continuous Integration

### GitHub Actions Workflow (Optional)

```yaml
name: Tests
on: [push, pull_request]
jobs:
  test:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v2
      - uses: r-lib/actions/setup-r@v2
      - name: Install dependencies
        run: Rscript -e 'devtools::install_deps(dependencies = TRUE)'
      - name: Run tests
        run: Rscript -e 'devtools::test()'
```

## Next Steps

1. **Implement Unit Tests**: Create test files in `tests/testthat/`
2. **Set Up CI/CD**: Configure GitHub Actions or similar
3. **Create Test Data**: Develop sample datasets for testing
4. **Document Issues**: Track bugs and improvements
5. **Performance Optimization**: Address any bottlenecks

## Contact

For testing questions or issues, contact the development team.

---

**Last Updated**: March 18, 2026
