# BioGeoBEARS Shiny Application - Installation Guide

## Quick Start

### Step 1: Install R and RStudio

- Download R from [https://cran.r-project.org/](https://cran.r-project.org/)
- Download RStudio from [https://posit.co/download/rstudio-desktop/](https://posit.co/download/rstudio-desktop/)

### Step 2: Install Required Packages

Open R or RStudio and run:

```R
# Install required dependencies
install.packages(c(
  "shiny",
  "shinydashboard",
  "shinyjs",
  "shinyWidgets",
  "shinythemes",
  "leaflet",
  "DT",
  "ape",
  "dplyr",
  "tidyr",
  "ggplot2",
  "plotly",
  "devtools"
))

# Install optional dependencies for full functionality
# (These are suggested but not required for basic operation)
install.packages(c(
  "BioGeoBEARS",
  "GenSA",
  "sf",
  "terra",
  "sp",
  "raster",
  "geosphere",
  "phytools",
  "markdown",
  "knitr"
))
```

### Step 3: Install the biogeoshiny Package

```R
# Install from local directory
devtools::install_local("path/to/biogeoshiny")

# Or install from GitHub (when available)
# devtools::install_github("your-username/biogeoshiny")
```

### Step 4: Run the Application

```R
library(biogeoshiny)
run_biogeoshiny()
```

The application will open in your default web browser.

## Detailed Installation

### For Linux Users

```bash
# Install system dependencies
sudo apt-get update
sudo apt-get install r-base r-base-dev
sudo apt-get install libssl-dev libcurl4-openssl-dev libxml2-dev

# Install RStudio (optional but recommended)
wget https://download1.rstudio.org/rstudio-xenial-1.4.1103-amd64.deb
sudo dpkg -i rstudio-xenial-1.4.1103-amd64.deb
```

### For macOS Users

```bash
# Using Homebrew
brew install r
brew install --cask rstudio

# Or download from official websites
```

### For Windows Users

1. Download R from [https://cran.r-project.org/bin/windows/base/](https://cran.r-project.org/bin/windows/base/)
2. Download RStudio from [https://posit.co/download/rstudio-desktop/](https://posit.co/download/rstudio-desktop/)
3. Run the installers and follow the prompts

## Troubleshooting

### Issue: "Package 'X' not found"

**Solution**: Install the missing package:
```R
install.packages("X")
```

### Issue: "Could not find shiny app directory"

**Solution**: Make sure the package is properly installed:
```R
# Reinstall the package
devtools::install_local("path/to/biogeoshiny", force = TRUE)

# Or check the installation
library(biogeoshiny)
```

### Issue: BioGeoBEARS installation fails

**Solution**: Install dependencies first:
```R
# Install dependencies
install.packages(c("GenSA", "FD", "snow", "parallel", "rexpokit", "cladoRcpp"))

# Then install BioGeoBEARS
install.packages("BioGeoBEARS")
```

### Issue: "Cannot open file" when uploading data

**Solution**: Make sure your file is in the correct format:
- Occurrence data: CSV with columns `species`, `longitude`, `latitude`
- Tree file: Newick or Nexus format

### Issue: Application runs slowly

**Solution**: 
- Close other applications to free up memory
- Use smaller datasets for testing
- Increase R memory allocation:
```R
memory.limit(size = 8000)  # Set to 8GB
```

## System Requirements

### Minimum Requirements
- **RAM**: 4 GB
- **Disk Space**: 2 GB
- **Processor**: Dual-core, 2 GHz

### Recommended Requirements
- **RAM**: 8 GB or more
- **Disk Space**: 5 GB or more
- **Processor**: Quad-core, 2.5 GHz or faster
- **Internet**: For downloading packages and data

## Performance Tips

1. **Update R and Packages**: Keep R and packages up to date
   ```R
   update.packages()
   ```

2. **Use RStudio**: Provides better memory management and performance

3. **Optimize Data**: 
   - Use smaller geographic areas
   - Reduce number of species
   - Increase grid resolution

4. **Parallel Processing**: Use multiple cores for BioGeoBEARS
   ```R
   # Set in BioGeoBEARS configuration
   num_cores <- parallel::detectCores() - 1
   ```

## Next Steps

After installation, refer to:
- **README.md** for user documentation
- **DEVELOPMENT.md** for development information
- In-app help for workflow guidance

## Getting Help

- Check the troubleshooting section above
- Review the in-app help documentation
- Consult BioGeoBEARS documentation: [http://phylo.wikidot.com/biogeobears](http://phylo.wikidot.com/biogeobears)
- Visit Wallace project: [https://wallaceecomod.github.io/wallace/](https://wallaceecomod.github.io/wallace/)

## Uninstallation

To remove the package:

```R
# Remove the package
remove.packages("biogeoshiny")

# Or reinstall if needed
devtools::install_local("path/to/biogeoshiny", force = TRUE)
```

## Version Information

- **Package Version**: 0.1.0
- **R Version Required**: >= 4.0.0
- **Last Updated**: 2025

## License

GPL (>= 3)
