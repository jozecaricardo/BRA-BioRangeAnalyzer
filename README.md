# BRA: BioRangeAnalyzer

![Logotipo do LEBIP](inst/images/logo_lebip_readme.png)

BRA is an **R Shiny application** designed to support biogeographic workflows from occurrence-data preparation to ancestral-range estimation and export. The application integrates data checking, range extrapolation, BioGeoBEARS analyses, map-based visualisation, and file export in a single interface, with an emphasis on usability for teaching and research in historical biogeography.

The project page is available at [biogeografiaAULAS](https://jozecaricardo.github.io/biogeografiaAULAS/).

BRA was developed in the context of the course **"Métodos em Biogeografia Histórica"** and is associated with collaborations involving the **Laboratório de Estudos da Biodiversidade Pampiana (LEBIP)** and the **Laboratório de Entomologia, Sistemática e Biogeografia (LESB/FURG)**. Public institutional information indicates that LESB/FURG was created in 2014 by **Prof. Augusto Ferrari**, within the **Instituto de Ciências Biológicas** of the **Universidade Federal do Rio Grande (FURG)**, as described on the public institutional page of the laboratory.

## Institutional identity

| Laboratory | Logo |
| --- | --- |
| **LEBIP** — Laboratório de Estudos da Biodiversidade Pampiana | ![Logotipo do LEBIP](inst/images/logo_lebip_readme.png) |
| **LESB/FURG** — Laboratório de Entomologia, Sistemática e Biogeografia | ![Logotipo do LESB/FURG](inst/images/logo_lesb_readme.png) |

## Authors

| Name | Affiliation |
| --- | --- |
| **José Ricardo Inacio Ribeiro** | Universidade Federal do Pampa (UNIPAMPA), campus São Gabriel, Rio Grande do Sul, Brazil |
| **Augusto Ferrari** | Universidade Federal do Rio Grande (FURG), Rio Grande, Rio Grande do Sul, Brazil |

## Main features

BRA provides an integrated interface for the main analytical stages commonly required in historical biogeography. In practical terms, the application supports occurrence cleaning and validation, phylogenetic input, alternative range-extrapolation routines, BioGeoBEARS model fitting, visual inspection of ancestral-range results, and export of downstream files for complementary analyses.

| Component | Scope |
| --- | --- |
| **Data input and validation** | Import of occurrence tables (`.csv` and `.txt`) and phylogenetic trees (`.nwk`, `.newick`, `.tre`, `.txt`, `.nex`, `.nexus`) |
| **Range extrapolation** | Buffers, convex hulls, minimum spanning tree (MST), and irregular-bin workflows |
| **BioGeoBEARS analyses** | DEC, DEC+J, DIVALIKE, DIVALIKE+J, BAYAREALIKE, and BAYAREALIKE+J |
| **Outputs and exports** | Export to BioGeoBEARS, TNT, NEXUS, and NDM-compatible formats, together with reports and visual outputs |
| **Visual interpretation** | Ancestral-range plots, uncertainty plots, regular-grid maps, and extrapolation layers |

## Installation

BRA requires **R 4.0.0 or later**. Use of **RStudio** is recommended for installation, package loading, and local execution.

### Install from GitHub

**Option 1: Basic installation**

```r
devtools::install_github("jozecaricardo/BRA-BioRangeAnalyzer")
```

**Option 2: Installation with all dependencies (recommended)**

If you encounter issues with missing packages (especially `BioGeoBEARS`, `raster`, `terra`, or `sp`), use the `dependencies = TRUE` flag to automatically install all required dependencies:

```r
devtools::install_github("jozecaricardo/BRA-BioRangeAnalyzer", dependencies = TRUE)
```

**Note on BioGeoBEARS:** If `BioGeoBEARS` installation fails during compilation, you may need to install it separately:

```r
install.packages('devtools', repos='https://cloud.r-project.org')
devtools::install_github('nmatzke/BioGeoBEARS', INSTALL_opts='--byte-compile', upgrade='never')
```

## Run the app

### Recommended method (to avoid compatibility issues)

**⚠️ Important:** Due to namespace conflicts between `raster` and `terra` packages, the most reliable way to run BRA is to download and execute it directly without installing as a package. This method avoids potential compatibility issues:

```r
# Download the repository
download.file(
  "https://github.com/jozecaricardo/BRA-BioRangeAnalyzer/archive/refs/heads/main.zip",
  "~/BRA-BioRangeAnalyzer-download.zip",
  mode = "wb"
)

# Extract
unzip("~/BRA-BioRangeAnalyzer-download.zip", exdir = "~")

# Run the app
shiny::runApp("~/BRA-BioRangeAnalyzer-main/inst/shiny")
```

### Alternative method (if you prefer package installation)

After installation, start the application with:

```r
library(biogeoshiny)
biogeoshiny::run_biogeoshiny()
```

**Note:** This method may encounter namespace conflicts on some systems. If you experience errors related to `as.data.frame` or raster objects, use the recommended method above.

## Input formats

BRA expects a simple occurrence table with species names and geographic coordinates. Column names are interpreted case-insensitively.

### Occurrence file

| Required column | Meaning |
| --- | --- |
| `spp` | Species name |
| `long` | Longitude |
| `lat` | Latitude |

Example:

```text
spp,long,lat
Species1,-50.5,-25.3
Species1,-51.2,-26.1
Species2,-48.3,-24.5
```

### Tree file

| Accepted format | Extensions |
| --- | --- |
| **Newick** | `.nwk`, `.newick`, `.tre`, `.txt` |
| **NEXUS** | `.nex`, `.nexus` |

## How to cite the app

If you use **BRA: BioRangeAnalyzer** in teaching material, course activities, reports, theses, dissertations, manuscripts, or other scientific products, it is advisable to cite the application explicitly as software. Because the repository does not yet expose a formal software DOI, the most transparent approach at the moment is to cite the **software title, version, authorship, platform, year, and repository URL**.

A suggested citation format is:

> Ribeiro, J. R. I., & Ferrari, A. (2026). **BRA: BioRangeAnalyzer** (Version 0.1.0) [R Shiny application]. GitHub. https://github.com/jozecaricardo/BRA-BioRangeAnalyzer

If you prefer a shorter in-text form, you may cite it as **BRA: BioRangeAnalyzer v0.1.0**.

## How to cite BioGeoBEARS

BRA integrates workflows based on **BioGeoBEARS**, and the public BioGeoBEARS repository explicitly recommends citing the package when it is used. For that reason, whenever BRA is used to run ancestral-range estimation through BioGeoBEARS, the software underlying that analysis should also be cited.

A practical citation form is:

> Matzke, N. J. (2013). **BioGeoBEARS: BioGeography with Bayesian (and likelihood) Evolutionary Analysis in R scripts** [R package/software].

## Project context

The present repository documents an application intended to make biogeographic analyses more accessible in teaching and research settings. Its design combines user-guided workflows with reproducible exports, so that users can move from raw occurrence data to ancestral-range interpretation and external analytical pipelines with fewer manual steps.

## License

BRA is distributed under the **GPL (>= 3)** license.
