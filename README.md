# lcmsPlot

[![R-CMD-check](https://github.com/computational-metabolomics/lcmsPlot/actions/workflows/R-CMD-check.yaml/badge.svg?branch=devel)](https://github.com/computational-metabolomics/lcmsPlot/actions/workflows/R-CMD-check.yaml)
[![Bioconductor devel build](https://bioconductor.org/shields/build/devel/bioc/lcmsPlot.svg)](https://bioconductor.org/checkResults/devel/bioc-LATEST/lcmsPlot/)
[![License: GPL-3](https://img.shields.io/badge/license-GPL--3-blue.svg)](LICENSE)

[lcmsPlot](https://bioconductor.org/packages/release/bioc/html/lcmsPlot.html) is an R/Bioconductor package designed to generate composable, reproducible, interoperable, and high-quality publication-ready visualisations of Liquid Chromatography-Mass Spectrometry (LC-MS) data and processing results. The package enables users to generate and
customise chromatograms, mass traces, spectra, and more with fine-tuned
aesthetics and annotation options.
By facilitating clear, reproducible, and customisable figures, _lcmsPlot_ enhances interpretation and communication within LC-MS-based metabolomics analyses.

## Overview

The motivation for developing a dedicated R package for LC–MS data visualization comes from the limitations of existing solutions such as base `xcms` plotting functions or generic plotting libraries.
While xcms and related tools are powerful for data processing and feature detection, their plotting capabilities are often designed primarily for *diagnostic purposes*, offering limited flexibility and aesthetics.

Researchers frequently need **publication-quality** chromatograms, spectra, and comparative visualisations that clearly highlight subtle differences across samples or conditions, which can be cumbersome to generate with the default methods.

The lcmsPlot package addresses this gap by providing:

- A consistent, intuitive, ggplot2-like interface for LC–MS visualizations,
  where plot layers are combined with `+`.
- Streamlined generation of high-quality, customisable plots: chromatograms
  (BPC, TIC, XIC), mass traces, spectra and mirror plots, intensity maps, peak
  density, retention-time alignment, and precursor purity.
- Improved reproducibility and reduced need for bespoke scripts: `lp_save()`
  records the provenance of a figure (software versions, pipeline, input
  files) next to it, and `lp_export_data()` writes the plotted data as
  CSV with a data dictionary.
- Interoperability with raw files, Bioconductor objects, and third-party LC-MS
  software (see [Interoperability](#interoperability)).
- Better performance at scale, with parallel extraction (`BPPARAM`) and
  batched plotting of large studies (`batch_size` + `iterate_plot_batches()`).
- An interactive Shiny app for exploring data without writing code
  (see [Interactive exploration](#interactive-exploration)).

### Example

Suppose you want to plot chromatograms for a specified m/z and retention time window across all samples, grouped by sample class, with an indicator line at a specific retention time.

With `xcms` base plotting, this requires several manual steps: extracting chromatograms, arranging plots by sample, handling group colors, and adding annotations:

```r
# Extract chromatograms
chr <- chromatogram(raw_data, mz = c(300, 305), rt = c(2500, 2550))

# Colors for groups
group_cols <- c("KO" = "red", "WT" = "blue")

# Plot setup: one panel per sample
n_samples <- length(chr)
n_col <- ceiling(sqrt(n_samples))
n_row <- ceiling(n_samples / n_col)

par(mfrow = c(n_row, n_col), mar = c(4, 4, 2, 1))

for (i in seq_along(chr)) {
  chri <- chr[[i]]
  df <- as.data.frame(chri)
  group <- pData(raw_data)$sample_group[i]
  
  # Plot
  plot(df$rtime, df$intensity,
       type = "l", col = group_cols[group],
       xlab = "Retention Time (s)", ylab = "Intensity",
       main = paste0(names(chr)[i], " (", group, ")"))
  
  # Add vertical line
  abline(v = 2520, lty = 2)
}
```

With `lcmsPlot`, the same task becomes much simpler and more reproducible:

```r
library(lcmsPlot)

lcmsPlot(raw_data, sample_id_column = "sample_name") +
  lp_chromatogram(features = rbind(c(
    mzmin = 300,
    mzmax = 305,
    rtmin = 2500,
    rtmax = 2550))) +
  lp_rt_line(intercept = 2520) +
  lp_arrange(group_by = "sample_group") +
  lp_facets(facets = "sample_id", ncol = 3)
```

## Installation

To install this package:

```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("lcmsPlot")
```

To install the development version:

```r
if (!require("remotes", quietly = TRUE))
    install.packages("remotes")

remotes::install_github("computational-metabolomics/lcmsPlot")
```

Optional features rely on suggested packages: `rawrr` for ThermoFisher `.raw`
files, `msPurity` for precursor purity plots, and `shiny`/`shinytoastr` for
the interactive app.

## Interoperability

lcmsPlot plots directly from the data you already have: raw files,
Bioconductor objects, and the results of commonly used LC-MS software. The
same `lp_*()` layers work for every input, so a figure can be reproduced
whichever tool produced the data.

### Supported inputs

| Source | Input | Versions | How to plot it |
|---|---|---|---|
| **Raw files** | mzML, mzXML, netCDF (via `mzR`) | Those read by mzR | `lcmsPlot(c("a.mzML", "b.mzML"), metadata = ...)` |
| | ThermoFisher `.raw` (via `rawrr`) | Those read by the Thermo RawFileReader | `lcmsPlot(c("a.raw", "b.raw"))` |
| **Bioconductor** | `MsExperiment` / `XcmsExperiment`, `XCMSnExp` (*xcms*) | xcms 4 and later | `lcmsPlot(xdata)` |
| | `xcmsRaw` (legacy *xcms*) | Any | `lcmsPlot(create_xcms_raw_list(files))` |
| | `XChromatograms`, `XChromatogram`, `MChromatograms` | xcms 4 and later | `lcmsPlot(chrs)` |
| | `purityA` (*msPurity*) | Any | `lcmsPlot(pa)` + `lp_purity_*()` / `lp_isolation_window()` |
| **Compound Discoverer** | `.cdResult` files | 3.x | `lcmsPlot("study.cdResult") + lp_compound_discoverer()` |
| | Scripting-node exports | 3.3, 3.4, 3.5 | `CompoundDiscovererNodeSource()` |
| **MZmine** | Feature lists | 2 and 3 or later | `MZmineFeatureListsSource()` |
| **MS-DIAL** | Peak tables | Any with the standard peak list columns | `MsDialPeaksSource()` |
| **LipidSearch** | Lipid result files | 4.2 and 5.2 | `LipidSearchSource()` + `lp_lipid_search()` |
| **mzTab-M** | Exports from MS-DIAL, MZmine, OpenMS, ... | 2.0.0-M | `MzTabMSource()` + `lp_mztab()` |
| **Spectral libraries** | MSP (via `MsBackendMsp`) or mzML reference spectra | Any | `lp_spectra(spectral_match_db = "library.msp")` for mirror plots |

See `?lcmsPlot-formats` for the columns each reader expects and the packages
each input needs.

Results from external software are paired with the raw files they were
computed from (`sample_paths`), so chromatograms are extracted from the
original data and the peaks reported by the software can be highlighted on
top of them. For example, plotting the six most abundant lipids from a
LipidSearch result:

```r
lipids <- LipidSearchSource(
    results_path = "path/to/LIPIDS_POS_5ppm.txt",
    sample_paths = c("path/to/sample_1.mzML", "path/to/sample_2.mzML"))

lcmsPlot(lipids) +
    lp_lipid_search(lipids_query = "lipid_rank <= 6") +
    lp_chromatogram(highlight_peaks = TRUE, ppm = 5) +
    lp_grid(rows = "sample_id", cols = "name", free_x = TRUE)
```

or only the compounds checked by hand in Compound Discoverer:

```r
lcmsPlot("study.cdResult") +
    lp_compound_discoverer(compounds_query = "checked") +
    lp_chromatogram(highlight_peaks = TRUE) +
    lp_grid(rows = "sample_id", cols = "name", free_x = TRUE)
```

### Compound Discoverer scripting node

lcmsPlot ships a ready-made Compound Discoverer post-processing node,
**Compounds Chromatograms**. It renders one chromatogram per compound from the
original `.raw` files and adds a clickable plot column to the Compounds table,
so each compound links to its chromatogram inside Compound Discoverer. The
node is installed with the package:

```r
system.file("compound-discoverer-nodes/CompoundsChromatograms",
            package = "lcmsPlot")
```

See its [README](inst/compound-discoverer-nodes/CompoundsChromatograms/README.md)
for installation and parameters.

## Interactive exploration

`lcmsPlotApp()` launches a Shiny dashboard for exploring data without writing
code: BPC/TIC, extracted ion chromatograms, peak density and spectra, with
options for faceting, grouping, legends and titles. It accepts mzML/CDF and
ThermoFisher `.raw` files, Compound Discoverer `.cdResult` files, and
`XCMSnExp`/`MsExperiment` objects saved as `.rds` or `.RData`.

```r
shiny::runApp(lcmsPlotApp())
```

The app needs the suggested packages `shiny` and `shinytoastr`.

## Resources

- Package-level help page: `?\`lcmsPlot-package\``
- Functions/classes help pages: `?<function>`
- Vignettes: `browseVignettes("lcmsPlot")`
- Changes between versions: [NEWS.md](NEWS.md)

## Contributing

Contributions are welcome. Please read [CONTRIBUTING.md](CONTRIBUTING.md) for
the branch naming, commit message, and pull request conventions, and report
bugs on the [issue tracker](https://github.com/computational-metabolomics/lcmsPlot/issues).
