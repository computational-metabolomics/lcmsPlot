# lcmsPlot

`lcmsPlot` is an R package designed for visualising
Liquid Chromatography-Mass Spectrometry (LC-MS) data with high-quality plots
suitable for publication. The package enables users to generate and
customise chromatograms, mass traces, spectra, and more with fine-tuned
aesthetics and annotation options.

## Overview

The motivation for developing a dedicated R package for LC–MS data visualization comes from the limitations of existing solutions such as base `xcms` plotting functions or generic plotting libraries.
While xcms and related tools are powerful for data processing and feature detection, their plotting capabilities are often designed primarily for *diagnostic purposes*, offering limited flexibility and aesthetics.

Researchers frequently need **publication-quality** chromatograms, spectra, and comparative visualisations that clearly highlight subtle differences across samples or conditions, which can be cumbersome to generate with the default methods.

The lcmsPlot package addresses this gap by providing:

- A consistent, intuitive interface for LC–MS visualizations.
- Streamlined generation of high-quality, customisable plots.
- Improved reproducibility and reduced need for bespoke scripts.
- Interoperability with the existing ecosystem (e.g., `MsExperiment`, `MSnbase`)
- Better performance at scale.

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
  chromatogram(features = rbind(c(
    mzmin = 300,
    mzmax = 305,
    rtmin = 2500,
    rtmax = 2550))) +
  rt_line(intercept = 2520) +
  arrange(group_by = "sample_group") +
  facets(facets = "sample_id", ncol = 3)
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
