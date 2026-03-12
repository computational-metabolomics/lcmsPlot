# lcmsPlot 0.99.18

- Added support for `xcmsRaw` objects via the new `XcmsRawList` S4 container
  and `create_xcms_raw_list()` convenience helper that reads files in parallel
  via `BiocParallel`.
- Added `XcmsRawReader`, a `MsRawReader` subclass that reads scan headers and
  peaks directly from in-memory `xcmsRaw` slots without opening a file
  connection.
- Added support for `XChromatograms`, `XChromatogram`, and `MChromatograms`
  as direct data inputs to `lcmsPlot()`.
- `XChromatograms` now automatically derive a `feature_id` from the m/z and
  RT ranges of each row, enabling faceting and gridding on `feature_id`.
- Refactored the chromatogram creator into a unified `create_chromatograms()`
  S4 generic replacing the four previous specialised generics.
- Added `na.rm` argument to `lp_chromatogram()` to remove data points with
  `NA` intensity before plotting.
- Added GitHub Actions CI workflow (`.github/workflows/R-CMD-check.yaml`) that
  runs `devtools::test()` and `rcmdcheck::rcmdcheck()` against the Bioconductor
  devel Docker image on every push and pull request to `devel`.

# lcmsPlot 0.99.16

- Fixed `<-` vs `=` assignment NOTE raised by `R CMD CHECK`.

# lcmsPlot 0.99.15

- Addressed second round of Bioconductor pre-acceptance review comments.

# lcmsPlot 0.99.14

- Added titles to all vignette code chunks for improved readability.

# lcmsPlot 0.99.13

- Added runnable `\dontrun{}` examples to `MZmineFeatureListsSource` and
  `MsDialPeaksSource` documentation.

# lcmsPlot 0.99.12

- Added support for MZmine (v2+) feature lists via `MZmineFeatureListsSource`.
- Added support for MS-DIAL peak tables via `MsDialPeaksSource`.
- Updated vignette with MZmine and MS-DIAL interoperability sections.

# lcmsPlot 0.99.11

- Version bump.

# lcmsPlot 0.99.10

- Added support for Compound Discoverer results files (`.cdResult`) via a new
  SQL-based data source backed by `DBI` / `RSQLite`.
- Added `lp_compound_discoverer()` for loading and querying Compound Discoverer
  results directly within the `lcmsPlot` pipeline.
- Updated `lcmsPlot` class documentation.

# lcmsPlot 0.99.9

- Removed redundant example from `lcmsPlotClass` documentation.

# lcmsPlot 0.99.8

- Major codebase refactoring for Bioconductor standards compliance.
- `get_metadata()` for `XCMSnExp` and `MsExperiment` now falls back to
  deriving `sample_id` from file basenames when no `sample_id_column` is
  present.

# lcmsPlot 0.99.7

- Added figure dimension options (`fig.height`, `fig.width`) to vignette
  chunks for consistent rendered output.

# lcmsPlot 0.99.6

- Version bump; standardised `@return` documentation across all exported
  functions.

# lcmsPlot 0.99.5

- Added `@return` documentation to exported functions.
- Fixed xcms namespace qualification in the test helper.

# lcmsPlot 0.99.4

- Fixed use of `AnnotatedDataFrame` (was incorrectly using
  `NAnnotatedDataFrame` in some code paths).

# lcmsPlot 0.99.3

- Fixed namespace qualification for `xcms` and `patchwork` calls.

# lcmsPlot 0.99.2

- Fixed `MulticoreParam` call in vignette.

# lcmsPlot 0.99.1

- Fixed `BiocParallel` namespace registration in vignette.

# lcmsPlot 0.99.0

- Initial Bioconductor submission.
- Added comprehensive roxygen2 documentation with examples for all exported
  functions.
- Added `testthat` (edition 3) test suite covering chromatogram, spectra, and
  intensity map creators.
- Updated README with package overview and usage example.

# lcmsPlot 0.1.0

- Initial implementation of `lcmsPlot`.
- `lcmsPlot()` entry point accepting `MsExperiment`, `XCMSnExp`, and raw file
  paths (`character`) as data sources.
- `lp_chromatogram()`: base peak chromatograms (BPC), total ion chromatograms
  (TIC), and extracted ion chromatograms (XIC) with configurable m/z and RT
  tolerances.
- `lp_spectra()`: mass spectra extraction with `"closest"`, `"closest_apex"`,
  and `"across_peak"` scan-selection modes; standalone and chromatogram-linked
  display; spectral mirror plots.
- `lp_intensity_map()`: two-dimensional m/z vs RT intensity maps with optional
  density smoothing.
- `lp_total_ion_current()`: TIC distributions as violin or box plots.
- `lp_mass_trace()`: mass trace overlays below chromatogram panels.
- `lp_rt_line()`: vertical reference lines at specified retention times.
- `lp_rt_diff_plot()`: retention-time correction diagnostic plot.
- `lp_facets()`, `lp_grid()`, `lp_arrange()`, `lp_layout()`: flexible panel
  layout helpers.
- `lp_labels()`, `lp_legend()`: axis and legend annotation helpers.
- `highlight_peaks` and `highlight_apices` options in `lp_chromatogram()`.
- `MsExperiment` support including retention-time corrected chromatograms
  (`rt_type = "corrected"` / `"both"`).
