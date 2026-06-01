# lcmsPlot 1.1.3

- Redesigned `lcmsPlotApp()` as a Tailwind-based dashboard: sticky topbar,
  left navigation rail with a card-style file uploader and sample picker,
  and a content area with the plot in a card. The Tailwind CSS is
  precompiled and shipped under `inst/www/lcmsPlot.css`; no runtime network
  or Node dependency. Dev sources live in `tools/tailwind/` (excluded from
  the installed package).
- Added a per-tab Options panel that exposes `lp_facets()`, `lp_arrange()`,
  `lp_legend()`, and `lp_labels()` directly from the UI. Users can now
  facet, group, reposition the legend, and add titles without writing R
  code.
- Dashboard CSS is now served via `htmltools::htmlDependency()` and the
  Tailwind content scan + safelist were hardened so utility classes
  referenced from R templates actually end up in the compiled stylesheet.
  The UI is now wrapped in `shiny::bootstrapPage()` so the hidden
  `tabsetPanel` switches correctly. Card corners reduced from
  `rounded-2xl` to `rounded-md`, sidebar gap bumped to `gap-6`, and the
  initial nav highlight is now set statically in markup.
- The Shiny uploader now accepts additional data-source types: Thermo
  `.raw` files (multi-file), Compound Discoverer `.cdResult` SQLite
  results (single file), serialised `XCMSnExp` / `MsExperiment` objects
  via `.rds` (single file), and saved workspaces containing one such
  object via `.RData` / `.rda` (single file). All existing mzML / CDF
  workflows are unchanged.

# lcmsPlot 1.1.2

- Added `lcmsPlotApp()`, an interactive Shiny app for exploring raw LC-MS
  files. Users can upload one or more mzML / CDF files and view base-peak
  and total-ion chromatograms, extract ion chromatograms by m/z and ppm
  (with an optional retention-time window), inspect peak density across
  samples, and plot spectra by scan index. The app follows Bioconductor's
  Shiny code-organisation guidelines: all UI/server code lives in `R/`
  and `lcmsPlotApp()` returns a `shinyApp` object rather than calling
  `runApp()` internally. `shiny`, `shinytoastr`, and `shinytest2` were
  added to `Suggests`.

# lcmsPlot 1.1.1

- Added support for `purityA` objects (msPurity) as the `data_obj` input to
  `lcmsPlot()`, enabling direct visualisation of precursor ion purity results.
- Added `lp_purity_overlay()`: overlays per-scan `inPurity` scores as coloured
  points on a chromatogram, with a diverging colour scale and an optional
  threshold midpoint.
- Added `lp_purity_timeline()`: scatter plot of `inPurity` versus retention
  time per sample, with an optional horizontal threshold line.
- Added `lp_purity_distribution()`: violin or boxplot distribution of
  `inPurity` scores grouped by sample, with an optional threshold line.
- Added `lp_isolation_window()`: MS1 spectrum plot annotated with the
  isolation window rectangle, precursor m/z dashed line, and `inPurity` score
  as a subtitle.
- Added `purity_scores` slot to `lcmsPlotDataContainer` for storing per-scan
  purity data extracted from `purityA@puritydf`.
- Added `msPurity` to `Suggests` in DESCRIPTION.

# lcmsPlot 0.99.20

- Added `lp_peak_density()` for peak density plots that mirror
  `xcms::plotChromPeakDensity()`: y-axis shows sample indices positioned
  within the kernel density range, x-axis shows RT, coloured points mark
  individual detected peaks per sample, and a density line is overlaid.
  When `min_fraction` and `sample_groups` are provided the density-descent
  grouping algorithm is simulated and semi-transparent rectangles highlight
  RT regions that would form feature groups.
  When used after `lp_chromatogram()`, `features` is inherited automatically.
  Supports `rt_unit` (`"second"` / `"minute"`), `bw`, `min_samples`,
  `max_features`, and multiple m/z windows with auto-faceting.
- Fixed `highlight_peaks_mode = "rectangle"` and `"point"` in
  `lp_chromatogram()`: geoms now use `inherit.aes = FALSE` to avoid
  evaluating the global `x = rt_plot` aesthetic against `detected_peaks`,
  which previously caused an `object 'rt_plot' not found` error.
  RT values for rectangle and point modes are now scaled correctly when
  `rt_unit = "minute"`.
- Added `line_type` parameter to `lp_chromatogram()` (passed through to
  `geom_line()`; any ggplot2 linetype string is accepted).
- Added `x_dim`, `y_dim`, and `fill_scale` parameters to
  `lp_intensity_map()`, allowing the m/z and RT axes to be swapped and
  the fill colour scale to be replaced with any ggplot2 scale object.
- Expanded the "Plot peak density" section of the `lcms_data_plotting`
  vignette with explanatory prose and a combined chromatogram + peak
  density example.

# lcmsPlot 0.99.19

- Migrated all internal data structures from `data.frame` to `tibble`;
  `tibble` is now a formal `Imports` dependency.
- Fixed `all_of()` call in the `DBIConnection` chromatogram creator to pass
  column names as a single character vector.
- Qualified `tibble()` calls in `MZmineFeatureListsSource` and
  `MsDialPeaksSource` examples as `tibble::tibble()` to avoid
  `could not find function "tibble"` errors during `R CMD check`.
- Expanded `get_metadata()` documentation with dedicated sections describing
  the behaviour for all nine dispatch methods: `character`, `XCMSnExp`,
  `MsExperiment`, `MChromatograms`, `XChromatograms`, `XChromatogram`,
  `XcmsRawList`, `ExternalDataSource`, and `DBIConnection`.

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
