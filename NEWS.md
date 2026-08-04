# lcmsPlot 1.1.8

- Closed six gaps between lcmsPlot and the plotting methods of *xcms*, so the
  `backend = "lcmsPlot"` variants of those methods have an equivalent to call.

- The *Easy-to-use, intuitive, and efficient LC-MS data plotting with lcmsPlot*
  vignette now covers all of this. A new *Quality control of peak detection*
  section documents the two new layers, and the four extended arguments are
  documented in the sections that already covered those layers.

- **New layer `lp_chrom_peak_rects()`** draws one rectangle per detected
  chromatographic peak, spanning `rtmin`-`rtmax` by `mzmin`-`mzmax`. It covers
  both xcms methods that draw that primitive. Called on its own it owns a panel,
  drawing the rectangles on an otherwise empty rt / m/z frame with one facet per
  sample, as `xcms::plotChromPeaks()` does; `rt_range` and `mz_range` restrict it
  the way that method's `xlim` and `ylim` do. Called after `lp_mass_trace()` or
  `lp_intensity_map()` it decorates that panel instead, as
  `xcms::plot(type = "XIC")` does, following the host's `x_dim`, window and
  `sample_ids` so the axes are not widened and no extra facets appear. Which mode
  applies is decided from the finished plot, so it does not depend on the order
  the layers were added in. Every box is given a minimum height so that peaks on
  m/z-binned data (where `mzmin == mzmax`) read as segments instead of vanishing.

- **New layer `lp_peak_count_image()`** reproduces
  `xcms::plotChromPeakImage()`: retention-time bins on x, samples on y, filled by
  the number of chromatographic peaks per bin. Samples are ordered by injection
  order, and bins with no peaks are kept at zero, so a sample that stopped
  producing peaks half-way through a run reads as an empty stretch rather than
  disappearing. Counts match the xcms method's binning exactly.

- **`lp_peak_density()` now accepts `XChromatograms` and `XChromatogram`
  objects**, which is what `xcms::plotChromPeakDensity()` takes. The type gate is
  now a capability check through `get_detected_peaks()` rather than a class
  check, the creator goes through that adapter instead of calling
  `xcms::chromPeaks()` directly (the sample column is named `column`, not
  `sample`, on chromatogram objects), and `features` is redundant for these
  inputs because each extracted ion chromatogram already carries its own m/z
  window. A new `simulate` argument mirrors the xcms method: `TRUE` descends the
  density curve for candidate features, `FALSE` draws the feature definitions the
  object already stores.

- **`lp_chromatogram()` gains `stacked` and `transform`**, matching
  `xcms::plotChromatogramsOverlay()`. `stacked` offsets each series up the y axis
  in m/z order by a fraction of the intensity range so co-eluting traces stop
  occluding each other, using the same offset formula as xcms, and suppresses the
  y axis because a stacked axis has no single meaning; the offsets are returned
  on the plot as a `stacked_offsets` attribute. `transform` (for example
  `log10`) is applied to the intensities and to the peak-highlight geometry
  alike, so shaded peaks stay attached to their traces.

- **`lp_chromatogram(aggregation_fun = "mean")`** produces the averaged ion
  chromatogram of `xcms::plotChrom(base = FALSE)`. The `ms_header()` contract has
  gained `peaksCount` and `meanIntensity`, letting each backend define the mean
  faithfully: `XcmsRawList` inputs carrying a profile matrix use the
  profile-matrix column mean and so reproduce the xcms trace exactly, while
  file-based inputs average over the measured peaks of each scan. Unknown
  aggregation functions now raise an error rather than silently returning a TIC.

- **`lp_intensity_map()` gains `geom`, `point_size`, `bin_rt`, `bin_mz` and
  `colour_scale`.** `geom = "point"` draws the individual centroids as
  `xcms::plotRaw()` does, skipping the binning so gaps in the mass traces stay
  visible instead of being implied away by a tile grid. The retention-time and
  m/z bin widths, previously hard-coded at 0.1, are now arguments.

## Breaking changes

- The `density` argument of `lp_intensity_map()` has been **removed**. It is
  replaced by `geom = "density"`, which selects the same rendering. Two separate
  switches could contradict each other, so there is now a single one. Replace
  `lp_intensity_map(..., density = TRUE)` with
  `lp_intensity_map(..., geom = "density")`; `density = FALSE` was the default
  and can simply be dropped.

## Bug fixes

- `lp_intensity_map()` no longer errors with *"argument is of length zero"* on
  narrow m/z windows. The peak matrix was subset without `drop = FALSE`, so any
  scan with exactly one peak in range collapsed to a vector and the following
  `nrow()` check failed. A narrow window is precisely what an XIC plot passes.

- `plot_chromatogram()` no longer joins every series in a facet into a single
  zig-zagging line when the colour aesthetic is something other than
  `sample_id`. The ggplot `group` was pinned to `sample_id` while `lp_arrange()`
  set only the colour; it is now a composite of the sample and the arrangement
  column, so colouring by `feature_id` splits the lines while coarser factors
  such as `sample_group` keep one line per sample as before.

- `lp_peak_density(simulate = TRUE)` no longer draws overlapping feature
  rectangles that all start at the same retention time. The internal descent
  used to find each feature's boundaries compared neighbouring density values
  with `<=`, so once the previous feature's range had been zeroed it walked
  straight across that plateau and re-collected the peaks already assigned. Each
  later rectangle therefore reached back to the first feature's leftmost peak,
  and passed the `min_fraction` test on the strength of those peaks. The
  comparison is now strict, matching the `DescendMin` routine of *xcms*, and the
  rectangles agree with `xcms::plotChromPeakDensity()`.

# lcmsPlot 1.1.7

- Every man page documenting an exported object now carries a runnable example.
  The four `\dontrun{}` blocks that referenced unavailable vendor data
  (`CompoundDiscovererNodeSource()`, `LipidSearchSource()`,
  `lp_compound_discoverer()`, `lp_lipid_search()`) have been replaced with
  self-contained code, and the four precursor-purity layers
  (`lp_purity_overlay()`, `lp_purity_timeline()`, `lp_purity_distribution()`,
  `lp_isolation_window()`) have gained examples where they previously had none.
  This satisfies the Bioconductor requirement that at least 80% of such pages
  have runnable examples; coverage is now 100% (33 of 33, up from 25 of 33).

- The purity examples use the DDA files already shipped in
  `inst/extdata/standards-mzml.zip`, whose MS2 scans carry real precursor
  selection and isolation-window metadata, so `msPurity::purityA()` yields six
  fragmentation events with informative `inPurity` scores. They are guarded with
  `@examplesIf requireNamespace("msPurity")` because `msPurity` is a suggested
  dependency. `lp_isolation_window()` uses `half_width = 0.6` to match the
  isolation width recorded in those files.

- The Compound Discoverer scripting-node and LipidSearch examples build a
  minimal vendor export in `tempdir()`, following the existing
  `MZmineFeatureListsSource()` and `MsDialPeaksSource()` examples, and point
  `sample_paths` at the shipped mzML files so the documented plots extract real
  chromatograms rather than only constructing a data source.

# lcmsPlot 1.1.6

- Added support for LipidSearch result files (both **4.2** and **5.2**; the
  version is auto-detected). The new `LipidSearchSource()` constructor parses the
  lipid table, reshapes the per-sample `Area` / `Height` / RT / observed-m/z
  columns into an xcms-style peak table, and builds a plottable data source.
  Rejected lipids are dropped unless `keep_rejected = TRUE`, and lipid ions
  reported at several retention times get distinct plot labels. Plotting uses the
  new `lp_lipid_search(lipids_query, rt_extend)` layer, whose query can reference
  the lipid annotations (`class`, `sub_class`, `grade`, `adduct`, `lipid_rank`,
  ...).

- LipidSearch **5.2** files carry no raw-file names, so their samples are keyed
  `s1`, `s2`, … (from the `OrgMeanArea[...]` columns). For 5.2, map raw files by
  naming `sample_paths` with those keys (`c("s1" = "a.mzML", "s2" = "b.mzML")`) —
  order-independent and covering any subset — or pass an unnamed vector matched
  positionally. 5.2 exports add a `sub_class` annotation and an explicit adduct,
  and their per-lipid `BaseRt` provides an extraction window even for lipids
  detected in no sample. 4.2's basename matching is unchanged. The `rej`
  annotation is now a `logical` for both versions.

- Samples plotted from a LipidSearch source are **not** limited to the ones
  declared in the result file: `sample_paths` determines the sample list and is
  matched to the declarations by file basename without extension (so a result
  file listing `.raw` files works with converted `.mzML` files). A supplied
  path that matches no declaration is still a full sample - every queried
  lipid's chromatogram is extracted there using the lipid's consensus m/z and
  retention-time window, it simply has no reported peak to highlight.

- The per-compound chromatogram extraction loop is now shared between the
  Compound Discoverer scripting-node and LipidSearch sources
  (`create_compound_chromatograms()`), along with the compound ranking and
  column-resolution helpers.

# lcmsPlot 1.1.5

- `lp_total_ion_current()` now works when raw files are passed directly to
  `lcmsPlot()` (a `character` vector of `.mzML`, `.mzXML`, `.CDF`, or `.raw`
  paths), in addition to `XCMSnExp` / `MsExperiment` objects. The per-scan TIC
  is read from the raw file headers via the shared raw-file reader interface.

- Added a runnable worked example for large, multi-sample LC-MS studies. Using
  a 50-sample mzML study (the published Sacurine dataset, MetaboLights
  MTBLS404), it demonstrates the batching API (`batch_size` +
  `iterate_plot_batches()`) together with `patchwork` to build one
  self-contained composite figure per sample (TIC and BPC, an EIC + mass trace
  for a shared compound, and a zoomed intensity map), written to a multi-page
  PDF report. It is distributed as a standalone script rather than a packaged
  vignette, since it downloads and processes ~50 raw files.

# lcmsPlot 1.1.4

- Added support for using `lcmsPlot` inside a Compound Discoverer custom
  Scripting Node. The new `CompoundDiscovererNodeSource()` constructor reads
  the `node_args.json` file passed by Compound Discoverer together with the
  tab-delimited table exports it references (Compounds, Compounds per File,
  Features, and their link tables), selects a representative ion per compound
  and study file (preferring `[M+H]+1` / `[M-H]-1`, otherwise the most
  abundant feature), and builds a plottable data source. Chromatograms are
  extracted from the raw (mzML) files supplied via `sample_paths`, which are
  matched to study files by file basename (falling back to positional
  matching by study-file ID). Plotting reuses the existing
  `lp_compound_discoverer(compounds_query, rt_extend)` interface. `jsonlite`
  was added to `Imports`.

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
