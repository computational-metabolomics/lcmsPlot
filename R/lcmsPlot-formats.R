#' Supported input formats
#'
#' lcmsPlot reads raw LC-MS data, Bioconductor objects, and the results of
#' commonly used LC-MS software. The table lists every input, the versions it
#' supports, the way to read it, and the packages it needs beyond the
#' imports of lcmsPlot.
#'
#' | Source | Input | Versions | Read with | Needs |
#' |---|---|---|---|---|
#' | Raw files | mzML, mzXML, netCDF (ANDI-MS `.CDF`) | Those read by mzR | `lcmsPlot(paths)` | |
#' | | ThermoFisher `.raw` | Those read by the Thermo RawFileReader | `lcmsPlot(paths)` | rawrr, MsBackendRawFileReader |
#' | xcms | `XcmsExperiment`, `MsExperiment`, `XCMSnExp` | xcms 4 and later | `lcmsPlot(object)` | |
#' | | `XChromatograms`, `XChromatogram`, `MChromatograms` | xcms 4 and later | `lcmsPlot(object)` | |
#' | | `xcmsRaw` (legacy xcms) | Any | [create_xcms_raw_list()] | |
#' | msPurity | `purityA` | Any | `lcmsPlot(object)` | msPurity |
#' | Compound Discoverer | `.cdResult` files (SQLite) | 3.x | `lcmsPlot(path)` + [lp_compound_discoverer()] | |
#' | | Scripting node exports | 3.3, 3.4, 3.5 | [CompoundDiscovererNodeSource()] | |
#' | MZmine | Feature list CSV exports | 2 (`row m/z` columns) and 3 or later (`datafile:` columns) | [MZmineFeatureListsSource()] | |
#' | MS-DIAL | Peak list CSV or TSV exports | Any with the `Precursor m/z`, `RT (min)`, `RT left(min)`, `RT right (min)`, `Area`, and `Height` columns | [MsDialPeaksSource()] | |
#' | LipidSearch | Result files | 4.2 and 5.2 | [LipidSearchSource()] + [lp_lipid_search()] | |
#' | mzTab-M | Small molecule features (SMF) | 2.0.0-M | [MzTabMSource()] + [lp_mztab()] | |
#' | Spectral libraries | MSP, or mzML reference spectra | Any | `lp_spectra(spectral_match_db = )` | |
#'
#' Results from external software are paired with the raw files they were
#' computed from (`sample_paths`), so that the chromatograms are extracted from
#' the original data and the peaks reported by the software are drawn on top
#' of them.
#'
#' @return None: this page documents the supported inputs.
#' @name lcmsPlot-formats
#' @seealso [lcmsPlot()] to read the data, and the "Interoperability" section
#' of `vignette("lcms_data_plotting", package = "lcmsPlot")` for examples.
NULL
