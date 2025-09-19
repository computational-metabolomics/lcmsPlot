default_options <- list(
  sample_id_column = "sample_id",
  parallel_param = NULL,
  bypass_plot_generation = FALSE,
  batch_size = NULL,
  batch_index = 1,
  chromatograms = list(
    show = FALSE,
    features = NULL,
    sample_ids = NULL,
    ppm = 10,
    rt_tol = 10,
    highlight_peaks = FALSE,
    highlight_peaks_color = NULL,
    highlight_peaks_factor = "sample_id",
    aggregation_fun = "max",
    rt_adjusted = FALSE,
    rt_unit = "second",
    fill_gaps = FALSE
  ),
  mass_traces = list(
    show = FALSE
  ),
  spectra = list(
    show = FALSE,
    sample_ids = NULL,
    mode = 'closest_apex', # One of: closest, closest_apex, across_peak
    ms_level = 1,
    rt = NULL, # For mode=closest
    scan_index = NULL,
    interval = 3, # mode=across_peak
    spectral_match_db = NULL,
    match_target_index = 1
  ),
  total_ion_current = list(
    show = FALSE,
    sample_ids = NULL,
    type = "boxplot"
  ),
  intensity_maps = list(
    show = FALSE,
    sample_ids = NULL,
    mz_range = NULL,
    rt_range = NULL,
    density = FALSE
  ),
  rt_diff = list(
    show = FALSE
  ),
  facets = list(
    facets = NULL,
    ncol = NULL,
    nrow = NULL,
    free_x = FALSE,
    free_y = FALSE
  ),
  grid = list(
    rows = NULL,
    cols = NULL,
    free_y = FALSE
  ),
  arrangement = list(
    group_by = NULL
  ),
  labels = list(
    title = NULL,
    legend = NULL
  ),
  legend = list(
    position = "right"
  ),
  rt_lines = list(),
  layout = list(
    design = NULL
  )
)

get_grouping_variables <- function(opts) {
  if (!is.null(opts$facets$facets)) {
    opts$facets$facets
  } else {
    c(opts$grid$rows, opts$grid$cols)
  }
}
