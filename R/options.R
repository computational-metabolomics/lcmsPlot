default_options <- list(
  sample_id_column = "sample_id",
  chromatograms = list(
    show = FALSE,
    features = NULL,
    sample_ids = NULL,
    ppm = 10,
    rt_tol = 10,
    highlight_peaks = FALSE,
    highlight_peaks_color = NULL,
    aggregation_fun = "max"
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
    interval = 3, # mode=across_peak
    spectral_match_db = NULL,
    match_target_index = 1
  ),
  intensity_maps = list(
    show = FALSE,
    mz_range = NULL,
    rt_range = NULL,
    density = FALSE
  ),
  facets = list(
    facets = NULL,
    ncol = NULL,
    nrow = NULL
  ),
  grid = list(
    rows = NULL,
    cols = NULL
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
