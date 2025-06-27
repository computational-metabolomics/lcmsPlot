default_options <- list(
  sample_id_column = "sample_id",
  chromatograms = list(
    show = TRUE,
    highlight_peaks = FALSE
  ),
  mass_traces = list(
    show = FALSE
  ),
  spectra = list(
    show = FALSE,
    mode = 'closest_apex', # One of: closest, closest_apex, across_peak
    ms_level = 1,
    rt = NULL, # For mode=closest
    interval = 3 # mode=across_peak
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
