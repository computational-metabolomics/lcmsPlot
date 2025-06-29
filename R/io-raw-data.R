io_utils <- list(
  get_raw_data = function(sample_paths) {
    raw_data <- list()

    for (sample_path in sample_paths) {
      raw_data[[sample_path]] <- mzR::openMSfile(sample_path)
    }

    return(raw_data)
  },

  close_raw_data = function(raw_data) {
    for (raw_obj in raw_data) {
      mzR::close(raw_obj)
    }
  }
)
