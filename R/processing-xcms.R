xcms_utils <- list(
  format_feature_identifiers = function(features, num_digits_rt = 0, num_digits_mz = 0) {
    features %>%
      # Extract decorations from 'name' (text after the first '_')
      mutate(idsDeco = stringr::str_extract(name, "_.*$")) %>%
      # Replace NA with empty string for idsDeco
      mutate(idsDeco = ifelse(is.na(idsDeco), "", idsDeco)) %>%
      # Create the custom name
      mutate(namecustom = make.unique(
        paste0("M", round(mz, num_digits_mz),
               "T", round(rt, num_digits_rt),
               idsDeco)
      )) %>%
      # Collapse peakidx column to a comma-separated string
      mutate(peakidx = sapply(peakidx, function(x) paste(x, collapse = ","))) %>%
      # Reorder columns
      relocate(name, namecustom) %>%
      select(-idsDeco)
  }

  # get_xraw_list_from_files = function(sample_paths) {
  #   xraw_list <- list()
  #
  #   for (sample_path in sample_paths) {
  #     xraw_list[[sample_path]] <- xcms::xcmsRaw(sample_path)
  #   }
  #
  #   return(xraw_list)
  # },
  #
  # get_xraw_list = function(directory) {
  #   sample_paths <- list.files(directory, full.names = TRUE)
  #   sample_paths <- sample_paths[!grepl('MS2', sample_paths)]
  #
  #   return(xcms_utils$get_xraw_list_from_files(sample_paths))
  # }
)
