#' Creates an lcmsPlotDataContainer from a TIC
#'
#' @param obj A lcmsPlotDataContainer object
#' @param options The lcmsPlot options
#' @returns A lcmsPlotDataContainer object
#' @export
setGeneric(
  "create_total_ion_current",
  function(obj, options) standardGeneric("create_total_ion_current")
)

#' @rdname create_total_ion_current
setMethod(
  f = "create_total_ion_current",
  signature = c("lcmsPlotDataContainer", "list"),
  definition = function(obj, options) {
    metadata <- obj@metadata %>% filter(.data$sample_id %in% options$total_ion_current$sample_ids)

    tc <- xcms::spectra(obj@data_obj) %>%
      Spectra::tic() %>%
      split(f = xcms::fromFile(obj@data_obj))
    tc <- lapply(seq_along(tc), function (i) {
      x <- tc[[i]] %>%
        as.data.frame() %>%
        mutate(metadata_index = i, additional_metadata_index = i)
      colnames(x)[1] <- "intensity"
      return(x)
    }) %>%
      do.call(rbind, .)

    obj@total_ion_current <- tc

    validObject(obj)

    return(obj)
  }
)
