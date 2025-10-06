#' Creates an lcmsPlotDataContainer from a TIC
#'
#' @param obj A lcmsPlotDataContainer object
#' @param options The lcmsPlot options
#' @return A lcmsPlotDataContainer object
#' @export
#' @examples
#' data_obj <- get_XCMSnExp_object_example(indices = c(1, 7))
#' data_container <- create_data_container_from_obj(data_obj, sample_id_column = "sample_name", metadata = NULL)
#' opts <- default_options()
#' data_container <- create_total_ion_current(data_container, opts)
setGeneric(
  "create_total_ion_current",
  function(obj, options) standardGeneric("create_total_ion_current")
)

#' @rdname create_total_ion_current
setMethod(
  f = "create_total_ion_current",
  signature = c("lcmsPlotDataContainer", "list"),
  definition = function(obj, options) {
    if (is(obj@data_obj, "MsExperiment")) {
      tc <- xcms::spectra(obj@data_obj) %>% Spectra::tic()
    } else {
      tc <- MSnbase::tic(obj@data_obj)
    }

    tc <- tc %>% split(f = xcms::fromFile(obj@data_obj))
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
