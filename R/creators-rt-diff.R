#' Creates an lcmsPlotDataContainer from an RT diff dataset
#'
#' @param obj A lcmsPlotDataContainer object
#' @param options The lcmsPlot options
#' @returns A lcmsPlotDataContainer object
#' @export
#' @examples
#' data_obj <- get_XCMSnExp_object_example(indices = c(1, 7), should_group_peaks = TRUE)
#' data_container <- create_data_container_from_obj(data_obj, sample_id_column = "sample_name", metadata = NULL)
#' opts <- default_options()
#' data_container <- create_rt_diff(data_container, opts)
setGeneric(
  "create_rt_diff",
  function(obj, options) standardGeneric("create_rt_diff")
)

#' @rdname create_rt_diff
setMethod(
  f = "create_rt_diff",
  signature = c("lcmsPlotDataContainer", "list"),
  definition = function(obj, options) {
    rt_raw <- xcms::rtime(obj@data_obj, adjusted = FALSE)
    rt_adj <- xcms::adjustedRtime(obj@data_obj)

    rt_diff <- data.frame(
      rt_raw = rt_raw,
      rt_adj = rt_adj,
      diff = rt_adj - rt_raw,
      metadata_index = xcms::fromFile(obj@data_obj),
      additional_metadata_index = xcms::fromFile(obj@data_obj)
    )

    obj@rt_diff <- rt_diff

    validObject(obj)

    return(obj)
  }
)
