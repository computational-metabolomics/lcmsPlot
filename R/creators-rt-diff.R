#' Create an instance of class `lcmsPlotDataContainer` from an RT adjustment
#' dataset
#'
#' This function creates an `lcmsPlotDataContainer` object from
#' a dataset that contains the data to generate RT raw vs adjusted plots.
#'
#' @param obj An instance of class `lcmsPlotDataContainer`.
#' @param options A `list` representing the plot object's options.
#' @return An instance of class `lcmsPlotDataContainer` with the created
#' RT adjustment dataset, a `data.frame`
#' with columns `rt_raw`, `rt_adj`, and `diff`.
#' @keywords internal
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
            feature_metadata_id = xcms::fromFile(obj@data_obj)
        )

        obj@rt_diff <- rt_diff

        validObject(obj)

        return(obj)
    }
)
