#' Create an instance of class `lcmsPlotDataContainer` from total ion current
#' (TIC) dataset
#'
#' This function creates an `lcmsPlotDataContainer` object from
#' a dataset that contains TIC values for each sample.
#'
#' @param obj An instance of class `lcmsPlotDataContainer`.
#' @param options A `list` representing the plot object's options.
#' @return An instance of class `lcmsPlotDataContainer` with the created
#' RT adjustment dataset, a `data.frame`
#' with columns `sample_id` and `intensity`.
#' @keywords internal
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
            tc <- xcms::spectra(obj@data_obj) |> Spectra::tic()
        } else {
            tc <- MSnbase::tic(obj@data_obj)
        }

        tc <- tc |> split(f = xcms::fromFile(obj@data_obj))
        tc <- do.call(
            rbind,
            lapply(seq_along(tc), function(i) {
                x <- tc[[i]] |>
                    as.data.frame() |>
                    mutate(metadata_index = i, additional_metadata_index = i)
                colnames(x)[1] <- "intensity"
                x
            })
        )

        obj@total_ion_current <- tc

        validObject(obj)

        return(obj)
    }
)
