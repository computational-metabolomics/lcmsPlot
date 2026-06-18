#' Create an instance of class `lcmsPlotDataContainer` from total ion current
#' (TIC) dataset
#'
#' This function creates an `lcmsPlotDataContainer` object from
#' a dataset that contains TIC values for each sample.
#'
#' @param obj An instance of class `lcmsPlotDataContainer`.
#' @param options A `list` representing the plot object's options.
#' @return An instance of class `lcmsPlotDataContainer` with the created
#' RT adjustment dataset, a `tibble`
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
        if (is.character(obj@data_obj)) {
            tc <- .tic_from_raw_files(obj@metadata, options)
        } else {
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
                        as_tibble() |>
                        mutate(metadata_index = i, feature_metadata_id = i)
                    colnames(x)[1] <- "intensity"
                    x
                })
            )
        }

        obj@total_ion_current <- tc

        validObject(obj)

        return(obj)
    }
)

#' Build the total ion current dataset from raw MS files
#'
#' Reads the per-scan total ion current (TIC) from the headers of the raw
#' files referenced in `metadata` and returns it in the layout expected by the
#' `total_ion_current` slot of an `lcmsPlotDataContainer`.
#'
#' @param metadata A `data.frame` of sample metadata with at least the columns
#' `sample_id`, `sample_index`, and `sample_path`.
#' @param options A `list` of plot options. `options$total_ion_current$sample_ids`
#' selects which samples to read; `options$parallel_param`, when set, is used to
#' read files in parallel.
#' @return A `tibble` with columns `intensity`, `metadata_index`, and
#' `feature_metadata_id` (one row per MS1 scan).
#' @keywords internal
.tic_from_raw_files <- function(metadata, options) {
    metadata <- metadata |>
        filter(.data$sample_id %in% options$total_ion_current$sample_ids)
    raw_data <- io_get_raw_data(metadata$sample_path)

    process_sample <- function(i) {
        sample_metadata <- metadata[i, ]
        raw_obj <- raw_data[[sample_metadata$sample_path]]

        hdr <- ms_header(raw_obj)
        hdr <- hdr[hdr$msLevel == 1, ]

        tibble(
            intensity = hdr$totIonCurrent,
            metadata_index = sample_metadata$sample_index,
            feature_metadata_id = sample_metadata$sample_index
        )
    }

    if (!is.null(options$parallel_param)) {
        tic_list <- BiocParallel::bplapply(
            seq_len(nrow(metadata)),
            process_sample,
            BPPARAM = options$parallel_param)
    } else {
        tic_list <- lapply(seq_len(nrow(metadata)), process_sample)
    }

    io_close_raw_data(raw_data)

    do.call(rbind, tic_list)
}
