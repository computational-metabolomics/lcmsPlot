#' Creates an lcmsPlotDataContainer from an intensity map
#'
#' @param obj A lcmsPlotDataContainer object
#' @param options The lcmsPlot options
#' @return A lcmsPlotDataContainer object
#' @export
#' @examples
#' data_obj <- get_XCMSnExp_object_example(indices = c(1, 7))
#' data_container <- create_data_container_from_obj(data_obj, sample_id_column = "sample_name", metadata = NULL)
#'
#' opts <- default_options()
#' opts$intensity_maps$rt_range <- c(2000, 2100)
#' opts$intensity_maps$mz_range <- c(200, 300)
#' opts$intensity_maps$sample_ids <- c('ko15', 'wt15')
#'
#' data_container <- create_intensity_map(data_container, opts)
setGeneric(
  "create_intensity_map",
  function(obj, options) standardGeneric("create_intensity_map")
)

#' @rdname create_intensity_map
setMethod(
  f = "create_intensity_map",
  signature = c("lcmsPlotDataContainer", "list"),
  definition = function(obj, options) {
    metadata <- obj@metadata %>%
      filter(.data$sample_id %in% options$intensity_maps$sample_ids)

    process_sample <- function(i) {
      sample_metadata <- metadata[i, ]
      ms <- mzR::openMSfile(sample_metadata$sample_path)
      hdr <- mzR::header(ms)

      rt_range <- options$intensity_maps$rt_range
      mz_range <- options$intensity_maps$mz_range

      idx <- which(hdr$retentionTime >= rt_range[1] &
                     hdr$retentionTime <= rt_range[2])

      scans <- lapply(idx, function(j) {
        pk <- mzR::peaks(ms, j)
        if (nrow(pk) > 0) {
          pk <- pk[pk[,1] >= mz_range[1] & pk[,1] <= mz_range[2], ]
          if (nrow(pk) > 0) {
            data.frame(
              rt = rep(hdr$retentionTime[j], nrow(pk)),
              mz = pk[,1],
              intensity = pk[,2]
            )
          } else {
            NULL
          }
        } else {
          NULL
        }
      })

      mzR::close(ms)

      df <- do.call(rbind, scans)

      if (!is.null(df) && nrow(df) > 0) {
        df %>%
          mutate(
            rt = round(.data$rt, 1),
            mz = round(.data$mz, 1)
          ) %>%
          group_by(.data$rt, .data$mz) %>%
          summarize(intensity = sum(.data$intensity), .groups = "drop") %>%
          mutate(
            metadata_index = sample_metadata$sample_index,
            additional_metadata_index = sample_metadata$sample_index
          )
      } else {
        NULL
      }
    }

    if (!is.null(options$parallel_param)) {
      results <- BiocParallel::bplapply(
        seq_len(nrow(metadata)),
        process_sample,
        BPPARAM = options$parallel_param
      )
    } else {
      results <- lapply(seq_len(nrow(metadata)), process_sample)
    }

    obj@intensity_maps <- bind_rows(results)

    validObject(obj)

    return(obj)
  }
)
