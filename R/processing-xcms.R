xcms_utils <- list(
    has_rt_alignment_been_performed = function(xcms_obj) {
        xcms::processHistory(xcms_obj) |>
            lapply(function(p) p@type == "Retention time correction") |>
            unlist() |>
            any()
    },

    # Feature identifiers come from the row names of the feature definitions,
    # i.e. whatever the user has set them to (or the "FT0001"-style defaults
    # xcms assigns during correspondence). Only when an object carries no row
    # names at all do we synthesise ids in the same "FT" style, so that a single
    # nomenclature is used throughout.
    feature_names = function(feature_definitions) {
        names <- rownames(feature_definitions)
        n <- nrow(feature_definitions)

        if (is.null(names) || length(names) != n || anyNA(names) ||
            !all(nzchar(names))) {
            width <- max(nchar(as.character(n)), 1L)
            names <- sprintf(paste0("FT%0", width, "d"), seq_len(n))
        }

        return(as.character(names))
    },

    collapse_peak_indices = function(features) {
        features |>
            # Collapse peakidx column to a comma-separated string
            mutate(
                peakidx = vapply(
                    .data$peakidx,
                    function(x) paste(x, collapse = ","),
                    FUN.VALUE = character(1))
            ) |>
            # Reorder columns
            relocate("name")
    }
)
