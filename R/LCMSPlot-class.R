#' Create a LCMSPlot object
#'
#' @param dataset An object of type XCMSnExp or MsExperiment
#' @param sample_id_column Which column should be used as the sample ID
#' @returns A LCMSPlot object
#' @export
lcmsPlot <- function(dataset, sample_id_column = 'sample_id') {
  if (inherits(dataset, "XCMSnExp")) {
    from_XCMSnExp(dataset, sample_id_column)
  } else {
    return(NULL) # TODO
  }
}

setOldClass(c("gg", "ggplot"))

#' LCMSPlot class
#'
#' @slot samples The samples metadata data frame
#' @slot features The features (grouped raw peaks) data frame
#' @slot peaks The raw peaks data frame
#' @slot options The object options
#' @slot .data The PeaklyChromatogramContainer object
#' @slot .plot The underlying ggplot2 object
#' @export
setClass(
  "LCMSPlot",
  slots = list(
    samples = "data.frame",
    features = "data.frame",
    peaks = "data.frame",
    options = "list",
    .data = "PeaklyChromatogramContainer",
    .plot = "ANY"
  ),
  prototype = list(
    samples = NULL,
    features = NULL,
    peaks = NULL,
    options = default_options,
    .data = NULL,
    .plot = NULL
  )
)

#' Apply a function to the LCMSPlot object
#'
#' @param e1 A LCMSPlot object
#' @param e2 A function that takes a LCMSPlot object and returns another
#' @export
setMethod(
  f = "+",
  signature = c("LCMSPlot", "function"),
  definition = function(e1, e2) { e2(e1) }
)

#' Create a LCMSPlot object from an XCMSnExp object
#'
#' @param xcmsnExp The XCMSnExp object
#' @param sample_id_column Which column should be used as the sample ID
#' @returns A LCMSPlot object
#' @export
setGeneric(
  "from_XCMSnExp",
  function(xcmsnExp, sample_id_column) standardGeneric("from_XCMSnExp")
)

#' @rdname from_XCMSnExp
#' @importFrom xcms fileNames phenoData featureDefinitions chromPeaks groupnames
setMethod(
  f = "from_XCMSnExp",
  signature = c("XCMSnExp", "character"),
  definition = function(xcmsnExp, sample_id_column) {
    sample_metadata <- phenoData(xcmsnExp)@data %>%
      mutate(
        sample_index = row_number(),
        sample_id = .data[[sample_id_column]],
        sample_path = fileNames(xcmsnExp)
      )

    features <- as.data.frame(featureDefinitions(xcmsnExp)) %>%
      rename(mz = mzmed, rt = rtmed) %>%
      mutate(name = groupnames(xcmsnExp)) %>%
      xcms_utils$format_feature_identifiers(num_digits_rt = 0, num_digits_mz = 4)

    new("LCMSPlot",
        samples = sample_metadata,
        features = features,
        peaks = as.data.frame(chromPeaks(xcmsnExp)),
        .data = new("PeaklyChromatogramContainer"),
        .plot = ggplot()
    )
  }
)

#' Create chromatograms from features
#'
#' @param obj The LCMSPlot object
#' @param feature_ids A vector of feature IDs to plot
#' @param sample_ids A vector of sample IDs to plot
#' @param ppm The ppm error for the chromatograms
#' @param rt_tol The RT tolerance for the chromatograms
#' @param highlight_peaks Whether to highlight the picked peaks
#' @returns The updated LCMSPlot object
#' @export
setGeneric(
  "lcmsPlot_chromatogram_from_features",
  function(obj, feature_ids, sample_ids, ppm, rt_tol, highlight_peaks) standardGeneric("lcmsPlot_chromatogram_from_features")
)

#' @rdname lcmsPlot_chromatogram_from_features
setMethod(
  f = "lcmsPlot_chromatogram_from_features",
  signature = c("LCMSPlot", "character", "character", "numeric", "numeric", "logical"),
  definition = function(obj, feature_ids, sample_ids, ppm, rt_tol, highlight_peaks) {
    samples <- obj@samples %>% filter(sample_id %in% sample_ids)
    paths <- samples %>% pull(sample_path)
    # xraw_list <- xcms_utils$get_xraw_list_from_files(paths) # TODO: cache
    raw_data <- io_utils$get_raw_data(paths)
    features_matrix <- obj@features %>% filter(name %in% feature_ids)

    obj@.data <- create_chromatograms_from_features(
      raw_data,
      features_matrix,
      obj@peaks,
      samples,
      ppm,
      rt_tol)

    obj@options$chromatograms$highlight_peaks <- highlight_peaks

    return(obj)
  }
)

#' Create chromatograms from features
#'
#' @param obj The LCMSPlot object
#' @param samples ...
#' @param rt_range ...
#' @param mz_range ...
#' @returns The updated LCMSPlot object
#' @export
setGeneric(
  "lcmsPlot_chromatogram_from_raw",
  function(obj, features, sample_ids, ppm, rt_tol) standardGeneric("lcmsPlot_chromatogram_from_raw")
)

#' @rdname lcmsPlot_chromatogram_from_raw
setMethod(
  f = "lcmsPlot_chromatogram_from_raw",
  signature = c("LCMSPlot", "matrix", "character", "numeric", "numeric"),
  definition = function(obj, features, sample_ids, ppm, rt_tol) {
    samples <- obj@samples %>% filter(sample_id %in% sample_ids)
    paths <- samples %>% pull(sample_path)
    raw_data <- io_utils$get_raw_data(paths)
    obj@.data <- create_chromatograms_from_raw(
      raw_data,
      features,
      samples,
      ppm,
      rt_tol)

    return(obj)
  }
)

#' Set mass traces
#'
#' @export
setGeneric(
  "lcmsPlot_mass_traces",
  function(obj) standardGeneric("lcmsPlot_mass_traces")
)

#' @rdname lcmsPlot_mass_traces
setMethod(
  f = "lcmsPlot_mass_traces",
  signature = c("LCMSPlot"),
  definition = function(obj) {
    obj@options$mass_traces$show <- TRUE
    return(obj)
  }
)

#' Set the arrangement of the chromatograms
#'
#' @param obj The LCMSPlot object
#' @param group_by ...
#' @returns The updated LCMSPlot object
#' @export
setGeneric(
  "lcmsPlot_arrange",
  function(obj, group_by) standardGeneric("lcmsPlot_arrange")
)

#' @rdname lcmsPlot_arrange
setMethod(
  f = "lcmsPlot_arrange",
  signature = c("LCMSPlot", "characterOrNULL"),
  definition = function(obj, group_by) {
    obj@options$arrangement <- list(
      group_by = group_by
    )
    return(obj)
  }
)

#' Set the arrangement of the chromatograms
#'
#' @param obj The LCMSPlot object
#' @param facets ...
#' @param ncol ...
#' @param nrow ...
#' @returns The updated LCMSPlot object
#' @export
setGeneric(
  "lcmsPlot_facets",
  function(obj, facets, ncol, nrow) standardGeneric("lcmsPlot_facets")
)

#' @rdname lcmsPlot_facets
setMethod(
  f = "lcmsPlot_facets",
  signature = c("LCMSPlot", "characterOrNULL", "numericOrNULL", "numericOrNULL"),
  definition = function(obj, facets, ncol, nrow) {
    obj@options$facets <- list(
      facets = facets,
      ncol = ncol,
      nrow = nrow
    )
    return(obj)
  }
)

#' Set the arrangement of the chromatograms
#'
#' @param obj The LCMSPlot object
#' @param rows ...
#' @param cols ...
#' @returns The updated LCMSPlot object
#' @export
setGeneric(
  "lcmsPlot_grid",
  function(obj, rows, cols) standardGeneric("lcmsPlot_grid")
)

#' @rdname lcmsPlot_grid
setMethod(
  f = "lcmsPlot_grid",
  signature = c("LCMSPlot", "characterOrNULL", "characterOrNULL"),
  definition = function(obj, rows, cols) {
    obj@options$grid <- list(
      rows = rows,
      cols = cols
    )
    return(obj)
  }
)

#' Set the labels of the plot
#'
#' @param obj The LCMSPlot object
#' @param title The title of the plot
#' @param legend The legend's title
#' @returns The updated LCMSPlot object
#' @export
setGeneric(
  "lcmsPlot_labels",
  function(obj, title, legend) standardGeneric("lcmsPlot_labels")
)

#' @rdname lcmsPlot_labels
setMethod(
  f = "lcmsPlot_labels",
  signature = c("LCMSPlot", "characterOrNULL", "characterOrNULL"),
  definition = function(obj, title, legend) {
    obj@options$labels <- list(
      title = title,
      legend = legend
    )
    return(obj)
  }
)

#' Set the legend options
#'
#' @param obj The LCMSPlot object
#' @param position The legend's position in the plot
#' @returns The updated LCMSPlot object
#' @export
setGeneric(
  "lcmsPlot_legend",
  function(obj, position) standardGeneric("lcmsPlot_legend")
)

#' @rdname lcmsPlot_legend
setMethod(
  f = "lcmsPlot_legend",
  signature = c("LCMSPlot", "characterOrNULL"),
  definition = function(obj, position) {
    obj@options$legend <- list(
      position = position
    )
    return(obj)
  }
)

#' Plot the LCMSPlot object
#'
#' @param object The LCMSPlot object
#'
#' @importFrom ggplot2 ggplot aes geom_line geom_point geom_ribbon facet_wrap
#' @export
setMethod(
  f = "show",
  signature = "LCMSPlot",
  function(object) {
    dataset_types <- c("chromatograms", "mass_traces")
    datasets <- lapply(dataset_types, function(dataset_name) {
      if (object@options[[dataset_name]]$show) {
        dataset <- slot(object@.data, dataset_name)
        dataset <- merge_by_index(dataset, object@.data@metadata, index_col = 'metadata_index')
        return(dataset)
      } else {
        return(NULL)
      }
    })
    names(datasets) <- dataset_types
    datasets <- remove_null_elements(datasets)

    object@.plot <- plot_data(datasets, object)

    print(object@.plot)
  }
)

#' Define the chromatograms to plot
#'
#' @param features A vector of feature IDs to plot (TODO: change)
#' @param sample_ids A vector of sample IDs to plot
#' @param ppm The ppm error for the chromatograms
#' @param rt_tol The RT tolerance for the chromatograms
#' @param highlight_peaks Whether to highlight the picked peaks
#' @returns A function that takes and returns a LCMSPlot object
#' @export
chromatogram <- function(
  features,
  sample_ids,
  ppm,
  rt_tol,
  highlight_peaks = FALSE
) {
  function(obj) {
    if (is.character(features)) {
      obj <- lcmsPlot_chromatogram_from_features(obj, features, sample_ids, ppm, rt_tol, highlight_peaks)
    } else {
      obj <- lcmsPlot_chromatogram_from_raw(obj, features, sample_ids, ppm, rt_tol)
    }

    return(obj)
  }
}

#' Define the mass trace to plot
#'
#' @returns A function that takes and returns a LCMSPlot object
#' @export
mass_trace <- function() {
  function(obj) {
    obj <- lcmsPlot_mass_traces(obj)
    return(obj)
  }
}

#' Define the arrangement of chromatograms
#'
#' @param group_by The column to group by (in the samples metadata)
#' @returns A function that takes and returns a LCMSPlot object
#' @export
arrange <- create_lcmsPlot_function("lcmsPlot_arrange", default_options$arrangement)

#' Define the arrangement of chromatograms
#'
#' @param rows ...
#' @param cols ...
#' @returns A function that takes and returns a LCMSPlot object
#' @export
grid <- create_lcmsPlot_function("lcmsPlot_grid", default_options$grid)

#' Define the arrangement of chromatograms
#'
#' @param facets ...
#' @param ncol ...
#' @param nrow ...
#' @returns A function that takes and returns a LCMSPlot object
#' @export
facets <- create_lcmsPlot_function("lcmsPlot_facets", default_options$facets)

#' Define the labels of the plot, such as title and legend
#'
#' @param title The plot title
#' @param legend The legend's title
#' @returns A function that takes and returns a LCMSPlot object
#' @export
labels <- create_lcmsPlot_function("lcmsPlot_labels", default_options$labels)

#' Define the legend layout
#'
#' @param position The legend position
#' @returns A function that takes and returns a LCMSPlot object
#' @export
legend <- create_lcmsPlot_function("lcmsPlot_legend", default_options$legend)
