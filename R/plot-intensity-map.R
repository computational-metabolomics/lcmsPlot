#' Plot an LC-MS intensity map
#'
#' `plot_intensity_map()` generates a ggplot2 from an intensity map which is a
#' point-cloud of m/z and RT points.
#'
#' @param datasets A named `list` of data frames containing the primary datasets
#' to plot. For a 2D intensity map the used key is `intensity_maps`.
#' @param supporting_datasets A `list` of supporting data frames, such as
#' detected peaks, used for highlighting features.
#' @param options A list of plot options, controlling units, faceting,
#' highlighting, and other visual parameters.
#' @param single A `logical` value that indicates whether it should treat
#' the plot as a single dataset variant, which can affect faceting
#' and layout behavior. Default is `FALSE`.
#' @return A `ggplot` object representing the 2D intensity map plot.
#' @keywords internal
plot_intensity_map <- function(
    datasets,
    supporting_datasets,
    options,
    single = FALSE
) {
    dataset <- datasets$intensity_maps

    extra_layers <- list(
        legend_title(options),
        faceting(options, single),
        grid_layout(options, single)
    )
    extra_layers <- remove_null_elements(extra_layers)

    x_dim <- options$intensity_maps$x_dim
    y_dim <- options$intensity_maps$y_dim
    fill_scale <- options$intensity_maps$fill_scale

    dim_label <- function(dim) if (dim == "rt") "RT (sec)" else "m/z"

    if (options$intensity_maps$density) {
        ggplot(dataset, aes(x = .data[[x_dim]], y = .data[[y_dim]])) +
            geom_density_2d_filled(contour_var = "ndensity") +
            (if (!is.null(fill_scale)) fill_scale else scale_fill_viridis_d()) +
            labs(x = dim_label(x_dim), y = dim_label(y_dim), fill = "Density") +
            theme_minimal() +
            extra_layers
    } else {
        ggplot(
            dataset,
            aes(
                x = .data[[x_dim]],
                y = .data[[y_dim]],
                fill = log1p(.data$intensity))) +
            geom_tile() +
            (if (!is.null(fill_scale)) fill_scale else scale_fill_viridis_c()) +
            labs(x = dim_label(x_dim), y = dim_label(y_dim), fill = "log(Intensity)") +
            theme_minimal() +
            extra_layers
    }
}
