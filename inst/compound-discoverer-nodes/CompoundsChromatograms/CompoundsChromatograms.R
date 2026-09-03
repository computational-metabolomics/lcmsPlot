#!/usr/bin/env Rscript
# =============================================================================
# lcmsPlot Compound Discoverer Scripting-Node script
# -----------------------------------------------------------------------------
# Runs inside a Compound Discoverer (CD) custom Scripting Node. CD passes the
# path to a node_args.json file as the 6th command-line argument. The script:
#   1. Loads lcmsPlot (the installed package, or a source tree if configured).
#   2. Reads the raw (.raw) file paths from CD's exported Input Files table.
#   3. Builds a CompoundDiscovererNodeSource, extracts chromatograms once, and
#      saves one plot per compound (faceted by sample) into a directory next to
#      the CD result file.
#   4. Writes a node_response.json that adds a plot-filename column to the
#      Compounds table (via CD's filename cell renderer).
#
# Configuration comes from the node's parameters (declared in node.json and
# delivered in node_args.json's "NodeParameters"); each has a default.
#
# Run inside CD, or standalone:
#   Rscript CompoundsChromatograms.R <node_args.json>
#
# Requested Tables and Columns (set on the node):
#   Compounds; Compounds per File; Input Files
# =============================================================================

# Directory of this script (from the --file= argument).
.this_file <- local({
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", args, value = TRUE)
    if (length(file_arg) == 1) {
        normalizePath(sub("^--file=", "", file_arg))
    } else {
        NA_character_
    }
})
.script_dir <- if (!is.na(.this_file)) dirname(.this_file) else getwd()

# -----------------------------------------------------------------------------
# Logging: tee progress to stdout (shown as CD "Info") and a flushed log file;
# problems go to stderr (CD "Warning"). Log path: LCMSPLOT_LOG env >
# <script_dir>/lcmsPlot-cd-node.log > <tempdir>. A missing log file after a run
# means R never executed.
# -----------------------------------------------------------------------------
.log_file <- local({
    candidates <- c(
        Sys.getenv("LCMSPLOT_LOG", unset = ""),
        file.path(.script_dir, "lcmsPlot-cd-node.log"),
        file.path(tempdir(), "lcmsPlot-cd-node.log")
    )
    candidates <- candidates[nzchar(candidates)]
    for (cand in candidates) {
        ok <- tryCatch({
            con <- file(cand, open = "at")
            close(con)
            TRUE
        }, error = function(e) FALSE)
        if (isTRUE(ok)) return(cand)
    }
    ""
})

# Append one timestamped line to the log file. Never throws.
.log_write <- function(text) {
    if (!nzchar(.log_file)) {
        return(invisible(NULL))
    }
    line <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%OS3"), "  ", text)
    tryCatch({
        con <- file(.log_file, open = "at")
        on.exit(close(con))
        writeLines(line, con)
    }, error = function(e) invisible(NULL))
    invisible(NULL)
}

# Progress -> stdout + log file (CD shows stdout as "Info").
log_msg <- function(...) {
    text <- paste0(..., collapse = "")
    cat(text, "\n", sep = "")
    .log_write(text)
    invisible(NULL)
}

# Problems -> stderr + log file (CD shows stderr as "Warning").
log_err <- function(...) {
    text <- paste0(..., collapse = "")
    message(text)
    .log_write(text)
    invisible(NULL)
}

options(warn = 1)

# -----------------------------------------------------------------------------
# Configuration is read from the node parameters delivered in node_args.json's
# "NodeParameters" (keyed by each parameter's Name; see node.json). Missing
# parameters fall back to the defaults set in main().
# -----------------------------------------------------------------------------

# Load lcmsPlot: from the installed package when `dir` is empty (deployment), or
# from a source tree via pkgload when `dir` is set (development). Load-time
# messages/warnings go to the log file only.
load_lcmsplot <- function(dir) {
    quiet <- function(expr) {
        withCallingHandlers(
            expr,
            message = function(m) {
                .log_write(paste0("[load msg] ",
                                  sub("\n$", "", conditionMessage(m))))
                invokeRestart("muffleMessage")
            },
            warning = function(w) {
                .log_write(paste0("[load warn] ", conditionMessage(w)))
                invokeRestart("muffleWarning")
            }
        )
    }

    if (nzchar(dir)) {
        if (!requireNamespace("pkgload", quietly = TRUE)) {
            stop(
                "'pkgload' is required to load lcmsPlot from a source tree. ",
                "Install it with install.packages(\"pkgload\"), or clear the ",
                "'lcmsPlot Source Directory' parameter to use the installed ",
                "package.",
                call. = FALSE
            )
        }
        if (!dir.exists(dir)) {
            stop("lcmsPlot Source Directory not found: ", dir, call. = FALSE)
        }
        log_msg("STEP: loading lcmsPlot from source: ", dir)
        quiet(pkgload::load_all(dir, quiet = TRUE))
    } else {
        log_msg("STEP: loading the installed lcmsPlot package")
        if (!requireNamespace("lcmsPlot", quietly = TRUE)) {
            stop(
                "The 'lcmsPlot' package is not installed in this R. Install it ",
                "(see the node README), or set the 'lcmsPlot Source Directory' ",
                "parameter to a source tree.",
                call. = FALSE
            )
        }
        quiet(suppressPackageStartupMessages(
            library("lcmsPlot", character.only = TRUE)))
    }
    log_msg("  -> lcmsPlot loaded")
}

# Node parameters from node_args.json (empty list if absent / unreadable).
.read_node_parameters <- function(node_args) {
    tryCatch({
        np <- jsonlite::fromJSON(node_args, simplifyVector = TRUE)$NodeParameters
        if (is.null(np)) list() else as.list(np)
    }, error = function(e) list())
}

# One node parameter as character, or `default` when absent/blank.
.param <- function(params, name, default = "") {
    v <- params[[name]]
    if (is.null(v) || length(v) != 1 || is.na(v) ||
        !nzchar(trimws(as.character(v)))) {
        default
    } else {
        trimws(as.character(v))
    }
}

# One numeric node parameter, or `default` when absent/invalid.
.num_param <- function(params, name, default) {
    v <- suppressWarnings(as.numeric(.param(params, name, as.character(default))))
    if (is.na(v)) default else v
}

# =============================================================================

main <- function() {
    log_msg("STEP: resolving node_args.json path")
    node_args <- resolve_node_args()
    log_msg("Using node_args.json: ", node_args)

    log_msg("STEP: reading node parameters")
    params <- .read_node_parameters(node_args)
    lcmsplot_dir <- .param(params, "lcmsPlot Source Directory", "")
    compounds_query <- .param(params, "Compounds Query", "compound_rank <= 20")
    output_dir <- .param(params, "Output Directory", "")
    plot_column_name <- .param(params, "Plot Column Name", "Plot")
    plot_renderer <- .param(params, "Plot Cell Renderer GUID",
                               "EB29D794-4F2E-4785-8B80-A24D8C0FB3E4")
    plot_position <- .param(params, "Plot Column Position After", "Name")
    plot_path_mode <- .param(params, "Plot Path Mode", "absolute")
    plot_width <- .num_param(params, "Plot Width", 8)
    plot_height <- .num_param(params, "Plot Height", 5)
    plot_dpi <- .num_param(params, "Plot DPI", 150)
    log_msg("  lcmsPlot source: ",
            if (nzchar(lcmsplot_dir)) lcmsplot_dir else "(installed package)")
    log_msg("  compounds query: ",
            if (nzchar(compounds_query)) compounds_query else "(all)")
    log_msg("  output dir: ",
            if (nzchar(output_dir)) output_dir else "(next to result file)")
    log_msg("  plot: ", plot_width, "x", plot_height, " in @ ", plot_dpi,
            " dpi; column '", plot_column_name, "'; path mode '",
            plot_path_mode, "'")

    load_lcmsplot(lcmsplot_dir)

    log_msg("STEP: registering BiocParallel::SerialParam()")
    BiocParallel::register(BiocParallel::SerialParam())

    log_msg("STEP: resolving raw file paths from the Input Files table")
    sample_paths <- build_sample_paths(node_args)
    log_msg("Resolved ", length(sample_paths), " raw file path(s):")
    for (sp in sample_paths) log_msg("  - ", sp)

    log_msg("STEP: constructing CompoundDiscovererNodeSource")
    ds <- CompoundDiscovererNodeSource(
        node_args = node_args,
        sample_paths = sample_paths
    )
    .log_write(paste(utils::capture.output(show(ds)), collapse = "\n"))
    show(ds)

    # `+ lp_chromatogram()` extracts the chromatograms from the raw files (once).
    log_msg("STEP: building the lcmsPlot plot object (extracts chromatograms)")
    query <- if (nzchar(compounds_query) &&
                 !identical(toupper(compounds_query), "NULL")) {
        compounds_query
    } else {
        NULL
    }
    p <- lcmsPlot(ds) +
        lp_compound_discoverer(compounds_query = query, rt_extend = 5) +
        lp_chromatogram(highlight_peaks = TRUE)

    plot_dir <- resolve_plot_dir(node_args, output_dir)
    dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
    log_msg("STEP: rendering per-compound plots into: ", plot_dir)

    # Map each compound label to its "Compounds ID" for the response column.
    id_by_name <- ds@compounds |>
        dplyr::distinct(name, compound_id)

    names_to_plot <- unique(p@data@feature_metadata$name)
    log_msg("  ", length(names_to_plot), " compound(s) to plot")

    id_to_value <- list()
    n_done <- 0L
    for (nm in names_to_plot) {
        cid <- id_by_name$compound_id[match(nm, id_by_name$name)]
        out_path <- file.path(
            plot_dir, paste0(cid, "_", .sanitize_filename(nm), ".png"))

        render_compound_plot(p, nm, out_path, title = nm,
                             width = plot_width, height = plot_height,
                             dpi = plot_dpi)

        id_to_value[[length(id_to_value) + 1]] <- data.frame(
            compound_id = cid,
            value = plot_cell_value(out_path, node_args, plot_path_mode),
            stringsAsFactors = FALSE
        )
        n_done <- n_done + 1L
        if (n_done %% 10L == 0L || n_done == length(names_to_plot)) {
            log_msg("  rendered ", n_done, "/", length(names_to_plot))
        }
    }
    id_to_value <- do.call(rbind, id_to_value)
    log_msg("Per-compound plots written to: ", plot_dir)

    if (is.null(id_to_value) || nrow(id_to_value) == 0) {
        log_msg("STEP: no compounds plotted; skipping node_response.json.")
        return(invisible(NULL))
    }
    log_msg("STEP: writing node_response.json (adds the '", plot_column_name,
            "' column)")
    response_path <- lcmsPlot:::.cd_node_write_plot_column(
        node_args = node_args,
        id_to_value = id_to_value,
        column_name = plot_column_name,
        renderer = plot_renderer,
        position_after = plot_position
    )
    if (is.null(response_path)) {
        log_msg("  no ExpectedResponsePath in node_args; skipped (standalone).")
    } else {
        log_msg("  node_response.json written to: ", response_path)
    }
}

# Render one compound's chromatograms (faceted by sample) to `out_path`, reusing
# the already-extracted data in `p` by filtering its stored datasets to one
# compound (by its unique `name`). Reaches into lcmsPlot's internal slots.
render_compound_plot <- function(p, nm, out_path, title, width, height, dpi) {
    fids <- p@data@feature_metadata$feature_metadata_id[
        p@data@feature_metadata$name == nm]

    pc <- p
    pc@data@chromatograms <- dplyr::filter(
        p@data@chromatograms, feature_metadata_id %in% fids)
    pc@data@feature_metadata <- dplyr::filter(
        p@data@feature_metadata, name == nm)
    pc@data@detected_peaks <- dplyr::filter(
        p@data@detected_peaks, name == nm)

    pc@options$grid <- NULL
    pc@options$facets <- list(
        facets = "sample_id", ncol = NULL, nrow = NULL,
        free_x = TRUE, free_y = FALSE)
    pc@options$labels <- list(title = title, legend = "Sample")

    ggplot2::ggsave(
        filename = out_path,
        plot = pc + lp_get_plot(),
        width = width,
        height = height,
        dpi = dpi
    )
}

# Make a compound label safe for use as a file name.
.sanitize_filename <- function(x) {
    x <- gsub("[^A-Za-z0-9._-]+", "_", x)
    x <- gsub("^_+|_+$", "", x)
    if (!nzchar(x)) "compound" else x
}

# Cell value for the Compounds plot column per `path_mode`. The plot directory
# is a sub-folder of the result-file directory.
plot_cell_value <- function(out_path, node_args, path_mode) {
    if (identical(path_mode, "filename")) {
        return(basename(out_path))
    }
    if (identical(path_mode, "relative")) {
        return(file.path(basename(dirname(out_path)), basename(out_path)))
    }
    normalizePath(out_path, winslash = "\\", mustWork = FALSE)
}

# node_args.json path: the 6th commandArgs() value in CD, or the first non-flag
# trailing argument / LCMSPLOT_NODE_ARGS when standalone.
resolve_node_args <- function() {
    full_args <- commandArgs(trailingOnly = FALSE)
    trailing <- commandArgs(trailingOnly = TRUE)

    candidate <- NA_character_
    if (length(full_args) >= 6 && file.exists(full_args[6]) &&
        grepl("\\.json$", full_args[6], ignore.case = TRUE)) {
        candidate <- full_args[6]
    } else if (any(!grepl("^--", trailing) & nzchar(trailing))) {
        candidate <- trailing[!grepl("^--", trailing) & nzchar(trailing)][1]
    } else if (nzchar(Sys.getenv("LCMSPLOT_NODE_ARGS"))) {
        candidate <- Sys.getenv("LCMSPLOT_NODE_ARGS")
    }

    if (is.na(candidate) || !nzchar(candidate)) {
        stop(
            "Could not determine the node_args.json path (the 6th command-line ",
            "argument in Compound Discoverer, or the first argument standalone).",
            call. = FALSE
        )
    }
    if (!file.exists(candidate)) {
        stop("node_args.json file not found: ", candidate, call. = FALSE)
    }
    normalizePath(candidate)
}

# Raw sample paths from CD's exported Input Files (study-files) table.
build_sample_paths <- function(node_args) {
    tables <- lcmsPlot:::.cd_node_read_tables(node_args)
    tabs <- lcmsPlot:::.cd_node_classify_tables(tables)
    study_files <- tabs$study_files

    if (is.null(study_files)) {
        stop(
            "No Input Files (study-files) table was exported, so raw file paths ",
            "cannot be resolved. Request the Input Files table with a File Name ",
            "column.",
            call. = FALSE
        )
    }

    fname_col <- lcmsPlot:::first_matching_column(
        study_files, lcmsPlot:::.cd_node_cols$file_name)
    if (is.na(fname_col)) {
        stop(
            "The Input Files table has no recognised file-name column. Expected ",
            "one of: ",
            paste(lcmsPlot:::.cd_node_cols$file_name, collapse = ", "),
            call. = FALSE
        )
    }

    sample_paths <- study_files[[fname_col]]
    missing <- !file.exists(sample_paths)
    if (any(missing)) {
        stop(
            "The following raw file(s) referenced by the Input Files table were ",
            "not found:\n", paste0("  - ", sample_paths[missing], collapse = "\n"),
            call. = FALSE
        )
    }
    sample_paths
}

# Compound Discoverer result-file directory (persistent), or NA if unavailable.
resolve_result_dir <- function(node_args) {
    result_path <- tryCatch({
        info <- jsonlite::fromJSON(node_args, simplifyVector = TRUE)
        rp <- info$ResultFilePath
        if (!is.null(rp) && length(rp) == 1 && nzchar(rp)) rp else NA_character_
    }, error = function(e) NA_character_)

    if (!is.na(result_path) && dir.exists(dirname(result_path))) {
        return(dirname(result_path))
    }
    NA_character_
}

# Directory for the per-compound PNGs. CD deletes its per-run scratch dir, so
# plots go next to the result file. Priority: `output_dir` >
# <result-dir>/<result>_lcmsPlot_plots > dir(node_args)/lcmsPlot_plots.
resolve_plot_dir <- function(node_args, output_dir) {
    if (nzchar(output_dir)) {
        return(output_dir)
    }

    result_dir <- resolve_result_dir(node_args)
    if (!is.na(result_dir)) {
        info <- jsonlite::fromJSON(node_args, simplifyVector = TRUE)
        base <- tools::file_path_sans_ext(basename(info$ResultFilePath))
        return(file.path(result_dir, paste0(base, "_lcmsPlot_plots")))
    }

    log_err(
        "WARNING: could not resolve a persistent output directory from ",
        "ResultFilePath; writing plots beside node_args.json (may be removed by ",
        "Compound Discoverer after the run). Set the 'Output Directory' ",
        "parameter for a fixed path."
    )
    file.path(dirname(node_args), "lcmsPlot_plots")
}

# -----------------------------------------------------------------------------
# Startup diagnostics + run.
# -----------------------------------------------------------------------------
log_msg("=== lcmsPlot CD node script started ===")
log_msg("R version: ", R.version.string)
log_msg("Working directory: ", getwd())
log_msg("Script file: ", .this_file)
log_msg("Log file: ",
        if (nzchar(.log_file)) .log_file else "(file logging disabled)")
log_msg("commandArgs(FALSE): ",
        paste(commandArgs(trailingOnly = FALSE), collapse = " "))

tryCatch(
    withCallingHandlers(
        main(),
        error = function(e) {
            log_err("Call stack at error:")
            calls <- sys.calls()
            for (i in seq_along(calls)) {
                .log_write(paste0(
                    "  [", i, "] ",
                    paste(deparse(calls[[i]]), collapse = " ")
                ))
            }
        }
    ),
    error = function(e) {
        log_err("FATAL: lcmsPlot CD node script failed: ", conditionMessage(e))
        log_err("=== script ended with error (exit 1) ===")
        quit(status = 1, save = "no")
    }
)
log_msg("=== script completed successfully (exit 0) ===")
