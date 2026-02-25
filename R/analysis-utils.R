#' Save a ggplot2 plot to one or more file formats with options for size and directory creation.
#'
#' This function is a wrapper around `ggsave()` that allows for saving the same plot in multiple formats and handles
#' directory creation if needed. It also includes an option to skip saving if the file already exists, based on a
#' global option.
#'
#' @param filename Character. The path and base file name for saving the plot. If a file extension is included, it will
#' override `file_formats` and save only in that format.
#' @param plot ggplot object. The plot to save. Defaults to the last plot produced (`last_plot()`).
#' @param width Numeric. The width of the plot in inches. Default is 10.
#' @param height Numeric. The height of the plot in inches. Default is 10.
#' @param create_dir Logical. Whether to create the directory if it does not exist. Default is TRUE.
#' @param file_formats Character vector. File formats to save to (e.g., "png", "pdf"). Default is c("png", "pdf").
#' @param overwrite Logical. Whether to overwrite existing files. Default is taken from the global option
#' "export_plot.overwrite" (default FALSE).
#'
#' @return Invisibly returns NULL. Used for its side effect of writing plot files.
#' @examples
#' \dontrun{
#' export_plot("results/myplot", plot = myplot, width = 8, height = 6)
#' }
#' # Change overwrite behavior globally
#' options(export_plot.overwrite = TRUE)
#'
#' # Now existing files will be overwritten without needing to specify overwrite in the function call
#' \dontrun{
#' export_plot("results/myplot", plot = myplot, width = 8, height = 6)
#' }
#'
#' @export
#'
export_plot <-
  function(
    filename,
    plot = ggplot2::last_plot(),
    width = 10,
    height = 10,
    create_dir = TRUE,
    file_formats = c("png", "pdf"),
    overwrite = getOption("export_plot.overwrite", FALSE)
  ) {
    assert_single_value(filename, type = "string")
    assert_class(plot, "ggplot")
    assert_single_value(width, type = "numeric")
    assert_single_value(height, type = "numeric")
    assert_single_value(create_dir, type = "bool")
    assert_vector(file_formats, type = "character", n = 1)
    assert_single_value(overwrite, type = "bool")

    ext <- fs::path_ext(filename)

    if (nchar(ext) > 0) {
      file_formats <- ext
      filename <- fs::path_ext_remove(filename)
    }

    for (format in file_formats) {
      out_file <- paste0(filename, ".", format)

      if (file.exists(out_file) && !overwrite) {
        rlang::warn(glue::glue("Skipping existing file: {out_file}"))
        next
      }

      ggsave(
        filename = out_file,
        plot = plot,
        width = width,
        height = height,
        create.dir = create_dir,
        limitsize = FALSE
      )
    }
    invisible(NULL)
  }
