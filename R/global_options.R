#' Global Options in pixelatorR
#'
#' Information about global options provided in pixelatorR.
#'
#' Available options:
#' \itemize{
#'    \item{\code{pixelatorR.arrow_outdir} The directory where the edgelists are stored.}
#'    \item{\code{pixelatorR.arrowdir_maxsize} The maximum allowed size of the edgelists directory.}
#'    \item{\code{pixelatorR.auto_cleanup} Whether to run \code{edgelists_directories_clean()}
#'    automatically when the size of the edgelists directory exceeds \code{pixelatorR.arrowdir_maxsize}.
#'    The cleanup is triggered by functions which create new edgelist copies, such as \code{subset},
#'    \code{merge} and \code{RenameCells}.}
#'    \item{\code{pixelatorR.verbose} Whether to print messages to the console.}
#' }
#'
#' @section Verbosity:
#' Verbosity can be turned off globally for all pixelatorR function by setting
#' \code{options(pixelatorR.verbose = FALSE)}.
#'
#' @section Edgelists stored on disk:
#' When creating Seurat objects with pixelatorR, the graph data will be stored on disk
#' in a directory that can be set with the global option \code{pixelatorR.arrow_outdir}.
#' By default, this directory is set to './edgelists' which will be created in the
#' current working directory.
#'
#' Certain operations, such as \code{subset}, \code{merge}
#' and \code{RenameCells}, will create copies of the graph data in new sub directories within
#' this folder to make sure that each Seurat object that is created has it's own
#' edgelist file(s) associated with it.
#'
#' A downside with this strategy is that
#' the size of the directory can grow very large if many Seurat objects are created
#' or if the functions mentioned above are called many times. When Seurat
#' objects are overwritten or removed, some directories may no longer be linked to any
#' global variable in the current R session, but these directories
#' will not be removed automatically. pixelatorR offers a utility function called
#' \code{edgelists_directories_clean()} to remove them.
#'
#' To make the user aware of this behavior, a message will be displayed once the
#' first time \code{ReadMPX_Seurat} that creates a new edgelist directory is called.
#'
#' Moreover, when the size of the directory exceeds the size limit set with
#' the global option \code{pixelatorR.arrowdir_maxsize} (which is set to 5GB by default),
#' \code{edgelists_directories_clean()} will be run automatically. Users will be given
#' a prompt to confirm the cleanup before it is executed. This behavior can be turned off
#' by setting \code{options(pixelatorR.auto_cleanup = TRUE)}. You can also run
#' \code{edgelist_directories_du()} to get the current size of the directory.
#'
#' @seealso [edgelist_directories_clean()] to clean up the edgelist directory,
#' [edgelist_directories_du()] to get the current size of the edgelist directory.
#'
#' @examples
#' library(dplyr)
#'
#' # Set edgelist directory to a temporary folder
#' options(pixelatorR.arrow_outdir = tempdir())
#'
#' # Check the current size of the edgelist directory
#' edgelist_directories_du()
#'
#' # List files in the edgelist directory
#' getOption("pixelatorR.arrow_outdir") %>% fs::dir_ls()
#'
#' # Increase maximum allowed size of edgelist directory to 5GB
#' # This will generate a prompt asking users to cleanup the edgelist directory
#' # if the total size exceeds 50GB
#' options(pixelatorR.arrowdir_maxsize = fs::fs_bytes(5 * 1024^3))
#'
#' # Turn of prompts for automatic cleanup
#' options(pixelatorR.auto_cleanup = TRUE)
#'
#' # Turn off verbosity
#' options(pixelatorR.verbose = FALSE)
#'
#' @name pixelatorR_options
#' @rdname pixelatorR_options
NULL
