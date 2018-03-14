# util.R
# Miscellaneous utility functions that don't fit elsewhere

#' A `dplyr`-friendly equivalent to the old `plyr::ldply()`
#' 
#' @param ll a list of objects to which to apply a function
#' @param fn the function to apply
#' @param ... extra arguments passed to \code{fn}
#' @param .id if not \code{NULL}, name of column to use for group IDs (inferred from \code{attr(ll, "names")})
#' @return a \code{tibble} with concatenated results of applying \code{fn()} to each element of \code{ll}
#' @details Hard work of combining results is done under the hood by \code{dplyr::bind_rows()}; consult that
#'     function's documentation to see how incomplete rows, missing columns, conflicting data types, etc. are handled.
#' 
#' @export
myldply <- function(ll, fn, ..., .id = NULL) {
	rez <- lapply(ll, fn, ...)
	dplyr::bind_rows(rez, .id = .id)
}