## graphics.R

#' Discrete scale of different point shapes that repeats
#' 
#' @param ... passed through to underlying \code{ggplot2::scale_...()} functions
#' 
#' @export
scale_shape_repeater <- function(...) {
	
	.shapes <- function(n) {
		allshapes <- c(0:10, 15:18)
		maxn <- length(allshapes)
		#function(n) {
		if (n > maxn)
			shapes <- allshapes[ rep(seq_len(maxn), length.out = n) ]
		else
			shapes <- allshapes[ seq_len(n) ]
		return(shapes)
		#}
	}
	
	ggplot2::discrete_scale("shape", "shape_repeater", .shapes, ...)
	
}

#' Combined shape and colour scale for uniquely identifying many groups of points
#' 
#' @param ... passed through to underlying \code{ggplot2::scale_...()} functions
#' @param pal name of an \code{RColorBrewer} palette to use 
#' @details This follows the usual convention on PCA plots of denoting populations by a shape-colour combination,
#'   so that similar colours are more easily discriminated when there are many (>10) groups of points on the plot.
#' 
#' @export
scale_combo <- function(..., pal = "Spectral") {
	
	.cols <- function(n) {
		#function(n) {
		maxn <- RColorBrewer::brewer.pal.info[pal,"maxcolors"]-2
		if (maxn < n)
			cols <- RColorBrewer::brewer.pal(maxn, pal)[ rep(seq_len(maxn), length.out = n) ]
		else
			cols <- RColorBrewer::brewer.pal(n, pal)
		return(cols)
		#}
	}
	
	
	list( ggplot2::discrete_scale("colour", "my_colour", .cols, ...),
		  scale_shape_repeater(...) )
	
}