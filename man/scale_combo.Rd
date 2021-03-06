% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/graphics.R
\name{scale_combo}
\alias{scale_combo}
\title{Combined shape and colour scale for uniquely identifying many groups of points}
\usage{
scale_combo(..., pal = "Spectral")
}
\arguments{
\item{...}{passed through to underlying \code{ggplot2::scale_...()} functions}

\item{pal}{name of an \code{RColorBrewer} palette to use}
}
\description{
Combined shape and colour scale for uniquely identifying many groups of points
}
\details{
This follows the usual convention on PCA plots of denoting populations by a shape-colour combination,
  so that similar colours are more easily discriminated when there are many (>10) groups of points on the plot.
}
