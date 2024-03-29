% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/treemix.R
\name{plot_treemix}
\alias{plot_treemix}
\title{Plot population tree}
\usage{
plot_treemix(
  obj,
  plot.nodes = TRUE,
  plot.migration = TRUE,
  branch.colour = "grey",
  branch.width = 0.5,
  label = TRUE,
  ...
)
}
\arguments{
\item{obj}{pieces of a TreeMix run from \code{read_treemix()}}

\item{plot.nodes}{show dots for internal and leaf nodes of population tree}

\item{plot.migration}{show migration edges as curved arrows}

\item{branch.colour}{colour to use when drawing edges of population tree}

\item{branch.width}{weight of drift edges of population tree}
}
\value{
A \code{ggplot2} plot, with dataframes fo tips and edges of tree returned as attributes \code{attr(,"tips")}
    and \code{attr(,"edges")} respectively
}
\description{
Plot population tree
}
