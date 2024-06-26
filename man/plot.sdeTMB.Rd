% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods.R
\name{plot.sdeTMB}
\alias{plot.sdeTMB}
\title{Basic summary of objects of class 'sdeTMB'}
\usage{
\method{plot}{sdeTMB}(
  object,
  plot.obs = 1,
  pacf = FALSE,
  extended = FALSE,
  ggtheme = getggplot2theme()
)
}
\arguments{
\item{plot.obs}{a vector to indicate which observations should be plotted for. If multiple
are chosen a list of plots for each observation is returned.}

\item{pacf}{logical to indicate whether or not the partial autocorrelations should be returned.
The default is FALSE in which case a histogram is returned instead.}

\item{extended}{logical. if TRUE additional information is printed}

\item{ggtheme}{ggplot2 theme to use for creating the ggplot.}
}
\value{
A huge amount of information
}
\description{
Basic summary of objects of class 'sdeTMB'
}
