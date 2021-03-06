% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/output_processing.R
\name{plot_model}
\alias{plot_model}
\title{Plot model fits}
\usage{
plot_model(model_output, type = "biphasic", detection_threshold = 20,
  textsize = 9, pointsize = 1, linesize = 0.5, facet_col = NULL)
}
\arguments{
\item{model_output}{output from model fitting using ushr().}

\item{type}{character string indicating whether the biphasic or single phase fits should be plotted. Must be either "biphasic", "single", or "triphasic". Defaults to "biphasic".}

\item{detection_threshold}{numeric value indicating the detection threshold of the assay used to measure viral load. Default value is 20.}

\item{textsize}{numeric value for base text size in ggplot. Default is 9.}

\item{pointsize}{numeric value for point size in ggplot. Default is 1.}

\item{linesize}{numeric value for line width in ggplot. Default is 0.5.}

\item{facet_col}{numeric value for number of columns to use when faceting subject panels. Defaults to NULL (i.e. ggplot default).}
}
\description{
This function plots the output from model fitting.
}
\examples{

set.seed(1234567)

simulated_data <- simulate_data(nsubjects = 20)

model_output <- ushr(data = simulated_data)

plot_model(model_output, type = "biphasic")

}
