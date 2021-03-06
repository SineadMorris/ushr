% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulate_data.R
\name{simulate_time}
\alias{simulate_time}
\title{Simulate timepoints for subjects}
\usage{
simulate_time(npoints, censortime, id, index, max_datapoints)
}
\arguments{
\item{npoints}{numeric value indicating the number of observations to be sampled.}

\item{censortime}{numeric value indicating the maximum time point to include in the analysis.}

\item{id}{subject id. Can be numeric or a character.}

\item{index}{numeric identifier for each subject/model combination.}

\item{max_datapoints}{numeric value indicating the maximum number of data points collected from any subject.}
}
\description{
This function chooses the correct function for sampling observation times.
}
