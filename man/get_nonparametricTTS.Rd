% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_TTS.R
\name{get_nonparametricTTS}
\alias{get_nonparametricTTS}
\title{Non-parametric TTS function}
\usage{
get_nonparametricTTS(vl, suppression_threshold, time, npoints)
}
\arguments{
\item{vl}{numeric vector of viral load measurements.}

\item{suppression_threshold}{numeric value for the suppression threshold: measurements below this value will be assumed to represent viral suppression. Typically this would be the detection threshold of the assay. Default value is 20.}

\item{time}{numeric vector indicating the time when vl measurements were taken.}

\item{npoints}{numeric value indicating the number of interpolation points to be considered.}
}
\description{
This function computes the non-parametric form of the time to suppression
}
