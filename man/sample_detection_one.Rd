% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sample_detection_one.R
\name{sample_detection_one}
\alias{sample_detection_one}
\title{Probability of detection one}
\usage{
sample_detection_one(J = 50, K = 8, theta = 0.06,
  p_detection = 0.3)
}
\arguments{
\item{J}{number of samples per sampleing event}

\item{K}{number of molecular replicates}

\item{theta}{sample-level capture probability}

\item{p_detection}{molecular-level detection probability}
}
\description{
This function calculates the probability
of detecting eDNA in at least one sample
assuming eDNA is present at a site.
}
\examples{
sample_detection_one()


}
