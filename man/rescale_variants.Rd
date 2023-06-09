% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bionickomics_utils.r
\name{rescale_variants}
\alias{rescale_variants}
\title{Rescale variants based upon new intronic sizes (fixed)}
\usage{
rescale_variants(variantInfo, transcript_grange, minIntronSize)
}
\arguments{
\item{variantInfo}{\code{dataframe} of variant information for gene of interest}

\item{transcript_grange}{\code{GRanges} object of the transcript model for the gene of interest}

\item{minIntronSize}{\code{numeric} vector to fix all \code{introns} to for representation}
}
\value{
rescaled GRanges object
}
\description{
Rescale variants based upon new intronic sizes (fixed)
}
\examples{
rescale_variants(variantdataframe, GRangeObject, minIntronSize)
}
