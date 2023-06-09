% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bionickomics_utils.r
\name{get_ClinVar_variants}
\alias{get_ClinVar_variants}
\title{Get ClinVar Variants}
\usage{
get_ClinVar_variants(gene_of_interest)
}
\arguments{
\item{gene_of_interest}{}
}
\value{
a dataframe of tidied NCBI ClinVar pathological variants specific to the gene of interest
}
\description{
\code{get_ClinVar_variants()} is used to obtain the latest pathological variants for a gene of interest from the NCBI ClinVar database.

This wrapper function downloads the latest ClinVar information, filters on the \code{gene_of_interest}, which is then subsequently formatted for lolliplotR processing.
}
\details{
Variants will be cached for future access.
}
\examples{
get_ClinVar_variants("CRB1")
get_ClinVar_variants("SMN1")
}
