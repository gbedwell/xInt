% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/classic_prop_tests.R
\name{classic_prop_tests}
\alias{classic_prop_tests}
\title{Classic Proportion Comparisons}
\usage{
classic_prop_tests(xint.obj, comparison = c("all", "pooled"), p.adj = "BH")
}
\arguments{
\item{xint.obj}{The xInt object containing the relevant data.}

\item{comparison}{The extent of the pairwise comparisons. \if{html}{\out{<code>}}all\if{html}{\out{</code>}} will compare every dataset to all other datasets. \if{html}{\out{<code>}}pooled\if{html}{\out{</code>}} pools the information by condition and performs pairwise comparisons across conditions.}

\item{p.adj}{The method by which to correct p-values for multiple comparisons. Can be NULL. Defaults to "BH". See p.adjust() documentation for more details.}
}
\value{
A matrix of p-values and effect sizes for each pairwise comparison performed.
}
\description{
Performs global pairwise comparisons on an xInt object using classic tests for equal proportions.
}
