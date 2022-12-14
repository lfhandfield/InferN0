% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wrap_functions.R
\name{infern0ZeroCorrectedWilcox2}
\alias{infern0ZeroCorrectedWilcox2}
\title{Compute Wilcox test adjusted for missing data (Hypergeometric-Normal mixture for rank sum test)}
\usage{
infern0ZeroCorrectedWilcox2(input, list.pos, list.neg, total_count,
  nb.threads = 4, list.partition = c(), nbpart = 4,
  use.clusterID = F, print.progress = T, do.quartile.average = T,
  gene.use = c(), do.hypergeo.correct = T, do.plotZscores = F,
  do.downsample = F, do.output.overwritematrix = F,
  do.fdr.correct = T)
}
\arguments{
\item{input.matrix:}{InferN0 Scope with transposed matrix...}

\item{list.pos:}{list of rows for Negative Partition, or T/F vector matching dimention}

\item{list.neg:}{list of rows for Positive Partition, or T/F vector matching dimention}

\item{total_count:}{cells read depth, used to partition into quartiles and normalize fold change}
}
\description{
Compute Wilcox test adjusted for missing data (Hypergeometric-Normal mixture for rank sum test)
}
