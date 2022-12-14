% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/INFRN_IO.R
\name{infern0ExportDataAsMtx}
\alias{infern0ExportDataAsMtx}
\title{Recover Raw Counts}
\usage{
infern0ExportDataAsMtx(infrn0, pathprefix, matname = c(),
  matrix.filename = "matrix.mtx", gene_filename = "genes.tsv",
  cell_filename = "barcodes.tsv", use.normalized = F, do.gzip = F)
}
\arguments{
\item{path}{or list of paths to tab delimited text table}

\item{pathprefix:}{output prefix for the 3 output files}
}
\description{
Recovers the 'raw.data' using 'data' assuming cells had their expression normalized (no transcript-length/gene specific normalization)
}
