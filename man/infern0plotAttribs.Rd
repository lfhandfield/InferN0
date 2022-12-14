% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/INFRN_Plot.R
\name{infern0plotAttribs}
\alias{infern0plotAttribs}
\title{Function used to easily populate ggplot related options stored in a list, that can be used by ploting function "Overlay" as plot.attrib argument be applied on a input ggplot.}
\usage{
infern0plotAttribs(..., ggplot.or.baseAttribs = c(), baseAttribs = c(),
  no.legend = F, no.jitter = F, no.xticks = F, no.xlabel = F,
  no.yticks = F, no.ylabel = F, title = c(), xtitle = c(),
  ytitle = c(), xlabel = c(), ylabel = c(), xnames.bold = F,
  xnames.color = c(), xnames.size = c(), ynames.bold = F,
  ynames.color = c(), ynames.size = c(), legnames.size = c(),
  no.xnamedticks = F, no.ynamedticks = F, nb.col.legend = c())
}
\arguments{
\item{ggplot.or.plot.attrib:}{Either a ggplot object or a default list of attributes that is to be modified (effectively equivalent to setting the "baseAttribs" instead). If a ggplot is provided, outputs an ggplot object with the applied modifications instead.}

\item{baseAttribs:}{list of attributes that is to be modified}
}
\description{
Function used to easily populate ggplot related options stored in a list, that can be used by ploting function "Overlay" as plot.attrib argument be applied on a input ggplot.
}
