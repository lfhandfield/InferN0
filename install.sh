#!/bin/bash

R -e 'install.packages(c("Rcpp", "RcppArmadillo","ggplot2"))
R CMD INSTALL ./
