% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pbartCV.R
\name{pbartCV}
\alias{pbartCV}
\title{cross-validated pbart}
\usage{
pbartCV(x, y, x.test, y.test, ntrees = c(50, 200), ks = c(0.5, 2, 5),
  powers = c(2, 3, 5, 10, nrow(x)), k_folds = 5, verbose = F)
}
\value{
predicted label and test AUC: ypred,Test_AUC,best_params
}
\description{
cross-validated pbart
}
