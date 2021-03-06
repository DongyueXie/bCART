% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/adaBARTr.R
\name{adaBARTr}
\alias{adaBARTr}
\title{Fit an adaptive bayesian additive regression tree(adaBART) model for classification}
\usage{
adaBARTr(X, y, x.test, cutoff = 0.5, k = 2, binaryOffset = NULL,
  power = 2, base = 0.95, p_split = "CGM", r = 2, ntree = 50,
  ndpost = 700, nskip = 300, Tmin = 2, printevery = 100,
  p_modify = c(0.5, 0.5, 0), save_trees = F, rule = "bart",
  pre_train = F, n_pre_train = 100)
}
\arguments{
\item{all}{see ?pBARTr}
}
\description{
Fit an adaptive bayesian additive regression tree(adaBART) model for classification
}
