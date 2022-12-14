\name{gammatheta_off}
\alias{gammatheta_off}
\title{
d^2/d_{gamma_l}d{theta_m}
}
\description{
Calculates the Information matrix term of I_{eta gamma} for outcomes l and m
}
\usage{
gammatheta_off(l,m, HL_array, HR_array, xAll, apply_diffs, temp_beta,
   A_i, no_l_all, no_two_all, gMat, a1, a2, k, d)
}
\arguments{
   \item{l}{Index of first outcome of interest.}
   \item{m}{Index of second outcome of interest.}
  \item{HL_array}{n x k matrix containing all the hazard values for the left times.}
  \item{HR_array}{n x k matrix containing all the hazard values for the right times.}
  \item{xAll}{List of design matrices and censoring terms.}
  \item{apply_diffs}{Matrix containing the differences between survival functions of the left and right time intervals.}
  \item{temp_beta}{vector of fitted coefficients.}
  \item{A_i}{Product of apply_diffs across all outcomes K summed over all quadrature nodes d.}
  \item{no_l_all}{n x (k - 1) matrix containing the product of apply_diffs across all outcomes k excluding the current outcome l.}
  \item{no_two_all}{n x (k - 2) matrix containing the product of apply_diffs across all outcomes k excluding the outcomes l and m.}
  \item{gMat}{n x q matrix of genetic information.}
  \item{a1}{First shape parameter of beta parameter.}
  \item{a2}{Second shape parameter of beta parameter.}
  \item{k}{Total number of outcomes.}
  \item{d}{Number of quadrature nodes.}
}
\value{
The output is a matrix containing the component of the information matrix of the gamma and theta parameters for outcomes l and m.
}
