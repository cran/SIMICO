\name{gammatheta}
\alias{gammatheta}
\title{
d^2/d_{gamma_k}d{theta_k}
}
\description{
Calculates the Information matrix term of I_{eta gamma} for outcome k.
}
\usage{
gammatheta(l, HL_array, HR_array, tpos_all, obs_all, apply_diffs,
   temp_beta, A_i, xDats, no_l_all, gMat, a1, a2, d)
}
\arguments{
   \item{l}{Index of first outcome of interest.}
  \item{HL_array}{n x k matrix containing all the hazard values for the left times.}
  \item{HR_array}{n x k matrix containing all the hazard values for the right times.}
  \item{tpos_all}{n x k matrix containing a indictor for whether that time is left-censored or not.}
  \item{obs_all}{n x k matrix containing a indictor for whether that time is right-censored or not.}
    \item{apply_diffs}{Matrix containing the differences between survival functions of the left and right time intervals.}
  \item{temp_beta}{vector of fitted coefficients.}
  \item{A_i}{Product of apply_diffs across all outcomes k summed over all quadrature nodes d.}
  \item{xDats}{List of design matrices for all outcomes.}
  \item{no_l_all}{n x (k - 1) matrix containing the product of apply_diffs across all outcomes k excluding the current outcome l.}
  \item{gMat}{n x q matrix of genetic information.}
  \item{a1}{First shape parameter of beta parameter.}
  \item{a2}{Second shape parameter of beta parameter.}
  \item{d}{Number of quadrature nodes.}

}
\value{
The output is a matrix containing the component of the information matrix of the gamma and theta parameters for outcome l.
}
