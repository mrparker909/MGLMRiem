% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mglm_spd.R
\name{mglm_spd}
\alias{mglm_spd}
\title{mglm_spd}
\usage{
mglm_spd(X, Y, maxiter = 500, pKarcher = F, enableCheckpoint = F,
  checkpointPath = "./", Memory = 0, tol = 1e-10)
}
\arguments{
\item{X}{X is a dimX by N matrix of column vectors, each row representing a covariate, each column representing an observation.}

\item{Y}{Y is a p by p by N array of p by p symmetric positive definite response matrices.}

\item{maxiter}{maxiter is the maximum number of iterations before the optimization algorithm stops searching for an optimum. If the algorithm stops before reaching maxiters, then the "converged" variable will be set to TRUE, otherwise it will be set to FALSE.}

\item{pKarcher}{if TRUE, the estimated base point p will start at the Karcher mean of the observed SPD matrices Y, otherwise the starting point will be the identity matrix.}

\item{enableCheckpoint}{if TRUE, will create a checkpoint file at the end of each iteration. The checkpoint file may be loaded into R using load(checkpoint.rda), and then mglm_spd_checkpoint(checkpoint) can be run to continue running where the MGLM algorithm left off.}

\item{checkpointPath}{path to write checkpoint.rda file (if enableCheckpoint=TRUE).}

\item{Memory}{Memory is the maximum length of the returned list. This can be useful when the number of iterations needed for convergence is very large. Memory==0 gives 'unlimited' list length.}

\item{tol}{tolerance for threshold below which gradient is determined to be zero (stopping condition for algorithm).}
}
\value{
returns a named list containing the following elements: p (the estimated base point on the manifold), V (the set of estimated covariate coefficient tangent vectors), E (the value of the objective function, which is the sum of squared geodesic error, at each iteration), Yhat (the fitted response values), gnorm (the norm of the gradient at each iteration), converged (a flag indicating whether the algorithm converged before maxiter was reached), MGLMsteps (number of iterations taken by the algorithm).
}
\description{
This is the MGLM algorithm for performing Riemann Regression on SPD manifolds, developed by: Kim, H. J., Adluru, N., Collins, M. D., Chung, M. K., Bendin, B. B., Johnson, S. C., … Singh, V. (2014). Multivariate General Linear Models (MGLM) on Riemannian Manifolds with Applications to Statistical Analysis of Diffusion Weighted Images. 2014 IEEE Conference on Computer Vision and Pattern Recognition, 2705–2712. https://doi.org/10.1109/CVPR.2014.352. This algorithm was adapted from the publicly available matlab code to R by Matthew RP Parker; blame for errors in the transcription go to Matthew RP Parker 2019. Any issues arising from additions and modifications to the code are also the fault of Matthew RP Parker.
}
