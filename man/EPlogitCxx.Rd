\name{EPlogitCxx}
\alias{EPlogitCxx}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
C++ internal function to compute the EP approximation (use EPlogit instead).
%%  ~~function to do ... ~~
}
\description{
The function is the C++ internal function used by EPlogit
}
\usage{
EPlogitCxx(...)
}
\arguments{
  \item{...}{
See arguments of EPlogit
%%     ~~Describe \code{X} here~~
}
}
%- maybe also 'usage' for other objects documented here.
\value{
\item{m}{ Mean of the Gaussian approximation}
\item{V}{ Variance matrix of the Gaussian approximation}
\item{Z}{ Approximated log marginal likelihood}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
N. Chopin and J. Ridgway. Leave Pima Indians alone: binary regression as a benchmark for Bayesian computation. arxiv:1506.08640.
}
\author{
James Ridgway
}
\note{
More priors and models should be available shortly.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~
\keyword{ Bayesian Logit }
\keyword{ Expectation Propagation }% __ONLY ONE__ keyword per line
