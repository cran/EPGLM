\name{EPlogit}
\alias{EPlogit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute the EP approximation of a logit regression.
%%  ~~function to do ... ~~
}
\description{
The function computes the EP approximation of a logit regression with Gaussian prior.
}
\usage{
EPlogit(X, Y, s)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
Design matrix. Should include a constant column if a bias is to be considered.
%%     ~~Describe \code{X} here~~
}
  \item{Y}{
Response vector, souhld take values 0 and 1.
%%     ~~Describe \code{Y} here~~
}
  \item{s}{
Prior variance
%%     ~~Describe \code{s} here~~
}
}
\details{
The implementation is based on the remarks of Chopin and Ridgway (2015) and computes a Gaussian approximation to the Bayesian logit model. The approximation can serve as a very efficient estimation or as the starting point to Monte Carlo algorithms. 
}
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
N. Chopin and J. Ridgway. Leave Pima Indians alone: binary regression as a benchmark for Bayesian computation. arxiv:1506.08640
}
\author{
James Ridgway
}
\note{
More priors and models should be available shortly.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~


\examples{
data(Pima.tr)
Y<-as.matrix(as.numeric(Pima.tr[,8]))-1% formating the response
X<-cbind(1,data.matrix(Pima.tr[,1:7]))% adding a constnat to the model
Sol<-EPlogit(X,Y,100)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Bayesian Logit }
\keyword{ Expectation Propagation }% __ONLY ONE__ keyword per line
