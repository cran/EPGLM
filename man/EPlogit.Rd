\name{EPlogit}
\alias{EPlogit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute the EP (expectation propagation) approximation of a logit regression.
}
\description{
The function computes the EP approximation of a logit regression with Gaussian prior. The user must specify the design matrix, the response vector
and the prior variance. For more information on the default prior variance see Chopin and Ridgway [2015].
}
\usage{
EPlogit(X, Y, s)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
Design matrix. The matrix should include a constant column if a bias is to be considered. 
}
  \item{Y}{
Response vector. The vector should take values in {0,1}. 
}
  \item{s}{
Prior variance. The prior is taken to be spherical Gaussian, the variance must therefore be specified in the form of a scalar.
For default choices see Chopin and Ridgway [2015]. 
}
}
\details{
The implementation is based on the remarks of Chopin and Ridgway (2015) and computes a Gaussian approximation to the Bayesian logit model. The approximation can serve as a very efficient estimation or as the starting point to Monte Carlo algorithms. The output value is given in the form 
of the parameters of the Gaussian approximation (mean and variance matrix) and an approximation to the log marginal likelihood. 
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

 \section{Warning }{
The current implementation does not include damping or the possibility to use fractional EP (hopefully it will in a future version).
This might results in poor performance for large datasets.  
}


\examples{
data(Pima.tr)
Y<-as.matrix(as.numeric(Pima.tr[,8]))-1% formating the response
X<-cbind(1,data.matrix(Pima.tr[,1:7]))% adding a constnat to the model
Sol<-EPlogit(X,Y,100)%EP approximation of the model Y~X with prior variance 100.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Bayesian Logit }
\keyword{ Expectation Propagation }% __ONLY ONE__ keyword per line
