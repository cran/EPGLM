#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(BH)]]
#define BOOST_DISABLE_ASSERTS
#include "header.h"


// [[Rcpp::export]]
Rcpp::List EPprobitCxx(mat X, mat Y, double s)
{
	EPprobit ep;
	ep(X,Y,1,s);
	mat Theta=ep.Get_m();
	mat V=ep.Get_V();
	double Z=ep.Get_Z();
	return Rcpp::List::create(Rcpp::Named("m") = Theta,
		                  Rcpp::Named("V") = V,
				  Rcpp::Named("Z") = Z);
}

RCPP_MODULE(EPprobitCxx) {
	function( "EPprobitCxx", &EPprobitCxx);
}
// [[Rcpp::export]]
Rcpp::List EPlogitCxx(mat X, mat Y, double s)
{
	EPlogit ep2;
	ep2(X,Y,1,s);
	mat Theta=ep2.Get_m();
	mat V=ep2.Get_V();
	double Z=ep2.Get_Z();
	return Rcpp::List::create(Rcpp::Named("m") = Theta,
		                  Rcpp::Named("V") = V,
				  Rcpp::Named("Z") = Z);
}
RCPP_MODULE(EPlogitCxx) {
	function( "EPlogitCxx", &EPlogitCxx);
}
