// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// calcNormS
arma::mat calcNormS(const List& D, int ncols);
RcppExport SEXP _hogsvdR_calcNormS(SEXP DSEXP, SEXP ncolsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const List& >::type D(DSEXP);
    Rcpp::traits::input_parameter< int >::type ncols(ncolsSEXP);
    rcpp_result_gen = Rcpp::wrap(calcNormS(D, ncols));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_hogsvdR_calcNormS", (DL_FUNC) &_hogsvdR_calcNormS, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_hogsvdR(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
