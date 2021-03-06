// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// geno_cor_new
RcppExport SEXP geno_cor_new(SEXP YY, SEXP GDD, SEXP ww, SEXP coorr, SEXP mm, SEXP KK);
RcppExport SEXP _HDGENE_geno_cor_new(SEXP YYSEXP, SEXP GDDSEXP, SEXP wwSEXP, SEXP coorrSEXP, SEXP mmSEXP, SEXP KKSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type YY(YYSEXP);
    Rcpp::traits::input_parameter< SEXP >::type GDD(GDDSEXP);
    Rcpp::traits::input_parameter< SEXP >::type ww(wwSEXP);
    Rcpp::traits::input_parameter< SEXP >::type coorr(coorrSEXP);
    Rcpp::traits::input_parameter< SEXP >::type mm(mmSEXP);
    Rcpp::traits::input_parameter< SEXP >::type KK(KKSEXP);
    rcpp_result_gen = Rcpp::wrap(geno_cor_new(YY, GDD, ww, coorr, mm, KK));
    return rcpp_result_gen;
END_RCPP
}
// geno_cor_new2
RcppExport SEXP geno_cor_new2(SEXP YY, SEXP GDD, SEXP GDDD, SEXP ww, SEXP coorr, SEXP mm, SEXP KK);
RcppExport SEXP _HDGENE_geno_cor_new2(SEXP YYSEXP, SEXP GDDSEXP, SEXP GDDDSEXP, SEXP wwSEXP, SEXP coorrSEXP, SEXP mmSEXP, SEXP KKSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type YY(YYSEXP);
    Rcpp::traits::input_parameter< SEXP >::type GDD(GDDSEXP);
    Rcpp::traits::input_parameter< SEXP >::type GDDD(GDDDSEXP);
    Rcpp::traits::input_parameter< SEXP >::type ww(wwSEXP);
    Rcpp::traits::input_parameter< SEXP >::type coorr(coorrSEXP);
    Rcpp::traits::input_parameter< SEXP >::type mm(mmSEXP);
    Rcpp::traits::input_parameter< SEXP >::type KK(KKSEXP);
    rcpp_result_gen = Rcpp::wrap(geno_cor_new2(YY, GDD, GDDD, ww, coorr, mm, KK));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_HDGENE_geno_cor_new", (DL_FUNC) &_HDGENE_geno_cor_new, 6},
    {"_HDGENE_geno_cor_new2", (DL_FUNC) &_HDGENE_geno_cor_new2, 7},
    {NULL, NULL, 0}
};

RcppExport void R_init_HDGENE(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
