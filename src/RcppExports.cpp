// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// thinCloud
LogicalVector thinCloud(NumericMatrix& las, double voxel);
RcppExport SEXP _TreeLS_thinCloud(SEXP lasSEXP, SEXP voxelSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix& >::type las(lasSEXP);
    Rcpp::traits::input_parameter< double >::type voxel(voxelSEXP);
    rcpp_result_gen = Rcpp::wrap(thinCloud(las, voxel));
    return rcpp_result_gen;
END_RCPP
}
// getCircle
List getCircle(NumericMatrix& las, double pixel, double rad_max, double min_den, unsigned int min_votes);
RcppExport SEXP _TreeLS_getCircle(SEXP lasSEXP, SEXP pixelSEXP, SEXP rad_maxSEXP, SEXP min_denSEXP, SEXP min_votesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix& >::type las(lasSEXP);
    Rcpp::traits::input_parameter< double >::type pixel(pixelSEXP);
    Rcpp::traits::input_parameter< double >::type rad_max(rad_maxSEXP);
    Rcpp::traits::input_parameter< double >::type min_den(min_denSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type min_votes(min_votesSEXP);
    rcpp_result_gen = Rcpp::wrap(getCircle(las, pixel, rad_max, min_den, min_votes));
    return rcpp_result_gen;
END_RCPP
}
// singleStack
List singleStack(NumericMatrix& las, double pixel, double rad_max, double min_den, unsigned int min_votes);
RcppExport SEXP _TreeLS_singleStack(SEXP lasSEXP, SEXP pixelSEXP, SEXP rad_maxSEXP, SEXP min_denSEXP, SEXP min_votesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix& >::type las(lasSEXP);
    Rcpp::traits::input_parameter< double >::type pixel(pixelSEXP);
    Rcpp::traits::input_parameter< double >::type rad_max(rad_maxSEXP);
    Rcpp::traits::input_parameter< double >::type min_den(min_denSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type min_votes(min_votesSEXP);
    rcpp_result_gen = Rcpp::wrap(singleStack(las, pixel, rad_max, min_den, min_votes));
    return rcpp_result_gen;
END_RCPP
}
// stackMap
List stackMap(NumericMatrix& las, double hmin, double hmax, double hstep, double pixel, double rad_max, double min_den, unsigned int min_votes);
RcppExport SEXP _TreeLS_stackMap(SEXP lasSEXP, SEXP hminSEXP, SEXP hmaxSEXP, SEXP hstepSEXP, SEXP pixelSEXP, SEXP rad_maxSEXP, SEXP min_denSEXP, SEXP min_votesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix& >::type las(lasSEXP);
    Rcpp::traits::input_parameter< double >::type hmin(hminSEXP);
    Rcpp::traits::input_parameter< double >::type hmax(hmaxSEXP);
    Rcpp::traits::input_parameter< double >::type hstep(hstepSEXP);
    Rcpp::traits::input_parameter< double >::type pixel(pixelSEXP);
    Rcpp::traits::input_parameter< double >::type rad_max(rad_maxSEXP);
    Rcpp::traits::input_parameter< double >::type min_den(min_denSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type min_votes(min_votesSEXP);
    rcpp_result_gen = Rcpp::wrap(stackMap(las, hmin, hmax, hstep, pixel, rad_max, min_den, min_votes));
    return rcpp_result_gen;
END_RCPP
}
// houghStemPoints
LogicalVector houghStemPoints(NumericMatrix& las, double h1, double h2, double hstep, double radius, double pixel, double density, unsigned int votes);
RcppExport SEXP _TreeLS_houghStemPoints(SEXP lasSEXP, SEXP h1SEXP, SEXP h2SEXP, SEXP hstepSEXP, SEXP radiusSEXP, SEXP pixelSEXP, SEXP densitySEXP, SEXP votesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix& >::type las(lasSEXP);
    Rcpp::traits::input_parameter< double >::type h1(h1SEXP);
    Rcpp::traits::input_parameter< double >::type h2(h2SEXP);
    Rcpp::traits::input_parameter< double >::type hstep(hstepSEXP);
    Rcpp::traits::input_parameter< double >::type radius(radiusSEXP);
    Rcpp::traits::input_parameter< double >::type pixel(pixelSEXP);
    Rcpp::traits::input_parameter< double >::type density(densitySEXP);
    Rcpp::traits::input_parameter< unsigned int >::type votes(votesSEXP);
    rcpp_result_gen = Rcpp::wrap(houghStemPoints(las, h1, h2, hstep, radius, pixel, density, votes));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_TreeLS_thinCloud", (DL_FUNC) &_TreeLS_thinCloud, 2},
    {"_TreeLS_getCircle", (DL_FUNC) &_TreeLS_getCircle, 5},
    {"_TreeLS_singleStack", (DL_FUNC) &_TreeLS_singleStack, 5},
    {"_TreeLS_stackMap", (DL_FUNC) &_TreeLS_stackMap, 8},
    {"_TreeLS_houghStemPoints", (DL_FUNC) &_TreeLS_houghStemPoints, 8},
    {NULL, NULL, 0}
};

RcppExport void R_init_TreeLS(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
