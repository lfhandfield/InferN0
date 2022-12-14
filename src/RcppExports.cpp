// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <RcppArmadillo.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// InferN0_identifyL0Network
Rcpp::List InferN0_identifyL0Network(const Rcpp::NumericMatrix& x, const double& y);
RcppExport SEXP _InferN0_InferN0_identifyL0Network(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const double& >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(InferN0_identifyL0Network(x, y));
    return rcpp_result_gen;
END_RCPP
}
// InferN0_identifyL0NetworkGold
Rcpp::List InferN0_identifyL0NetworkGold(const Rcpp::NumericMatrix& x, const arma::mat& y, const double& z);
RcppExport SEXP _InferN0_InferN0_identifyL0NetworkGold(SEXP xSEXP, SEXP ySEXP, SEXP zSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const double& >::type z(zSEXP);
    rcpp_result_gen = Rcpp::wrap(InferN0_identifyL0NetworkGold(x, y, z));
    return rcpp_result_gen;
END_RCPP
}
// InferN0_printstringlist
Rcpp::List InferN0_printstringlist(SEXP x);
RcppExport SEXP _InferN0_InferN0_printstringlist(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(InferN0_printstringlist(x));
    return rcpp_result_gen;
END_RCPP
}
// InferN0_identifyL0NetworkGoldResurrected
Rcpp::List InferN0_identifyL0NetworkGoldResurrected(Rcpp::S4 input, Rcpp::List list, Rcpp::List traincov, Rcpp::List trainmean);
RcppExport SEXP _InferN0_InferN0_identifyL0NetworkGoldResurrected(SEXP inputSEXP, SEXP listSEXP, SEXP traincovSEXP, SEXP trainmeanSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::S4 >::type input(inputSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type list(listSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type traincov(traincovSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type trainmean(trainmeanSEXP);
    rcpp_result_gen = Rcpp::wrap(InferN0_identifyL0NetworkGoldResurrected(input, list, traincov, trainmean));
    return rcpp_result_gen;
END_RCPP
}
// infernal_wilcox_geneset
Rcpp::List infernal_wilcox_geneset(SEXP scope, Rcpp::List list);
RcppExport SEXP _InferN0_infernal_wilcox_geneset(SEXP scopeSEXP, SEXP listSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type scope(scopeSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type list(listSEXP);
    rcpp_result_gen = Rcpp::wrap(infernal_wilcox_geneset(scope, list));
    return rcpp_result_gen;
END_RCPP
}
// infernal_wilcox_scope
Rcpp::NumericMatrix infernal_wilcox_scope(SEXP scope, Rcpp::List list);
RcppExport SEXP _InferN0_infernal_wilcox_scope(SEXP scopeSEXP, SEXP listSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type scope(scopeSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type list(listSEXP);
    rcpp_result_gen = Rcpp::wrap(infernal_wilcox_scope(scope, list));
    return rcpp_result_gen;
END_RCPP
}
// infernal_wilcox_matrix
Rcpp::NumericMatrix infernal_wilcox_matrix(Rcpp::S4 input, Rcpp::List list);
RcppExport SEXP _InferN0_infernal_wilcox_matrix(SEXP inputSEXP, SEXP listSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::S4 >::type input(inputSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type list(listSEXP);
    rcpp_result_gen = Rcpp::wrap(infernal_wilcox_matrix(input, list));
    return rcpp_result_gen;
END_RCPP
}
// infernal_find_DE
Rcpp::List infernal_find_DE(SEXP scope, Rcpp::List list);
RcppExport SEXP _InferN0_infernal_find_DE(SEXP scopeSEXP, SEXP listSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type scope(scopeSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type list(listSEXP);
    rcpp_result_gen = Rcpp::wrap(infernal_find_DE(scope, list));
    return rcpp_result_gen;
END_RCPP
}
// infernal_find_markers
Rcpp::List infernal_find_markers(SEXP scope, Rcpp::List rlist);
RcppExport SEXP _InferN0_infernal_find_markers(SEXP scopeSEXP, SEXP rlistSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type scope(scopeSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type rlist(rlistSEXP);
    rcpp_result_gen = Rcpp::wrap(infernal_find_markers(scope, rlist));
    return rcpp_result_gen;
END_RCPP
}
// infernal_sharedDE
Rcpp::List infernal_sharedDE(Rcpp::S4 input);
RcppExport SEXP _InferN0_infernal_sharedDE(SEXP inputSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::S4 >::type input(inputSEXP);
    rcpp_result_gen = Rcpp::wrap(infernal_sharedDE(input));
    return rcpp_result_gen;
END_RCPP
}
// infernal_cmpInverses
Rcpp::List infernal_cmpInverses(Rcpp::NumericMatrix input);
RcppExport SEXP _InferN0_infernal_cmpInverses(SEXP inputSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type input(inputSEXP);
    rcpp_result_gen = Rcpp::wrap(infernal_cmpInverses(input));
    return rcpp_result_gen;
END_RCPP
}
// testSlot
void testSlot(Rcpp::S4 obj, SEXP slot);
RcppExport SEXP _InferN0_testSlot(SEXP objSEXP, SEXP slotSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::S4 >::type obj(objSEXP);
    Rcpp::traits::input_parameter< SEXP >::type slot(slotSEXP);
    testSlot(obj, slot);
    return R_NilValue;
END_RCPP
}
// infernal_serialize
Rcpp::RawVector infernal_serialize(Rcpp::S4 obj);
RcppExport SEXP _InferN0_infernal_serialize(SEXP objSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::S4 >::type obj(objSEXP);
    rcpp_result_gen = Rcpp::wrap(infernal_serialize(obj));
    return rcpp_result_gen;
END_RCPP
}
// infernal_deserialize
void infernal_deserialize(Rcpp::S4 obj, Rcpp::RawVector src);
RcppExport SEXP _InferN0_infernal_deserialize(SEXP objSEXP, SEXP srcSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::S4 >::type obj(objSEXP);
    Rcpp::traits::input_parameter< Rcpp::RawVector >::type src(srcSEXP);
    infernal_deserialize(obj, src);
    return R_NilValue;
END_RCPP
}
// infernal_createScope
SEXP infernal_createScope();
RcppExport SEXP _InferN0_infernal_createScope() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(infernal_createScope());
    return rcpp_result_gen;
END_RCPP
}
// infernal_readScope
SEXP infernal_readScope(SEXP path);
RcppExport SEXP _InferN0_infernal_readScope(SEXP pathSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type path(pathSEXP);
    rcpp_result_gen = Rcpp::wrap(infernal_readScope(path));
    return rcpp_result_gen;
END_RCPP
}
// infernal_readScope2
SEXP infernal_readScope2(SEXP path);
RcppExport SEXP _InferN0_infernal_readScope2(SEXP pathSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type path(pathSEXP);
    rcpp_result_gen = Rcpp::wrap(infernal_readScope2(path));
    return rcpp_result_gen;
END_RCPP
}
// infernal_writeScope
void infernal_writeScope(SEXP infscope, SEXP path);
RcppExport SEXP _InferN0_infernal_writeScope(SEXP infscopeSEXP, SEXP pathSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type infscope(infscopeSEXP);
    Rcpp::traits::input_parameter< SEXP >::type path(pathSEXP);
    infernal_writeScope(infscope, path);
    return R_NilValue;
END_RCPP
}
// infernal_writeScope2
void infernal_writeScope2(SEXP infscope, Rcpp::RawVector serialized, SEXP path);
RcppExport SEXP _InferN0_infernal_writeScope2(SEXP infscopeSEXP, SEXP serializedSEXP, SEXP pathSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type infscope(infscopeSEXP);
    Rcpp::traits::input_parameter< Rcpp::RawVector >::type serialized(serializedSEXP);
    Rcpp::traits::input_parameter< SEXP >::type path(pathSEXP);
    infernal_writeScope2(infscope, serialized, path);
    return R_NilValue;
END_RCPP
}
// infernal_readMatrixFolder
Rcpp::List infernal_readMatrixFolder(Rcpp::List list);
RcppExport SEXP _InferN0_infernal_readMatrixFolder(SEXP listSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type list(listSEXP);
    rcpp_result_gen = Rcpp::wrap(infernal_readMatrixFolder(list));
    return rcpp_result_gen;
END_RCPP
}
// infernal_convertDenseTextMatrixToSparseFormat
void infernal_convertDenseTextMatrixToSparseFormat(Rcpp::List list);
RcppExport SEXP _InferN0_infernal_convertDenseTextMatrixToSparseFormat(SEXP listSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type list(listSEXP);
    infernal_convertDenseTextMatrixToSparseFormat(list);
    return R_NilValue;
END_RCPP
}
// infernal_exportDataAsMtx
void infernal_exportDataAsMtx(SEXP infscope, Rcpp::List list);
RcppExport SEXP _InferN0_infernal_exportDataAsMtx(SEXP infscopeSEXP, SEXP listSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type infscope(infscopeSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type list(listSEXP);
    infernal_exportDataAsMtx(infscope, list);
    return R_NilValue;
END_RCPP
}
// infernal_reverseTCM
void infernal_reverseTCM(SEXP rscp, Rcpp::List list);
RcppExport SEXP _InferN0_infernal_reverseTCM(SEXP rscpSEXP, SEXP listSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type rscp(rscpSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type list(listSEXP);
    infernal_reverseTCM(rscp, list);
    return R_NilValue;
END_RCPP
}
// infernal_reverseTCM_inMTX
void infernal_reverseTCM_inMTX(Rcpp::List list);
RcppExport SEXP _InferN0_infernal_reverseTCM_inMTX(SEXP listSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type list(listSEXP);
    infernal_reverseTCM_inMTX(list);
    return R_NilValue;
END_RCPP
}
// infernal_readSparseMatrix
void infernal_readSparseMatrix(SEXP scope, Rcpp::S4 x, bool isRaw);
RcppExport SEXP _InferN0_infernal_readSparseMatrix(SEXP scopeSEXP, SEXP xSEXP, SEXP isRawSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type scope(scopeSEXP);
    Rcpp::traits::input_parameter< Rcpp::S4 >::type x(xSEXP);
    Rcpp::traits::input_parameter< bool >::type isRaw(isRawSEXP);
    infernal_readSparseMatrix(scope, x, isRaw);
    return R_NilValue;
END_RCPP
}
// infernal_readDenseMatrix
void infernal_readDenseMatrix(SEXP scope, Rcpp::NumericMatrix x, bool isRaw);
RcppExport SEXP _InferN0_infernal_readDenseMatrix(SEXP scopeSEXP, SEXP xSEXP, SEXP isRawSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type scope(scopeSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< bool >::type isRaw(isRawSEXP);
    infernal_readDenseMatrix(scope, x, isRaw);
    return R_NilValue;
END_RCPP
}
// infernal_createOutput
Rcpp::List infernal_createOutput(SEXP data);
RcppExport SEXP _InferN0_infernal_createOutput(SEXP dataSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type data(dataSEXP);
    rcpp_result_gen = Rcpp::wrap(infernal_createOutput(data));
    return rcpp_result_gen;
END_RCPP
}
// infernal_fitNegativeBinomial
Rcpp::List infernal_fitNegativeBinomial(SEXP data);
RcppExport SEXP _InferN0_infernal_fitNegativeBinomial(SEXP dataSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type data(dataSEXP);
    rcpp_result_gen = Rcpp::wrap(infernal_fitNegativeBinomial(data));
    return rcpp_result_gen;
END_RCPP
}
// infernal_getParameterList
Rcpp::CharacterVector infernal_getParameterList(SEXP rscp);
RcppExport SEXP _InferN0_infernal_getParameterList(SEXP rscpSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type rscp(rscpSEXP);
    rcpp_result_gen = Rcpp::wrap(infernal_getParameterList(rscp));
    return rcpp_result_gen;
END_RCPP
}
// infernal_setParameter
uint32_t infernal_setParameter(SEXP rscp, uint32_t which, SEXP value);
RcppExport SEXP _InferN0_infernal_setParameter(SEXP rscpSEXP, SEXP whichSEXP, SEXP valueSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type rscp(rscpSEXP);
    Rcpp::traits::input_parameter< uint32_t >::type which(whichSEXP);
    Rcpp::traits::input_parameter< SEXP >::type value(valueSEXP);
    rcpp_result_gen = Rcpp::wrap(infernal_setParameter(rscp, which, value));
    return rcpp_result_gen;
END_RCPP
}
// infernal_getParameter
SEXP infernal_getParameter(SEXP rscp, uint32_t which);
RcppExport SEXP _InferN0_infernal_getParameter(SEXP rscpSEXP, SEXP whichSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type rscp(rscpSEXP);
    Rcpp::traits::input_parameter< uint32_t >::type which(whichSEXP);
    rcpp_result_gen = Rcpp::wrap(infernal_getParameter(rscp, which));
    return rcpp_result_gen;
END_RCPP
}
// infernal_getDataRow
SEXP infernal_getDataRow(SEXP rscp, Rcpp::List list);
RcppExport SEXP _InferN0_infernal_getDataRow(SEXP rscpSEXP, SEXP listSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type rscp(rscpSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type list(listSEXP);
    rcpp_result_gen = Rcpp::wrap(infernal_getDataRow(rscp, list));
    return rcpp_result_gen;
END_RCPP
}
// infernal_runTestFunction
void infernal_runTestFunction(SEXP data);
RcppExport SEXP _InferN0_infernal_runTestFunction(SEXP dataSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type data(dataSEXP);
    infernal_runTestFunction(data);
    return R_NilValue;
END_RCPP
}
// infernal_HierarchicalClustering
void infernal_HierarchicalClustering(SEXP infscope, Rcpp::List list);
RcppExport SEXP _InferN0_infernal_HierarchicalClustering(SEXP infscopeSEXP, SEXP listSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type infscope(infscopeSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type list(listSEXP);
    infernal_HierarchicalClustering(infscope, list);
    return R_NilValue;
END_RCPP
}
// infernal_computePartialCorrelation
Rcpp::List infernal_computePartialCorrelation(SEXP infscope, Rcpp::List list);
RcppExport SEXP _InferN0_infernal_computePartialCorrelation(SEXP infscopeSEXP, SEXP listSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type infscope(infscopeSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type list(listSEXP);
    rcpp_result_gen = Rcpp::wrap(infernal_computePartialCorrelation(infscope, list));
    return rcpp_result_gen;
END_RCPP
}
// infernal_IdentifyNetwork
Rcpp::List infernal_IdentifyNetwork(Rcpp::S4 input, Rcpp::List list, Rcpp::List traincov, Rcpp::List trainmean);
RcppExport SEXP _InferN0_infernal_IdentifyNetwork(SEXP inputSEXP, SEXP listSEXP, SEXP traincovSEXP, SEXP trainmeanSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::S4 >::type input(inputSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type list(listSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type traincov(traincovSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type trainmean(trainmeanSEXP);
    rcpp_result_gen = Rcpp::wrap(infernal_IdentifyNetwork(input, list, traincov, trainmean));
    return rcpp_result_gen;
END_RCPP
}
// infernal_IdentifyConstrainedCovar
Rcpp::List infernal_IdentifyConstrainedCovar(Rcpp::S4 input, Rcpp::List list, Rcpp::List traincov, Rcpp::List trainmean);
RcppExport SEXP _InferN0_infernal_IdentifyConstrainedCovar(SEXP inputSEXP, SEXP listSEXP, SEXP traincovSEXP, SEXP trainmeanSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::S4 >::type input(inputSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type list(listSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type traincov(traincovSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type trainmean(trainmeanSEXP);
    rcpp_result_gen = Rcpp::wrap(infernal_IdentifyConstrainedCovar(input, list, traincov, trainmean));
    return rcpp_result_gen;
END_RCPP
}
// infernal_saveHierarchical
void infernal_saveHierarchical(SEXP data, SEXP path);
RcppExport SEXP _InferN0_infernal_saveHierarchical(SEXP dataSEXP, SEXP pathSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type data(dataSEXP);
    Rcpp::traits::input_parameter< SEXP >::type path(pathSEXP);
    infernal_saveHierarchical(data, path);
    return R_NilValue;
END_RCPP
}
// infernal_show
void infernal_show(SEXP data);
RcppExport SEXP _InferN0_infernal_show(SEXP dataSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type data(dataSEXP);
    infernal_show(data);
    return R_NilValue;
END_RCPP
}
// infernal_computeCovar
Rcpp::List infernal_computeCovar(SEXP scptrexp, Rcpp::List list);
RcppExport SEXP _InferN0_infernal_computeCovar(SEXP scptrexpSEXP, SEXP listSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type scptrexp(scptrexpSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type list(listSEXP);
    rcpp_result_gen = Rcpp::wrap(infernal_computeCovar(scptrexp, list));
    return rcpp_result_gen;
END_RCPP
}
// infernal_getFullDeviations
Rcpp::NumericMatrix infernal_getFullDeviations(SEXP InferN0scp, Rcpp::List list);
RcppExport SEXP _InferN0_infernal_getFullDeviations(SEXP InferN0scpSEXP, SEXP listSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type InferN0scp(InferN0scpSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type list(listSEXP);
    rcpp_result_gen = Rcpp::wrap(infernal_getFullDeviations(InferN0scp, list));
    return rcpp_result_gen;
END_RCPP
}
// infernal_cmpModeledVariance
Rcpp::List infernal_cmpModeledVariance(SEXP infscope, Rcpp::List list);
RcppExport SEXP _InferN0_infernal_cmpModeledVariance(SEXP infscopeSEXP, SEXP listSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type infscope(infscopeSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type list(listSEXP);
    rcpp_result_gen = Rcpp::wrap(infernal_cmpModeledVariance(infscope, list));
    return rcpp_result_gen;
END_RCPP
}
// infernal_cmpVariance
Rcpp::List infernal_cmpVariance(SEXP infscope, Rcpp::List list);
RcppExport SEXP _InferN0_infernal_cmpVariance(SEXP infscopeSEXP, SEXP listSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type infscope(infscopeSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type list(listSEXP);
    rcpp_result_gen = Rcpp::wrap(infernal_cmpVariance(infscope, list));
    return rcpp_result_gen;
END_RCPP
}
// internal_Infern0_hierarchicalClustering
Rcpp::List internal_Infern0_hierarchicalClustering(SEXP data, SEXP cellgroup, int nbclusters);
RcppExport SEXP _InferN0_internal_Infern0_hierarchicalClustering(SEXP dataSEXP, SEXP cellgroupSEXP, SEXP nbclustersSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type data(dataSEXP);
    Rcpp::traits::input_parameter< SEXP >::type cellgroup(cellgroupSEXP);
    Rcpp::traits::input_parameter< int >::type nbclusters(nbclustersSEXP);
    rcpp_result_gen = Rcpp::wrap(internal_Infern0_hierarchicalClustering(data, cellgroup, nbclusters));
    return rcpp_result_gen;
END_RCPP
}
// infernal_mergeSparseMatrices
Rcpp::S4 infernal_mergeSparseMatrices(Rcpp::S4 matrixA, Rcpp::S4 matrixB);
RcppExport SEXP _InferN0_infernal_mergeSparseMatrices(SEXP matrixASEXP, SEXP matrixBSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::S4 >::type matrixA(matrixASEXP);
    Rcpp::traits::input_parameter< Rcpp::S4 >::type matrixB(matrixBSEXP);
    rcpp_result_gen = Rcpp::wrap(infernal_mergeSparseMatrices(matrixA, matrixB));
    return rcpp_result_gen;
END_RCPP
}
// infernal_genSynthetic
Rcpp::List infernal_genSynthetic(SEXP infscope, Rcpp::List data);
RcppExport SEXP _InferN0_infernal_genSynthetic(SEXP infscopeSEXP, SEXP dataSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type infscope(infscopeSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type data(dataSEXP);
    rcpp_result_gen = Rcpp::wrap(infernal_genSynthetic(infscope, data));
    return rcpp_result_gen;
END_RCPP
}
// infernal_initModelWithClustering
void infernal_initModelWithClustering(SEXP infscope, Rcpp::List data);
RcppExport SEXP _InferN0_infernal_initModelWithClustering(SEXP infscopeSEXP, SEXP dataSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type infscope(infscopeSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type data(dataSEXP);
    infernal_initModelWithClustering(infscope, data);
    return R_NilValue;
END_RCPP
}
// infernal_modelhidden
Rcpp::List infernal_modelhidden(SEXP infscope, Rcpp::List data);
RcppExport SEXP _InferN0_infernal_modelhidden(SEXP infscopeSEXP, SEXP dataSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type infscope(infscopeSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type data(dataSEXP);
    rcpp_result_gen = Rcpp::wrap(infernal_modelhidden(infscope, data));
    return rcpp_result_gen;
END_RCPP
}
// infernal_modelhidden_v2
Rcpp::List infernal_modelhidden_v2(SEXP infscope, Rcpp::List data);
RcppExport SEXP _InferN0_infernal_modelhidden_v2(SEXP infscopeSEXP, SEXP dataSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type infscope(infscopeSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type data(dataSEXP);
    rcpp_result_gen = Rcpp::wrap(infernal_modelhidden_v2(infscope, data));
    return rcpp_result_gen;
END_RCPP
}
// infernal_exportpcs
Rcpp::List infernal_exportpcs(SEXP infscope, Rcpp::List rlist);
RcppExport SEXP _InferN0_infernal_exportpcs(SEXP infscopeSEXP, SEXP rlistSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type infscope(infscopeSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type rlist(rlistSEXP);
    rcpp_result_gen = Rcpp::wrap(infernal_exportpcs(infscope, rlist));
    return rcpp_result_gen;
END_RCPP
}
// infernalDensity
Rcpp::List infernalDensity(Rcpp::NumericMatrix data, Rcpp::List inlist);
RcppExport SEXP _InferN0_infernalDensity(SEXP dataSEXP, SEXP inlistSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type data(dataSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type inlist(inlistSEXP);
    rcpp_result_gen = Rcpp::wrap(infernalDensity(data, inlist));
    return rcpp_result_gen;
END_RCPP
}
// infernalTest
Rcpp::List infernalTest();
RcppExport SEXP _InferN0_infernalTest() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(infernalTest());
    return rcpp_result_gen;
END_RCPP
}
// infernalKendallCols
Rcpp::List infernalKendallCols(Rcpp::S4 x);
RcppExport SEXP _InferN0_infernalKendallCols(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::S4 >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(infernalKendallCols(x));
    return rcpp_result_gen;
END_RCPP
}
// infernalCumulant
Rcpp::List infernalCumulant(Rcpp::List inlist);
RcppExport SEXP _InferN0_infernalCumulant(SEXP inlistSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type inlist(inlistSEXP);
    rcpp_result_gen = Rcpp::wrap(infernalCumulant(inlist));
    return rcpp_result_gen;
END_RCPP
}
// infernalDeconvolve2D
Rcpp::List infernalDeconvolve2D(Rcpp::S4 scdata_in, Rcpp::S4 bulk_in, Rcpp::List inlist);
RcppExport SEXP _InferN0_infernalDeconvolve2D(SEXP scdata_inSEXP, SEXP bulk_inSEXP, SEXP inlistSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::S4 >::type scdata_in(scdata_inSEXP);
    Rcpp::traits::input_parameter< Rcpp::S4 >::type bulk_in(bulk_inSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type inlist(inlistSEXP);
    rcpp_result_gen = Rcpp::wrap(infernalDeconvolve2D(scdata_in, bulk_in, inlist));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_InferN0_InferN0_identifyL0Network", (DL_FUNC) &_InferN0_InferN0_identifyL0Network, 2},
    {"_InferN0_InferN0_identifyL0NetworkGold", (DL_FUNC) &_InferN0_InferN0_identifyL0NetworkGold, 3},
    {"_InferN0_InferN0_printstringlist", (DL_FUNC) &_InferN0_InferN0_printstringlist, 1},
    {"_InferN0_InferN0_identifyL0NetworkGoldResurrected", (DL_FUNC) &_InferN0_InferN0_identifyL0NetworkGoldResurrected, 4},
    {"_InferN0_infernal_wilcox_geneset", (DL_FUNC) &_InferN0_infernal_wilcox_geneset, 2},
    {"_InferN0_infernal_wilcox_scope", (DL_FUNC) &_InferN0_infernal_wilcox_scope, 2},
    {"_InferN0_infernal_wilcox_matrix", (DL_FUNC) &_InferN0_infernal_wilcox_matrix, 2},
    {"_InferN0_infernal_find_DE", (DL_FUNC) &_InferN0_infernal_find_DE, 2},
    {"_InferN0_infernal_find_markers", (DL_FUNC) &_InferN0_infernal_find_markers, 2},
    {"_InferN0_infernal_sharedDE", (DL_FUNC) &_InferN0_infernal_sharedDE, 1},
    {"_InferN0_infernal_cmpInverses", (DL_FUNC) &_InferN0_infernal_cmpInverses, 1},
    {"_InferN0_testSlot", (DL_FUNC) &_InferN0_testSlot, 2},
    {"_InferN0_infernal_serialize", (DL_FUNC) &_InferN0_infernal_serialize, 1},
    {"_InferN0_infernal_deserialize", (DL_FUNC) &_InferN0_infernal_deserialize, 2},
    {"_InferN0_infernal_createScope", (DL_FUNC) &_InferN0_infernal_createScope, 0},
    {"_InferN0_infernal_readScope", (DL_FUNC) &_InferN0_infernal_readScope, 1},
    {"_InferN0_infernal_readScope2", (DL_FUNC) &_InferN0_infernal_readScope2, 1},
    {"_InferN0_infernal_writeScope", (DL_FUNC) &_InferN0_infernal_writeScope, 2},
    {"_InferN0_infernal_writeScope2", (DL_FUNC) &_InferN0_infernal_writeScope2, 3},
    {"_InferN0_infernal_readMatrixFolder", (DL_FUNC) &_InferN0_infernal_readMatrixFolder, 1},
    {"_InferN0_infernal_convertDenseTextMatrixToSparseFormat", (DL_FUNC) &_InferN0_infernal_convertDenseTextMatrixToSparseFormat, 1},
    {"_InferN0_infernal_exportDataAsMtx", (DL_FUNC) &_InferN0_infernal_exportDataAsMtx, 2},
    {"_InferN0_infernal_reverseTCM", (DL_FUNC) &_InferN0_infernal_reverseTCM, 2},
    {"_InferN0_infernal_reverseTCM_inMTX", (DL_FUNC) &_InferN0_infernal_reverseTCM_inMTX, 1},
    {"_InferN0_infernal_readSparseMatrix", (DL_FUNC) &_InferN0_infernal_readSparseMatrix, 3},
    {"_InferN0_infernal_readDenseMatrix", (DL_FUNC) &_InferN0_infernal_readDenseMatrix, 3},
    {"_InferN0_infernal_createOutput", (DL_FUNC) &_InferN0_infernal_createOutput, 1},
    {"_InferN0_infernal_fitNegativeBinomial", (DL_FUNC) &_InferN0_infernal_fitNegativeBinomial, 1},
    {"_InferN0_infernal_getParameterList", (DL_FUNC) &_InferN0_infernal_getParameterList, 1},
    {"_InferN0_infernal_setParameter", (DL_FUNC) &_InferN0_infernal_setParameter, 3},
    {"_InferN0_infernal_getParameter", (DL_FUNC) &_InferN0_infernal_getParameter, 2},
    {"_InferN0_infernal_getDataRow", (DL_FUNC) &_InferN0_infernal_getDataRow, 2},
    {"_InferN0_infernal_runTestFunction", (DL_FUNC) &_InferN0_infernal_runTestFunction, 1},
    {"_InferN0_infernal_HierarchicalClustering", (DL_FUNC) &_InferN0_infernal_HierarchicalClustering, 2},
    {"_InferN0_infernal_computePartialCorrelation", (DL_FUNC) &_InferN0_infernal_computePartialCorrelation, 2},
    {"_InferN0_infernal_IdentifyNetwork", (DL_FUNC) &_InferN0_infernal_IdentifyNetwork, 4},
    {"_InferN0_infernal_IdentifyConstrainedCovar", (DL_FUNC) &_InferN0_infernal_IdentifyConstrainedCovar, 4},
    {"_InferN0_infernal_saveHierarchical", (DL_FUNC) &_InferN0_infernal_saveHierarchical, 2},
    {"_InferN0_infernal_show", (DL_FUNC) &_InferN0_infernal_show, 1},
    {"_InferN0_infernal_computeCovar", (DL_FUNC) &_InferN0_infernal_computeCovar, 2},
    {"_InferN0_infernal_getFullDeviations", (DL_FUNC) &_InferN0_infernal_getFullDeviations, 2},
    {"_InferN0_infernal_cmpModeledVariance", (DL_FUNC) &_InferN0_infernal_cmpModeledVariance, 2},
    {"_InferN0_infernal_cmpVariance", (DL_FUNC) &_InferN0_infernal_cmpVariance, 2},
    {"_InferN0_internal_Infern0_hierarchicalClustering", (DL_FUNC) &_InferN0_internal_Infern0_hierarchicalClustering, 3},
    {"_InferN0_infernal_mergeSparseMatrices", (DL_FUNC) &_InferN0_infernal_mergeSparseMatrices, 2},
    {"_InferN0_infernal_genSynthetic", (DL_FUNC) &_InferN0_infernal_genSynthetic, 2},
    {"_InferN0_infernal_initModelWithClustering", (DL_FUNC) &_InferN0_infernal_initModelWithClustering, 2},
    {"_InferN0_infernal_modelhidden", (DL_FUNC) &_InferN0_infernal_modelhidden, 2},
    {"_InferN0_infernal_modelhidden_v2", (DL_FUNC) &_InferN0_infernal_modelhidden_v2, 2},
    {"_InferN0_infernal_exportpcs", (DL_FUNC) &_InferN0_infernal_exportpcs, 2},
    {"_InferN0_infernalDensity", (DL_FUNC) &_InferN0_infernalDensity, 2},
    {"_InferN0_infernalTest", (DL_FUNC) &_InferN0_infernalTest, 0},
    {"_InferN0_infernalKendallCols", (DL_FUNC) &_InferN0_infernalKendallCols, 1},
    {"_InferN0_infernalCumulant", (DL_FUNC) &_InferN0_infernalCumulant, 1},
    {"_InferN0_infernalDeconvolve2D", (DL_FUNC) &_InferN0_infernalDeconvolve2D, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_InferN0(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
