#include <Rcpp.h>
using namespace Rcpp;

//' @importFrom Rcpp evalCpp
//' @useDynLib OptSurvCutR, .registration = TRUE
// [[Rcpp::export]]
IntegerVector cpp_get_group_assignments(NumericVector predictor, NumericVector cuts) {
    int n = predictor.size();
    int num_cuts = cuts.size();
    IntegerVector groups(n);
    
    // Process allocations down down inside continuous memory blocks
    for (int i = 0; i < n; i++) {
        int g = 1; // Instantiate Base Cohort Tier
        while (g <= num_cuts && predictor[i] > cuts[g - 1]) {
            g++;
        }
        groups[i] = g;
    }
    return groups;
}
