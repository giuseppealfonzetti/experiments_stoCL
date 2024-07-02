#ifndef utils_H
#define utils_H

// Fisher transformation for correlation
// [[Rcpp::export]]
double zofr_cpp(const double r){
    double z = .5*log((1+r)/(1-r));
    return z;
}

// Inverse of Fisher transformation
// [[Rcpp::export]]
double rofz_cpp(const double z){
    double r = (exp(2*z)-1)/(exp(2*z)+1);
    return r;
}

// Derivative of inverse of Fisher transformation
// [[Rcpp::export]]
double drofz_cpp(const double z){
    double r = rofz_cpp(z);
    double dr = 2*(1-r)*exp(2*z)/(exp(2*z)+1);
    return dr;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix rmultinom_wrapper(const double prob, const unsigned int classes, const unsigned int batch, const unsigned int K) {

    Rcpp::NumericVector probs(classes, prob);
    Rcpp::IntegerVector outcome(classes);
    R::rmultinom(batch, probs.begin(), classes, outcome.begin());


    Rcpp::NumericMatrix out(classes, K);
    for(unsigned int j = 0; j < K; j++){
        out(Rcpp::_,j) = outcome;
    }

    return out;
}

 // [[Rcpp::export]]
 std::vector<int> hyper_sampling(const unsigned int K, const unsigned int N, const unsigned int SEED){
     std::mt19937 randomizer(SEED);
     std::vector<int> pool(N*K);
     std::iota (std::begin(pool), std::end(pool), 0);
     std::shuffle(pool.begin(), pool.end(), randomizer);
     return(pool);
 }

// [[Rcpp::export]]
 std::vector<int> unit_sampling(const unsigned int N, const unsigned int SEED){
     std::mt19937 randomizer(SEED);
     std::vector<int> pool(N);
     std::iota (std::begin(pool), std::end(pool), 0);
     std::shuffle(pool.begin(), pool.end(), randomizer);
     return(pool);
}

//'@export
// [[Rcpp::export]]
 std::vector<int> components_given_unit(const unsigned int UNIT, const unsigned int K){
     std::vector<int> pool(K);
     std::iota (std::begin(pool), std::end(pool), UNIT*K);
     return(pool);
}

// [[Rcpp::export]]
std::vector<int> bernoulli_sampling(const unsigned int K, const unsigned int N, const double PROB){
     std::vector<int> pool;
     for( int iterator = 0; iterator < N*K; iterator++){
         if(R::runif(0,1) < PROB ) pool.push_back(iterator);
     }
     return(pool);
}

// [[Rcpp::export]]
std::vector<int> index_to_component(const unsigned int P, const unsigned int N, const unsigned int INDEX){
    int K = P*(P-1)/2;
    if(INDEX >= N*K) Rcpp::stop("'INDEX' must be less than N*K (count starts from 0)");
    int i = INDEX / K;
    int pair = INDEX % K;
    int j = P - 2 - floor(sqrt(-8*pair + 4*P*(P-1)-7)/2.0 - 0.5);
    int jp = pair + j + 1 - P*(P-1)/2 + (P-j)*((P-j)-1)/2;
    std::vector<int> component{i, j, jp};
    return(component);
}
#endif
