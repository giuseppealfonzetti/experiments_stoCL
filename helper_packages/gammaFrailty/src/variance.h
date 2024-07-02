#ifndef variance_H
#define variance_H

#include "pairClass.h"

//' @export
// [[Rcpp::export]]
Eigen::MatrixXd sampleJ(
        Eigen::VectorXd THETA,
        Eigen::MatrixXd DATA,
        Eigen::MatrixXd X,
        unsigned int STRUCT,
        const bool PRINTFLAG = false,
        const int PAIRS_RANGE = 100
){
    unsigned int d = THETA.size();
    unsigned int n = DATA.rows();
    unsigned int p = DATA.cols();
    unsigned int r = X.cols();

    double lxi  = THETA(0);
    double artanhrho = THETA(1);
    Eigen::VectorXd beta = THETA.segment(2, r);
    Eigen::MatrixXd J = Eigen::MatrixXd::Zero(d,d);



    pair_class pair;
    for(unsigned int i = 0; i < n; i++){
        const Eigen::VectorXd x_i  = X.row(i);
        Eigen::VectorXd gradient_i = Eigen::VectorXd::Zero(d);

        for(unsigned int j = 1; j < p; j++){
            unsigned int n_j = DATA(i, j);
            double alpha_j = THETA(r+2+j);

            for( unsigned int jp = std::max(0, int(j - PAIRS_RANGE)); jp < j; jp++){
                const unsigned int n_jp = DATA(i, jp);
                const double alpha_jp = THETA(r+2+jp);

                pair.setup_(j, jp, n_j, n_jp, x_i, beta, alpha_j, alpha_jp, lxi, artanhrho, p, STRUCT);
                pair.compute_intermediate_( );
                pair.compute_dintermediate_();
                const double ll = pair.compute_ll_()/n;
                gradient_i += pair.compute_gradient_();

                if(PRINTFLAG &!(ll==ll)) Rcpp::Rcout << "(i, j, j') = ("<<i<<", "<<j<<", "<< jp << ") -> "<< ll << "\n";

            }
        }

        J += gradient_i*gradient_i.transpose();
    }


    J /= n;
    return J;
}

// Sample estimate of H (or -hessian)
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd sampleH(
        Eigen::VectorXd THETA,
        Eigen::MatrixXd DATA,
        Eigen::MatrixXd X,
        unsigned int STRUCT,
        const bool PRINTFLAG = false,
        const bool INVERTFLAG = false,
        const int PAIRS_RANGE = 100

){
    unsigned int d = THETA.size();
    unsigned int n = DATA.rows();
    unsigned int p = DATA.cols();
    unsigned int r = X.cols();

    double lxi  = THETA(0);
    double artanhrho = THETA(1);
    Eigen::VectorXd beta = THETA.segment(2, r);
    Eigen::MatrixXd H = Eigen::MatrixXd::Zero(d,d);



    pair_class pair;
    for(unsigned int i = 0; i < n; i++){
        Eigen::VectorXd x_i = X.row(i);

        for(unsigned int j = 1; j < p; j++){
            unsigned int n_j = DATA(i, j);
            double alpha_j = THETA(r+2+j);

            for( unsigned int jp = std::max(0, int(j - PAIRS_RANGE)); jp < j; jp++){
                const unsigned int n_jp = DATA(i, jp);
                const double alpha_jp = THETA(r+2+jp);

                pair.setup_(j, jp, n_j, n_jp, x_i, beta, alpha_j, alpha_jp, lxi, artanhrho, p, STRUCT);
                pair.compute_intermediate_( );
                pair.compute_dintermediate_();
                const double ll = pair.compute_ll_()/n;
                const Eigen::VectorXd gradient_k = pair.compute_gradient_();

                H += gradient_k*gradient_k.transpose();
                if(PRINTFLAG &!(ll==ll)) Rcpp::Rcout << "(i, j, j') = ("<<i<<", "<<j<<", "<< jp << ") -> "<< ll << "\n";

            }
        }


    }


    Eigen::MatrixXd out = H/n;
    if(INVERTFLAG){
        const Eigen::MatrixXd Id = Eigen::MatrixXd::Identity(d,d);
        Eigen::LLT<Eigen::MatrixXd> llt;
        llt.compute(out);
        out = llt.solve(Id);
    }
    return out;
}

//' @export
// [[Rcpp::export]]
Rcpp::List sampleVar(
        Eigen::VectorXd THETA,
        Eigen::MatrixXd DATA,
        Eigen::MatrixXd X,
        unsigned int STRUCT,
        const unsigned int NU,
        const unsigned int METHOD,
        const unsigned int RANGE,
        const bool TOTFLAG,
        const bool PRINTFLAG,
        const int PAIRS_RANGE = 100
){

    // Identify dimensions
    const unsigned int d = THETA.size();
    const unsigned int n = DATA.rows();

    // Initialise variance matrices and std errors vectors
    Eigen::MatrixXd tot, jmat, sandwich, cond;
    Eigen::VectorXd stochastic_se, statistical_se, tot_se;


    Rcpp::List var;
    Rcpp::List se;

    // Compute inverse of negative hessian
    Eigen::MatrixXd inv_hmat = sampleH(THETA, DATA, X, STRUCT, PRINTFLAG, true, PAIRS_RANGE);

    if(METHOD == 0){
        jmat = sampleJ(THETA, DATA, X, STRUCT, PRINTFLAG, PAIRS_RANGE);
        sandwich = inv_hmat*jmat*inv_hmat;
        tot = sandwich/n;

        var = Rcpp::List::create( Rcpp::Named("var_tot") =  sandwich/static_cast<double>(n));
        se  = Rcpp::List::create( Rcpp::Named("se_tot" ) = (sandwich/static_cast<double>(n)).diagonal().array().sqrt() );

    }

    if(METHOD == 1){
        jmat = sampleJ(THETA, DATA, X, STRUCT, PRINTFLAG, PAIRS_RANGE);
        sandwich = inv_hmat*jmat*inv_hmat;
        double st_scale = NU*RANGE;

        if(TOTFLAG) {

            var = Rcpp::List::create(
                Rcpp::Named("var_stoc") = sandwich/st_scale,
                Rcpp::Named("var_stat") = sandwich/n,
                Rcpp::Named("var_tot")  = sandwich.array()*( (1/st_scale) + (1/static_cast<double>(n)) )
                );

            se = Rcpp::List::create(
                Rcpp::Named("se_stoc") = (sandwich/st_scale).diagonal().array().sqrt(),
                Rcpp::Named("se_stat") = (sandwich/n).diagonal().array().sqrt(),
                Rcpp::Named("se_tot")  = (sandwich*( (1/st_scale) + (1/static_cast<double>(n)) )).diagonal().array().sqrt()
            );
        }else{
            var = Rcpp::List::create(Rcpp::Named("var_stoc") = sandwich/st_scale);
            se  = Rcpp::List::create(Rcpp::Named("se_stoc") = (sandwich/st_scale).diagonal().array().sqrt());
        }
    }

    if(METHOD == 2){

        double st_scale = NU*RANGE;

        if(TOTFLAG) {

            jmat = sampleJ(THETA, DATA, X, STRUCT, PRINTFLAG, PAIRS_RANGE);
            sandwich = inv_hmat*jmat*inv_hmat;
            var = Rcpp::List::create(
                Rcpp::Named("var_stoc") = inv_hmat/st_scale,
                Rcpp::Named("var_stat") = sandwich/static_cast<double>(n),
                Rcpp::Named("var_tot")  = (sandwich/static_cast<double>(n))+(inv_hmat/st_scale)
            );

            se = Rcpp::List::create(
                Rcpp::Named("se_stoc") = (inv_hmat/st_scale).diagonal().array().sqrt(),
                Rcpp::Named("se_stat") = (sandwich/static_cast<double>(n)).diagonal().array().sqrt(),
                Rcpp::Named("se_tot")  = ((sandwich/static_cast<double>(n))+(inv_hmat/st_scale)).diagonal().array().sqrt()
            );
        }else{
            var = Rcpp::List::create(Rcpp::Named("var_stoc") = inv_hmat/st_scale);
            se  = Rcpp::List::create(Rcpp::Named("se_stoc") = (inv_hmat/st_scale).diagonal().array().sqrt());
        }
    }


    Rcpp::List out = Rcpp::List::create(Rcpp::Named("var") = var, Rcpp::Named("se") = se);
    return out;
}

#endif

