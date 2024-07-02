#include <Rcpp.h>
#define EIGEN_PERMANENTLY_DISABLE_STUPID_WARNINGS
#define EIGEN_DONT_PARALLELIZE
#include <RcppEigen.h>
#include <RcppClock.h>
#include <math.h>
#include <random>
#include "utils.h"
#include "variance.h"
#include "binarynodeClass.h"

// [[Rcpp::depends(RcppEigen, RcppClock)]]

//' @export
// [[Rcpp::export]]
Rcpp::List ncl(
        Eigen::MatrixXd &DATA,
        Eigen::VectorXd &THETA,
        std::vector<bool> &CONSTRAINTS,
        const bool VERBOSEFLAG = false
){
    // Identify dimensions
    const unsigned int n = DATA.rows();
    const unsigned int d = THETA.size();
    const unsigned int p = (- 1 + sqrt(1 + 8*d))/2;

    binarynodeClass Node;
    double cl = 0;
    Eigen::VectorXd gradient = Eigen::VectorXd::Zero(d);
    for(unsigned int i = 0; i < n; i++){

        Eigen::VectorXd data_i = DATA.row(i);
        for(unsigned int node = 0; node < p; node++){

            Node.setup_(data_i, THETA, CONSTRAINTS, p, node);
            cl -= Node.ll_();
            gradient -= Node.gradient_();
        }
    }

    Rcpp::List output =
        Rcpp::List::create(
            Rcpp::Named("nll") = cl,
            Rcpp::Named("ngradient") = gradient
        );

    return output;

}



//' @export
// [[Rcpp::export]]
Rcpp::List isingGraph3(
        const Eigen::MatrixXd &DATA,
        const Eigen::MatrixXd &HOLDOUT,
        const Eigen::VectorXd &THETA_INIT,
        const std::vector<bool> &CONSTRAINTS,
        const unsigned int MAXT,
        const unsigned int BURN,
        const double STEPSIZE,
        const double NU,
        const int METHODFLAG,
        Eigen::VectorXd SCALEVEC,
        const unsigned int SEED = 123,
        const bool VERBOSEFLAG = false,
        const bool HOLDOUTFLAG = false,
        const double PAR1 = 1,
        const double PAR2 = 1,
        const double PAR3 = .75,
        const unsigned int SAMPLING_WINDOW = 1,
        const unsigned int EACH = 1,
        const unsigned int EACHCLOCK = 100,
        const int STEPSIZEFLAG = 0,
        const int T_INIT = 1
){
    // Set up clock monitor to export to R session trough RcppClock
    Rcpp::Clock clock;
    clock.tick("main");

    // Identify dimensions
    const unsigned int n = DATA.rows();
    const unsigned int d = THETA_INIT.size();
    const unsigned int p = (- 1 + sqrt(1 + 8*d))/2;
    const unsigned int kk= p;

    // scaling constant according to sampling_scheme
    // Compute scaling constant
    double scale;
    switch(METHODFLAG){
    case 0:
        scale = 1/static_cast<double>(n) ;
        break;
    case 1:
        scale = 1/static_cast<double>(NU);
        break;
    case 2:
        scale = 1/static_cast<double>(NU);
        break;
    case 3:
        scale = 1/static_cast<double>(NU);
        break;
    case 4:
        scale = 1/static_cast<double>(NU);
        break;
    case 5:
        scale = 1/static_cast<double>(NU);
        break;
    }
    // Rcpp::Rcout << " Final scale = " << scale<< "\n";

    // Initialize storage for iterations quantities
    // Eigen::MatrixXd path_theta    = Eigen::MatrixXd::Zero(MAXT + 1, d); path_theta.row(0)    = THETA_INIT;
    // Eigen::MatrixXd path_av_theta = Eigen::MatrixXd::Zero(MAXT + 1, d); path_av_theta.row(0) = THETA_INIT;
    // Eigen::MatrixXd path_grad     = Eigen::MatrixXd::Zero(MAXT,     d);
    std::vector<Eigen::VectorXd> path_theta; path_theta.push_back(THETA_INIT);
    std::vector<Eigen::VectorXd> path_av_theta; path_av_theta.push_back(THETA_INIT);
    std::vector<Eigen::VectorXd> path_grad;

    std::vector<double> path_nll;
    std::vector<double> iter_idx; iter_idx.push_back(T_INIT-1);
    int last_stored_iter = T_INIT-1;


    ///////////////////////
    /* OPTIMISATION LOOP */
    ///////////////////////

    binarynodeClass Node;
    std::vector<int> outloop_pool;
    Eigen::VectorXd theta_t = THETA_INIT;
    Eigen::VectorXd av_theta_t = THETA_INIT;
    Eigen::VectorXd av_theta_prev = THETA_INIT;

    unsigned int sampling_window_iterator = 0;
    unsigned int last_iter = MAXT;
    for(unsigned int t = T_INIT; t <= MAXT; t++){
        // check user interruption
        Rcpp::checkUserInterrupt();
        if(VERBOSEFLAG)Rcpp::Rcout << "\r Iter:" << t << " ";
        std::vector<int> pool_t;

        /////////////////////
        // SAMPLING SCHEME //
        /////////////////////
        if(t % EACHCLOCK == 0) clock.tick("sampling_step");


        switch(METHODFLAG){
        case 0:
            {pool_t.resize(n*kk);
            std::iota (std::begin(pool_t), std::end(pool_t), 0);
            break;}
        case 1:
            {std::vector<int> pool_units = unit_sampling(n, SEED + t);
                for(unsigned int i = 0; i < NU; i++){
                    std::vector<int> tmp = components_given_unit(int(pool_units[i]), kk);
                    pool_t.insert(pool_t.end(), std::begin(tmp), std::end(tmp));
                }
            break;}
        case 2:
            {double prob = static_cast<double>(NU)/static_cast<double>(n);
            // need for external seeding
            pool_t = bernoulli_sampling(kk, n, prob);
            break;}
        case 3:
            {std::vector<int> tmp = hyper_sampling(kk, n, SEED + t);
            pool_t = {tmp.begin(), tmp.begin()+int(NU*kk)};
            break;}
        case 4:
            {
                if(sampling_window_iterator == 0){
                    outloop_pool = unit_sampling(n, SEED + t);
                }
                for(unsigned int i = 0; i < NU; i++){
                    int tmp_ind = sampling_window_iterator*NU + i;
                    std::vector<int> tmp = components_given_unit(int(outloop_pool[tmp_ind]), kk);
                    pool_t.insert(pool_t.end(), std::begin(tmp), std::end(tmp));
                }
                sampling_window_iterator++;
                if(sampling_window_iterator==SAMPLING_WINDOW) sampling_window_iterator = 0;
                break;
            }
        case 5:
            {
                if(sampling_window_iterator == 0){
                outloop_pool = hyper_sampling(kk, n, SEED + t);
                }

                int tmp_ind = sampling_window_iterator*NU*kk;
                pool_t = {outloop_pool.begin() + tmp_ind, outloop_pool.begin() + tmp_ind + int(NU * kk)};

                sampling_window_iterator++;
                if(sampling_window_iterator==SAMPLING_WINDOW) sampling_window_iterator = 0;
                break;
            }

        }

        if(t % EACHCLOCK == 0) clock.tock("sampling_step");

        //initialize iteration quantities
        Eigen::VectorXd ngradient_t = Eigen::VectorXd::Zero(d);
        double ncl = 0;

        ///////////////////////////
        /* GRADIENT COMPUTATION  */
        ///////////////////////////
        if(t % EACHCLOCK == 0) clock.tick("stochastic_gradient");
        for(unsigned int index = 0; index < pool_t.size(); index++){
            std::vector<int> component = index_to_component(kk, n, int(pool_t[index]));
            Node.setup_(DATA.row(int(component[0])), theta_t, CONSTRAINTS, p, int(component[1]));
            ngradient_t -= Node.gradient_();
            ncl -= Node.ll_();
        }
        if(t % EACHCLOCK == 0) clock.tock("stochastic_gradient");

        ncl *= scale;
        ngradient_t *= scale;

        ///////////////////////////
        /*    PARAMETERS UPDATE  */
        ///////////////////////////
        if(t % EACHCLOCK == 0) clock.tick("update");
        double stepsize_t = STEPSIZE;
        switch(STEPSIZEFLAG){
        case 0:
            stepsize_t *= pow(t, -PAR3);
            break;
        case 1:
            stepsize_t *= PAR1 * pow(1 + PAR2*STEPSIZE*t, -PAR3);
            break;
        }
        theta_t -= Eigen::VectorXd(stepsize_t * SCALEVEC.array() * ngradient_t.array());
        if(t % EACHCLOCK == 0) clock.tock("update");


        /////////////////////////////////
        /*    ONLINE AVERAGE UPDATE    */
        /////////////////////////////////

        if(t <= BURN){
            av_theta_t = theta_t;
        }else{
            av_theta_t = ( (t - BURN - 1) * av_theta_prev + theta_t ) / (t - BURN);
            // path_av_theta.row(t) = ( (t - BURN - 1) * path_av_theta.row(t - 1) + path_theta.row(t) ) / (t - BURN);
        }
        av_theta_prev = av_theta_t;

        /////////////////////////////////
        /* STORE ITERATION QUANTITIES  */
        /////////////////////////////////

        if(t == 1 | t == BURN | (t-last_stored_iter) == EACH | t == last_iter){
            iter_idx.push_back(t);
            path_theta.push_back(theta_t);
            path_grad.push_back(ngradient_t);
            path_av_theta.push_back(av_theta_t);
            last_stored_iter = t;
        }

        // path_nll.push_back(ncl);
    }

    clock.tock("main");
    clock.stop("clock");




    Rcpp::List output = Rcpp::List::create(
        Rcpp::Named("path_theta") = path_theta,
        Rcpp::Named("path_av_theta") = path_av_theta,
        Rcpp::Named("path_grad") = path_grad,
        Rcpp::Named("scale") = scale,
        Rcpp::Named("n") = n,
        Rcpp::Named("iter_idx") = iter_idx,
        Rcpp::Named("methodflag") = METHODFLAG
    );

    return output;

}
