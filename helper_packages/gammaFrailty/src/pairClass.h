#ifndef pairClass_H
#define pairClass_H

#include <boost/math/special_functions.hpp>
#include "utils.h"

class pair_class{
    public:
        // Inputs
        unsigned int _j, _jp;        // periods indices
        unsigned int _n_j, _n_jp;    // count values
        double _alpha_j, _alpha_jp;  // intercepts
        Eigen::VectorXd _x;         // external regressors
        Eigen::VectorXd _beta;      // regression coefficients
        double _xi;                 // inverse of gamma parameetrs
        double _rho;                // correlation adjacent periods
        double _lxi;
        double _artanhrho;
        unsigned int _struct;

        // Derived helpers
        unsigned int _p;            // number of periods
        unsigned int _r;            // number of regressors
        unsigned int _d;            // number of parameters
        unsigned int _adiff;        // absolute difference indices
        double _drho_jjp;           // derivative of rho^{|j-j'|}

        // Intermediate quantities
        double _m_1, _m_2;
        double _u_j, _u_jp;
        double _log_u_j, _log_u_jp;
        double _rho_jjp;
        double _D_j, _D_jp;
        double _Delta;
        double _f;
        double _S;
        double _logS;

    //////////////////////////
    //////////////////////////
    /// INTERNALFUNCTIONS  ///
    //////////////////////////
    //////////////////////////

        // SETUP input quantities
        void setup_(unsigned int J,unsigned int JP,unsigned int N_J,unsigned int N_JP,Eigen::VectorXd X,
                    Eigen::VectorXd BETA,double ALPHA_J,double ALPHA_JP,double LXI,double artanhRHO,
                    unsigned int P, unsigned int STRUCT = 0);

        // public intermediate computing (only for debugging)
        double compute_Q_(unsigned int ind, bool verboseFLAG = false);
        double compute_Sind_(unsigned int ind, bool verboseFLAG = false);

        double compute_logQ0_(bool verboseFLAG = false);
        double compute_stdQ_(unsigned int ind, bool verboseFLAG = false);

        // compute intermediate quantities and derivatives
        void compute_intermediate_();
        void compute_dintermediate_(bool verboseS = false);

        // return log-likelihood and gradient contribution
        double compute_ll_();
        double compute_ll_stable_();

        Eigen::VectorXd compute_gradient_();

        // return lists with intermediate quantities (only for debugging)
        Rcpp::List get_intermediate_();
        Rcpp::List get_dintermediate_();

    private:

        // Compute intermediate quantities
        void compute_m_();
        void compute_rho_();
        void compute_u_();
        void compute_Delta_();
        void compute_D_();
        void compute_f_();
        void compute_S_(bool verboseFLAG = false);
        void compute_logS_stable_(bool verboseFLAG = false);
        double compute_Sind_stable_(unsigned int ind, bool verboseFLAG = false);

        /////////////////
        // DERIVATIVES //
        /////////////////
        // (xi, rho,   beta   , alpha_1,..., alpha_p)
        // (0 , 1  , 2,...,r+1, r+2    ,..., r+p+1  )

        Eigen::VectorXd _du_j; Eigen::VectorXd _du_jp;
        void compute_du_();

        double _dDelta_u_j, _dDelta_u_jp; Eigen::VectorXd _dDelta;
        void compute_dDelta_();

        double _dD_u; Eigen::VectorXd _dD_j, _dD_jp;
        void compute_dD_();

        double _df_D_j, _df_D_jp; double _df_Delta; Eigen::VectorXd _df;
        void compute_df_();

        double compute_dQ_xi_dividedbyQ_(unsigned int ind, bool verboseFLAG = false);

        Eigen::VectorXd compute_dSind_(unsigned int ind, bool verboseFLAG = false);

        Eigen::VectorXd _dlogS;
        void compute_dlogS_(bool verboseFLAG = false);
        void compute_dlogS_stable_(bool verboseFLAG = false);


};

// SETUP input quantities
void pair_class::setup_(
        unsigned int J,
        unsigned int JP,
        unsigned int N_J,
        unsigned int N_JP,
        Eigen::VectorXd X,
        Eigen::VectorXd BETA,
        double ALPHA_J,
        double ALPHA_JP,
        double LXI,
        double artanhRHO,
        unsigned int P,
        unsigned int STRUCT
){
    _j = J;
    _jp = JP;
    _n_j = N_J;
    _n_jp = N_JP;
    _x = X;
    _beta = BETA;
    _alpha_j = ALPHA_J;
    _alpha_jp = ALPHA_JP;
    // xi = exp(LXI);
    _xi = exp(-LXI);
    _lxi = LXI;
    //xi = rofz_cpp(LXI) + 1.1;
    //xi = pow(LXI, 2);
    _artanhrho = artanhRHO;
    _rho = rofz_cpp(artanhRHO);
    _p = P;
    _r = BETA.size();
    _d = _p + _r + 2;
    _struct = STRUCT;
}

// public intermediate computing (only for debugging)
double pair_class::compute_Q_(unsigned int ind, bool verboseFLAG){
    double Q = 1;
    for(unsigned int ind2 = _m_2; ind2 <= _m_1 + _m_2 - ind - 1; ind2++){
        Q *= 1 + ind2 * _xi;
    }

    return Q;
}
double pair_class::compute_Sind_(unsigned int ind, bool verboseFLAG){
    double Q = compute_Q_(ind);
    double S = pow(-_xi*_f, ind) * Rf_choose(_m_1, ind) * Rf_choose(_m_2, ind) * boost::math::factorial<double>(ind) * Q;

    if(verboseFLAG)Rcpp::Rcout << "pow(-xi*f, ind):"<<pow(-_xi*_f, ind)<<"\nRf_choose(m_1, ind):"<<Rf_choose(_m_1, ind)
                               <<"\nRf_choose(m_2, ind):"<<Rf_choose(_m_2, ind)<<"\nboost::math::factorial<double>(ind):"
                               << boost::math::factorial<double>(ind) << "\nQ:"<<Q<<"\n";
    return S;
}

// compute intermediate quantities and derivatives
void pair_class::compute_intermediate_(){
    compute_m_();
    compute_rho_();
    compute_u_();
    compute_Delta_();
    compute_D_();
    compute_f_();
    //compute_S_();
}
void pair_class::compute_dintermediate_(bool verboseS){

    compute_du_();
    compute_dDelta_();
    compute_dD_();
    compute_df_();
    compute_dlogS_(verboseS);
}

// return log-likelihood and gradient contribution
double pair_class::compute_ll_(){

    double out = 0;
    double tmp0 = 0;
    for(unsigned int ind = 0; ind < _m_2; ind++){
        tmp0 += log(1 + ind*_xi);
    }
    double tmp1 = _n_j*_log_u_j + _n_jp*_log_u_jp +
        _n_j*log(_D_j) + _n_jp*log(_D_jp) - (_n_j + _n_jp + 1/_xi)*log(_Delta);

    out = tmp0 + tmp1 + log(_S);

    return out;
}
Eigen::VectorXd pair_class::compute_gradient_(){
    Eigen::VectorXd out = Eigen::VectorXd::Zero(_d);

    double tmp_xi= 0;
    for(unsigned int ind = 0; ind < _m_2; ind++){
        tmp_xi += ind/(1 + ind*_xi);
    }

    Eigen::VectorXd tmp0 = -pow(_xi*_Delta,-1)*_dDelta; tmp0(0) += pow(_xi, -2)*log(_Delta);
    Eigen::VectorXd tmp1 = (_n_j/_u_j)*_du_j + (_n_jp/_u_jp)*_du_jp + (_n_j/_D_j)*_dD_j + (_n_jp/_D_jp)*_dD_jp - ((_n_j+_n_jp)/_Delta)*_dDelta;
    out = tmp0 + tmp1 + _dlogS;
    out[0] += tmp_xi;

    // out[0] *= xi;
    out[0] *= -_xi;
    //out[0] *= 2*lxi;
    //out[0] *= drofz_cpp(lxi);
    out[1] *= drofz_cpp(_artanhrho);

    return out;
}

// return lists with intermediate quantities (only for debugging)
Rcpp::List pair_class::get_intermediate_(){

    Rcpp::List output =
        Rcpp::List::create(
            Rcpp::Named("m_1"    ) = _m_1,
            Rcpp::Named("m_2"    ) = _m_2,
            Rcpp::Named("u_j"    ) = _u_j,
            Rcpp::Named("u_jp"   ) = _u_jp,
            Rcpp::Named("rho_jjp") = _rho_jjp,
            Rcpp::Named("D_j"    ) = _D_j,
            Rcpp::Named("D_jp"   ) = _D_jp,
            Rcpp::Named("Delta"  ) = _Delta,
            Rcpp::Named("f"      ) = _f,
            Rcpp::Named("S"      ) = _S
        );

    return output;
}
Rcpp::List pair_class::get_dintermediate_(){

    Rcpp::List output =
        Rcpp::List::create(
            Rcpp::Named("du_j"    ) = _du_j,
            Rcpp::Named("du_jp"   ) = _du_jp,
            Rcpp::Named("drho_jjp") = _drho_jjp,
            Rcpp::Named("dD_j"    ) = _dD_j,
            Rcpp::Named("dD_jp"   ) = _dD_jp,
            Rcpp::Named("dDelta"  ) = _dDelta,
            Rcpp::Named("df"      ) = _df,
            Rcpp::Named("dlogS"   ) = _dlogS
        );

    return output;
}

// Compute intermediate quantities
void pair_class::compute_m_(){
    _m_1 = std::min(_n_j, _n_jp);
    _m_2 = std::max(_n_j, _n_jp);

}
void pair_class::compute_rho_(){
    int diff;
    // if(_struct==0){
    //     diff = _j-_jp;
    // }elseif(_struct==1){
    //     diff = 1;
    // }
    switch(_struct){
    case 0:
        diff = _j-_jp;
        break;
    case 1:
        diff = 1;
        break;
    }
    _adiff = abs(diff);
    _rho_jjp = pow(_rho, _adiff);
    _drho_jjp = pow(_rho, _adiff-1)*_adiff;
}
void pair_class::compute_u_(){
    _log_u_j = _alpha_j + _beta.dot(_x);
    _u_j = exp(_log_u_j);
    _log_u_jp = _alpha_jp + _beta.dot(_x);
    _u_jp = exp(_log_u_jp);
}
void pair_class::compute_Delta_(){
    _Delta = 1+ _xi*(_u_j+_u_jp) + pow(_xi,2)*_u_j*_u_jp*(1-_rho_jjp);
}
void pair_class::compute_D_(){
    _D_j = 1 + _xi*_u_jp*(1-_rho_jjp);
    _D_jp = 1 + _xi*_u_j*(1-_rho_jjp);
}
void pair_class::compute_f_(){
    _f = (_Delta*(1-_rho_jjp))/(_D_j*_D_jp);
}
void pair_class::compute_S_(bool verboseFLAG){
    double tmp = 0;
    for(unsigned int ind = 0; ind <= _m_1; ind++){
        tmp += compute_Sind_(ind);
    }

    _S = tmp;
}
//// DERIVATIVES ////
// _du_j, _du_jp;
void pair_class::compute_du_(){
    _du_j.resize(_d); _du_jp.resize(_d);

    // wrt xi and rho
    _du_j.setZero(); _du_jp.setZero();

    // wrt beta
    _du_j.segment(2,_r) = _u_j*_x; _du_jp.segment(2,_r) = _u_jp*_x;

    // wrt alpha
    _du_j(_r+2+_j) = _u_j; _du_jp(_r+2+_jp) = _u_jp;
}
// _dDelta_u_j, _dDelta_u_jp, _dDelta;
void pair_class::compute_dDelta_(){
    _dDelta.resize(_d);
    _dDelta.setZero();

    // partial: wrt u_j and u_jp
    _dDelta_u_j = _xi*_D_j; _dDelta_u_jp = _xi*_D_jp;

    // wrt xi
    _dDelta(0) = _u_j + _u_jp + 2*_xi*_u_j*_u_jp*(1-_rho_jjp);

    // wrt rho
    _dDelta(1) = -pow(_xi,2)*_u_j*_u_jp*_drho_jjp;

    // wrt beta
    _dDelta.segment(2,_r) = _dDelta_u_j * _du_j.segment(2,_r) + _dDelta_u_jp * _du_jp.segment(2,_r);

    // wrt alpha
    _dDelta(_r+2+_j ) = _dDelta_u_j  * _du_j( _r+2+_j);
    _dDelta(_r+2+_jp) = _dDelta_u_jp * _du_jp(_r+2+_jp);
}
// _dD_u _dD_j, _dD_jp;
void pair_class::compute_dD_(){
    _dD_j.resize(_d); _dD_jp.resize(_d);
    _dD_j.setZero(); _dD_jp.setZero();

    // partial: D_j wrt u_jp and D_jp wrt u_j
    _dD_u = _xi*(1-_rho_jjp);

    // wrt xi
    _dD_j(0)  = _u_jp*(1-_rho_jjp);
    _dD_jp(0) = _u_j *(1-_rho_jjp);

    // wrt rho
    _dD_j(1)  = -_xi*_u_jp*_drho_jjp;
    _dD_jp(1) = -_xi*_u_j *_drho_jjp;

    // wrt beta
    _dD_j.segment(2,_r)  = _dD_u*_du_jp.segment(2,_r);
    _dD_jp.segment(2,_r) = _dD_u* _du_j.segment(2,_r);

    // wrt alpha
    _dD_j(_r+2+_jp) = _dD_u*_du_jp(_r+2+_jp);
    _dD_jp(_r+2+_j) = _dD_u*_du_j( _r+2+_j );
}
// _df_D_j, _df_D_jp;  _df_Delta, _df;
void pair_class::compute_df_(){
    _df.resize(_d);
    _df.setZero();

    // partial: wrt D_j and D_jp
    _df_D_j  = -_Delta * ( 1 - _rho_jjp )/( pow(_D_j, 2) * _D_jp);
    _df_D_jp = -_Delta * ( 1 - _rho_jjp )/( pow(_D_jp,2) * _D_j );

    // partial: wrt Delta
    _df_Delta = (1 - _rho_jjp)/(_D_j*_D_jp);

    // wrt xi
    _df(0) =  _df_Delta * _dDelta(0) + _df_D_j * _dD_j(0) + _df_D_jp * _dD_jp(0);

    // wrt rho
    _df(1) = _df_Delta * _dDelta(1) + _df_D_j * _dD_j(1) + _df_D_jp * _dD_jp(1) - (_Delta/(_D_j*_D_jp))*_drho_jjp;

    // wrt beta
    _df.segment(2,_r) = _df_Delta*_dDelta.segment(2,_r) + _df_D_j*_dD_j.segment(2,_r) + _df_D_jp*_dD_jp.segment(2,_r);

    // wrt alpha
    _df(_r+2+_j ) = _df_Delta*_dDelta(_r+2+_j ) + _df_D_jp*_dD_jp(_r+2+_j);
    _df(_r+2+_jp) = _df_Delta*_dDelta(_r+2+_jp) + _df_D_j *_dD_j(_r+2+_jp);
}

double pair_class::compute_dQ_xi_dividedbyQ_(unsigned int ind, bool verboseFLAG){
    double tmp_sum = 0;
    for(unsigned int ind2 = _m_2; ind2 <= _m_1 + _m_2 - ind - 1; ind2++){
        double tmp = 1 + ind2*_xi;
        tmp_sum += ind2/tmp;
    }

    return tmp_sum;
}

Eigen::VectorXd pair_class::compute_dSind_(unsigned int ind, bool verboseFLAG){
    Eigen::VectorXd out = Eigen::VectorXd::Zero(_d);

    double Sind = compute_Sind_(ind);

    // wrt xi
    out(0) = Sind*(compute_dQ_xi_dividedbyQ_(ind) + ind*pow(_xi,-1) + ind*pow(_f,-1)*_df(0));

    // wrt rho
    out(1) = Sind*(ind/_f)*_df(1);

    // wrt beta
    out.segment(2,_r) = Sind*(ind/_f)*_df.segment(2,_r);

    // wrt alpha
    out(_r+2+_jp) = Sind*(ind/_f)*_df(_r+2+_jp);
    out(_r+2+_j)  = Sind*(ind/_f)*_df(_r+2+_j );

    return out;
}

// _dlogS;
void pair_class::compute_dlogS_(bool verboseFLAG){
    _dlogS.resize(_d);
    _dlogS.setZero();

    double tmp_sum = 0;
    Eigen::VectorXd tmp_dSum = Eigen::VectorXd::Zero(_d);

    for(unsigned int ind = 0; ind <= _m_1; ind++){
        double Sind = compute_Sind_(ind);
        if(verboseFLAG) Rcpp::Rcout << "index:" << ind << " Sind:" << Sind << "\n";
        tmp_sum  += Sind;

        Eigen::VectorXd out = Eigen::VectorXd::Zero(_d);
        {
            // wrt xi
            out(0) = Sind*(compute_dQ_xi_dividedbyQ_(ind) + ind*pow(_xi,-1) + ind*pow(_f,-1)*_df(0));

            // wrt rho
            out(1) = Sind*(ind/_f)*_df(1);

            // wrt beta
            out.segment(2,_r) = Sind*(ind/_f)*_df.segment(2,_r);

            // wrt alpha
            out(_r+2+_jp) = Sind*(ind/_f)*_df(_r+2+_jp);
            out(_r+2+_j)  = Sind*(ind/_f)*_df(_r+2+_j );
        }

        tmp_dSum += out;
    }


    _S = tmp_sum;
    _dlogS = tmp_dSum/_S;
}

// functions for numerical stability
double pair_class::compute_logQ0_(bool verboseFLAG){
    double out = 0;
    for(unsigned int ind2 = _m_2; ind2 <= _m_1 + _m_2 - 1; ind2++){
        out += log(1 + ind2 * _xi);
    }

    return out;
}
double pair_class::compute_stdQ_(unsigned int ind, bool verboseFLAG){
    double Q = 1;
    for(unsigned int ind2 = _m_1 + _m_2 - ind; ind2 <= _m_1 + _m_2 - 1; ind2++){
        Q *= 1 + ind2 * _xi;
    }

    return Q;
}
double pair_class::compute_Sind_stable_(unsigned int ind, bool verboseFLAG){
    double std_fct = boost::math::factorial<double>(ind)/compute_stdQ_(ind);
    double out = pow(-_xi*_f, ind) * Rf_choose(_m_1, ind) * Rf_choose(_m_2, ind) * std_fct;

    if(verboseFLAG)Rcpp::Rcout << "pow(-xi*f, ind):"<<pow(-_xi*_f, ind)<<"\nRf_choose(m_1, ind):"<<Rf_choose(_m_1, ind)
                               <<"\nRf_choose(m_2, ind):"<<Rf_choose(_m_2, ind)<<"\nboost::math::factorial<double>(ind):"
                               << boost::math::factorial<double>(ind) << "\n std fct:"<<std_fct<<"\n";
    return out;
}
void pair_class::compute_logS_stable_(bool verboseFLAG){
    double out = compute_logQ0_();
    double tmp = 0;
    for(unsigned int ind = 1; ind <= _m_1; ind++){
        tmp += compute_Sind_stable_(ind);
    }

    out += log1p(tmp);
    _logS = out;
    _S = exp(out);
}
double pair_class::compute_ll_stable_(){

    double out = 0;
    double tmp0 = 0;
    for(unsigned int ind = 0; ind < _m_2; ind++){
        tmp0 += log(1 + ind*_xi);
    }
    double tmp1 = _n_j*_log_u_j + _n_jp*_log_u_jp +
        _n_j*log(_D_j) + _n_jp*log(_D_jp) - (_n_j + _n_jp + 1/_xi)*log(_Delta);

    out = tmp0 + tmp1 + _logS;

    return out;
}
void pair_class::compute_dlogS_stable_(bool verboseFLAG){
    _dlogS.resize(_d);
    _dlogS.setZero();

    double tmp_sum = 0;
    Eigen::VectorXd tmp_dSum = Eigen::VectorXd::Zero(_d);

    for(unsigned int ind = 0; ind <= _m_1; ind++){
        double Sind = compute_Sind_(ind);
        if(verboseFLAG) Rcpp::Rcout << "index:" << ind << " Sind:" << Sind << "\n";
        tmp_sum  += Sind;

        Eigen::VectorXd out = Eigen::VectorXd::Zero(_d);
        {
            // wrt xi
            out(0) = Sind*(compute_dQ_xi_dividedbyQ_(ind) + ind*pow(_xi,-1) + ind*pow(_f,-1)*_df(0));

            // wrt rho
            out(1) = Sind*(ind/_f)*_df(1);

            // wrt beta
            out.segment(2,_r) = Sind*(ind/_f)*_df.segment(2,_r);

            // wrt alpha
            out(_r+2+_jp) = Sind*(ind/_f)*_df(_r+2+_jp);
            out(_r+2+_j)  = Sind*(ind/_f)*_df(_r+2+_j );
        }

        tmp_dSum += out;
    }


    _S = tmp_sum;
    _dlogS = tmp_dSum/_S;
}

#endif
