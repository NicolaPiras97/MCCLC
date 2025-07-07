// [[Rcpp::depends(BH)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppEigen)]]
#include <iostream>
#include </opt/homebrew/Cellar/boost/1.81.0_1/include/boost/random.hpp>
#include <cmath>
#include <vector>
#include <string>
#include <armadillo>
#include</opt/homebrew/Cellar/eigen/3.4.0_1/include/eigen3/unsupported/Eigen/CXX11/Tensor>
#define _USE_MATH_DEFINES
#include <chrono>


using namespace Eigen;
using namespace std;
using namespace arma;


tuple <mat,int,int> vector_to_matrix(const std::vector< std::vector<int> > &student_vector)
{
  mat y(student_vector.size(),student_vector[0].size());
  
  for(int i = 0; i < student_vector.size(); i++)
  {
    for(int j = 0; j < student_vector[i].size(); j++)
    {
      y(i,j) = student_vector[i][j];
    }
  }
  
  int _n=student_vector.size();
  int _G=student_vector[0].size();
  
  return {y,_n,_G};
}


// [[Rcpp::export]]
int main2(vector<vector<int>> yy){
  mat y;
  int _n;
  int _G;
  tie(y,_n,_G)=vector_to_matrix(yy);
  
    auto start = chrono::steady_clock::now();
    std::vector< std::vector<int> > student_vector;


    int _L,_H,_R,_K,_Q,_iter,_iter0,_iter3,_starting,_burn,_thinn;
    _L=4;
    _H=3;
    _R=2;
    _iter=800;
    _iter0=2000;
    _iter3=500;
    _starting=16;
    _burn=200;
    _thinn=3;
    double eps0=0.000000001;
    
    
    _K = y(0, 1);
    for (int j = 0; j < _n; j++) {
        if (y(j, 1) > _K) {
            _K = y(j, 1);
        }
    }
    _Q = y(0, 2);
    for (int j = 0; j < _n; j++) {
        if (y(j, 2) > _Q) {
            _Q = y(j, 2);
        }
    }

  
    // thinning
    int nn=0;
    int nnn=_iter-_burn;
    while(_thinn<=nnn){
        nn=nn+1;
        nnn=nnn-_thinn;
    }
    vector<int> it(nn);
    it[0]=_burn+_thinn-1;
    for(int j=1;j<nn;j++){
        it[j]=it[j-1]+_thinn;
    }
    int mu=0;

    
// _vettoreqk contains the n observations ordered respect to Q
    vector<int> _vettoreqk(_n);
    int xx=0;
    for(int q=0;q<_Q;q++){
        for(int j=0;j<_n;j++){
            if(y(j,2)==q+1){
                _vettoreqk[xx]=y(j,0);
                xx=xx+1;
            }
        }
    }
   

// _indvecqk contains the positions of vettoreqk according to the previous order
    vector<int> _indvecqk(_n);
    int xxx=0;
    for(int j=0;j<_n;j++){
        if(y(j,0)==_vettoreqk[xxx]){
            _indvecqk[xxx]=j;
            xxx=xxx+1;
            if(xxx<=_n-1){
                j=-1;
            }
            else if(xxx==_n){
                j=_n-1;
            }
        }
    }
  
   

    // _nk counts how many observations fall in each possible k
    vector<int> _nk(_K);
    // _nq counts how many observations fall in each possible q
    vector<int> _nq(_Q);
    //_nkq contains the number of units in each combination of k and q
    mat _nkq(_K,_Q);

    for(int k=0;k<_K;k++){
        _nk[k]=0;
    }
    for(int k=0;k<_K;k++){
        for(int j=0;j<_n;j++){
            if(y(j,1)==k+1){
                _nk[k]+=1;
            }
        }
    }

    for(int q=0;q<_Q;q++){
        _nq[q]=0;
    }
    for(int q=0;q<_Q;q++){
        for(int j=0;j<_n;j++){
            if(y(j,2)==q+1){
                _nq[q]+=1;
            }
        }
    }

    for(int k=0;k<_K;k++){
        for(int q=0;q<_Q;q++){
            _nkq(k,q)=0;
        }
    }

    for(int k=0;k<_K;k++){
        for(int q=0;q<_Q;q++){
            for(int j=0;j<_n;j++){
                if(y(j,1)==k+1 && y(j,2)==q+1){
                    _nkq(k,q)+=1;
                }
            }
        }
    }

    int ckq=0;
    for(int k=0;k<_K;k++) {
        for (int q = 0; q < _Q; q++) {
            if(_nkq(k,q)!=0){
                ckq+=1;
            }
        }
    }

   
    //initialization of _pi_l 1/L
    arma::cube _pi_l(_L, _H, _R);
    for (int h = 0; h < _H; h++) {
        for (int r = 0; r < _R; r++) {
            for (int l = 0; l < _L; l++) {
                _pi_l(l, h, r) = (double) 1 / _L;
            }
        }
    }

    //initialization of p
    int _I = _G - 3;
    vector<int> vC(_I);
    int ma;
    for(int i=0;i<_I;i++){
        ma = y(0,i+3);
        for (int j = 0; j < _n; j++) {
            if (y(j,i+3) > ma) {
                ma = y(j,i+3);
            }
        }
        vC[i]=ma+1;
    }

    int _C; //max between the values that the variables assume +1
    int maa;
    maa = vC[0];
    for (int i = 0; i < _I; i++) {
        if (vC[i] > maa) {
            maa = vC[i];
        }
    }
    _C=maa;
   

    
    mat _Prob1(_K, _H);
    mat _Prob2(_Q, _R);
    //boost::mt19937 generator2(1234);
    std::random_device rdtest;
    std::mt19937 generator(rdtest());
    mat _w(_K, _Q);
    mat _z(_Q, _R);
    vector<int> _numh(_H);
    vector<double> _p_h(_H);
    vector<int> _numr(_R);
    vector<double> _p_r(_R);

    vector<double> _p_rc_hat(_R);
    std::fill(_p_rc_hat.begin(), _p_rc_hat.end(), 0);
    
    vector<double> _p_hc_hat(_H);
    std::fill(_p_hc_hat.begin(), _p_hc_hat.end(), 0);
    
    mat _w_hat(_K, _H);
    mat _z_hat(_Q, _R);
    for (int k = 0; k < _K; k++) {
        for (int h = 0; h < _H; h++) {
            _w_hat(k, h) = 0;
        }
    }

    for (int q = 0; q < _Q; q++) {
        for (int r = 0; r < _R; r++) {
            _z_hat(q, r) = 0;
        }
    }

    //definitions

    Tensor<double, 4> _f(_n, _L, _H, _R);
    //Tensor<double,4> _f(_n,_L,_H,_R);
    Tensor<double, 4> _res(_K, _Q, _H, _R);
    mat _ff(_K, _H);
    mat _num(_n, _L);
    vector<double> _den(_n);
    vector<double> _denom(_n);
    mat _Prob3(_n, _L);
    mat _num2(_K, _H);
    vector<double> _den2(_K);
    mat _ff2(_Q, _R);
    Tensor<double, 4> _res2(_Q, _K, _H, _R);
    mat _num3(_Q, _R);
    vector<double> _den3(_Q);
    mat _wjk(_n, _H);//extension of w in which the k-th row of w is repeted _nk[k] times
    mat _zjq(_n, _R);//extension of z in which the q-th row of z is repeted _nq[q] times
    arma::cube _B(_n, _H, _R);
    mat _den4(_H, _R);
    //Tensor<double,4> _aa(_n,_L,_H,_R);
    arma::cube _num4(_L, _H, _R);
    vector<double> _den5(_L);
    mat _den6(_L,_R);
    mat _den7(_L,_H);
    //arma::cube _den5(_L,_H,_R);
    //mat _num5(_I,_L);
    Tensor<double,5> _num5(_I, _C, _L, _H, _R);
    //Tensor<double,5> _aaa(_I,_n,_L,_H,_R);
    Tensor<double,4> _P(_starting,_I,_C,_L);
    Tensor<double,2> _PL(_starting,_L);
    vector<int> Seed(_starting);
    vector<double> logLik(_starting);
    Tensor<double, 4> _pf(_n, _L, _H, _R);
    arma::cube _fykq(_n, _H, _R);
    Tensor<double, 4> _res3(_n, _L, _H, _R);
    mat _xjkq(_n, _L);
    mat _tt(_K,_Q);

    mat _x_hat(_n, _L);
    _x_hat.zeros(_n, _L);

    mat _x_hat_int(_n, _L);
    _x_hat_int.zeros(_n, _L);

    //classification at the second level
    vector<int> _indh(_K);
    vector<int> _indr(_Q);
    mat _classif2(_n, 3);

    //classification at the first level
    vector<int> _indl(_n);

    mat _classif(_n, 4);

    double _maxL;

    //final estimates
    arma::cube _pi_l_hat(_L, _H, _R);
    for (int l = 0; l < _L; l++) {
        for (int h = 0; h < _H; h++) {
            for (int r = 0; r < _R; r++) {
                _pi_l_hat(l, h, r) = 0;
            }
        }
    }

    Tensor<double, 5> _p_hat(_I, _C, _L, _H, _R);
    for (int i = 0; i < _I; i++) {
        for (int c = 0; c < _C; c++) {
            for (int l = 0; l < _L; l++) {
                for (int h = 0; h < _H; h++) {
                    for (int r = 0; r < _R; r++) {
                        _p_hat(i, c, l, h, r) = 0;
                    }
                }
            }
        }
    }

    vector<double> _p_h_hat(_H);
    std::fill(_p_h_hat.begin(), _p_h_hat.end(), 0);
    vector<double> _p_r_hat(_R);
    std::fill(_p_r_hat.begin(), _p_r_hat.end(), 0);

    //I save the parameters for each iteratio and each starting point
    mat _p_h_all(_iter, _H);
    mat _p_r_all(_iter, _R);
    //vector<vector<vector<double>>> _p_h_all;
    mat _logL_all(_starting, _iter);

    arma_rng::set_seed_random();
    arma::cube _p1(_I, _C, _L, fill::randu);
    Tensor<double, 5> _p(_I, _C, _L, _H, _R);
    //_p.setRandom();
    for (int i = 0; i < _I; i++) {
        for (int c = 0; c < _C; c++) {
            for (int l = 0; l < _L; l++) {
                for (int h = 0; h < _H; h++) {
                    for (int r = 0; r < _R; r++) {
                        _p(i, c, l, h, r) = 0;
                    }
                }
            }
        }
    }

    vector<double> _pl(_L);
    for (int l = 0; l < _L; l++) {
        _pl[l] = (double) 1 / _L;
    }

    for (int h = 0; h < _H; h++) {
        for (int r = 0; r < _R; r++) {
            for (int l = 0; l < _L; l++) {
                _pi_l(l, h, r) = (double) 1 / _L;
            }
        }
    }

    
    vector<double> _de(_n);
    vector<double> _de2(_L);
    mat _tau(_n, _L);

    for (int st = 0; st < _starting; st++) {

        arma::cube _pb(_I, _C, _L);
        double eqprob;
        vector<double> vec(_C * _L * _I);
        eqprob = (double) 1 / (_C * _L * _I);
        std::fill(vec.begin(), vec.end(), eqprob);
        std::random_device rdtest;
        std::mt19937 generator2(rdtest());
        boost::random::discrete_distribution<int> distribution3(vec.begin(), vec.end());
        for (int i = 0; i < _I; i++) {
            for (int c = 0; c < _C; c++) {
                for (int l = 0; l < _L; l++) {
                    int sample = distribution3(generator2);
                    _pb(i, c, l) = sample + 1;
                }
            }
        }

        mat sumc(_I, _L);
        for (int i = 0; i < _I; i++) {
            for (int l = 0; l < _L; l++) {
                sumc(i, l) = 0;
            }
        }
        for (int i = 0; i < _I; i++) {
            for (int l = 0; l < _L; l++) {
                for (int c = 0; c < vC[i]; c++) {
                    sumc(i, l) += _pb(i, c, l);
                }
            }
        }
        for (int i = 0; i < _I; i++) {
            for (int l = 0; l < _L; l++) {
                for (int c = 0; c < vC[i]; c++) {
                    _p1(i, c, l) = (double) _pb(i, c, l) / sumc(i, l);
                }
            }
        }

        
        vector<int> _z0(_L);
        double eqprob0;
        vector<double> vec0(_L);
        eqprob0 = (double) 1 / (_L);
        std::fill(vec0.begin(), vec0.end(), eqprob0);
        boost::random::discrete_distribution<int> distribution0(vec0.begin(), vec0.end());
        for (int l = 0; l < _L; l++) {
            int sample = distribution0(generator2);
            _z0[l] = sample + 1;
        }

        int sumc0=0;
        for (int l = 0; l < _L; l++) {
            sumc0+= _z0[l];
        }

        for (int l = 0; l < _L; l++) {
            _pl[l]= (double) _z0[l]/ sumc0;
        }

        mat _f1(_n, _L);
        for (int j = 0; j < _n; j++) {
            for (int l = 0; l < _L; l++) {
                _f1(j, l) = 1;
            }
        }

        for (int j = 0; j < _n; j++) {
            for (int l = 0; l < _L; l++) {
                for (int i = 0; i < _I; i++) {
                    _f1(j, l) = _f1(j, l) * _p1(i, y(j, i + 3), l);
                }
            }
        }

        double L, L1;
        vector<double> _flog(_n);
        for (int j = 0; j < _n; j++) {
            for (int l = 0; l < _L; l++) {
                _flog[j] += _pl[l] * _f1(j, l);
            }
        }
        for (int j = 0; j < _n; j++) {
            _flog[j] = log(_flog[j]);
        }
        L1 = 0;
        for (int j = 0; j < _n; j++) {
            L1 += _flog[j];
        }

        for (int u = 0; u < _iter0; u++) {

          mat num0(_n, _L);
          for (int j = 0; j < _n; j++) {
            for (int l = 0; l < _L; l++) {
              num0(j, l) = (double) log(_f1(j, l)) + log(_pl[l]);
            }
          }
          vector<double> _max0(_n);
          for (int j = 0; j < _n; j++) {
            _max0[j] = num0(j, 0);
            for (int l = 0; l < _L; l++) {
              if (num0(j, l) > _max0[j]) {
                _max0[j] = num0(j, l);
              }
            }
          }
          
          for (int j = 0; j < _n; j++) {
            for (int l = 0; l < _L; l++) {
              num0(j, l) = num0(j, l) - _max0[j];
            }
          }
          
          for (int j = 0; j < _n; j++) {
            for (int l = 0; l < _L; l++) {
              num0(j, l) = exp(num0(j, l));
            }
          }
          
          vector<double> _den0(_n);
          for (int j = 0; j < _n; j++) {
            _den0[j] = 0;
          }
          
          for (int j = 0; j < _n; j++) {
            for (int l = 0; l < _L; l++) {
              _den0[j] += num0(j, l);
            }
          }
          
          for (int j = 0; j < _n; j++) {
            for (int l = 0; l < _L; l++) {
              _tau(j, l) = (double) num0(j, l) / _den0[j];
            }
          }

            for (int l = 0; l < _L; l++) {
                _pl[l] = 0;
            }

            for (int l = 0; l < _L; l++) {
                for (int j = 0; j < _n; j++) {
                    _pl[l] += _tau(j, l);
                }
            }

            for (int l = 0; l < _L; l++) {
                _pl[l] = _pl[l] / _n;
            }


            for (int l = 0; l < _L; l++) {
                _de2[l] = 0;
            }

            for (int l = 0; l < _L; l++) {
                for (int j = 0; j < _n; j++) {
                    _de2[l] += _tau(j, l);
                }
            }

            for (int i = 0; i < _I; i++) {
                for (int c = 0; c < _C; c++) {
                    for (int l = 0; l < _L; l++) {
                        _p1(i, c, l) = 0;
                    }
                }
            }

            for (int i = 0; i < _I; i++) {
                for (int l = 0; l < _L; l++) {
                    for (int j = 0; j < _n; j++) {
                        _p1(i, y(j, i + 3), l) += _tau(j, l);
                    }
                }
            }

            for (int i = 0; i < _I; i++) {
                for (int c = 0; c < _C; c++) {
                    for (int l = 0; l < _L; l++) {
                        _p1(i, c, l) = (double) _p1(i, c, l) / _de2[l];
                    }
                }
            }

            for (int i = 0; i < _I; i++) {
              for (int l = 0; l < _L; l++) {
                for (int c = 0; c < _C; c++) {
                  if (_p1(i, c, l) > 0.9999) {
                    _p1(i, c, l) = 0.999;
                  }
                }
              }
            }
            for (int i = 0; i < _I; i++) {
              for (int l = 0; l < _L; l++) {
                for (int c = 0; c < _C; c++) {
                  if (_p1(i, c, l) < 0.0001) {
                    _p1(i, c, l) = 0.001;
                  }
                }
              }
            }
            for (int j = 0; j < _n; j++) {
                for (int l = 0; l < _L; l++) {
                    _f1(j, l) = 1;
                }
            }

            for (int j = 0; j < _n; j++) {
                for (int l = 0; l < _L; l++) {
                    for (int i = 0; i < _I; i++) {
                        _f1(j, l) = _f1(j, l) * _p1(i, y(j, i + 3), l);
                    }
                }
            }

            for (int j = 0; j < _n; j++) {
                _flog[j] = 0;
            }
            for (int j = 0; j < _n; j++) {
                for (int l = 0; l < _L; l++) {
                    _flog[j] += _pl[l] * _f1(j, l);
                }
            }
            for (int j = 0; j < _n; j++) {
                _flog[j] = log(_flog[j]);
            }
            L = 0;
            for (int j = 0; j < _n; j++) {
                L += _flog[j];
            }
            if (abs(L1 - L) < eps0) {
                u = _iter0 - 1;
            }
            L1 = L;

        }
        
        double min1,temp1;
        vector<double> tempp(_C);
        for(int l=0;l<_L-1;l++){
          min1=l;
          for(int j=l+1;j<_L;j++){
            if(_pl[j]<_pl[min1]){
              min1=j;
            }
          }
          temp1=_pl[l];
          _pl[l]=_pl[min1];
          _pl[min1]=temp1;
          for (int i = 0; i < _I; i++) {
            for (int c = 0; c < _C; c++) {
              tempp[c] = _p1(i, c, l);
              _p1(i,c, l) = _p1(i,c, min1);
              _p1(i,c, min1) = tempp[c];
            }
          }
        }

        Seed[st] = rdtest();
        logLik[st] = L;
        for (int i = 0; i < _I; i++) {
            for (int c = 0; c < _C; c++) {
                for (int l = 0; l < _L; l++) {
                    _P(st,i, c, l) = _p1(i, c, l);
                }
            }
        }

        for (int l = 0; l < _L; l++) {
            _PL(st,l)=_pl[l];
        }

    }

    double _maxlogl;
    int indicestarting;
    _maxlogl = logLik[0];
    for (int st = 0; st < _starting; st++) {
        if (logLik[st] > _maxlogl) {
            _maxlogl = logLik[st];
            indicestarting = st;
        }
    }

    for (int i = 0; i < _I; i++) {
        for (int c = 0; c < _C; c++) {
            for (int l = 0; l < _L; l++) {
                for (int h = 0; h < _H; h++) {
                    for (int r = 0; r < _R; r++) {
                        _p(i, c, l, h, r) = _P(indicestarting, i, c, l);
                    }
                }
            }
        }
    }

   
   
   
   
    //check about empty clusters
    vector<int> _cv(_H);
    int _pcv=0;
    vector<int> _cv2(_R);
    int _pcv2=0;

        while(_pcv==0) {

            _w.zeros(_K, _H);
            vector<double> vecc(_H);
            double eqprobb = (double) 1 / _H;
            std::fill(vecc.begin(), vecc.end(), eqprobb);
            //boost::mt19937 generator(1234);
            std::random_device rdtest;
            std::mt19937 generator(rdtest());
            boost::random::discrete_distribution<int> distribution(vecc.begin(), vecc.end());
            for (int k = 0; k < _K; k++) {
                int sample = distribution(generator);
                _w(k, sample) = 1;
            }

            for (int h = 0; h < _H; h++) {
                _cv[h] = 0;
            }
            _pcv = 1;
            for(int h=0;h<_H;h++){
                for(int k=0;k<_K;k++){
                    _cv[h]+=_w(k,h);
                }
            }
            for(int h=0;h<_H;h++){
                _pcv=_pcv*_cv[h];
            }
        }

        while(_pcv2==0){
            for (int r = 0; r < _R; r++) {
                _cv2[r] = 0;
            }
            _pcv2 = 1;
        _z.zeros(_Q, _R);
        vector<double> vec2(_R);
        double eqprob2 = (double) 1 / _R;
        std::fill(vec2.begin(), vec2.end(), eqprob2);
        //boost::mt19937 generator2(1234);
        boost::random::discrete_distribution<int> distribution2(vec2.begin(), vec2.end());
        for (int q = 0; q < _Q; q++) {
            int sample = distribution2(generator);
            _z(q, sample) = 1;
        }

            for(int r=0;r<_R;r++){
                for(int q=0;q<_Q;q++){
                    _cv2[r]+=_z(q,r);
                }
            }
            for(int r=0;r<_R;r++){
                _pcv2=_pcv2*_cv2[r];
            }
        }

    std::fill(_numh.begin(), _numh.end(), 0);

    for (int h = 0; h < _H; h++) {
        for (int k = 0; k < _K; k++) {
            _numh[h] += _w(k, h);
        }
    }


    for (int h = 0; h < _H; h++) {
        _p_h[h] = (double) _numh[h] / _K;
    }
    
    double min1,temp1;
    vector<int> temp2(_K);
   

    std::fill(_numr.begin(), _numr.end(), 0);

    for (int r = 0; r < _R; r++) {
        for (int q = 0; q < _Q; q++) {
            _numr[r] += _z(q, r);
        }
    }

  

    for (int r = 0; r < _R; r++) {
        _p_r[r] = (double) _numr[r] / _Q;
    }
    
    double min2,temp3;
    vector<int> temp4(_Q);
    
    
    
   
    for(int u1=0;u1<_iter;u1++) {

        for (int j = 0; j < _n; j++) {
            for (int l = 0; l < _L; l++) {
                for (int h = 0; h < _H; h++) {
                    for (int r = 0; r < _R; r++) {
                        _f(j, l, h, r) = 1;
                    }
                }
            }
        }

        for (int j = 0; j < _n; j++) {
            for (int l = 0; l < _L; l++) {
                for (int h = 0; h < _H; h++) {
                    for (int r = 0; r < _R; r++) {
                        for (int i = 0; i < _I; i++) {
                            _f(j, l, h, r) = _f(j, l, h, r) * _p(i, y(j, i + 3), l, h, r);
                        }
                    }
                }
            }
        }

        for (int j = 0; j < _n; j++) {
            for (int h = 0; h < _H; h++) {
                for (int r = 0; r < _R; r++) {
                    for (int l = 0; l < _L; l++) {
                        _pf(j, l, h, r) = (double) _pi_l(l, h, r) * _f(j, l, h, r);
                    }
                }
            }
        }

        for (int j = 0; j < _n; j++) {
            for (int h = 0; h < _H; h++) {
                for (int r = 0; r < _R; r++) {
                    _fykq(j, h, r) = 0;
                }
            }
        }

        for (int j = 0; j < _n; j++) {
            for (int h = 0; h < _H; h++) {
                for (int r = 0; r < _R; r++) {
                    for (int l = 0; l < _L; l++) {
                        _fykq(j, h, r) += _pf(j, l, h, r);
                    }
                }
            }
        }

        //for(int gi=0;gi<_iter2;gi++){

        for (int k = 0; k < _K; k++) {
            for (int h = 0; h < _H; h++) {
                _ff(k, h) = 0;
            }
        }

        for (int k = 0; k < _K; k++) {
            for (int q = 0; q < _Q; q++) {
                for (int h = 0; h < _H; h++) {
                    for (int r = 0; r < _R; r++) {
                        _res(k, q, h, r) = 0;
                    }
                }
            }
        }

        //we wrote f(y_k|z,w_k=h) in logarithm form
        int m = 0;
        for (int k = 0; k < _K; k++) {
            for (int q = 0; q < _Q; q++) {
                for (int h = 0; h < _H; h++) {
                    for (int r = 0; r < _R; r++) {
                        for (int j = m; j < _nkq(k, q) + m; j++) {
                            _res(k, q, h, r) += log(_fykq(j, h, r));
                        }
                    }
                }
                m = m + _nkq(k, q);
            }
        }

        for (int k = 0; k < _K; k++) {
            for (int h = 0; h < _H; h++) {
                for (int q = 0; q < _Q; q++) {
                    for (int r = 0; r < _R; r++) {
                        _res(k, q, h, r) = _z(q, r) * _res(k, q, h, r);
                    }
                }
            }
        }

        for (int k = 0; k < _K; k++) {
            for (int h = 0; h < _H; h++) {
                for (int r = 0; r < _R; r++) {
                    for (int q = 0; q < _Q; q++) {
                        _ff(k, h) += _res(k, q, h, r);
                    }
                }
            }
        }

        for (int k = 0; k < _K; k++) {
            for (int h = 0; h < _H; h++) {
                _num2(k, h) = log(_p_h[h]) + _ff(k, h);
            }
        }

        vector<double> _maxh(_K);
        for (int k = 0; k < _K; k++) {
            _maxh[k] = _num2(k, 0);
            for (int h = 0; h < _H; h++) {
                if (_num2(k, h) > _maxh[k]) {
                    _maxh[k] = _num2(k, h);
                }
            }
        }

        for (int k = 0; k < _K; k++) {
            for (int h = 0; h < _H; h++) {
                _num2(k, h) = _num2(k, h) - _maxh[k];
            }
        }

        for (int k = 0; k < _K; k++) {
            for (int h = 0; h < _H; h++) {
                _num2(k, h) = exp(_num2(k, h));
            }
        }

        for (int k = 0; k < _K; k++) {
            _den2[k] = 0;
        }

        for (int k = 0; k < _K; k++) {
            for (int h = 0; h < _H; h++) {
                _den2[k] += _num2(k, h);
            }
        }

        
        for (int k = 0; k < _K; k++) {
            for (int h = 0; h < _H; h++) {
                _Prob1(k, h) = (double) _num2(k, h) / _den2[k];
            }
        }

        _pcv=0;
        while(_pcv==0) {
            _pcv = 1;
            for (int h = 0; h < _H; h++) {
                _cv[h]=0;
            }
             _w.zeros(_K, _H);
            //boost::mt19937 generator2(1234);
            for (int k = 0; k < _K; k++) {
                boost::random::discrete_distribution<int> distribution4(_Prob1.row(k).begin(), _Prob1.row(k).end());
                int sample = distribution4(generator);
                _w(k, sample) = 1;
            }

            for (int h = 0; h < _H; h++) {
                for (int k = 0; k < _K; k++) {
                    _cv[h] += _w(k, h);
                }
            }
            for (int h = 0; h < _H; h++) {
                _pcv = _pcv * _cv[h];
            }
        }

        std::fill(_numh.begin(), _numh.end(), 0);

        for (int h = 0; h < _H; h++) {
            for (int k = 0; k < _K; k++) {
                _numh[h] += _w(k, h);
            }
        }

        //update of p_h

        for (int h = 0; h < _H; h++) {
            _p_h[h] = (double) _numh[h] / _K;
        }
        
        vector<double> _p_hc(_H);
        for (int h = 0; h < _H; h++) {
          _p_hc[h]=_p_h[h];
        }
        
        double min1, temp1;
        vector<int> temp2(_K);
        mat temp3(_n,_R);
        arma::cube temp4(_n,_L,_R);
        mat temp5(_L,_R);
        for (int i = 0; i < _H - 1; i++) {
          min1 = i;
          for (int j = i + 1; j < _H; j++) {
            if (_p_h[j] < _p_h[min1]) {
              min1 = j;
            }
          }
          temp1 = _p_h[i];
          _p_h[i] = _p_h[min1];
          _p_h[min1] = temp1;
          
        }
        
        
        for (int q = 0; q < _Q; q++) {
            for (int r = 0; r < _R; r++) {
                _ff2(q, r) = 0;
            }
        }

        for (int q = 0; q < _Q; q++) {
            for (int k = 0; k < _K; k++) {
                for (int h = 0; h < _H; h++) {
                    for (int r = 0; r < _R; r++) {
                        _res2(q, k, h, r) = 0;
                    }
                }
            }
        }

        int mm = 0;
        for (int q = 0; q < _Q; q++) {
            for (int k = 0; k < _K; k++) {
                for (int h = 0; h < _H; h++) {
                    for (int r = 0; r < _R; r++) {
                        for (int j = mm; j < _nkq(k, q) + mm; j++) {
                            _res2(q, k, h, r) += log(_fykq(_indvecqk[j], h, r));
                        }
                    }
                }
                mm = mm + _nkq(k, q);
            }
        }

        for (int q = 0; q < _Q; q++) {
            for (int r = 0; r < _R; r++) {
                for (int k = 0; k < _K; k++) {
                    for (int h = 0; h < _H; h++) {
                        _res2(q, k, h, r) = _w(k, h) * _res2(q, k, h, r);
                    }
                }
            }
        }

        for (int q = 0; q < _Q; q++) {
            for (int r = 0; r < _R; r++) {
                for (int k = 0; k < _K; k++) {
                    for (int h = 0; h < _H; h++) {
                        _ff2(q, r) += _res2(q, k, h, r);
                    }
                }
            }
        }

        for (int q = 0; q < _Q; q++) {
            for (int r = 0; r < _R; r++) {
                _num3(q, r) = log(_p_r[r]) + _ff2(q, r);
            }
        }

        vector<double> _maxr(_Q);
        for (int q = 0; q < _Q; q++) {
            _maxr[q] = _num3(q, 0);
            for (int r = 0; r < _R; r++) {
                if (_num3(q, r) > _maxr[q]) {
                    _maxr[q] = _num3(q, r);
                }
            }
        }

        for (int q = 0; q < _Q; q++) {
            for (int r = 0; r < _R; r++) {
                _num3(q, r) = _num3(q, r) - _maxr[q];
            }
        }

        for (int q = 0; q < _Q; q++) {
            for (int r = 0; r < _R; r++) {
                _num3(q, r) = exp(_num3(q, r));
            }
        }


        for (int q = 0; q < _Q; q++) {
            _den3[q] = 0;
        }

        for (int q = 0; q < _Q; q++) {
            for (int r = 0; r < _R; r++) {
                _den3[q] += _num3(q, r);
            }
        }

        for (int q = 0; q < _Q; q++) {
            for (int r = 0; r < _R; r++) {
                _Prob2(q, r) = (double) _num3(q, r) / _den3[q];
            }
        }

        _pcv2=0;
        while(_pcv2==0) {
            _pcv2 = 1;
            for (int r = 0; r < _R; r++) {
                _cv2[r]=0;
            }
            _z.zeros(_Q, _R);
            for (int q = 0; q < _Q; q++) {
                boost::random::discrete_distribution<int> distribution5(_Prob2.row(q).begin(), _Prob2.row(q).end());
                int sample = distribution5(generator);
                _z(q, sample) = 1;
            }

            for (int r = 0; r < _R; r++) {
                for (int q = 0; q < _Q; q++) {
                    _cv2[r] += _z(q, r);
                }
            }
            for (int r = 0; r < _R; r++) {
                _pcv2 = _pcv2 * _cv2[r];
            }
        }

        std::fill(_numr.begin(), _numr.end(), 0);

        for (int r = 0; r < _R; r++) {
            for (int q = 0; q < _Q; q++) {
                _numr[r] += _z(q, r);
            }
        }

        //update of p_r

        for (int r = 0; r < _R; r++) {
            _p_r[r] = (double) _numr[r] / _Q;
        }
        
        vector<double> _p_rc(_R);
        for (int r = 0; r < _R; r++) {
          _p_rc[r]=_p_r[r];
        }
        
        double min2, temp6;
        vector<int> temp7(_Q);
        mat temp8(_n,_H);
        arma::cube temp9(_n,_L,_H);
        mat temp10(_L,_H);
        for (int i = 0; i < _R - 1; i++) {
          min2 = i;
          for (int j = i + 1; j < _R; j++) {
            if (_p_r[j] < _p_r[min2]) {
              min2 = j;
            }
          }
          temp6 = _p_r[i];
          _p_r[i] = _p_r[min2];
          _p_r[min2] = temp6;
          
          
        }
        
        for (int j = 0; j < _n; j++) {
          for (int h = 0; h < _H; h++) {
            _wjk(j,h)=0;
          }
        }
        for (int j = 0; j < _n; j++) {
          for (int r = 0; r < _R; r++) {
            _zjq(j,r)=0;
          }
        }
        for(int j=0;j<_n;j++) {
          for (int h = 0; h < _H; h++) {
            for (int r = 0; r < _R; r++) {
              _B(j,h,r)=0;
            }
          }
        }
        
        for (int j = 0; j < _n; j++) {
            for (int h = 0; h < _H; h++) {
                _wjk(j, h) = _w(y(j, 1) - 1, h);
            }
        }

        for (int j = 0; j < _n; j++) {
            for (int r = 0; r < _R; r++) {
                _zjq(j, r) = _z(y(j, 2) - 1, r);
            }
        }

        for (int j = 0; j < _n; j++) {
            for (int h = 0; h < _H; h++) {
                for (int r = 0; r < _R; r++) {
                    _B(j, h, r) = _wjk(j, h) * _zjq(j, r);
                }
            }
        }

        for (int j = 0; j < _n; j++) {
            for (int l = 0; l < _L; l++) {
                for (int h = 0; h < _H; h++) {
                    for (int r = 0; r < _R; r++) {
                        if (_B(j, h, r) == 1) {
                            _Prob3(j, l) = log(_pi_l(l, h, r)) + log(_f(j, l,h,r));
                        }
                    }
                }
            }
        }

        vector<double> _maxl(_n);
        for (int j = 0; j < _n; j++) {
            _maxl[j] = _Prob3(j, 0);
            for (int l = 0; l < _L; l++) {
                if (_Prob3(j, l) > _maxl[j]) {
                    _maxl[j] = _Prob3(j, l);
                }
            }
        }

        for (int j = 0; j < _n; j++) {
            for (int l = 0; l < _L; l++) {
                _Prob3(j, l) = _Prob3(j, l) - _maxl[j];
            }
        }

        for (int j = 0; j < _n; j++) {
            for (int l = 0; l < _L; l++) {
                _Prob3(j, l) = exp(_Prob3(j, l));
            }
        }

        for (int j = 0; j < _n; j++) {
            _denom[j] = 0;
        }

        for (int j = 0; j < _n; j++) {
            for (int l = 0; l < _L; l++) {
                _denom[j] += _Prob3(j, l);
            }
        }

        for (int j = 0; j < _n; j++) {
            for (int l = 0; l < _L; l++) {
                _Prob3(j, l) = (double) _Prob3(j, l) / _denom[j];
            }
        }

                  for (int l = 0; l < _L; l++) {
            for (int h = 0; h < _H; h++) {
                for (int r = 0; r < _R; r++) {
                    _num4(l, h, r) = 0;
                }
            }
        }

        for (int l = 0; l < _L; l++) {
            for (int h = 0; h < _H; h++) {
                for (int r = 0; r < _R; r++) {
                    for (int j = 0; j < _n; j++) {
                        if(_B(j,h,r)==1) {
                            _num4(l, h, r) += _Prob3(j, l);
                        }
                    }
                }
            }
        }

        for (int h = 0; h < _H; h++) {
            for (int r = 0; r < _R; r++) {
                _den4(h, r) = 0;
            }
        }
        for (int h = 0; h < _H; h++) {
            for (int r = 0; r < _R; r++) {
                for (int l = 0; l < _L; l++) {
                    _den4(h, r) += _num4(l, h, r);
                }
            }
        }

        for (int l = 0; l < _L; l++) {
            for (int h = 0; h < _H; h++) {
                for (int r = 0; r < _R; r++) {
                    if (_den4(h, r) != 0) {
                        _pi_l(l, h, r) = (double) _num4(l, h, r) / _den4(h, r);
                    }
                }
            }
        }

        for (int h = 0; h < _H; h++) {
          for (int r = 0; r < _R; r++) {
            for (int l = 0; l < _L; l++) {
              if (_pi_l(l, h, r) == 0) {
                _pi_l(l,h,r)=0.001;
              }
              if (_pi_l(l, h, r) == 1) {
                _pi_l(l,h,r)=0.999;
              }
            }
          }
        }

        //update of p

        vector<int> _dephr(_I);
        for (int i = 0; i < _I; i++) {
            _dephr[i]=0;
        }

        for (int i = 0; i < _I; i++) {
            for (int c = 0; c < _C; c++) {
                for (int l = 0; l < _L; l++) {
                    for (int h = 0; h < _H; h++) {
                        for (int r = 0; r < _R; r++) {
                            _num5(i, c, l, h, r) = 0;
                        }
                    }
                }
            }
        }
        for (int i = 0; i < _I; i++) {
            for (int c = 0; c < _C; c++) {
                for (int l = 0; l < _L; l++) {
                    for (int h = 0; h < _H; h++) {
                        for (int r = 0; r < _R; r++) {
                            for (int j = 0; j < _n; j++) {
                                if(_dephr[i]==3) {
                                    if (y(j, i + 3) == c && _B(j, h, r) == 1) {
                                        _num5(i, c, l, h, r) += _Prob3(j,l);
                                    }
                                }
                                else if(_dephr[i]==2){
                                    if (y(j, i + 3) == c  && _zjq(j,r)==1) {
                                        _num5(i, c, l, h, r) += _Prob3(j,l);
                                    }
                                }
                                else if(_dephr[i]==1){
                                    if (y(j, i + 3) == c &&  _wjk(j,h)==1) {
                                        _num5(i, c, l, h, r) += _Prob3(j,l);
                                    }
                                }
                                else if(_dephr[i]==0){
                                    if (y(j, i + 3) == c) {
                                        _num5(i, c, l, h, r) += _Prob3(j,l);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }


        for (int l = 0; l < _L; l++) {
            _den5[l] = 0;
        }
        for (int l = 0; l < _L; l++) {
            for (int j = 0; j < _n; j++) {
                _den5[l] += _Prob3(j, l);
            }
        }

        for (int l = 0; l < _L; l++) {
            for (int r = 0; r < _R; r++) {
                _den6(l, r) = 0;
            }
        }

        for (int l = 0; l < _L; l++) {
            for (int r = 0; r < _R; r++) {
                for (int j = 0; j < _n; j++) {
                    if(_zjq(j,r)==1) {
                        _den6(l, r) +=_Prob3(j, l);
                    }
                }
            }
        }

        for (int l = 0; l < _L; l++) {
            for (int h = 0; h < _H; h++) {
                _den7(l, h) = 0;
            }
        }

        for (int l = 0; l < _L; l++) {
            for (int h = 0; h < _H; h++) {
                for (int j = 0; j < _n; j++) {
                    if(_wjk(j,h)==1) {
                        _den7(l, h) += _Prob3(j, l);
                    }
                }
            }
        }

      for (int i = 0; i < _I; i++) {
            for (int c = 0; c < _C; c++) {
                for (int l = 0; l < _L; l++) {
                    for (int h = 0; h < _H; h++) {
                        for (int r = 0; r < _R; r++) {
                            if(_dephr[i]==3) {
                                if (_num4(l, h, r) != 0) {
                                    _p(i, c, l, h, r) = (double) _num5(i, c, l, h, r) / _num4(l, h, r);
                                }
                            }
                            else if(_dephr[i]==2) {
                                if (_den6(l,r) != 0) {
                                    _p(i, c, l, h, r) = (double) _num5(i, c, l, h, r) / _den6(l,r);
                                }
                            }
                            else if(_dephr[i]==1) {
                                if (_den7(l,h) != 0) {
                                    _p(i, c, l, h, r) = (double) _num5(i, c, l, h, r) / _den7(l,h);
                                }
                            }
                            else if(_dephr[i]==0) {
                                if (_den5[l] != 0) {
                                    _p(i, c, l, h, r) = (double) _num5(i, c, l, h, r) / _den5[l];
                                }
                            }
                        }
                    }
                }
            }
        }

        for (int i = 0; i < _I; i++) {
          for (int h = 0; h < _H; h++) {
            for (int r = 0; r < _R; r++) {
              for (int l = 0; l < _L; l++) {
                for (int c = 0; c < _C; c++) {
                  if (_p(i, c, l, h, r) > 0.9999) {
                    _p(i, c, l, h, r) = 0.999;
                  }
                }
              }
            }
          }
        }
        for (int i = 0; i < _I; i++) {
          for (int h = 0; h < _H; h++) {
            for (int r = 0; r < _R; r++) {
              for (int l = 0; l < _L; l++) {
                for (int c = 0; c < _C; c++) {
                  if (_p(i, c, l, h, r) < 0.0001) {
                    _p(i, c, l, h, r) = 0.001;
                  }
                }
              }
            }
          }
        }
        
        //mean of parameters

        if(u1>=_burn && u1==it[mu]){
            
            mu=mu+1;

            
            for(int h=0;h<_H;h++){
                _p_h_hat[h] +=_p_h[h];
            }

            
            for(int r=0;r<_R;r++){
                _p_r_hat[r] +=_p_r[r];
            }
            
            for (int h = 0; h < _H; h++) {
              _p_hc_hat[h] += _p_hc[h];
            }
            
            for (int r = 0; r < _R; r++) {
              _p_rc_hat[r] += _p_rc[r];
            }

            
            for (int l = 0; l < _L; l++) {
                for (int h = 0; h < _H; h++) {
                    for (int r = 0; r < _R; r++) {
                        _pi_l_hat(l,h,r) += _pi_l(l,h,r);
                    }
                }
            }

            
            for (int i = 0; i < _I; i++) {
                for (int c = 0; c < _C; c++) {
                    for(int l=0;l<_L;l++){
                        for (int h = 0; h < _H; h++) {
                            for (int r = 0; r < _R; r++) {
                                _p_hat(i, c, l, h, r) += _p(i, c, l, h, r);
                            }
                        }
                    }
                }
            }

        }

        if(u1==_iter-1){

            //_p_h_hat
            for(int h=0;h<_H;h++){
                _p_h_hat[h]=(double)_p_h_hat[h]/(nn);
            }

            cout<<"estimated ph"<<endl;
          
            for(int h=0;h<_H;h++){
              
              cout<< _p_h_hat[h]<<",";
            }
            cout<<"\n";


            for (int h = 0; h < _H; h++) {
              _p_hc_hat[h] = (double) _p_hc_hat[h] / (nn);
            }
           
            //_p_r_hat
            for(int r=0;r<_R;r++){
                _p_r_hat[r]=(double)_p_r_hat[r]/(nn);
            }

            cout<<"estimated pr"<<endl;
            
            for(int r=0;r<_R;r++){
              
              cout<< _p_r_hat[r]<<",";
            }
            cout<<"\n";

            
            for (int r = 0; r < _R; r++) {
              _p_rc_hat[r] = (double) _p_rc_hat[r] / (nn);
            }
            
          
            //_pi_l_hat
            for (int l = 0; l < _L; l++) {
                for (int h = 0; h < _H; h++) {
                    for (int r = 0; r < _R; r++) {
                        _pi_l_hat(l,h,r) = (double)_pi_l_hat(l,h,r)/(nn);
                    }
                }
            }
            
            
            for (int i = 0; i < (_H - 1); i++) {
              min1 = i;
              for (int j = (i + 1); j < _H; j++) {
                if (_p_hc_hat[j] < _p_hc_hat[min1]) {
                  min1 = j;
                }
              }
              
              temp1 = _p_hc_hat[i];
              _p_hc_hat[i] = _p_hc_hat[min1];
              _p_hc_hat[min1] = temp1;
              
              for (int l = 0; l < _L; l++) {
                for (int r = 0; r < _R; r++) {
                  temp5(l, r) =(double) _pi_l_hat(l, i, r);
                  _pi_l_hat(l, i, r) = _pi_l_hat(l, min1, r);
                  _pi_l_hat(l, min1, r) = temp5(l, r);
                }
              }
              
            }
            
            for (int i = 0; i < (_R - 1); i++) {
              min2 = i;
              for (int j = (i + 1); j < _R; j++) {
                if (_p_rc_hat[j] < _p_rc_hat[min2]) {
                  min2 = j;
                }
              }
              
              temp6 = _p_rc_hat[i];
              _p_rc_hat[i] = _p_rc_hat[min2];
              _p_rc_hat[min2] = temp6;
              
              for (int l = 0; l < _L; l++) {
                for (int h = 0; h < _H; h++) {
                  temp10(l, h) =(double) _pi_l_hat(l,h,i);
                  _pi_l_hat(l,h,i) = _pi_l_hat(l,h,min2);
                  _pi_l_hat(l,h,min2) = temp10(l, h);
                }
              }
              
            }

            cout << "_pi_l_hat(l,0,0)";
            cout << "\n";
            for (int l = 0; l < _L; l++) {
                cout << _pi_l_hat(l, 0, 0) << ",";
            }
            cout << "\n";
            cout << "_pi_l_hat(l,0,1)";
            cout << "\n";
            for (int l = 0; l < _L; l++) {
                 cout << _pi_l_hat(l, 0, 1) << ",";
            }
            cout << "\n";
            cout << "_pi_l_hat(l,1,0)";
            cout << "\n";
            for (int l = 0; l < _L; l++) {
                cout << _pi_l_hat(l, 1, 0) << ",";
            }
            cout << "\n";
            cout << "_pi_l_hat(l,1,1)";
            cout << "\n";
            for (int l = 0; l < _L; l++) {
                cout << _pi_l_hat(l, 1, 1) << ",";
            }
            cout << "\n";

           
            //_p_hat
            for (int i = 0; i < _I; i++) {
                for (int c = 0; c < _C; c++) {
                    for(int l=0;l<_L;l++){
                        for (int h = 0; h < _H; h++) {
                            for (int r = 0; r < _R; r++) {
                                _p_hat(i, c, l, h, r) = (double) _p_hat(i, c, l, h, r) / (nn);
                            }
                        }
                    }
                }
            }

            cout << "\n";
            cout << "_p_hat(0,c,l,0,0)";
            cout << "\n";   
            for (int c = 0; c < _C; c++) {
                for (int l = 0; l < _L; l++) {
                    cout << _p_hat(0,c,l, 0, 0) << ",";
               }
               cout << "\n";
            }
            cout << "\n";

            cout << "\n";
            cout << "_p_hat(1,c,l,0,0)";
            cout << "\n";   
            for (int c = 0; c < _C; c++) {
                for (int l = 0; l < _L; l++) {
                    cout << _p_hat(1,c,l, 0, 0) << ",";
               }
               cout << "\n";
            }
            cout << "\n";

           
        }

          
          for (int h = 0; h < _H; h++) {
            _p_h[h]=_p_hc[h];
          }
          for (int r = 0; r < _R; r++) {
            _p_r[r]=_p_rc[r];
          }
        
        
    }
    
    //-------------------------------------------------------------------
    //classification

    
    for (int j = 0; j < _n; j++) {
        for (int l = 0; l < _L; l++) {
            for (int h = 0; h < _H; h++) {
                for (int r = 0; r < _R; r++) {
                    _f(j, l, h, r) = 1;
                }
            }
        }
    }

    for (int j = 0; j < _n; j++) {
        for (int l = 0; l < _L; l++) {
            for (int h = 0; h < _H; h++) {
                for (int r = 0; r < _R; r++) {
                    for (int i = 0; i < _I; i++) {
                        _f(j, l, h, r) = _f(j, l, h, r) * _p_hat(i, y(j, i + 3), l, h, r);
                    }
                }
            }
        }
    }

    for (int j = 0; j < _n; j++) {
        for (int h = 0; h < _H; h++) {
            for (int r = 0; r < _R; r++) {
                for (int l = 0; l < _L; l++) {
                    _pf(j,l,h,r)= (double)_pi_l_hat(l,h,r) * _f(j,l,h,r);
                }
            }
        }
    }


    for (int j = 0; j < _n; j++) {
        for (int h = 0; h < _H; h++) {
            for (int r = 0; r < _R; r++) {
                _fykq(j,h,r)=0;
            }
        }
    }

    for (int j = 0; j < _n; j++) {
        for (int h = 0; h < _H; h++) {
            for (int r = 0; r < _R; r++) {
                for (int l = 0; l < _L; l++) {
                    _fykq(j,h,r) += _pf(j,l,h,r);
                }
            }
        }
    }

    _w_hat.zeros(_K,_H);
    _z_hat.zeros(_Q,_R);
    _x_hat.zeros(_n,_L);

    for(int gi=0;gi<_iter3;gi++){

        for (int k = 0; k < _K; k++) {
            for (int h = 0; h < _H; h++) {
                _ff(k,h)=0;
            }
        }

        for (int k = 0; k < _K; k++) {
            for (int q = 0; q < _Q; q++) {
                for (int h = 0; h < _H; h++) {
                    for (int r = 0; r < _R; r++) {
                        _res(k,q,h,r)=0;
                    }
                }
            }
        }


        int m=0;
        for (int k = 0; k < _K; k++) {
            for (int q = 0; q < _Q; q++) {
                for (int h = 0; h < _H; h++) {
                    for (int r = 0; r < _R; r++) {
                        for (int j = m; j < _nkq(k,q)+m; j++) {
                            _res(k,q,h,r) += log(_fykq(j,h,r));
                        }
                    }
                }
                m=m+_nkq(k,q);
            }
        }


        for (int k = 0; k < _K; k++) {
            for (int h = 0; h < _H; h++) {
                for (int q = 0; q < _Q; q++) {
                    for (int r = 0; r < _R; r++) {
                        _res(k,q,h,r) =_z(q,r) * _res(k,q,h,r);
                    }
                }
            }
        }


        for (int k = 0; k < _K; k++) {
            for (int h = 0; h < _H; h++) {
                for (int r = 0; r < _R; r++) {
                    for (int q = 0; q < _Q; q++) {
                        _ff(k,h) += _res(k,q,h,r);
                    }
                }
            }
        }

        for(int k=0;k<_K;k++){
            for(int h=0;h<_H;h++){
                _num2(k,h)=log(_p_h[h]) + _ff(k,h);
            }
        }

        vector<double>_maxh(_K);
        for(int k=0;k<_K;k++){
            _maxh[k]=_num2(k,0);
            for(int h=0;h<_H;h++){
                if(_num2(k,h)>_maxh[k]){
                    _maxh[k]=_num2(k,h);
                }
            }
        }

        for(int k=0;k<_K;k++){
            for(int h=0;h<_H;h++){
                _num2(k,h)=_num2(k,h)-_maxh[k];
            }
        }

        for(int k=0;k<_K;k++){
            for(int h=0;h<_H;h++){
                _num2(k,h)=exp(_num2(k,h));
            }
        }


        for(int k=0;k<_K;k++){
            _den2[k]=0;
        }

        for(int k=0;k<_K;k++){
            for(int h=0;h<_H;h++){
                _den2[k] +=_num2(k,h);
            }
        }


        for(int k=0;k<_K;k++){
            for(int h=0;h<_H;h++){
                _Prob1(k,h)=(double)_num2(k,h) / _den2[k];
            }
        }


        _w.zeros(_K,_H);
//boost::mt19937 generator2(1234);
        for(int k=0;k<_K;k++){
            boost::random::discrete_distribution<int> distribution7 (_Prob1.row(k).begin(),_Prob1.row(k).end());
            int sample = distribution7(generator);
            _w(k,sample) = 1;
        }

        for(int k = 0; k < _K; k++){
            for(int h = 0; h < _H; h++){
                _w_hat(k,h)+=_w(k,h);
            }
        }


        for (int q = 0; q < _Q; q++) {
            for (int r = 0; r < _R; r++) {
                _ff2(q,r)=0;
            }
        }

        for (int q = 0; q < _Q; q++) {
            for (int k = 0; k < _K; k++) {
                for (int h = 0; h < _H; h++) {
                    for (int r = 0; r < _R; r++) {
                        _res2(q,k,h,r)=0;
                    }
                }
            }
        }

        int mm=0;
        for (int q = 0; q < _Q; q++) {
            for (int k = 0; k < _K; k++) {
                for (int h = 0; h < _H; h++) {
                    for (int r = 0; r < _R; r++) {
                        for (int j = mm; j < _nkq(k,q)+mm; j++) {
                            _res2(q,k,h,r) += log(_fykq(_indvecqk[j],h,r));
                        }
                    }
                }
                mm=mm+_nkq(k,q);
            }
        }

        for (int q = 0; q < _Q; q++) {
            for (int r = 0; r < _R; r++) {
                for (int k = 0; k < _K; k++) {
                    for (int h = 0; h < _H; h++) {
                        _res2(q,k,h,r) =_w(k,h) * _res2(q,k,h,r);
                    }
                }
            }
        }


        for (int q = 0; q < _Q; q++) {
            for (int r = 0; r < _R; r++) {
                for (int k = 0; k < _K; k++) {
                    for (int h = 0; h < _H; h++) {
                        _ff2(q,r) += _res2(q,k,h,r);
                    }
                }
            }
        }

        for(int q=0;q<_Q;q++){
            for(int r=0;r<_R;r++){
                _num3(q,r)=log(_p_r[r]) + _ff2(q,r);
            }
        }

        vector<double>_maxr(_Q);
        for(int q=0;q<_Q;q++){
            _maxr[q]=_num3(q,0);
            for(int r=0;r<_R;r++){
                if(_num3(q,r)>_maxr[q]){
                    _maxr[q]=_num3(q,r);
                }
            }
        }

        for(int q=0;q<_Q;q++){
            for(int r=0;r<_R;r++){
                _num3(q,r)=_num3(q,r)-_maxr[q];
            }
        }

        for(int q=0;q<_Q;q++){
            for(int r=0;r<_R;r++){
                _num3(q,r)=exp(_num3(q,r));
            }
        }

        for(int q=0;q<_Q;q++){
            _den3[q]=0;
        }

        for(int q=0;q<_Q;q++){
            for(int r=0;r<_R;r++){
                _den3[q] +=_num3(q,r);
            }
        }

        for(int q=0;q<_Q;q++){
            for(int r=0;r<_R;r++){
                _Prob2(q,r)=(double)_num3(q,r) / _den3[q];
            }
        }

        _z.zeros(_Q,_R);
        for(int q=0;q<_Q;q++){
            boost::random::discrete_distribution<int> distribution8 (_Prob2.row(q).begin(),_Prob2.row(q).end());
            int sample = distribution8(generator);
            _z(q,sample) = 1;
        }

        for(int q = 0; q < _Q; q++){
            for(int r = 0; r < _R; r++){
                _z_hat(q,r)+=_z(q,r);
            }
        }

    }

    std::fill(_indh.begin(),_indh.end(),0);
    for(int k=0;k<_K;k++){
        int _max=_w_hat(k,0);
        for(int h=0;h<_H;h++){
            if(_w_hat(k,h)>_max){
                _max=_w_hat(k,h);
                _indh[k]=h;
            }
        }
    }

    
    
    
    std::fill(_indr.begin(),_indr.end(),0);
    for(int q=0;q<_Q;q++){
        int _max2=_z_hat(q,0);
        for(int r=0;r<_R;r++){
            if(_z_hat(q,r)>_max2){
                _max2=_z_hat(q,r);
                _indr[q]=r;
            }
        }
    }

   

    _wjk.zeros(_n,_H);

    for(int j=0;j<_n;j++){
        for(int h=0;h<_H;h++){
            if(_indh[y(j,1)-1]==h){
                _wjk(j,h)=1;
            }
        }
    }

    _zjq.zeros(_n,_R);

    for(int j=0;j<_n;j++){
        for(int r=0;r<_R;r++){
            if(_indr[y(j,2)-1]==r){
                _zjq(j,r)=1;
            }
        }
    }
    
    for(int j=0;j<_n;j++) {
      for (int h = 0; h < _H; h++) {
        for (int r = 0; r < _R; r++) {
          _B(j,h,r)=0;
        }
      }
    }


    for(int j=0;j<_n;j++){
        for(int h=0;h<_H;h++){
            for(int r=0;r<_R;r++){
                _B(j,h,r)=_wjk(j,h)*_zjq(j,r);
            }
        }
    }

    for (int j = 0; j < _n; j++) {
        for (int l = 0; l < _L; l++) {
            for (int h = 0; h < _H; h++) {
                for (int r = 0; r < _R; r++) {
                    if (_B(j, h, r) == 1) {
                        _Prob3(j, l) = log(_pi_l_hat(l, h, r)) + log(_f(j, l,h,r));
                    }
                }
            }
        }
    }

    vector<double> _maxl(_n);
    for (int j = 0; j < _n; j++) {
        _maxl[j] = _Prob3(j, 0);
        for (int l = 0; l < _L; l++) {
            if (_Prob3(j, l) > _maxl[j]) {
                _maxl[j] = _Prob3(j, l);
            }
        }
    }

    for (int j = 0; j < _n; j++) {
        for (int l = 0; l < _L; l++) {
            _Prob3(j, l) = _Prob3(j, l) - _maxl[j];
        }
    }

    for (int j = 0; j < _n; j++) {
        for (int l = 0; l < _L; l++) {
            _Prob3(j, l) = exp(_Prob3(j, l));
        }
    }

    for (int j = 0; j < _n; j++) {
        _denom[j] = 0;
    }

    for (int j = 0; j < _n; j++) {
        for (int l = 0; l < _L; l++) {
            _denom[j] += _Prob3(j, l);
        }
    }

    for (int j = 0; j < _n; j++) {
        for (int l = 0; l < _L; l++) {
            _Prob3(j, l) = (double) _Prob3(j, l) / _denom[j];
        }
    }

    for(int vv=0;vv<_iter3;vv++){

        for (int j = 0; j < _n; j++) {
            for (int l = 0; l < _L; l++){
                _xjkq(j,l)=0;
            }
        }


        for(int j=0;j<_n;j++){
            boost::random::discrete_distribution<int> distribution9 (_Prob3.row(j).begin(),_Prob3.row(j).end());
            int sample = distribution9(generator);
            _xjkq(j,sample) = 1;
        }

        for(int j = 0; j < _n; j++){
            for(int l = 0; l < _L; l++){
                _x_hat(j,l)+=_xjkq(j,l);
            }
        }
    }


    int _b=1;
    for(int n=0;n<_n;n++){
        if(y(n,0)==_b){
            _classif2(_b-1,0)=_b;
            for(int k=0;k<_K;k++){
                if(y(n,1)==k+1){
                    _classif2(_b-1,1)=_indh[k]+1;
                }
            }
            for(int q=0;q<_Q;q++){
                if(y(n,2)==q+1){
                    _classif2(_b-1,2)=_indr[q]+1;
                }
            }
            _b=_b+1;
            n=-1;
        }
    }

    //classification at the first level

    mat _x_hat_ord(_n,_L);
    int _bb=1;
    for(int j=0;j<_n;j++){
        if(y(j,0)==_bb){
            for(int l=0;l<_L;l++){
                _x_hat_ord(_bb-1,l)=_x_hat(j,l);
            }
            _bb=_bb+1;
            j=-1;
        }
    }

    std::fill(_indl.begin(),_indl.end(),0);
    for(int j=0;j<_n;j++){
        double _max3 = _x_hat_ord(j,0);
        for(int l=0;l<_L;l++){
            if(_x_hat_ord(j,l) > _max3){
                _max3 = _x_hat_ord(j,l);
                _indl[j]=l;
            }
        }
    }



    //final classification
    for(int j=0;j<_n;j++){
        _classif(j,0)=j+1;
        _classif(j,1)=_indl[j]+1;
        _classif(j,2)=_classif2(j,1);
        _classif(j,3)=_classif2(j,2);
    }

   

    auto end = chrono::steady_clock::now();
    int elapsed_time = chrono::duration_cast<chrono::seconds>(end - start).count();
    std::cout << "Time: " << elapsed_time << " sec" << std::endl;

}
