// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]

#include <RcppArmadillo.h>
#include <iostream>
#include <Rcpp.h>
#include <RcppNumerical.h>

using namespace Rcpp;
using namespace std;
using namespace arma;
using namespace stats;
using namespace Numer;
#define Rcpp__stats__random_rnorm_h



// [[Rcpp::export]]
arma::mat RCPPdiffX(arma::mat X, int choice){
/* X is a mat whose dimension is J*nAttr 
 choice is a mat whose dimension is a index number */
	double Alts = X.n_rows;
	arma::mat L = ones<mat>(Alts-1,Alts);
	arma::mat L1 = -eye<mat>(Alts-1,Alts-1);
	int k = 1;
	for (int i = 1; i <= Alts; ++i) {
		if (choice!=i) {
			L.col(i-1)=L1.col(k-1);
      			k=k+1;
    		}



  }
	arma::mat diffX = L*X;
	return diffX;

}

// [[Rcpp::export]]
arma::mat RCPPweight(arma::mat diffX, arma::vec beta, double lam1, double lam2){
	arma::mat temp = diagmat(beta);
	arma::mat temp2 = diffX * temp;
	
	arma::mat W = exp(-lam1*temp2%(temp2>=0)+lam2*temp2%(temp2<0));
	//arma::mat W = exp(-lam1*temp2%(diffX>=0)+lam2*temp2%(diffX<0));
	return W;
}


// [[Rcpp::export]]
double RCPPdriftmean(arma::mat X, arma::vec beta, arma::vec zeta, double lam1, double lam2, double I0, int choice ){
	arma::mat diffX = RCPPdiffX(X, choice);


	arma::mat W = RCPPweight(diffX, beta, lam1, lam2);
	//double d = accu(W%diffX*diagmat(beta))/accu(W)+I0+zeta(choice-1);
	double d = accu(W%diffX*diagmat(beta))+I0+zeta(choice-1); 
	//double d = accu(W%diffX*diagmat(beta))+zeta(choice-1);
	d = d*(d>0);
	//d= exp(d);
	return d;	
}


// [[Rcpp::export]]
double RCPPcdf_rt(arma::mat X, arma::vec beta, arma::vec zeta, double lam1, double lam2, double I0,double b, double s, double A, int choice, double rt){
	double d = RCPPdriftmean(X, beta, zeta, lam1, lam2, I0, choice);
	/* make sure the d is bigger than 0 */
  
  double cdf = (1 + (b-A-d*rt)/A* normcdf((b-A-d*rt)/(s*rt)) - (b-d*rt)/A* normcdf((b-d*rt)/(s*rt)) + s*rt/A* normpdf((b-A-d*rt)/(s*rt))-s*rt/A* normpdf((b-d*rt)/(s*rt)))/normcdf(d/s);

	if (A<1e-10) cdf = normcdf(b/rt,d,s)/normcdf(d/s);
	if(cdf<1e-10|cdf>1.1) {cdf = 0;}
	return cdf;
	
}





// [[Rcpp::export]]
double RCPPpdf_rt(arma::mat X, arma::vec beta, arma::vec zeta, double lam1, double lam2, double I0,double b, double s, double A, double choice, double rt){
	double d = RCPPdriftmean(X, beta, zeta, lam1, lam2, I0, choice);
  	/* make sure the d is bigger than 0 */
  double pdf_imp = (-d*normcdf((b-A-d*rt)/(s*rt))+d*normcdf((b-d*rt)/(s*rt))+s*normpdf((b-A-d*rt)/(s*rt))-s*normpdf((b-d*rt)/(s*rt)))/A/normcdf(d/s);
	
	//pdf_imp = pdf_imp;
	if (A<1e-10) pdf_imp = normpdf(b/rt,d,s)*b/rt/rt/normcdf(d/s);
	if(pdf_imp<1e-10) {pdf_imp = 0;}
	return pdf_imp;
	
}


// [[Rcpp::export]]
double RCPPpdf_MLBA(arma::mat X, arma::vec beta, arma::vec zeta, double lam1, double lam2, double I0,double b, double s, double A,double choice, double rt){
  double temp = 1;
	double Alts = X.n_rows;

	for (int i = 1; i <= Alts; ++i) {
		if(i!=choice){
			temp = temp*(1- RCPPcdf_rt(X, beta, zeta, lam1, lam2, I0, b, s, A, i, rt));	
		}
	}
	
	double pdf_MLBA_imp = temp*RCPPpdf_rt(X, beta, zeta, lam1, lam2, I0, b, s, A, choice, rt);
	if (pdf_MLBA_imp<1e-10) {pdf_MLBA_imp = 0;}
	return pdf_MLBA_imp;
}


// [[Rcpp::export]]
double RCPPpdf_MLBA_rtknown(arma::mat X, arma::vec beta, arma::vec zeta, double lam1, double lam2, double I0, double b, double s, double A, double choice, double rt){

	//double Alts = X.n_rows;
	//double temp = 0;
  	//double test = 1;
	//for (int i = 1; i <= Alts; ++i) {
	//	temp = temp + RCPPpdf_MLBA(X, beta, zeta, lam1, lam2, I0, b, s, A, i, rt);
	// if(RCPPpdf_MLBA(X, beta, zeta, lam1, lam2, I0, b, s, A, i, rt)<0){
	//   test = 0;
	//  }
	//} 
	//temp = temp*test;
	double cdf_MLBA;
 	cdf_MLBA = RCPPpdf_MLBA(X, beta, zeta, lam1, lam2, I0, b, s, A, choice, rt);
 	//if (cdf_MLBA<1e-10) cdf_MLBA = 1e-10;
 	// Evan's assumption
	//if(temp==0) {	cdf_MLBA = 0;}

	// Scott's  assumption
	//if(temp==0) {cdf_MLBA = 1/Alts;}
	//else {
	//cdf_MLBA = RCPPpdf_MLBA(X, beta, zeta, lam1, lam2, I0, b, s, A, choice, rt)/temp;
	//} 

	return cdf_MLBA;
}

// [[Rcpp::export]]
arma::vec MLBA_rtknown_all(arma::mat X, arma::vec  beta, arma::vec zeta, double lam1, double lam2, double I0, double b, double s, double A, arma::vec choice, arma::vec rt){
	double nAlt = zeta.n_elem;
	double nObs = X.n_rows/nAlt;
	arma::vec MLBA = -ones(nObs);
	for (int i = 1; i <= nObs; ++i) {
		
		MLBA(i-1) = RCPPpdf_MLBA_rtknown(X.rows(nAlt*(i-1), (nAlt*i-1)), beta, zeta, lam1, lam2, I0, b, s, A, choice(i-1),rt(i-1));
	}  

	return MLBA;
}

// [[Rcpp::export]]

double RCPPMLBA_Lik_rtknown(arma::mat X, arma::vec  beta, arma::vec zeta, double lam1, double lam2, double I0, double b, double s, double A, arma::vec choice, arma::vec rt){
	arma::vec MLBA = MLBA_rtknown_all( X, beta, zeta, lam1, lam2, I0, b, s, A, choice, rt);
	double log_Lik =  accu(log(MLBA));
	return log_Lik;
}



// [[Rcpp::export]]
double RCPP_brier_CRT(arma::mat X, arma::vec  beta, arma::vec zeta, double lam1, double lam2, double I0, double b, double s, double A, arma::vec choice, arma::vec rt){
	double nAlt = zeta.n_elem;
	double nObs = X.n_rows/nAlt;
	double temp = 0;
	arma::vec choice_o = -ones(nObs);
	
	for(int i =1;i<=nAlt; ++i){
		arma::vec choice_s = ones(nObs)*i;
		
	  
	  for (int j=0;j<nObs;++j){
	    int temp2= choice(j)-choice_s(j);
	    if(temp2==0){
	      choice_o(j) = 1;
	    }
	    else{
	      choice_o(j) = 0;
	    }
	    
	  }
 
 	arma::vec MLBA_s = MLBA_rtknown_all( X, beta, zeta, lam1, lam2, I0, b, s, A, choice_s, rt);
	
	//temp = temp+accu(MLBA_s-choice_o);
	temp = temp+accu(pow(MLBA_s-choice_o,2));
}
	temp = temp/nObs;

	return temp;
}

/////////////////////////////////////////////////////////////

class MLBA_CDF: public Func
{
private:
	arma::mat X;
	arma::vec beta;
	arma::vec zeta;
	double lam1;
	double lam2;
 	double I0;
	double b;
	double s;
 	double A;
	double choice;

public:
    MLBA_CDF(arma::mat X_, arma::vec beta_, arma::vec zeta_, double lam1_, double lam2_, double I0_,double b_, double s_, double A_,double choice_) : X(X_), beta(beta_), zeta(zeta_), lam1(lam1_), lam2(lam2_), I0(I0_), b(b_),s(s_),A(A_),choice(choice_){}

    double operator()(const double& x) const
    {
        return RCPPpdf_MLBA(X, beta, zeta, lam1, lam2, I0, b ,s, A,choice,x);
    }
};


// [[Rcpp::export]]

double  RCPPcdf_MLBA_CO(arma::mat X, arma::vec beta, arma::vec zeta, double lam1, double lam2, double I0, double b, double s, double A, double choice){
	const double lower = 0, upper = 10;
	MLBA_CDF f(X, beta,zeta,lam1, lam2,I0, b,s, A,  choice) ;
	double err_est;
  int err_code;
	double res = integrate(f, lower, upper, err_est, err_code);
	if(err_est>0.01) res = 0;
	if(err_code!=0) res = 0;
	return res;
}





// [[Rcpp::export]]

double  RCPPcdf_MLBA_RTI(arma::mat X, arma::vec beta, arma::vec zeta, double lam1, double lam2, double I0, double b, double s, double A, double choice, double rt){
	const double lower = 0, upper = rt;
	MLBA_CDF f(X, beta,zeta,lam1, lam2,I0, b,s, A,  choice) ;
	double err_est;
    	int err_code;
	double res = integrate(f, lower, upper, err_est, err_code);
	//if(err_est>0.01) res = 0;
	if(err_code!=0) res = 0;
	return res;
}

////////////////////////////////////////////////////////


// [[Rcpp::export]]

arma::vec RCPPcdf_MLBA_CO_all(arma::mat X, arma::vec beta, arma::vec zeta, double lam1, double lam2, double I0, double b, double s, double A, arma::vec choice){
	double nAlt = zeta.n_elem;
	double nObs = X.n_rows/nAlt;
	arma::vec MLBA = -ones(nObs);
	for (int i = 1; i <= nObs; ++i) {
		
		MLBA(i-1) = RCPPcdf_MLBA_CO(X.rows(nAlt*(i-1), (nAlt*i-1)), beta, zeta, lam1, lam2, I0, b, s, A, choice(i-1));
	}  

	return MLBA;
}


// [[Rcpp::export]]

double RCPPMLBA_Lik_CO(arma::mat X, arma::vec beta, arma::vec zeta, double lam1, double lam2, double I0, double b, double s, double A, arma::vec choice){
	arma::vec MLBA = RCPPcdf_MLBA_CO_all(X, beta, zeta, lam1, lam2, I0, b, s, A, choice);
	double log_Lik = accu(log(MLBA+1e-12));
	return log_Lik;
}


// [[Rcpp::export]]
double RCPP_brier_CO(arma::mat X, arma::vec  beta, arma::vec zeta, double lam1, double lam2, double I0, double b, double s, double A, arma::vec choice){
	double nAlt = zeta.n_elem;
	double nObs = X.n_rows/nAlt;
	double temp = 0;
	arma::vec choice_o = -ones(nObs);
	
	for(int i =1;i<=nAlt; ++i){
	  arma::vec choice_s = ones(nObs)*i;

	  
	  for (int j=0;j<nObs;++j){
	    double temp2= choice(j)-choice_s(i);
	    if(temp2==0){
	      choice_o(j) = 1;
	    }
	    else{
	      choice_o(j) = 0;
	    }
	  }
	  arma::vec MLBA_s = RCPPcdf_MLBA_CO_all( X, beta, zeta, lam1, lam2, I0, b, s, A, choice_s);
	  temp = temp+accu(pow(MLBA_s-choice_o,2));
	}
	temp = temp/nObs;
	return temp;
}


//////////////////////////////////////////////

// [[Rcpp::export]]

arma::vec RCPPcdf_MLBA_RTI_all(arma::mat X, arma::vec beta, arma::vec zeta, double lam1, double lam2, double I0, double b, double s, double A, arma::vec choice,arma::vec rt){
	double nAlt = zeta.n_elem;
	double nObs = X.n_rows/nAlt;
	arma::vec MLBA = -ones(nObs);
	for (int i = 1; i <= nObs; ++i) {
		
		MLBA(i-1) = RCPPcdf_MLBA_RTI(X.rows(nAlt*(i-1), (nAlt*i-1)), beta, zeta, lam1, lam2, I0, b, s, A, choice(i-1),rt(i-1));
	}  

	return MLBA;
}


// [[Rcpp::export]]

double RCPPMLBA_Lik_RTI(arma::mat X, arma::vec beta, arma::vec zeta, double lam1, double lam2, double I0, double b, double s, double A, arma::vec choice,arma::vec rt){
	arma::vec MLBA = RCPPcdf_MLBA_RTI_all(X, beta, zeta, lam1, lam2, I0, b, s, A, choice,rt);
	double log_Lik = accu(log(MLBA));
	return log_Lik;
}

/////////////////////////////////////////////////////////////////// for indentification


// [[Rcpp::export]]
double RCPPcdf_rt_v2(double drift,double b, double s, double A, int choice, double rt){
	double d = drift;
	double cdf = (1 + (b-A-d*rt)/A* normcdf((b-A-d*rt)/(s*rt)) - (b-d*rt)/A*normcdf((b-d*rt)/(s*rt)) + s*rt/A* normpdf((b-A-d*rt)/(s*rt))-s*rt/A* normpdf((b-d*rt)/(s*rt)))/normcdf(d/s);
	if (A<1e-30) {cdf = normcdf(b/rt,d,s)/normcdf(d/s);}
	if(cdf<1e-30|cdf>1.1) {cdf =0;}
	
	return cdf;
	
	
}





// [[Rcpp::export]]
double RCPPpdf_rt_v2(double drift,double b, double s, double A, double choice, double rt){
	double d = drift;
  	
  double pdf_imp = (-d*normcdf((b-A-d*rt)/(s*rt))+d*normcdf((b-d*rt)/(s*rt))+s*normpdf((b-A-d*rt)/(s*rt))-s*normpdf((b-d*rt)/(s*rt)))/A/normcdf(d/s);
	
	
	if (A<1e-10) pdf_imp = normpdf(b/rt,d,s)*b/rt/rt/normcdf(d/s);
	if(pdf_imp<1e-30) {pdf_imp = 0;}
	if(pdf_imp>1e10) {pdf_imp = 1e10; printf("I am here.");}
	return pdf_imp;
	
}


// [[Rcpp::export]]
double RCPPpdf_MLBA_v2(arma::vec drift,double b, double s, double A,double choice, double rt){
  double temp = 1;
	double Alts = drift.n_elem;

	for (int i = 1; i <= Alts; ++i) {
		if(i!=choice){
			temp = temp*(1- RCPPcdf_rt_v2(drift(i-1), b, s, A, i, rt));	
		}
	}
	
	double pdf_MLBA_imp = temp*RCPPpdf_rt_v2(drift(choice-1), b, s, A, choice, rt);
	
	if (pdf_MLBA_imp<1e-30) {pdf_MLBA_imp = 0;}
	//if (pdf_MLBA_imp<1e-10) {pdf_MLBA_imp = 0;}
	if(pdf_MLBA_imp>1e10) {pdf_MLBA_imp = 1e10;}
	//if(pdf_MLBA_imp>100) {pdf_MLBA_imp = 100;}
	return pdf_MLBA_imp;
}


// [[Rcpp::export]]
double RCPPpdf_MLBA_rtknown_v2(arma::vec drift, double b, double s, double A, double choice, double rt){

	double Alts = drift.n_elem;
	double temp = 0;
  	//double test = 1;
	for (int i = 1; i <= Alts; ++i) {
		temp = temp + RCPPpdf_MLBA_v2(drift, b, s, A, i, rt);
	 //if(RCPPpdf_MLBA(X, beta, zeta, lam1, lam2, I0, b, s, A, i, rt)<0){
	 //  test = 0;
	 //}
	} 
	//temp = temp*test;
	double cdf_MLBA;
 	cdf_MLBA = RCPPpdf_MLBA_v2(drift, b, s, A, choice, rt);
 	//if (cdf_MLBA<1e-10) cdf_MLBA = 1e-10;
 	// Evan's assumption
	if(temp==0) {	cdf_MLBA = 0;}

	// Scott's  assumption
	//if(temp==0) {cdf_MLBA = 1/Alts;}

	else {
	cdf_MLBA = RCPPpdf_MLBA_v2(drift, b, s, A, choice, rt)/temp;
	} 

	return cdf_MLBA;
}



// [[Rcpp::export]]
double RCPPpdf_MLBA_crt_v2(arma::vec drift, double b, double s, double A, double choice, double rt){

	
	double cdf_MLBA;
 	cdf_MLBA = RCPPpdf_MLBA_v2(drift, b, s, A, choice, rt);

	return cdf_MLBA;
}

// [[Rcpp::export]]
arma::vec MLBA_rtknown_all_v2(arma::mat drift, double b, double s, double A, arma::vec choice, arma::vec rt){
	double nAlt = drift.n_cols;
	double nObs = drift.n_rows;
	arma::vec MLBA = -ones(nObs);
	for (int i = 1; i <= nObs; ++i) {
		
		MLBA(i-1) = RCPPpdf_MLBA_rtknown_v2( drift.row(i-1), b, s, A, choice(i-1),rt(i-1));
	}  

	return MLBA;
}

// [[Rcpp::export]]

double RCPPMLBA_Lik_rtknown_v2(arma::mat drift, double b, double s, double A, arma::vec choice, arma::vec rt){
	arma::vec MLBA = MLBA_rtknown_all_v2( drift, b, s, A, choice, rt);
	double log_Lik =  accu(log(MLBA+1e-12));
	return log_Lik;
}



// [[Rcpp::export]]
arma::vec MLBA_crt_all_v2(arma::mat drift, double b, double s, double A, arma::vec choice, arma::vec rt){
  double nAlt = drift.n_cols;
  double nObs = drift.n_rows;
  arma::vec MLBA = -ones(nObs);
  for (int i = 1; i <= nObs; ++i) {
    
    MLBA(i-1) = RCPPpdf_MLBA_crt_v2( drift.row(i-1), b, s, A, choice(i-1),rt(i-1));
  }  
  
  return MLBA;
}

// [[Rcpp::export]]

double RCPPMLBA_Lik_crt_v2(arma::mat drift, double b, double s, double A, arma::vec choice, arma::vec rt){
	arma::vec MLBA = MLBA_crt_all_v2( drift, b, s, A, choice, rt);
	double log_Lik =  accu(log(MLBA+1e-12));
	return log_Lik;
}



class MLBA_CDF_v2: public Func
{
private:
	arma::vec drift;
	double b;
	double s;
 	double A;
	double choice;

public:
    MLBA_CDF_v2(arma::vec drift_,double b_, double s_, double A_,double choice_) : drift(drift_), b(b_),s(s_),A(A_),choice(choice_){}
  
    double operator()(const double& x) const
    {
      //return RCPPpdf_rt_v2(drift(choice-1), b ,s, A,choice,x);
      return RCPPpdf_MLBA_v2(drift, b ,s, A,choice,x);
      
    }
};


// [[Rcpp::export]]

double  RCPPcdf_MLBA_CO_v2(arma::vec drift, double b, double s, double A, double choice){
	const double lower = 0, upper = 100;
	MLBA_CDF_v2 f(drift, b,s, A,choice);
	double err_est;
	int err_code;
	double res = integrate(f, lower, upper, err_est, err_code,100,1e-5);//,1000,0.001);
	//if(err_est>0.01) res = -1;
	if(err_code!=0) res = -166;//err_code;
	return res;
}

// [[Rcpp::export]]

arma::vec RCPPcdf_MLBA_CO_all_v2(arma::mat drift, double b, double s, double A, arma::vec choice){
  int nAlt = drift.n_cols;

  int nObs = drift.n_rows;
  arma::vec MLBA = -ones(nObs);
  arma::vec tmp = ones(2);
  for (int i = 1; i <= nObs; ++i) {
    for(int j =1;j<=nAlt;++j){
     tmp(j-1) =  drift.row(i-1)(j-1);
    }
    //arma::vec(drift.row(i-1))
    MLBA(i-1) = RCPPcdf_MLBA_CO_v2(tmp, b, s, A, choice(i-1));
  }  
  
  return MLBA;
  //return drift.row(0);
  
}


// [[Rcpp::export]]

double RCPPMLBA_Lik_CO_v2(arma::mat drift, double b, double s, double A, arma::vec choice){
  arma::vec MLBA = RCPPcdf_MLBA_CO_all_v2(drift, b, s, A, choice);
  double log_Lik = accu(log(MLBA+1e-12));
  return log_Lik;
}
