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
  // beta is non-negative here
	arma::mat  temp2= diffX * temp;
	
	arma::mat W = exp(-lam1*abs(temp2)%(temp2>=0)-lam2*abs(temp2)%(temp2<0));
	return W;
}


// [[Rcpp::export]]
double RCPPdriftmean(arma::mat X, arma::vec beta, arma::vec zeta, double lam1, double lam2, int choice ){
	arma::mat diffX = RCPPdiffX(X, choice);


	arma::mat W = RCPPweight(diffX, beta, lam1, lam2);
	//double d = accu(W%diffX*diagmat(beta))/accu(W)+I0+zeta(choice-1);
	
	double d = accu(W%diffX*diagmat(beta))+zeta(choice-1);//
	//double d = accu(W%diffX*diagmat(beta))+zeta(choice-1)+accu(X.row(choice-1)*diagmat(beta));
	// d can be negative
	//d = d*(d>0);
	//d= exp(d);
	return d;	
}




// [[Rcpp::export]]
double RCPPdriftmean_reverse(arma::vec drift_o, arma::vec drift_a, double drift_t){
//listed from smallest to largest;
	int index = accu(drift_a<drift_t);
	double res;
	if(index>=drift_a.n_elem){
		res = drift_t;
		}
	else{
	  if(drift_t >= 0.15){
	    
	    res = drift_o(index);
	  }
		else{
		  res = -10;
		}
	}

	return res;
}




// [[Rcpp::export]]
double RCPPdriftmean_last(arma::mat X, arma::vec beta, arma::vec zeta, double lam1, double lam2,double b, double s, double A, arma::vec choice,double rt_mean, arma::vec drift_o, arma::vec drift_a){
//This is a function to transfer drift rate mean to actual dirft rate mean
	double nAttr = X.n_rows;
	double nAlt = zeta.n_elem;
  	double nObs = X.n_rows/nAlt;
	arma::vec drift_mean = zeros(nObs-1);

	double temp;

	for (int i =1;i<=nObs-1;i++){	
		temp =  RCPPdriftmean(X.rows((i-1)*nAlt,i*nAlt-1), beta, zeta, lam1, lam2, choice(i-1));
		drift_mean(i-1) = 1/(temp+normpdf(temp/s)/normcdf(temp/s)*s);
	}

	temp = rt_mean*nObs/(b-A/2)-accu(drift_mean);

	if(temp>0) {
	//temp = RCPPdriftmean_reverse(drift_o,drift_a, temp*count/num(nAlt-1));
	temp = RCPPdriftmean_reverse(drift_o,drift_a, 1/temp);
	}
	else {
	temp = -10;
}

	return temp;
}





// [[Rcpp::export]]
arma::vec RCPPdriftmean_last_test(arma::mat X, arma::vec beta, arma::vec zeta, double lam1, double lam2,double b, double s, double A, arma::vec choice,double rt_mean, arma::vec drift_o, arma::vec drift_a){
//This is a function to transfer drift rate mean to actual dirft rate mean
	double nAttr = X.n_rows;
	double nAlt = zeta.n_elem;
  	double nObs = X.n_rows/nAlt;
	arma::vec drift_mean =zeros(nObs-1);

	double temp;

	for (int i =1;i<=nObs-1;i++){	
		temp =  RCPPdriftmean(X.rows((i-1)*nAlt,i*nAlt-1), beta, zeta, lam1, lam2, choice(i-1));
		drift_mean(i-1) = 1/(temp+normpdf(temp/s)/normcdf(temp/s)*s);
	}

	return drift_mean ;
}

// [[Rcpp::export]]
double RCPPcdf_rt(arma::mat X, arma::vec beta, arma::vec zeta, double lam1, double lam2, double b, double s, double A, int choice, double rt){
	int nAlt =  X.n_rows;
	double d;
	double cdf;
	d = RCPPdriftmean(X, beta, zeta, lam1, lam2, choice);
	if(d<(-6)) cdf = 1e-10;
	else{
	
	cdf = (1 + (b-A-d*rt)/A* normcdf((b-A-d*rt)/(s*rt)) - (b-d*rt)/A* normcdf((b-d*rt)/(s*rt)) + s*rt/A* normpdf((b-A-d*rt)/(s*rt))-s*rt/A* normpdf((b-d*rt)/(s*rt)))/normcdf(d/s);

	
	if (A<1e-10) cdf= normcdf(b/rt,d,s)/normcdf(d/s);
	if(cdf<1e-10) {cdf = 1e-10;}//{cdf = 1e-10;}
	}
	return cdf;
	
}

// [[Rcpp::export]]
double RCPPcdf_rt_last(arma::mat X, arma::vec beta, arma::vec zeta, double lam1, double lam2, double b, double s, double A, arma::vec choice,double rt_mean, arma::vec drift_o, arma::vec drift_a,double rt_last){
  	int nAlt = zeta.n_elem;
  	int nObs = X.n_rows/nAlt;
	  double cdf;
  	double d = RCPPdriftmean_last(X, beta, zeta, lam1, lam2, b, s, A, choice, rt_mean, drift_o, drift_a);
  	if(d<(-6)){cdf = 1e-10;}
	else{
		cdf = (1 + (b-A-d*rt_last)/A* normcdf((b-A-d*rt_last)/(s*rt_last)) - (b-d*rt_last)/A* normcdf((b-d*rt_last)/(s*rt_last)) + s*rt_last/A* normpdf((b-A-d*rt_last)/(s*rt_last))-s*rt_last/A* normpdf((b-d*rt_last)/(s*rt_last)))/normcdf(d/s);

  		if (A<1e-10) cdf= normcdf((b-d*rt_last)/(rt_last*s))/normcdf(d/s);
  		if(cdf<1e-10|d<(-6)) {cdf = 1e-10;}//{cdf = 1e-10;}
}
  
  return cdf;
  
}


// [[Rcpp::export]]
double RCPPpdf_rt(arma::mat X, arma::vec beta, arma::vec zeta, double lam1, double lam2, double b, double s, double A, double choice,double rt){
	int nAlt =  X.n_rows;
	double d;
	d = RCPPdriftmean(X, beta, zeta, lam1, lam2, choice);
	double pdf_tmp;
	if(d<(-6)) {pdf_tmp = 1e-10;}//{pdf_imp = 1e-10;}
	else{
	pdf_tmp = (-d*normcdf((b-A-d*rt)/(s*rt))+d*normcdf((b-d*rt)/(s*rt))+s*normpdf((b-A-d*rt)/(s*rt))-s*normpdf((b-d*rt)/(s*rt)))/A/normcdf(d/s);
	
	if (A<1e-10) pdf_tmp = normpdf(b/rt,d,s)*b/rt/rt/normcdf(d/s);
	if(pdf_tmp<1e-10) {pdf_tmp = 1e-10;}//{pdf_imp = 1e-10;}

}
		
	return pdf_tmp;
	
}


// [[Rcpp::export]]
double RCPPpdf_rt_last(arma::mat X, arma::vec beta, arma::vec zeta, double lam1, double lam2, double b, double s, double A, arma::vec choice,double rt_mean,arma::vec drift_o, arma::vec drift_a,double rt_last){
  	int nAlt = zeta.n_elem;
  	int nObs = X.n_rows/nAlt;
  	double d = RCPPdriftmean_last(X, beta, zeta, lam1, lam2, b, s, A, choice, rt_mean, drift_o, drift_a);
	double pdf;
  	if(d<-6) {pdf = 1e-10;}
	else{
	pdf = (-d*normcdf((b-A-d*rt_last)/(s*rt_last))+d*normcdf((b-d*rt_last)/(s*rt_last))+s*normpdf((b-A-d*rt_last)/(s*rt_last))-s*normpdf((b-d*rt_last)/(s*rt_last)))/A/normcdf(d/s);
  
  
  if (A<1e-10) pdf = normpdf(b/rt_last,d,s)*b/rt_last/rt_last/normcdf(d/s);
  if(pdf<1e-10) {pdf = 1e-10;}
}
  
  
  return pdf;
  
}


// [[Rcpp::export]]
double RCPPpdf_MLBA(arma::mat X, arma::vec beta, arma::vec zeta, double lam1, double lam2, double b, double s, double A,double choice,double rt){
	double temp = 1;
	double Alts = X.n_rows;

	for (int i = 1; i <= Alts; ++i) {
		if(i!=choice){
			temp = temp*(1- RCPPcdf_rt(X, beta, zeta, lam1, lam2,  b, s, A, i, rt));	
		}
	}
	
	double pdf_MLBA_imp = temp*RCPPpdf_rt(X, beta, zeta, lam1, lam2,  b, s, A, choice, rt);
	if (pdf_MLBA_imp<1e-10) {pdf_MLBA_imp = 1e-10;}//{pdf_MLBA_imp = 1e-10;}
	return pdf_MLBA_imp;
}


// [[Rcpp::export]]
double RCPPpdf_MLBA_last(arma::mat X, arma::vec beta, arma::vec zeta, double lam1, double lam2, double b, double s, double A, arma::vec choice,double rt_mean,arma::vec drift_o, arma::vec drift_a, double rt_last){
  double temp = 1;
  int Alts = zeta.n_elem;
  int nObs = X.n_rows/Alts;
  
  
  for (int i = 1; i <= Alts; ++i) {
    if(i!=choice(nObs-1)){
      temp = temp*(1- RCPPcdf_rt(X.rows((nObs-1)*3,(nObs*3-1)), beta, zeta, lam1, lam2,  b, s, A,i,rt_last));	
    }
  }
  
  double pdf_MLBA_imp = temp*RCPPpdf_rt_last(X, beta, zeta, lam1, lam2,  b, s, A, choice,rt_mean,drift_o,drift_a,rt_last);
  if (pdf_MLBA_imp<1e-10) {pdf_MLBA_imp = 1e-10;}//{pdf_MLBA_imp = 1e-10;}
  return pdf_MLBA_imp;
  
}


// [[Rcpp::export]]
double RCPPpdf_MLBA_rtg(arma::mat X, arma::vec beta, arma::vec zeta, double lam1, double lam2,  double b, double s, double A, double choice, double rt){
  
  double Alts = X.n_rows;
  double temp = 0;

  for (int i = 1; i <= Alts; ++i) {
  	temp = temp + RCPPpdf_MLBA(X, beta, zeta, lam1, lam2, b, s, A, i, rt);
  
  } 
  
  double cdf_MLBA;

  // Evan's assumption
  if(temp<1e-10) {	cdf_MLBA = 1e-10;}
  
  
  	else {
  	cdf_MLBA = RCPPpdf_MLBA(X, beta, zeta, lam1, lam2,  b, s, A, choice, rt)/temp;
  } 
  
  return cdf_MLBA;
}

// [[Rcpp::export]]
double RCPPpdf_MLBA_rtg_last(arma::mat X, arma::vec beta, arma::vec zeta, double lam1, double lam2, double b, double s, double A, arma::vec choice,double rt_mean,arma::vec drift_o, arma::vec drift_a, double rt_last){
  
  double Alts = zeta.n_elem;
  double nAttr = X.n_cols;
  double nObs = X.n_rows/Alts;
  double temp = 0;
  double temp1 = 0;
  double cdf_MLBA = 0;
  for (int i = 1; i <= Alts; ++i) {
	temp1 = RCPPpdf_MLBA_last(X, beta, zeta, lam1, lam2, b, s, A, choice,rt_mean,drift_o,drift_a,rt_last);
	if(i!=choice(nObs-1)){
		temp = temp + RCPPpdf_MLBA(X.rows((nObs-1)*Alts,(nObs*Alts-1)), beta, zeta, lam1, lam2, b, s, A, i, rt_last);
	}
	else{
		temp = temp + temp1;
	}
  	
  
  } 
  
  

  // Evan's assumption
  if(temp<1e-10) {	cdf_MLBA = 1e-10;}
  
  
  	else {
  	cdf_MLBA = temp1/temp;
  } 
  
  return cdf_MLBA;
}


// [[Rcpp::export]]
arma::vec MLBA_rtg_all(arma::mat X, arma::vec  beta, arma::vec zeta, double lam1, double lam2,  double b, double s, double A, arma::vec choice, double rt_mean,arma::vec drift_o, arma::vec drift_a, arma::vec rt){
  double nAlt = zeta.n_elem;
  double nObs = X.n_rows/nAlt;
  arma::vec MLBA = -ones(nObs);
  for (int i = 1; i <= nObs-1; ++i) {
    
    MLBA(i-1) = RCPPpdf_MLBA_rtg(X.rows(nAlt*(i-1), (nAlt*i-1)), beta, zeta, lam1, lam2, b, s, A, choice(i-1),rt(i-1));
  }  
    MLBA(nObs-1) = RCPPpdf_MLBA_rtg_last(X, beta, zeta, lam1, lam2, b, s, A, choice,rt_mean,drift_o,drift_a,rt(nObs-1));
  
  return MLBA;
}


// [[Rcpp::export]]

double RCPPMLBA_Lik_rtg(arma::mat X, arma::vec  beta, arma::vec zeta, double lam1, double lam2, double b, double s, double A, arma::vec choice, double rt_mean,arma::vec drift_o, arma::vec drift_a, arma::vec rt){
  arma::vec MLBA = MLBA_rtg_all( X, beta, zeta, lam1, lam2, b, s, A, choice,rt_mean,drift_o,drift_a, rt);
  double log_Lik =  accu(log(MLBA));
  return log_Lik;
}



// [[Rcpp::export]]
double RCPP_brier_RTG(arma::mat X, arma::vec  beta, arma::vec zeta, double lam1, double lam2,  double b, double s, double A, arma::vec choice, double rt_mean,arma::vec drift_o, arma::vec drift_a, arma::vec rt){
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
 
 	arma::vec MLBA_s = MLBA_rtg_all( X, beta, zeta, lam1, lam2, b, s, A, choice_s,rt_mean,drift_o,drift_a, rt);
	
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
	double b;
	double s;
 	double A;
	double choice;
	
	
public:
    MLBA_CDF(arma::mat X_, arma::vec beta_, arma::vec zeta_, double lam1_, double lam2_, double b_, double s_, double A_,double choice_) : X(X_), beta(beta_), zeta(zeta_), lam1(lam1_), lam2(lam2_), b(b_),s(s_),A(A_),choice(choice_){}

    double operator()(const double& x) const
    {
        return RCPPpdf_MLBA(X, beta, zeta, lam1, lam2, b ,s, A,choice,x);
    
}
};


class MLBA_CDF_last: public Func
{
private:
  arma::mat X;
  arma::vec beta;
  arma::vec zeta;
  double lam1;
  double lam2;
  double b;
  double s;
  double A;
  arma::vec choice;
  double rt_mean;
  arma::vec drift_o;
  arma::vec drift_a;
  
public:
  MLBA_CDF_last(arma::mat X_, arma::vec beta_, arma::vec zeta_, double lam1_, double lam2_, double b_, double s_, double A_,arma::vec choice_,double rt_mean_, arma::vec drift_o_, arma::vec drift_a_) : X(X_), beta(beta_), zeta(zeta_), lam1(lam1_), lam2(lam2_), b(b_),s(s_),A(A_),choice(choice_),rt_mean(rt_mean_),drift_o(drift_o_),drift_a(drift_a_){}
  
  double operator()(const double& x) const
  {
    return RCPPpdf_MLBA_last(X, beta, zeta, lam1, lam2, b ,s, A,choice,rt_mean,drift_o,drift_a,x);
    
  }
};


// [[Rcpp::export]]
double  RCPPcdf_MLBA_CO(arma::mat X, arma::vec beta, arma::vec zeta, double lam1, double lam2, double b, double s, double A, double choice){
	const double lower = 0, upper = 120;
	MLBA_CDF f(X, beta,zeta,lam1, lam2, b, s, A, choice);
	double err_est;
    	int err_code;
	double res = integrate(f, lower, upper, err_est, err_code,100,0.001);
	if(err_est>0.01) res = 1e-10;
	if(err_code!=0) res = 1e-10;
	return res;
}

// [[Rcpp::export]]
double  RCPPcdf_MLBA_CO_last(arma::mat X, arma::vec beta, arma::vec zeta, double lam1, double lam2, double b, double s, double A, arma::vec choice,double rt_mean,arma::vec drift_o, arma::vec drift_a){
  const double lower = 0, upper = 120;
  MLBA_CDF_last f(X, beta,zeta,lam1, lam2, b, s, A, choice, rt_mean, drift_o, drift_a);
  double err_est;
  int err_code;
  double res = integrate(f, lower, upper, err_est, err_code,100,0.001);
  if(err_est>0.01) res = 1e-10;
  if(err_code!=0) res = 1e-10;
  return res;
}




////////////////////////////////////////////////////////


// [[Rcpp::export]]

arma::vec RCPPcdf_MLBA_CO_all(arma::mat X, arma::vec beta, arma::vec zeta, double lam1, double lam2, double b, double s, double A, arma::vec choice,double rt_mean,arma::vec drift_o, arma::vec drift_a){
	double nAlt = zeta.n_elem;
	double nObs = X.n_rows/nAlt;
	arma::vec MLBA = -ones(nObs);
	
	for (int i = 1; i < nObs; ++i) {
	
		MLBA(i-1) = RCPPcdf_MLBA_CO(X.rows(nAlt*(i-1), (nAlt*i-1)), beta, zeta, lam1, lam2, b, s, A, choice(i-1));
	}  
  MLBA(nObs-1) = RCPPcdf_MLBA_CO_last(X, beta, zeta, lam1, lam2, b, s, A, choice,rt_mean, drift_o, drift_a);
	return MLBA;
}


// [[Rcpp::export]]

double RCPPMLBA_Lik_CO(arma::mat X, arma::vec beta, arma::vec zeta, double lam1, double lam2, double b, double s, double A, arma::vec choice,double rt_mean,arma::vec drift_o, arma::vec drift_a){
	arma::vec MLBA = RCPPcdf_MLBA_CO_all(X, beta, zeta, lam1, lam2, b, s, A, choice,rt_mean, drift_o, drift_a);
	double log_Lik = accu(log(MLBA));
	return log_Lik;
}


// [[Rcpp::export]]
double RCPP_brier_CO(arma::mat X, arma::vec  beta, arma::vec zeta, double lam1, double lam2,  double b, double s, double A, arma::vec choice,double rt_mean,arma::vec drift_o, arma::vec drift_a ){
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
	  arma::vec MLBA_s = RCPPcdf_MLBA_CO_all( X, beta, zeta, lam1, lam2,  b, s, A, choice_s,rt_mean,drift_o, drift_a);
	  
	  
	  temp = temp+accu(pow(MLBA_s-choice_o,2));
	}
	temp = temp/nObs;
	return temp;
}
//////////////////////////////////////////////  derivative calculation


//[[Rcpp::export]]
arma::vec derivative_di_beta(arma::mat X, arma::vec beta, arma::vec zeta, double lam1, double lam2, int choice){
  int nAlt = X.n_cols;
  int nAttr = X.n_rows;
  arma::mat diff_tmp =  RCPPdiffX(X, choice);
  arma::mat weight =   RCPPweight(diff_tmp,beta, lam1,lam2);
  arma::vec res = ones(nAttr);
  for (int i = 1;i<=nAttr;i++){
    // derivate of d_i in beta_k
    res(i-1) =  accu(weight.col(i-1)%diff_tmp.col(i-1)%(1+log(weight.col(i-1))));
  }
  return res;
  
}


//[[Rcpp::export]]
arma::vec derivative_di_lam(arma::mat X, arma::vec beta, arma::vec zeta, double lam1, double lam2, int choice){
  int nAlt = X.n_cols;
  int nAttr = X.n_rows;
  arma::mat diff_tmp =  RCPPdiffX(X, choice);
  arma::mat weight =   RCPPweight(diff_tmp,beta, lam1,lam2); //(nAlt-1)*nattr
  arma::vec res = ones(2);
  res(0) = -accu((weight%pow(diff_tmp,2)%(diff_tmp>=0))*pow(beta,2));
  res(1) = -accu((weight%pow(diff_tmp,2)%(diff_tmp<0))*pow(beta,2));
  
  return res;
  
}



//[[Rcpp::export]]
arma::mat derivative_eta_theta(arma::mat X, arma::vec beta, arma::vec zeta, double lam1, double lam2, int choice){
  //npar_max*(nAlt+1)
  int nAlt = X.n_cols;
  int nAttr = X.n_rows;
  int npar_max = nAttr+nAlt-1+2+1;
  arma::mat res = zeros(npar_max,nAlt+1);
  // input eta_b
  
  // input etai_beta & etai_lambda
  arma::vec tmp = zeros(npar_max);
  for(int i =1;i<=nAlt;i++){
    tmp = zeros(npar_max);
    tmp.subvec(0, nAttr-1) = derivative_di_beta(X, beta, zeta,lam1, lam2, i);
    if( i !=nAlt){tmp(nAttr-1+i)=1;}
    tmp.subvec(0, nAttr-1) = derivative_di_beta(X, beta, zeta,lam1, lam2, i);
    tmp.subvec(nAttr+nAlt-1, nAttr+nAlt) = derivative_di_lam(X, beta, zeta, lam1, lam2, i);
    res.col(i-1) = tmp;
  }
  
  res.at((npar_max-1),nAlt) = 1;
  return res;
}



////////////////////////////////////////
//[[Rcpp::export]]
double derivative_pdf_di(arma::mat X, arma::vec beta, arma::vec zeta, double lam1, double lam2, double b, double s, double A, int choice, double rt){
  
  double d = RCPPdriftmean(X, beta, zeta, lam1, lam2, choice);
  double pdf_imp1 = -normpdf(d/s)*(-d*normcdf((b-A-d*rt)/(s*rt))+d*normcdf((b-d*rt)/(s*rt))+s*normpdf((b-A-d*rt)/(s*rt))-s*normpdf((b-d*rt)/(s*rt)))/(A*normcdf(d/s)*normcdf(d/s)*s);
  double pdf_imp2 = (-normcdf((b-A-d*rt)/(s*rt))+((b-A)/(s*rt))*normpdf((b-A-d*rt)/(s*rt))+normcdf((b-d*rt)/(s*rt))-(b/(s*rt))*normpdf((b-d*rt)/(s*rt)))/(A*normcdf(d/s));
  
  double res = pdf_imp1+pdf_imp2;
  if (A<1e-10) res = -normpdf((b-d*rt)/(s*rt))*b*normpdf(d/s)/(rt*rt*s*normcdf(d/s)*normcdf(d/s))+normpdf((b-d*rt)/(s*rt))*b*(b-d*rt)/(rt*rt*rt*s*s*normcdf(d/s));
  
  return res;
}


//[[Rcpp::export]]
double derivative_cdf_di(arma::mat X, arma::vec beta, arma::vec zeta, double lam1, double lam2, double b, double s, double A, int choice, double rt){
  
  double d = RCPPdriftmean(X, beta, zeta, lam1, lam2, choice);
  
  
  
  
  double pdf_imp1 = -(1 + (b-A-d*rt)/A* normcdf((b-A-d*rt)/(s*rt)) - (b-d*rt)/A* normcdf((b-d*rt)/(s*rt)) + s*rt/A* normpdf((b-A-d*rt)/(s*rt))-s*rt/A* normpdf((b-d*rt)/(s*rt)))*normpdf(d/s)/(s*normcdf(d/s)*normcdf(d/s));
  double pdf_imp2 = rt*(normcdf((b-d*rt)/(s*rt))-normcdf((b-A-d*rt)/(s*rt)))/(A*normcdf(d/s));
  
  double res = pdf_imp1+pdf_imp2;
  normcdf(b/rt,d,s)/normcdf(d/s);
  if (A<1e-10) res = -normcdf(b/rt,d,s)/normcdf(d/s)/s*(normpdf(b/rt,d,s)/normcdf(b/rt,d,s)+normpdf(d/s)/normcdf(d/s));
  
  return res;
}

//[[Rcpp::export]]
double derivative_pdf_b(arma::mat X, arma::vec beta, arma::vec zeta, double lam1, double lam2, double b, double s, double A, int choice, double rt){
  
  double d = RCPPdriftmean(X, beta, zeta, lam1, lam2, choice);
  
  double res;
  res = (-(b-A)*normpdf((b-A-d*rt)/(s*rt))+b*normpdf((b-d*rt)/(s*rt)))/(rt*rt*s*A*normcdf(d/s));
  
  
  if (A<1e-10) res = normpdf((b-d*rt)/(s*rt))/(s*rt*normcdf(d/s));
  return res;
}

//[[Rcpp::export]]
double derivative_cdf_b(arma::mat X, arma::vec beta, arma::vec zeta, double lam1, double lam2, double b, double s, double A, int choice, double rt){
  
  double d = RCPPdriftmean(X, beta, zeta, lam1, lam2, choice);
  
  double res;
  res = (-normcdf((b-d*rt)/(s*rt))+normcdf((b-A-d*rt)/(s*rt)))/(A*normcdf(d/s));
  
  if (A<1e-10) res = normpdf((b-d*rt)/(s*rt))*(1-(b-d*rt)*b/(rt*rt*s*s))/(rt*rt*normcdf(d/s));
  return res;
}



//[[Rcpp::export]]
double derivative_mlba_crt_di(arma::mat X, arma::vec beta, arma::vec zeta, double lam1, double lam2, double b, double s, double A, int choice,int di, double rt){
  
  double temp = 1;
  double Alts = X.n_rows;
  double res;
  if(di==choice){
    for (int i = 1; i <= Alts; i++) {
      if(i!= di){
        temp = temp*(1- RCPPcdf_rt(X, beta, zeta, lam1, lam2,  b, s, A, i, rt));	
      }
    }
    res = temp*derivative_pdf_di(X,beta, zeta, lam1,lam2, b, s,A,choice,rt);
    
  }
  else{
    
    for (int i = 1; i <= Alts; i++) {
      if((i!= di)&(i!=choice)){
        temp = temp*(1- RCPPcdf_rt(X, beta, zeta, lam1, lam2,  b, s, A, i, rt));	
      }
    }
    res = -temp*derivative_cdf_di(X,beta, zeta, lam1,lam2, b, s,A,di,rt)*RCPPpdf_rt(X, beta, zeta, lam1, lam2,  b, s, A, choice, rt);
    
    
  }
  
  return res;
}





//[[Rcpp::export]]
double derivative_mlba_crt_b(arma::mat X, arma::vec beta, arma::vec zeta, double lam1, double lam2, double b, double s, double A, int choice, double rt){
  
  double temp = RCPPpdf_MLBA(X, beta, zeta, lam1, lam2, b, s, A, choice, rt);
  int Alts = X.n_rows;
  double res = 0;
  
  for (int i = 1;i<=Alts;i++){
    if(i!=choice){
      res = res-temp*derivative_cdf_b(X,beta, zeta, lam1,lam2, b, s,A,i,rt)/(1- RCPPcdf_rt(X, beta, zeta, lam1, lam2,  b, s, A, i, rt));
    }
  }
  res = res+temp*derivative_pdf_b(X,beta, zeta, lam1,lam2, b, s,A,choice,rt)/RCPPpdf_rt(X, beta, zeta, lam1, lam2,  b, s, A, choice, rt);
  
  
  return res;
}



class MLBA_deri_di: public Func
{
private:
  arma::mat X;
  arma::vec beta;
  arma::vec zeta;
  double lam1;
  double lam2;
  double b;
  double s;
  double A;
  double choice;
  double di;
public:
  MLBA_deri_di(arma::mat X_, arma::vec beta_, arma::vec zeta_, double lam1_, double lam2_, double b_, double s_, double A_,double choice_,double di_) : X(X_), beta(beta_), zeta(zeta_), lam1(lam1_), lam2(lam2_), b(b_),s(s_),A(A_),choice(choice_),di(di_){}
  
  double operator()(const double& x) const
  {
    return derivative_mlba_crt_di(X, beta, zeta, lam1, lam2, b ,s, A,choice,di,x);
  }
};

class MLBA_deri_b: public Func
{
private:
  arma::mat X;
  arma::vec beta;
  arma::vec zeta;
  double lam1;
  double lam2;
  double b;
  double s;
  double A;
  double choice;
public:
  MLBA_deri_b(arma::mat X_, arma::vec beta_, arma::vec zeta_, double lam1_, double lam2_, double b_, double s_, double A_,double choice_) : X(X_), beta(beta_), zeta(zeta_), lam1(lam1_), lam2(lam2_), b(b_),s(s_),A(A_),choice(choice_){}
  
  double operator()(const double& x) const
  {
    return derivative_mlba_crt_b(X, beta, zeta, lam1, lam2, b ,s, A,choice,x);
  }
};


// [[Rcpp::export]]
double  derivative_mlba_co_di(arma::mat X, arma::vec beta, arma::vec zeta, double lam1, double lam2, double b, double s, double A, double choice,double di){
  const double lower = 0, upper = 120;
  MLBA_deri_di f(X, beta,zeta,lam1, lam2, b, s, A,  choice,di) ;
  double err_est;
  int err_code;
  double res = integrate(f, lower, upper, err_est, err_code,100,0.001);
  if(err_est>0.01) res = 0;
  if(err_code!=0) res = 0;
  return res;
}

// [[Rcpp::export]]
double  derivative_mlba_co_b(arma::mat X, arma::vec beta, arma::vec zeta, double lam1, double lam2, double b, double s, double A, double choice){
  const double lower = 0, upper = 120;
  MLBA_deri_b f(X, beta,zeta,lam1, lam2, b, s, A,  choice) ;
  double err_est;
  int err_code;
  double res = integrate(f, lower, upper, err_est, err_code,100,0.001);
  if(err_est>0.01) res = 0;
  if(err_code!=0) res = 0;
  return res;
}




////////////////////////////////////////////////////

//[[Rcpp::export]]
arma::vec derivative_mlba_co_eta(arma::mat X, arma::vec beta, arma::vec zeta, double lam1, double lam2, double b, double s, double A, int choice){
  
  double nAlt = X.n_rows;
  arma::vec res = zeros(nAlt+1);
  // input crt_b
  res(nAlt) = derivative_mlba_co_b(X, beta, zeta, lam1, lam2, b, s, A, choice);
  
  // input crt_etai
  for(int i =1;i<=nAlt;i++){
    res(i-1) = derivative_mlba_co_di( X, beta, zeta, lam1, lam2, b, s, A, choice,i);}
  
  return res;
}


///////////////////////////////////////////////////// single fisher information matrix

//[[Rcpp::export]]
arma::mat derivative_mlba_co_per(arma::mat X, arma::vec beta, arma::vec zeta, double lam1, double lam2, double b, double s, double A, int choice){
  
  int nAlt = X.n_rows;
  int nAttr = X.n_cols;
  int maxpar = nAlt+nAttr-1+2+1;
  arma::mat res = zeros(maxpar,maxpar);
  
  double tmp = RCPPcdf_MLBA_CO(X, beta, zeta, lam1, lam2, b, s, A, choice);
  //if(tmp<= 1e-10) {tmp= 1e-10;}
  res = -derivative_eta_theta(X, beta, zeta, lam1, lam2, choice)*derivative_mlba_co_eta(X, beta, zeta, lam1, lam2, b, s, A, choice)*derivative_mlba_co_eta(X, beta, zeta, lam1, lam2, b, s, A, choice).t()*derivative_eta_theta(X, beta, zeta, lam1, lam2, choice).t()/(tmp*tmp);
  //res = res%(res>-1e2)%(res<1e2)-1e2*(res<-1e2)+1e2*(res>1e2);
  
  return res;
}

///////////////////////////////////////////////////// empirical fisher information matrix for all
// note only the first n-1 is calculated


//[[Rcpp::export]]
arma::mat derivative_mlba_co_all(arma::mat X, arma::vec beta, arma::vec zeta, double lam1, double lam2, double b, double s, double A, arma::vec choice){
  
  int nAlt = zeta.n_elem;
  int nAttrs = X.n_cols;
  int nObs = X.n_rows/nAlt;
  int maxpar = nAlt+nAttrs-1+2+1;
  arma::mat em_n_fisher = zeros(maxpar,maxpar);
  for (int i = 1; i <= nObs-1; ++i) {
    
    em_n_fisher = em_n_fisher+derivative_mlba_co_per(X.rows(nAlt*(i-1), (nAlt*i-1)), beta, zeta, lam1, lam2, b, s, A, choice(i-1));
  }  
  
  return em_n_fisher;
}

////////////////////////////////////
///////////////////////////////////////////////////// single score function


//[[Rcpp::export]]
arma::vec score_mlba_co_per(arma::mat X, arma::vec beta, arma::vec zeta, double lam1, double lam2, double b, double s, double A, int choice){
  
  int nAlt = X.n_rows;
  int nAttr = X.n_cols;
  int maxpar = nAlt+nAttr-1+2+1;
  arma::vec res = zeros(maxpar);
  double tmp = RCPPcdf_MLBA_CO(X, beta, zeta, lam1, lam2, b, s, A, choice);
  //if(tmp<= 1e-10) {tmp = 1e-10;}
  res = derivative_eta_theta(X, beta, zeta, lam1, lam2, choice)*derivative_mlba_co_eta(X, beta, zeta, lam1, lam2, b, s, A, choice)/tmp;
  //res = res%(res>(-1e2))%(res<1e2)-1e2*(res<=(-1e2))+1e2*(res>=1e2);
  
  
  return res;
}



///////////////////////////////////////////////////// empirical score information matrix for all

//[[Rcpp::export]]
arma::vec score_mlba_co_all(arma::mat X, arma::vec beta, arma::vec zeta, double lam1, double lam2, double b, double s, double A, arma::vec choice){
  
  int nAlt = zeta.n_elem;
  int nAttrs = X.n_cols;
  int nObs = X.n_rows/nAlt;
  int maxpar = nAlt+nAttrs-1+2+1;
  arma::vec em_score = zeros(maxpar);
  for (int i = 1; i <= nObs-1; ++i) {
    
    em_score = em_score+score_mlba_co_per(X.rows(nAlt*(i-1), (nAlt*i-1)), beta, zeta, lam1, lam2, b, s, A, choice(i-1));
  }  
  
  return em_score;
}



