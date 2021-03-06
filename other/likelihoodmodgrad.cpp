#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

//[[Rcpp::export]]
NumericVector likelihoodmodgrad(NumericVector theta, NumericVector y, NumericMatrix X,int J,double lambda1){
  
  int N = y.size();
  int K= X.ncol();
  NumericVector alpha(J+1);
  NumericVector palpha(J+1);
  NumericVector lambda(J);
  NumericVector beta(K);
  NumericVector g(2*J+K-2);
  alpha[0]=-9223372036854775807;
  alpha[J]=-alpha[0];
  for (int i=0 ;i<J-1;i++)
  {
    alpha[i+1]=theta[i];
  }
  for (int i=0 ;i<K;i++)
  {
    beta[i]=theta[i+J-1];
  }
 
  
  lambda[0]=lambda1;
  // lambda[0]=n[0]/N;
  for (int i=0 ;i<J-1;i++)
  {
    lambda[i+1]=theta[i+J+K-1];
  }

  for(int i=0;i<N;i++)
  {   
    double betax=0;
    for(int k=0;k<K;k++){
      betax= betax+beta[k]*X(i,k);
    }
    double dev=0;
    for(int k=0;k<J;k++){
      dev=dev+lambda[k]*(1/(1+exp(-alpha[k+1]+betax))-1/(1+exp(-alpha[k]+betax)));
    }
    double a=0;
    for(int k=0;k<J;k++){
      double jk=1/(1+exp(-alpha[k+1]+betax));
      double jkm=1/(1+exp(-alpha[k]+betax)); 
      a=a+lambda[k]*(jk*(1-jk)-jkm*(1-jkm));
    }
    
    if (y[i]==1)
    { 
      
      
      double pj=1/(1+exp(-alpha[1]+betax));
      double pjm=1/(1+exp(-alpha[0]+betax));  
      //beta_m
      
      for(int k=0;k<K;k++){
        g[J-1+k]=g[J-1+k]- X(i,k)*(1-pj-pjm)+a*X(i,k)/dev;
      }
      
      g[0]=g[0]+1-pj-(lambda[0]-lambda[1])*(1-pj)*pj/dev;
      
      for(int t=1;t<J-1;t++){
        double pt=1/(1+exp(-alpha[t+1]+betax));
        g[t]=g[t]-(lambda[t]-lambda[t+1])*(1-pt)*pt/dev;
      }
      for(int t=0;t<J-1;t++){
        g[t+J+K-1]=g[t+J+K-1]-(1/(1+exp(-alpha[t+2]+betax))-1/(1+exp(-alpha[t+1]+betax)))/dev;
      }
      
    }
    else if(y[i]<J){
      int l=y[i];
      
      
      
      double pj=1/(1+exp(-alpha[l]+betax));
      double pjm=1/(1+exp(-alpha[l-1]+betax));  
      
      for(int k=0;k<K;k++){
        g[J-1+k]=g[J-1+k]- X(i,k)*(1-pj-pjm)+a*X(i,k)/dev;
      }
      for (int t=0;t<J-1;t++){
        double jt=1/(1+exp(-alpha[t+1]+betax));
        g[t]=g[t]-(lambda[t]-lambda[t+1])*jt*(1-jt)/dev;;
      }
      
      g[l-2]=g[l-2]-pjm*(1-pjm)/(pj-pjm);
      g[l-1]=g[l-1]+pj*(1-pj)/(pj-pjm);
      
      for(int t=0;t<J-1;t++){
        g[J+t+K-1]=g[J+t+K-1]-(1/(1+exp(-alpha[t+2]+betax))-1/(1+exp(-alpha[t+1]+betax)))/dev;
      }
      g[J+l+K-3]=g[J+l+K-3]+1/lambda[l-1];
      
    }else{
      
      
      double pj=1/(1+exp(-alpha[J]+betax));
      double pjm=1/(1+exp(-alpha[J-1]+betax));  
      
      for(int k=0;k<K;k++){
        g[J-1+k]=g[J-1+k]- X(i,k)*(1-pj-pjm)+a*X(i,k)/dev;
      }
      
      g[J-2]=g[J-2]-pjm*(1-pjm)/(pj-pjm)-(lambda[J-2]-lambda[J-1])*pjm*(1-pjm)/dev;
      for(int t=0;t<J-2;t++){
        double pj = 1/(1+exp(-alpha[t+1]+betax));
        g[t]=g[t]-(lambda[t]-lambda[t+1])*pj*(1-pj)/dev;
      }
      
      for(int t=0;t<J-1;t++){
        g[J+t+K-1]=g[J+t+K-1]-(1/(1+exp(-alpha[t+2]+betax))-1/(1+exp(-alpha[t+1]+betax)))/dev;
      }
      g[2*J+K-3]=g[2*J+K-3]+1/lambda[J-1];
    }
    
    
  }
  
  return g;
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//




