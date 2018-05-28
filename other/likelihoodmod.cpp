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
double likelihoodmod(NumericVector theta, NumericVector y, NumericMatrix X,int J,double lambda1){
  
  int N = y.size();
  int K= X.ncol();
  NumericVector alpha(J+1);
  NumericVector lambda(J);
  NumericVector beta(K);
  alpha[0]=-9223372036854775807;
  alpha[J]=-alpha[0];
  for (int i=0 ;i<J-1;i++)
  {
    alpha[i+1]=theta[i];
  }
  lambda[0]=lambda1;
  // lambda[0]=n[0]/N;
  for (int i=0 ;i<J-1;i++)
  {
    lambda[i+1]=theta[i+J+K-1];
  }
  for (int i=0 ;i<K;i++)
  {
    beta[i]=theta[i+J-1];
  }
  double ln=0;

  for(int i=0;i<N;i++)
  {
    double betax=0;
    for(int k=0;k<K;k++){
      betax+= beta[k]*X(i,k);
    }
    
    for (int  j=0 ;j<J;j++)
    { 
      if (y[i]==j+1)
      {
        double a=lambda[j]*(1/(1+exp(-alpha[j+1]+betax)) - 1/(1+exp(-alpha[j]+betax)) );
        double b=0;
        for (int k=0;k<J;k++)
        {
          b+=lambda[k]*( 1/(1+exp(-alpha[k+1]+betax)) - 1/(1+exp(-alpha[k]+betax)) );
        }
        ln=ln+log(a/b);
      }
    }
  }

  return ln;
  
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//







