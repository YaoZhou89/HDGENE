#include <Rcpp.h>
using namespace Rcpp;
using namespace std;
// This is the C++ code for HDGENE.cor
//
// NumericMatrix  geno_cor_new (NumericVector & Y, NumericMatrix & GD, NumericMatrix & w,NumericMatrix & corr,int m, int K) {
//   for(int i = 1; i < m; i++){
//     for(int j =0; j < i; j++){
//       double pc = 0;
//       for(int k = 0; k < K; k++){
//         pc +=  GD(i,k)*GD(j,k)*Y(k);
//       }
//       corr(i,j) = pc;
//     }
//   }
//   return corr;
// }
// [[Rcpp::export]]
RcppExport SEXP  geno_cor_new ( SEXP  YY,SEXP GDD,SEXP ww,SEXP coorr,SEXP mm,SEXP KK) {
  NumericVector Y  = Rcpp::as<NumericVector>(YY);
  NumericMatrix GD = Rcpp::as<NumericMatrix>(GDD);
  NumericMatrix w = Rcpp::as<NumericMatrix>(ww);
  NumericMatrix corr = Rcpp::as<NumericMatrix>(coorr);
  int m =  Rcpp::as<int>(mm);
  int K =  Rcpp::as<int>(KK);
  for(int i = 1; i < m; i++){
    for(int j =0; j < i; j++){
      double pc = 0;
      for(int k = 0; k < K; k++){
        pc +=  GD(i,k)*GD(j,k)*Y(k);
      }
      corr(i,j) = pc;
    }
  }
  return wrap(corr);
}

// NumericVector product(NumericMatrix X, NumericVector Y, int P,int K){
//   NumericVector prod(P);
//   for(int p = 0; p <P; p++){
//     for(int k = 0; k < K; k++){
//       prod(p) += X(p,k)*Y(k) ;
//     }
//   }
//   return prod;
// }
// NumericVector crossprod(NumericMatrix X, NumericVector Y, int P,int K){
//   NumericVector prod(P);
//   for(int p = 0; p <P; p++){
//     for(int k = 0; k < K; k++){
//       prod(p) += X(k,p)*Y(k) ;
//     }
//   }
//   return prod;
// }
// // [[Rcpp::export]]
// NumericMatrix  geno_cor (NumericVector & Y, NumericMatrix & GD, NumericMatrix & tw,NumericMatrix & corr,int m, int K,int P) {
//   NumericVector GDi(K);
//   for (int k =0; k <K; k++){
//
//   }
//   GDi = product(GD,GD,1,K);
//   NumericVector GDj(P);
//   GDj = product(tw,GDi,P,K);
//
//   for(int i = 0; i < m; i++){
//     for(int j =0; j < i+1; j++){
//       double pc = 0;
//       NumericVector GDi(K);
//       NumericVector GDj(K);
//       for(int k = 0; k < K; k++){
//         GDi(k)=  GD(i,k)*GD(j,k);
//         for (int p=0;p<P;p++){
//           GDj(k) -= tw(p,k)*tw(p,k)*GDi(k);
//         }
//       }
//
//       for(int k = 0; k < K; k++){
//         pc +=  GD(i,k)*GD(j,k)*Y(k);
//       }
//       corr(i,j) = pc;
//     }
//   }
//   return corr;
// }
