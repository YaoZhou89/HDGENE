#include <Rcpp.h>
#include <cmath>
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
  for(int i = 0; i < m ; i++){
    for(int j =0; j < i ; j++){
      double pc = 0;
      for(int k = 0; k < K; k++){
        pc +=  GD(i,k)*GD(j,k)*Y(k);
      }
      corr(i,j) = pc;
    }
  }
  return wrap(corr);
}

// NumericVector product(NumericMatrix X, NumericVector Y, int n,int P){
//   // X %*% Y
//   //X: n*p; Y: p*1; prod: n*1
//   NumericVector prod(n);
//   for(int p = 0; p < n; p++){
//     for(int k = 0; k < P; k++){
//       prod(p) += X(p,k)*Y(k) ;
//     }
//   }
//   return prod;
// }
//
// NumericMatrix Mproduct(NumericMatrix X, NumericMatrix Y, int a,int b1,int b2,int c ){
//   // X * Y
//   NumericMatrix prod(a,c);
//   for(int p = 0; p < a; p++){
//     for(int k = 0; k < c; k++){
//       for( int t = 0; t < b1; t++){
//         prod(p,k) += X(p,t)*Y(t,k) ;
//       }
//     }
//   }
//   return prod;
// }
//
// NumericMatrix Mcrossproduct(NumericMatrix X, NumericMatrix Y, int a,int b,int a2,int c ){
//   // T(X) * Y
//   NumericMatrix prod(b,c);
//   for(int p = 0; p < b; p++){
//     for(int k = 0; k < c; k++){
//       for( int t = 0; t < a; t++){
//         prod(p,k) += X(t,p)*Y(t,k) ;
//       }
//     }
//   }
//   return prod;
// }
//
// NumericVector crossprod(NumericMatrix X, NumericVector Y, int P,int K){
//   // t(X) %*% Y
//   //X: p*n; Y: p*1, prod: p*1
//   NumericVector prod(P);
//   for(int p = 0; p < P; p++){
//     for(int k = 0; k < K; k++){
//       prod(p) += X(k,p)*Y(k) ;
//     }
//   }
//   return prod;
// }
//
// NumericVector getInteract(NumericVector X, NumericVector Y,int m){
//   NumericVector prod(m);
//   for(int p = 0; p < m; p++){
//       prod(p) = X(p)*Y(p) ;
//   }
//   return prod;
// }
// double colSumsq(NumericVector X, int n){
//   double sum = 0;
//   for ( int i = 0; i < n; i++){
//     sum += X(i)*X(i);
//   }
//   return sum;
// }
// double getAbs (NumericVector X, NumericVector Y,int n){
//   double r = 0;
//   for (int i = 0; i < n; i++){
//     r += X(i)*Y(i);
//   }
//   if (r < 0) r = -r;
//   return r;
// }
// // // NumericVector getMinus(NumericVector X, )
// // // inline static double sqrt_double( double x ){ return ::sqrt( x ); }
// // [[Rcpp::export]]
// NumericVector  geno_cor (NumericVector & Y, NumericMatrix & GD, NumericMatrix & w,int m, int n,int p) {
//   NumericMatrix corr(m,m);
//   NumericVector GDi(n);
//   NumericVector GDt(n),GDtt(n),GDttt(n);
//   double div,r,colsq;
//   for (int i = 0 ; i < m; i++){
//     for(int j = 0; j < i + 1; j++){
//       r = 0;
//       GDi = getInteract(GD(i,_),GD(j,_),n);
//       GDt = product(w,crossprod(w,GDi,p,n),n,p);
//       GDtt = GDi - GDt;
//       colsq = colSumsq(GDtt,n);
//       div = sqrt(colsq);
//       GDttt = GDtt/div;
//       r = getAbs(GDttt,Y,n);
//       corr(i,j) = r;
//     }
//   }
//   return corr;
// }
