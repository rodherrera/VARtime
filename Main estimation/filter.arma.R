library(Rcpp)
library(RcppArmadillo)

# filtro ARMA(1,1)
cppFunction(depends = "RcppArmadillo", ' 
arma::mat filter_cpp_arma(arma::mat ft,
                          arma::mat a,
                          arma::mat x) {
 
   int nrow_ft = ft.n_rows-1;
 
     for ( int i=0; i<nrow_ft; i++ ) {
     ft.row(i+1) =  a.row(0) + a.row(1) % ft.row(i) - a.row(2) % x.row(i);
   }
   return ft ;
}
')

cppFunction(depends = "RcppArmadillo", ' 
arma::mat rcppAR(arma::mat ft,
                          arma::mat a,
                          arma::mat x) {
 
   int nrow_ft = ft.n_rows-1;
 
     for ( int i=0; i<nrow_ft; i++ ) {
     ft.row(i+1) =   a.row(0) % ft.row(i) - a.row(1) % x.row(i);
   }
   return ft ;
}
')



cppFunction(depends = "RcppArmadillo", '
arma::mat rcppVAR1(arma::mat A, arma::mat B, arma::mat C, arma::mat zt) {
  int m = zt.n_rows; int n = zt.n_cols;
  arma::mat dt(m,n);
  dt.row(0) = arma::zeros<arma::mat>(1,n);
  for (int row=1; row<m; row++) {
    dt.row(row) = C.row(0) + dt.row(row-1)*trans(A) + B.row(0) % zt.row(row-1);
  }
  return dt;
}')

cppFunction(depends = "RcppArmadillo", '
arma::mat rcppVAR0(arma::mat A, arma::mat B, arma::mat C,arma::mat zt) {
  int m = zt.n_rows; int n = zt.n_cols;
  arma::mat dt(m,n);
  dt.row(0) = arma::zeros<arma::mat>(1,n);
  for (int row=1; row<m; row++) {
    dt.row(row) = C.row(0) + dt.row(row-1)*trans(A) + B.row(0) % zt.row(row-1);
  }
  return dt;
}')



cppFunction(depends = "RcppArmadillo", '
arma::mat rcppVAR2(arma::mat A, arma::mat B, arma::mat C, arma::mat D, arma::mat rv, arma::mat zt) {
  int m = zt.n_rows; int n = zt.n_cols;
  arma::mat dt(m,n);
  dt.row(0) = arma::zeros<arma::mat>(1,n);
  for (int row=1; row<m; row++) {
    dt.row(row) = C.row(0) + dt.row(row-1)*trans(A) + D.row(0) % rv.row(row-1) + B.row(0) % zt.row(row-1);
  }
  return dt;
}')



cppFunction(depends = "RcppArmadillo", ' 
arma::mat ewma(arma::mat ft,
                          arma::mat a,
                          arma::mat st) {
 
   int nrow_ft = ft.n_rows-1;
 
     for ( int i=0; i<nrow_ft; i++ ) {
     ft.row(i+1) =   (1-a.row(0)) % st.row(i+1) + a.row(0) % ft.row(i);
   }
   return ft ;
}
')

cppFunction(depends = "RcppArmadillo", '
arma::mat EWMA(arma::mat A, arma::mat zt) {
  int m = zt.n_rows; int n = zt.n_cols;
  arma::mat dt(m,n);
  dt.row(0) = arma::zeros<arma::mat>(1,n);
  for (int row=1; row<m; row++) {
    dt.row(row) = (1-A.row(0)) % zt.row(row) + A.row(0)%dt.row(row-1);
  }
  return dt;
}')


sourceCpp(
  code = 
    "
     #include <Rcpp.h>
     // [[Rcpp::export]]
     Rcpp::NumericVector ewmaRcpp(Rcpp::NumericVector x, double a){
       int n = x.length();
       Rcpp::NumericVector s(n);
       s[0] = x[0];
       if (n > 1) {
         for (int i = 1; i < n; i++) {
           s[i] = a * x[i] + (1 - a) * s[i-1];
         }
       }
       return s;
     }

    ")



cppFunction(depends = "RcppArmadillo", '
arma::mat rcppVAR3(arma::mat A, arma::mat B, arma::mat C, arma::mat D,  arma::mat E,arma::mat rv, arma::mat jumps, arma::mat zt) {
  int m = zt.n_rows; int n = zt.n_cols;
  arma::mat dt(m,n);
  dt.row(0) = arma::zeros<arma::mat>(1,n);
  for (int row=1; row<m; row++) {
    dt.row(row) =  C.row(0) + dt.row(row-1)*trans(A) + D.row(0) % rv.row(row-1) + E.row(0) % jumps.row(row-1) + B.row(0) % zt.row(row-1);
  }
  return dt;
}')



cppFunction(depends = "RcppArmadillo", '
arma::mat rcppVAR4(arma::mat A, arma::mat B, arma::mat C, arma::mat D, arma::mat F, arma::mat rv, arma::mat jumps, arma::mat zt) {
  int m = zt.n_rows; int n = zt.n_cols;
  arma::mat dt(m,n);
  dt.row(0) = arma::zeros<arma::mat>(1,n);
  for (int row=1; row<m; row++) {
    dt.row(row) =  C.row(0) + dt.row(row-1)*trans(A) + D.row(0) % rv.row(row-1) + F.row(0) % jumps.row(row-1) + B.row(0) % zt.row(row-1);
  }
  return dt;
}')
