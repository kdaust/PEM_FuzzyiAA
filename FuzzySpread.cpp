#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
void rSpread(NumericMatrix mat, double diag, double near){
  int nrow = mat.nrow();
  int ncol = mat.ncol();
  
  for(int i = 1; i < nrow-1; i++){
    for(int j = 1; j < ncol-1; j++){
      if(mat(i,j) == 1){
        //Rcout << i << j << "!NA \n";
        for(int a = -1; a <= 1; a++){
          for(int b = -1; b <= 1; b++){
            if(NumericMatrix::is_na(mat(i+a,j+b)) || (near > mat(i+a,j+b))){
              if(a != 0 && b != 0){
                mat(i+a,j+b) = diag;
              }else{
                mat(i+a,j+b) = near;
              }
            }
          }
        }
      }
    }
  }
}


