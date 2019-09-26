#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector vectorresize(NumericVector x,int size) {
      int nrow = x.size();
      int i,j,row;
      NumericVector out(size);
      NumericVector count(size);
      float binsize=float(nrow)/float(size);
      for(i=0;i<nrow;i++){
        row=int(i/binsize);
	if(row>size-1) row=size-1;
        out(row)=out(row)+x(i);
	count(row)+=count(row);
      }
      for(i=0;i<size;i++){
        if(count(i)>0) out(i)=out(i)/count(i);
      }
      return out;
   }
