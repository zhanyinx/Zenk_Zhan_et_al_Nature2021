#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix cresize(NumericMatrix x,int size) {
      int nrow = x.nrow();
      int i,j,row,column;
      NumericMatrix out(size,size);
      NumericMatrix count(size,size);
      float binsize=float(nrow)/float(size);
      for(i=0;i<nrow;i++){
        row=int(i/binsize);
	if(row>size-1) row=size-1;
        for(j=0;j<nrow;j++){
          column=int(j/binsize);
	  if(column>size-1) column=size-1;
          out(row,column)=out(row,column)+x(i,j);
          count(row,column)=count(row,column)+1;
        }
      }
      for(i=0;i<size;i++){
        for(j=0;j<size;j++){
          if(count(i,j)>0) out(i,j)=out(i,j)/count(i,j);
        }
      }
      return out;
   }

// [[Rcpp::export]]
NumericMatrix cresize_nozeros(NumericMatrix x,int size) {
      int nrow = x.nrow();
      int i,j,row,column;
      NumericMatrix out(size,size);
      NumericMatrix count(size,size);
      float binsize=float(nrow)/float(size);
      for(i=0;i<nrow;i++){
        row=int(i/binsize);
        if(row>size-1) row=size-1;
        for(j=0;j<nrow;j++){
          column=int(j/binsize);
          if(x(i,j)!=-1){
            if(column>size-1) column=size-1;
            out(row,column)=out(row,column)+x(i,j);
            count(row,column)=count(row,column)+1;
          }
        }
      }
      for(i=0;i<size;i++){
        for(j=0;j<size;j++){
          if(count(i,j)>0) out(i,j)=out(i,j)/count(i,j);
        }
      }
      return out;
   }
