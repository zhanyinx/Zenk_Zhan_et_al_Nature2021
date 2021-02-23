/*
 * memory.c
 *
 *  Created on: Nov 18, 2009
 *      Author: guido
 */
#include <stdlib.h>
#include "mhistogram.h"

/*****************************************************************************
 Allocate double matrix in the form x[l][m]
 *****************************************************************************/
double **AlloDoubleMat(int l, int m)
{
  double **x;
  int i,j;

  x = (double **) calloc(l,sizeof(double *));
  if (!x) Error("Cannot allocate double matrix");

  for (i=0;i<l;i++)
    {
      *(x+i) = calloc(m,sizeof(double));
      if (!(*(x+i))) Error("Cannot allocate double matrix");
      for (j=0;j<m;j++) *(*(x+i)+j) = 0.;
    }

  return x;
}

/*****************************************************************************
 Allocate double vector
 *****************************************************************************/
double *AlloDouble(int l)
{
  double *x;
  int i;

  x = (double *) calloc(l,sizeof(double));
  if (!x) Error("Cannot allocate double matrix");

  for (i=0;i<l;i++) *(x+i) =0;

  return x;
}

/*****************************************************************************
 Allocate int vector
 *****************************************************************************/
int *AlloInt(int l)
{
  int *x;
  int i;

  x = (int *) calloc(l,sizeof(int));
  if (!x) Error("Cannot allocate int");

  for (i=0;i<l;i++) *(x+i) =0;

  return x;
}
