/*
 * misc.c
 *
 *  Created on: Nov 6, 2009
 *      Author: guido
 */

#include "mhistogram.h"

/*****************************************************************************
 Allocate int matrix in the form x[l][m]
 *****************************************************************************/
int **AlloIntMat(int l, int m)
{
  int **x;
  int i,j;

  x = (int **) calloc(l,sizeof(int *));
  if (!x) Error("Cannot allocate int matrix");

  for (i=0;i<l;i++)
    {
      *(x+i) = calloc(m,sizeof(int));
      if (!(*(x+i))) Error("Cannot allocate int matrix");
      for (j=0;j<m;j++) *(*(x+i)+j) = 0;
    }

  return x;
}

/*****************************************************************************
 Allocate histogram structure
 *****************************************************************************/
struct histo_s *AlloHisto(int nhisto, int nbin)
{
	int i,j;
	struct histo_s *h;

	h = (struct histo_s *) calloc(nhisto,sizeof(struct histo_s));
	if (!h) Error("Cannot allocate histo_s");

	for (j=0;j<nhisto;j++)
	{
		(h+j)->histo = (double *) calloc(nbin,sizeof(double));
		if (!(h+j)->histo) Error("Cannot allocate histo_s");
		for (i=0;i<nbin;i++) ((h+j)->histo)[i] = 0;
	}

	return h;
}


/*****************************************************************************
 Rebin histogram, copying it to another one
 it frees the old histogram (!!!)
 *****************************************************************************/
struct histo_s *RebinHistogram(struct histo_s *old, struct param_s *p)
{
  double ebintot=0,emintot=9999.,emaxtot=-9999.,e;
  struct histo_s *new;
  int j,i,jnew,nbintot;

  if (p->debug>0) fprintf(stderr,"\nRebin histograms\n");

  // finds global maximum, minimum, ebin
  for (i=0;i<p->ntemp;i++)
  {
	  // find new min, max, bin
	    if ( (old+i)->emin < emintot ) emintot = (old+i)->emin;		// min of mins
	    if ( (old+i)->emax > emaxtot ) emaxtot = (old+i)->emax;		// max of maxs
	    if (p->ebin != 0) ebintot = p->ebin;						// set ebin from input...
	    else if ( (old+i)->ebin > ebintot) ebintot = (old+i)->ebin; // ...or choose coarser graining

	    if (p->debug>1) fprintf(stderr," t=%d\tmin=%lf\tmax=%lf\tbin=%lf\tnbin=%d\n",
	    		i,(old+i)->emin,(old+i)->emax,(old+i)->ebin,(old+i)->nbin);
  }

  // force an energy minimum
  if (p->femin==1) emintot = p->emin;

  // allocate new structures
  nbintot = (int) ((emaxtot-emintot)/ebintot) + 1;
  if (p->debug>0) fprintf(stderr," global\tmin=%lf\tmax=%lf\tbin=%lf\tnbin=%d\n",
											  emintot,emaxtot,ebintot,nbintot);
  new = AlloHisto(p->ntemp,nbintot+1);

  // copy histograms
  for (i=0;i<p->ntemp;i++)
   for (j=0;j<(old+i)->nbin;j++)					// run over old binning
   {
	 e = (old+i)->emin + (double) (old+i)->ebin * j;
	 jnew = (int)(( e - emintot ) / ebintot + 0.5);  									// new bin
	 if (jnew>=0 && jnew<=nbintot)
		 ((new+i)->histo)[jnew] += ((old+i)->histo)[j];							 // fill bin
   }

  // now all histograms contain the same ebin, emax, etc...
  for (i=0;i<p->ntemp;i++)
  {
	  (new+i)->ebin = ebintot;
	  (new+i)->emin = emintot;
	  (new+i)->emax = emaxtot;
	  (new+i)->nbin = nbintot;

	  if (p->debug>2)
	  {
		  fprintf(stderr,"HISTOGRAM %d:\n",i);
		  for (j=0;j<nbintot;j++)
		  {
			  fprintf(stderr," %d\t%lf\t%lf\t|",j,emintot+ebintot*j,((new+i)->histo[j]));
			  //for (k=0;k<(int)(((new+i)->histo[j])*20.);k++) fprintf(stderr,"#");
			  fprintf(stderr,"\n");
		  }
	  }
  }

  free(old);
  return new;
}

/*****************************************************************************
 Reorder histograms starting from highest temperatures
 returns the permutation vector *n; then one uses n[itemp]
 *****************************************************************************/
void OrderTemperatures(struct histo_s *h, struct param_s *p, int *n)
{
  int i,j,d;

  for (i=0;i<p->ntemp;i++) n[i]=i;

  for (i=0;i<p->ntemp-1;i++)
    for (j=i+1;j<p->ntemp;j++)
      if ( p->temp[n[i]] < p->temp[n[j]] )
	  {
		 d = n[i];
		 n[i] = n[j];
		 n[j] = d;
	  }

  if (p->debug>0) fprintf(stderr,"\nReoredering temperatures\n new order: ");
  if (p->debug>0)
	  for (i=0;i<p->ntemp;i++) fprintf(stderr,"%lf ",p->temp[n[i]]);
  if (p->debug>0) fprintf(stderr,"\n");
}

/*****************************************************************************
 Put to zero bins below a threshold
 *****************************************************************************/
void Threshold(struct histo_s *h, struct param_s *p)
{
   int i,j;
   double z=0.;

   for (i=0;i<p->ntemp;i++)
   {

	   // put to zero bins below threshold
	   for (j=0;j<(h+i)->nbin;j++)
	   {
		   if ( ((h+i)->histo)[j] < p->thresh  ) ((h+i)->histo)[j] =0;
		   else z += ((h+i)->histo)[j];
	   }

	   if (z<EPSILON) fprintf(stderr,"\nWARNING: histogram at T=%lf contains no data after thresholding\n\n",p->temp[i]);

   }
}


int irand(int r)
{
  int i;
  double q;


  q=(double) rand()/(RAND_MAX+1.0);
  q *= r;
  i= (int) q;
  return i;
}

double frand()
{
  double q;

  q=(double) rand()/(RAND_MAX+1.0);
  return q;
}
