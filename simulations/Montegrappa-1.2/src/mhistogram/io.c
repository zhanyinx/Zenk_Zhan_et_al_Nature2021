/*
 * io.c
 *
 *  Created on: Nov 6, 2009
 *      Author: guido
 */

#include "mhistogram.h"

void Error(char *text)
{
  time_t rawtime;

  time(&rawtime);
  fprintf(stderr,"\n* Error:\n");
  fprintf(stderr,"  %s\n",text);
  fprintf(stdout,"  after %3.1fs at %s\n",(double) ( clock() / CLOCKS_PER_SEC),
                                      asctime(localtime(&rawtime)));
  exit(1);
}

void ReadHistogram(struct histo_s *h, int itemp, struct param_s *p)
{
	int i,k,kmin=0,kmax=0;
	char form[5000],aux[5000];
	double daux1=0,daux2=0,emin=9999.,emax=-9999.,ebin,z=0,daux1old,histotmp[NBINMAX];
	FILE *fp;

	if (p->debug>0)
		fprintf(stderr,"Reading histogram from %s (columns %d, %d)\n",p->nfile[itemp],p->ncol1[itemp],p->ncol2[itemp]);

	// open file
	fp = fopen(p->nfile[itemp],"r");
	if (!fp) Error("Cannot open file for reading");

	// which columns
	strcpy(form,"");
	for (i=0;i<p->ncol1[itemp]-1;i++)
	  sprintf(form,"%s%%*s ",form);
	sprintf(form,"%s%%lf ",form);
	for (i=p->ncol1[itemp];i<p->ncol2[itemp]-1;i++)
	  sprintf(form,"%s%%*s ",form);
	sprintf(form,"%s%%lf",form);
	if (p->debug>0)
		fprintf(stderr," debug: format %s\n",form);

	// read file; finds max, min, deltae and normalization constant
	k=0;
	ebin=daux1old=0;
	while(fgets(aux,5000,fp)!=NULL)
		if (sscanf(aux,form,&daux1,&daux2)==2)
		{
			histotmp[k] = daux2;
	
											// record histogram
			if (k==1) ebin = daux1 - daux1old;								// find bin
			else if (k>1 && fabs(daux1-daux1old-ebin)>EPSILON)
				Error("Inhomogeneous binning in histogram");
			if (daux2>0 && daux1<emin) { emin=daux1; kmin=k; }				// find min
			if (daux2>0 && daux1>emax) { emax=daux1; kmax=k; }				// find max
			z += daux2;														// find normaliz. constant
			daux1old = daux1;												// to determine bin
			k++;
			if (k>=NBINMAX) Error("NBINMAX too small.");
		}

	if (z<EPSILON) Error("Empty histogram");

	for (i=kmin;i<=kmax;i++)
		((h+itemp)->histo)[i-kmin] = histotmp[i];

	if (p->debug>0)
	{
		fprintf(stderr," %d bins read.\n",kmax-kmin+1);
		fprintf(stderr," emin=%lf emax=%lf ebin=%lf\n",emin,emax,ebin);
	}

	(h+itemp)->nbin = kmax-kmin+1;
	(h+itemp)->emin = emin;
	(h+itemp)->emax = emax;
	(h+itemp)->ebin = ebin;

}
/*****************************************************************************
 Read directly data instead of histograms
 *****************************************************************************/
int ReadData(struct param_s *p, double *d, int itemp)
{
	int k=0,i;
	char form[1000],aux[5000];
	double daux;
	FILE *fp;

	if (p->debug>0) fprintf(stderr,"Read row data from file %s\n",p->nfile[itemp]);

	// open file
	fp = fopen(p->nfile[itemp],"r");
	if (!fp) Error("Cannot open file for reading");

	// format to be read
	strcpy(form,"");
	for (i=0;i<p->ncol1[itemp]-1;i++)
	  sprintf(form,"%s%%*s ",form);
	sprintf(form,"%s%%lf",form);

	if (p->debug>0)
		fprintf(stderr,"debug: format %s\n",form);

	// read file
	while(fgets(aux,5000,fp)!=NULL)
		if (sscanf(aux,form,&daux)==1)
		{
			d[k] = daux;
			k++;
			if (k>=NDATAMAX) Error("NDATAMAX too small");
		}

	if (p->debug>0)
		fprintf(stderr," %d data read.\n",k);

	return k;
}

/*****************************************************************************
 Make histogram out of row data
 *****************************************************************************/
void MakeHistogram(struct histo_s *h, double *d, int itemp, int ndata, struct param_s *p)
{
	int i,ibin,z=0;
	double min=999999.9,max=-9999999.9,ebin=1.;

	// find min and max
	for (i=0;i<ndata;i++)
	{
		if (i==0 || d[i]>max) max=d[i];
		if (i==0 || d[i]<min) min=d[i];
	}

	// set bin
	if (p->ebin != 0) ebin = p->ebin;		// given as input
	else ebin = (max-min)/100.;				// otherwise default

	// bin data
	for (i=0;i<ndata;i++)
	{
		ibin = (int) ( (double) (d[i]-min)/ebin + 0.5 );
		if (ibin>NBINMAX)
			{
				fprintf(stderr,"Try to bin energy %lf (line %d; emin=%lf bin=%d)\n",d[i],i,min,ibin);
				Error("NBINMAX too small");
			}
		(h+itemp)->histo[ibin] ++;
		z++;
	}


	(h+itemp)->emin = min;
	(h+itemp)->emax = max;
	(h+itemp)->ebin = ebin;
	(h+itemp)->nbin = (int)((max-min)/ebin) + 1;

	if (p->debug>0)
		fprintf(stderr,"Make histogram %d:\n min=%lf\tmax=%lf\tbin=%lf\t(nbin=%d)\n",itemp,min,max,ebin,(h+itemp)->nbin);

}

/*****************************************************************************
 Calculates the average energy from a straightforward average of the histogram
 *****************************************************************************/
void PrintAverageEnergies(FILE *fout, struct histo_s *h, double *t, int ntemp)
{
	int it,ie;
	double em,z;

	for (it=0;it<ntemp;it++)
	{
		em = 0;
		z = 0;

		for (ie=0;ie<(h+it)->nbin;ie++)
		{
			em += ( (h+it)->emin + (double) (h+it)->ebin * ie ) * (h+it)->histo[ie];
			z += (h+it)->histo[ie];
		}

		fprintf(fout,"%lf\t%lf\n",t[it],em/z);
	}
}

/*****************************************************************************
 Carries out a test of the algorithm, ignoring the actual data and creating
 a set of test files, containing the p(E) of a random energy model.
 It needs the temperatures (please specify some dumb file names) and ebin.
 *****************************************************************************/
void Test(struct param_s *p)
{
	int it,ie;
	double h,e,em,z;
	char nfile[20];
	FILE *fout;
	double ebin=1.;
	double media=0.;
	double sigma=20.;

	printf("Average energies of the test:\n\tT\t<E>\n");

	for (it=0;it<p->ntemp;it++)
	{
		sprintf(nfile,"test_%d.dat",it);
		fout = fopen(nfile,"w");
		em=0; z=0;
		for (ie=-100;ie<100;ie++)
		{
			// use as test the random energy model with some parameters....
			e = ebin * ie;
			h = exp(-(e-media)*(e-media)/(2.*sigma*sigma)) * exp(-e / p->temp[it]);
			fprintf(fout,"%lf\t%lf\n",p->ebin * ie,h);

			em += e*h;
			z += h;
		}
		fclose(fout);
		strcpy(p->nfile[it],nfile);
		p->ncol1[it]=1;
		p->ncol2[it]=2;

		printf("\t%lf\t\t%lf\n",p->temp[it],em/z);
	}

}
