/*
 * mhistogram.c
 *
 *  Created on: Nov 6, 2009
 *      Author: guido
 */

#include "mhistogram.h"

int main(int argc, char *argv[])
{
	struct param_s *p;
    int i,j,ndata;
    int *n,*binok;
    double *lg, **hh, *t, *data_tmp, **therm;
    struct histo_s *h,*h_tmp;
    FILE *fout,*ff=NULL;

	// Read parameters
	p = (struct param_s *) calloc(1,sizeof(struct param_s));
	if (argc != 2 || !strcmp(argv[1],"-h") || !strcmp(argv[1],"--help")) Help(stderr);
	Parser(argv[1],p);

	// Allocate stuff
	h_tmp = AlloHisto(p->ntemp,NBINMAX);

	// Do a test of the program
	if (p->test) Test(p);

	// Read files
	if (!p->nohisto)							// if you already have histograms
	   for (i=0;i<p->ntemp;i++)
        ReadHistogram(h_tmp,i,p);
	else										// if you have to make histograms from raw data
	{
		data_tmp = (double *) calloc(NDATAMAX,sizeof(double));
		for (i=0;i<p->ntemp;i++)
		 {
			ndata = ReadData(p,data_tmp,i);
			MakeHistogram(h_tmp,data_tmp,i,ndata,p);
		 }
		free(data_tmp);
	}

	// Rebin histograms (the new histogram is h, h_tmp is freed)
	h = RebinHistogram(h_tmp,p);

	// Threshold the histograms
	if (p->thresh>0)
		Threshold(h,p);

	// Print dumb averages obtained summing the histograms
	if (strcmp(p->nfene,""))
		{
			fout = fopen(p->nfene,"w");
			PrintAverageEnergies(fout,h,p->temp,p->ntemp);
			fclose(fout);
		}

	// Prepare data for multiple histogram
	n = AlloInt(p->ntemp);									// vector for permutation of temperatures
	OrderTemperatures(h,p,n);								// order temperatures from highest
	hh = AlloDoubleMat(p->ntemp,h->nbin);					// MHistogram wants double **
	for (i=0;i<p->ntemp;i++)								// copy histogram to double ** and reorder temperatures
		for (j=0;j<h->nbin;j++)
			hh[i][j] = ((h+n[i])->histo)[j];
	t = AlloDouble(p->ntemp);								// temperature vector
	for (i=0;i<p->ntemp;i++)
		t[n[i]] = p->temp[i];
	binok = AlloInt(h->nbin);

	// print out -log(Z)
	if (strcmp(p->nfz,""))
	{
		if (p->debug>0) fprintf(stderr,"Wants to print -log(Z) to %s.\n",p->nfz);
		ff = fopen(p->nfz,"w");
		if (!ff) FatalError("Cannot open file to write -log(Z)");
	}

	// Do multiple histogram
	lg = MHistogram(hh,t,p->ntemp,h->nbin,h->emin,h->ebin,p->debug,p->kb,NULL,NULL,p->paranoid,p->ignoreb,ff,p->deltat,binok);

	// Calculate thermodynamics
	therm = AlloDoubleMat(4,p->ntout);
	CalculateThermodynamics(lg,h->emin,h->ebin,h->nbin,p->tminout,p->tbinout,p->ntout,therm,p->kb,NULL,NULL,binok);

	// Print out thermodynamics
	if (!strcmp(p->nfout,"STDOUT"))
		fout = stdout;
	else
		fout = fopen(p->nfout,"w");
	if (p->debug>0) fprintf(stderr,"Printing thermodynamics to %s [columns: E <E>  sigma_E  Cv F]\n",p->nfout);
	PrintThermodynamics(fout,therm,p->tminout,p->tbinout,p->ntout);
	if (strcmp(p->nfout,"STDOUT")) fclose(fout);

	// Print out the density of states
	if (strcmp(p->nfg,""))
	{
		if (p->debug>0) fprintf(stderr,"Printing density of state to %s [columns: E log(g(E))]\n",p->nfg);
		fout = fopen(p->nfg,"w");
		for (i=0;i<h->nbin;i++)
				fprintf(fout,"%lf\t%lf\n",h->emin+h->ebin*i,lg[i]);
		fclose(fout);
	}

        fprintf(stderr,"Please cite: G. Tiana and L. Sutto, Equilibrium properties of realistic random heteropolymers and their relevance for globular and naturally unfolded proteins Phys. Rev. E 84 (2011) 061910\n");

	exit(0);
}

void Parser(char *fname, struct param_s *p)
{
	char aux[500];
	int i;
	FILE *fin;
	strcpy(p->title,"");
	strcpy(p->nfout,"STDOUT");
	strcpy(p->nfene,"");
	strcpy(p->nfene,"dumb_e.dat");
	strcpy(p->nfg,"");
	p->ntemp=0;
	p->thresh=0;
	p->nohisto=0;
	p->ebin=0;
	p->tminout=0;
	p->debug=0;
	p->tbinout=1;
	p->ntout=200;
	p->test=0;
	p->ncol2[0]=0;
	p->kb=1;
	p->paranoid=0;
	p->ignoreb=0;
	p->deltat=0;
	p->emin = 0;
	p->femin = 0;

	fin = fopen(fname,"r");
	if (!fin) Error("Input file missing.");

	while(fgets(aux,500,fin)!=NULL)
    {
       if (!strncmp(aux,"title",5))
        if ( !sscanf(aux,"title %s",p->title) ) Error("Cannot read title in inputfile");

       if (!strncmp(aux,"ntemp",5))
        if ( !sscanf(aux,"ntemp %d",&p->ntemp) ) Error("Cannot read ntemp in inputfile");

       if (!strncmp(aux,"ebin",4))
         if ( !sscanf(aux,"ebin %lf",&p->ebin) ) Error("Cannot read ebin in inputfile");

       if (!strncmp(aux,"emin",4))
       {
         if ( !sscanf(aux,"emin %lf",&p->emin) ) Error("Cannot read ebin in inputfile");
         p->femin=1;
       }
       if (!strncmp(aux,"nohisto",7)) p->nohisto=1;

       if (!strncmp(aux,"ignoreb",7))
    	   if ( !sscanf(aux,"ignoreb %d",&p->ignoreb) ) Error("Cannot read ignoreb in inputfile");

       if (!strncmp(aux,"files",5))
       {
    	   if (!p->ntemp) Error("ntemp must be defined before the list of files.");
    	   for (i=0;i<p->ntemp;i++)
    	     if ( fscanf(fin,"%s %lf %d %d",p->nfile[i],&p->temp[i],&p->ncol1[i],&p->ncol2[i]) < 3)
    	    	 {
					 fprintf(stderr,"file=%s temp=%lf col1=%d col2=%d\n",p->nfile[i],p->temp[i],p->ncol1[i],p->ncol2[i]);
					 Error("Cannot read file in inputfile");
    	    	 }
       }

       if (!strncmp(aux,"threshold",9))
               if ( !sscanf(aux,"threshold %lf",&p->thresh) ) Error("Cannot read threshold in inputfile");

       if (!strncmp(aux,"debug",5))
           if (sscanf(aux,"debug %d",&p->debug)!=1) p->debug = 1;

       if (!strncmp(aux,"outfile",7))
               if ( !sscanf(aux,"outfile %s",p->nfout) ) Error("Cannot read outfile in inputfile");

       if (!strncmp(aux,"enefile",7))
               if ( !sscanf(aux,"enefile %s",p->nfene) ) Error("Cannot read enefile in inputfile");

       if (!strncmp(aux,"gfile",5))
                     if ( !sscanf(aux,"gfile %s",p->nfg) ) Error("Cannot read gfile in inputfile");

       if (!strncmp(aux,"zfile",5))
                      if ( !sscanf(aux,"zfile %s",p->nfz) ) Error("Cannot read zfile in inputfile");

       if (!strncmp(aux,"tmin",4))
          if ( !sscanf(aux,"tmin %lf",&p->tminout) ) Error("Cannot read tmin in inputfile");

       if (!strncmp(aux,"tbin",4))
           if ( !sscanf(aux,"tbin %lf",&p->tbinout) ) Error("Cannot read tmin in inputfile");

       if (!strncmp(aux,"ntbin",4))
           if ( !sscanf(aux,"ntbin %d",&p->ntout) ) Error("Cannot read ntbin in inputfile");

       if (!strncmp(aux,"kb",2))
            if ( !sscanf(aux,"kb %lf",&p->kb) ) Error("Cannot read kb in inputfile");

       if (!strncmp(aux,"test",4)) p->test=1;

       if (!strncmp(aux,"paranoid",8)) p->paranoid=1;

       if (!strncmp(aux,"deltat",6))
              if ( !sscanf(aux,"deltat %lf",&p->deltat) ) Error("Cannot read deltat in inputfile");

    }

	if (p->ntemp==0) Error("Number of temperatures is set to zero");
	if (p->nohisto==1 && p->ncol2[0]!=0) Error("You specified two columns to read row data. You need only one.");

	if (p->debug>0) fprintf(stderr,"Debug level =\t\t%d\n",p->debug);
	if (p->debug>0)
	{
		fprintf(stderr,"INPUT:\n");
		fprintf(stderr,"title = \t%s\n",p->title);
		fprintf(stderr,"ntemp = \t%d\n",p->ntemp);
		fprintf(stderr,"ebin = \t\t%lf\n",p->ebin);
		if (p->femin==1) fprintf(stderr,"emin = \t\t%lf\n",p->emin);
		fprintf(stderr,"thresh = \t%lf\n",p->thresh);
		fprintf(stderr,"tmin = \t\t%lf\n",p->tminout);
		fprintf(stderr,"tbin = \t\t%lf\n",p->tbinout);
		fprintf(stderr,"#tbin = \t%d\n",p->ntout);
		fprintf(stderr,"Kb = \t\t%lf\n",p->kb);
		if (p->ignoreb) fprintf(stderr,"Ignore temperatures which prevent equation solving\n");
		fprintf(stderr,"\n");
	}
}

void Help(FILE *fp)
{
   fprintf(fp,"\n* usage: mhistogram <inputfile>\n");
   fprintf(fp,"G. Tiana, 2009\n");
   fprintf(fp,"Please cite: G. Tiana and L. Sutto, Equilibrium properties of realistic random heteropolymers and their relevance for globular and naturally unfolded proteins Phys. Rev. E 84 (2011) 061910\n");
   fprintf(fp,"NEEDED INPUT:\n");
   fprintf(fp," ntemp <INT>\t\t\tnumber of histograms to process.\n");
   fprintf(fp," files\n");
   fprintf(fp," <filename> <FLOAT> <INT> <INT>\tthe file containing the histogram, the temperature,\n");
   fprintf(fp,"\t\t\t\tthe columns of the containing the energy\n");
   fprintf(fp," ..or..\t\t\t\tand the counts.\n");
   fprintf(fp," <filename> <FLOAT> <INT>\tthe file containg the energy series, the temperature\n");
   fprintf(fp,"\t\t\t\tand the column containing the energy.\n\n");
   fprintf(fp,"OPTIONS:\n");
   fprintf(fp," title\t\t\t\ta name for the data set\n");
   fprintf(fp," nohisto\t\t\treads energy series instead of histograms\n");
   fprintf(fp," ebin <FLOAT>\t\t\tthe bin size for energies\n");
   fprintf(fp," emin <FLOAT>\t\t\tforce a specified energy minimum for the histograms\n");
   fprintf(fp," threshold <FLOAT>\t\tthreshold on the (normalized) count probabilit\n");
   fprintf(fp,"\t\t\t\tto neglect bins\n");
   fprintf(fp," debug <INT>\t\t\t0 = be silent\n");
   fprintf(fp,"\t\t\t\t1 = say what it is doing\n");
   fprintf(fp,"\t\t\t\t2 = give the details\n");
   fprintf(fp,"\t\t\t\t3 = print every cough\n");
   fprintf(fp," outfile <filename>\t\toutput file for canonical ensemble data (T, <E>, Cv, F)\n");
   fprintf(fp," enefile <filename>\t\toutput file for averages calculated dumbly from each T\n");
   fprintf(fp," gfile <filename>\t\toutput file for density of states\n");
   fprintf(fp," zfile <filename>\t\toutput file for -log(Z)\n");
   fprintf(fp," tmin <float>\t\t\ttemperature to start writing canonical ensemble data\n");
   fprintf(fp," tbin <float>\t\t\ttemperature bin for canonical ensemble data\n");
   fprintf(fp," kb <float>\t\t\tBoltzmann's constant (default 1)\n");
   fprintf(fp," ntbin <int>\t\t\t# of temperature bin for canonical ensemble data\n");
   fprintf(fp," test\t\t\t\ttest the algorithm with a random energy model\n");
   fprintf(fp," paranoid\t\t\tstops instead of give warning\n");
   fprintf(fp," ignoreb {0,1,2}\t\tignore temperatures which prevent equation solving\n");
   fprintf(fp," \t\t\t0=off, 1=equations have to converge (tougher), 2=a local minimum is enough\n");
   fprintf(fp," deltat <FLOAT>\t\tignore temperatures which are closer than it\n");

   exit(0);
}
