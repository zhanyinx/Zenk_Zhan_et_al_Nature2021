/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

 
    Copyright (C) 2014  G. Tiana, F. Villa, Y. Zhan, R. Capelli, C.Paissoni, P. Sormanni, R. Meloni

    MonteGrappa is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.


    This program is free. We just ask you to cite in the publications you get 
    making use of    MonteGrappa the following article:

    G. Tiana, F. Villa, Y. Zhan, R. Capelli, C.Paissoni, P. Sormanni, E. Heard,
    L. Giorgetti and R. Meloni, "MonteGrappa: an iterative Monte Carlo program
    to optimize biomolecular potentials in simplified models" Comp. Phys. Comm.
    Volume 186, January 2015, Pages 93–104.     DOI: 10.1016/j.cpc.2014.09.012

 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/*
 * do_mhistogram.c
 *
 *  Created on: Nov 10, 2009
 *      Author: guido
 *
 *  it needs also memory.c
 */

#include "do_mhistogram.h"

/*****************************************************************************
 Returns the logarithm of density of states from the histograms
	 h[temp][energy] is the normalized histograms ,
	 t[temp] is the list of temperatures
	 ntemp is the number of temperatures
	 nbin is the number of energy bins
	 emin is the energy corresponding to the first bin
	 ebin is the size of an energy bin
	 debug goes from 0 (silent) to 4 (every cough)
	 ignoreb !=0 means that it ignores temperatures which prevent equation solving (1=has to solve seriously, 2=it's sufficient to get a local minimum
	 f is for internal use
	 if it cannot solve the equations, returns NULL
 *****************************************************************************/
double *MHistogram(double **h, double *t, int ntemp, int nbin, double emin, double ebin, int debug, double kb,
		double *lg, double *f, int paranoid, int ignoreb, FILE *ff, double deltaT, int *binok)
{
	double e,fmax,z,good_f[NTEMPMAX],lastT=99999.;
    int i,ie,n,j,chk,r,useit[NTEMPMAX];
    FILE *fh;

    if (debug>0) fprintf(stderr,"Performing MHistogram\n");
    if (debug>1) fprintf(stderr," ntemp =\t%d\n nbin =\t%d\n emin =\t%lf\n ebin =\t%lf kb=%lf\n",ntemp,nbin,emin,ebin,kb);
    if (debug>1) fprintf(stderr," EPSILON=%e MHE_RESIDUAL=%e MHE_RESIDUAL2=%e NITERHIS=%d\n",EPSILON,MHE_RESIDUAL,MHE_RESIDUAL2,NITERHIS);
    if (debug>1) fprintf(stderr," NTEMPMAX=%d NEBINMAX=%d\n",NTEMPMAX,NEBINMAX);
    if (debug>1) { fh = fopen("histocheck.dat","w"); for (ie=0;ie<nbin;ie++) { fprintf(fh,"%d\t",ie); for (i=0;i<ntemp;i++) fprintf(fh,"%lf ",h[i][ie]);
    				fprintf(fh,"\n"); } fclose(fh); }
	if (debug>1) { for (j=0;j<ntemp;j++) fprintf(stderr," T[%d] = %lf  ",j,t[j]); fprintf(stderr,"\n"); }

	for (ie=0;ie<nbin;ie++) binok[ie] = 0;			// which bin to use (set by DensityofStates

	// Check if  there are empty histograms
	for (i=0;i<ntemp;i++)
	{
		z = 0;
		useit[i]=1;

		// if option is active, discard histograms which are too close in temperature
		if (deltaT>0)
		{
			if (lastT-t[i] < deltaT)
			{
				useit[i] = 0;
				if (debug>0) fprintf(stderr,"Discard T=%lf because DT=%lf\n",t[i],lastT-t[i]);
			}
			else lastT = t[i];
		}

		// check histograms
		for (ie=0;ie<nbin;ie++)  z += h[i][ie];
		if (z<EPSILON)
		{
			fprintf(stderr,"Histogram %d (T=%lf) is empty\n",i,t[i]);
			//MODIFICHE
			//FatalError("Empty histogram in MHistogram");
			//FINEMODIFICHE
			fprintf(stderr,"WARNING: cannot print thermodynamics, there is empty Histogram!");
		}
	}

	// Allocate stuff
	if (lg==NULL) lg = AlloDouble(nbin);				// logarithm of density of states
    if (f==NULL)  f = AlloDouble(ntemp);				// f=-log(Z)

    // One of the partition functions is arbitrary
    f[0] = 0;

	// Solve iteratively the MH equations, starting from the highest temperature
	for (i=1;i<ntemp;i++)					// Start with 2 equations
		if (useit[i]==1)
		{
			if (debug>2) {
				n=0;
				for (j=0;j<=i;j++) n+=useit[j];
				fprintf(stderr," Solve the %d equations for :\n",n);
				for (j=0;j<=i;j++) if (useit[j]==1) fprintf(stderr," T[%d] = %lf  ",j,t[j]);
				fprintf(stderr,"\n"); 
			}

			// Initial guess for i-th temperature
			// f_i = -log(\sum_E exp[-E/T + log(g(E)) - fmax]) - fmax

			DensityOfStates(h,t,f,emin,ebin,i,nbin,lg,debug,kb,useit,binok);	// g obtained from higher temperatures

			fmax = -9E19;
			chk = 0;
			for (ie=0;ie<nbin;ie++)								// find maximum exponential to avoid overflows
			{
				e = emin + (double) ebin * ie;
				if ( (-e/(kb*t[i]) + lg[ie] > fmax) && binok[ie] == 1 ) {fmax = -e/(kb*t[i]) + lg[ie]; chk = 1;}
			}

			if (chk == 0) FatalError("Cannot find fmax in MHistogram");

			z = 0;
			for (ie=0;ie<nbin;ie++)								// sum partition function
			{
				e = emin + (double) ebin * ie;
				if ( binok[ie] == 1 )
					z += exp(-e/(kb*t[i]) + lg[ie] -fmax);		// shift argument by fmax to prevent overflows
			}

			if (z<EPSILON)
			{
					fprintf(stderr,"Empty histogram in MHistogram when calculating initial guess for T=%lf",t[i]);
					exit(1);
			}

			f[i] = -log(z) - fmax;
			if (debug>2) fprintf(stderr," Initial guess: f[%d]=%lf (fmax=%lf)\n",i,f[i],fmax);

			// Solve the equation up to i-th temperature
			r = SolveEquations(h,t,f,emin,ebin,i+1,nbin,debug,kb,paranoid,useit,binok);
			if (debug>2) fprintf(stderr," Result :  f[%d]=%lf fmax=%e\n",i,f[i],fmax);


			// If option active, if cannot solve equation, redo without last temperature
			if (ignoreb!=0)
			{
				if ( (r != 1 && ignoreb==1) || (r==-1 && ignoreb==2))		// solver failed
				{
					if (debug>4) {
						//STOP RUN AT FIRST EQUATIONS NOT SOLVED
						fprintf(stderr,"\n\nntemp=%d\nnbin=%d\nemin=%lf\nebin=%lf\ndebug=%d\nkb=%lf\nignoreb=%d\ndeltaT=%lf\n",ntemp,nbin,emin,ebin,debug,kb,ignoreb,deltaT);
					//	PrintHistogram("dbg_histo.dat",h,t,ntemp,nbin,emin,ebin);
						fflush(NULL);
						exit(1);	
					}
					useit[i]=0;									// switch off current temperature
					for (j=0;j<ntemp;j++) f[j] = good_f[j];		// return to reliable f
					if (debug>1) fprintf(stderr,"MHistogram: cannot solve equation with T[%d]=%lf, switch it off.\n",i,t[i]);
					i--;										// repeat without it
				}
				else
					for (j=0;j<ntemp;j++) good_f[j] = f[j];		// record working f
			}
		}


	if (debug>0)
		if (ignoreb>0)
			{
				j=0;
				for (i=0;i<ntemp;i++) j += useit[i];
				fprintf(stderr,"Used %d/%d histograms in MHistogram",j,ntemp);
				if (j!=ntemp)
				{
					fprintf(stderr," (excluding ");
					for (i=0;i<ntemp;i++) if (useit[i]==0) fprintf(stderr,"T[%d]=%lf ",i,t[i]);
					fprintf(stderr,")\n");
				}
				else fprintf(stderr,"\n");

				if (j<(int)(ntemp/2)) fprintf(stderr,"WARNING: less than half of the histograms used in MHistograms");
			}

	// calculate the density of states with the f found above
	if (debug>1) fprintf(stderr,"\nCalculate the resulting density of states...\n");
	DensityOfStates(h,t,f,emin,ebin,ntemp,nbin,lg,debug,kb,useit,binok);

	if (ff != NULL)
		for (i=0;i<ntemp;i++)
			fprintf(ff,"%lf\t%lf\n",t[i],f[i]);



	return lg;
}

/*****************************************************************************
 Finds the logarithm of the density of states g from the free energies f
 it uses the highest ntemp temperatures
 *****************************************************************************/
void DensityOfStates(double **h, double *t, double *f, double emin, double ebin, int ntemp, int nbin, double *lg,
		int debug, double kb, int *useit, int *binok)
{
	int i,ie,chk;
	double num,den,e,fmax,n[NTEMPMAX];

	// Equation A.6 of Jesper's thesis:
	//
	//                              \sum_T h(T)
	//  g(E) = ---------------------------------------------------------------
	//          exp[Fmax] \sum_T ( chi_T(E) n_T exp[-E/T - log(Z_T) - Fmax] )
	//
	// Fmax is such that, for each energy E, the maximum argument of the exponential
	// at the denominator is 0.
	// f(T)=-log(Z)

	if (debug>2) fprintf(stderr," Density of states\n");

	// set binok
	for (ie=0;ie<nbin;ie++) binok[ie] = 0;
	for (i=0;i<ntemp;i++)
		for (ie=0;ie<nbin;ie++)
			if (h[i][ie]>0.) binok[ie] = 1;


	// normalization constant
	for (i=0;i<ntemp;i++)
		if (useit[i]==1)
		{
			n[i] = 0;
			for (ie=0;ie<nbin;ie++)
				if (binok[ie] == 1) n[i] += h[i][ie];
		}

	// calculate density of states
	chk=0;
	for (ie=0;ie<nbin;ie++)					// each energy bin is calculated separately
		if (binok[ie]==1)
		{
			num=0.;
			den=0.;

			e = emin + (double) ebin * ie;		// energy from energy bin

			// find the maximum exponent in the denominator, to prevent overflows
			fmax = -9E19;
			for (i=0;i<ntemp;i++)
				if (useit[i]==1)
						if ( (h[i][ie]>0) && (-e/(kb*t[i])+f[i]>fmax) ) fmax = -e/(kb*t[i])+f[i];

			// sum numerator and denominator
			for (i=0;i<ntemp;i++)
				if (useit[i]==1)
				{
					num += h[i][ie];								// numerator
					if (h[i][ie]>0.)								// sum only on the support of the histogram (Jesper's chi)
						den +=  n[i] * exp(-e/(kb*t[i])+f[i]-fmax);	// sum the denominator, shifting the exponent to prevent overflows
				}

			if (den>0. && num>0.)
			{
				lg[ie] = -fmax + log(num) - log(den);
				chk=1;
			}
			else lg[ie] = -9E19;

			if (debug>2 && binok[ie]==1) fprintf(stderr," %d\te=%lf\tlog(num)=%e\tlog(den)=%e\tfmax=%lf\tlog(g)=%lf\n",
					ie,e,log(num),log(den),fmax,lg[ie]);
		}

	if (chk==0)
	{
		for (i=0;i<ntemp;i++)
			if (binok[ie]==1)
				for (ie=0;ie<nbin;ie++) fprintf(stderr,"it=%d\tie=%d\th=%lf\n",i,ie,h[i][ie]);
		fprintf(stderr,"All histograms empty in DensityofStates\n\n");
		exit(1);
	}
}

/*****************************************************************************
 Calculates average energy, fluctuations, Cv
	 input:	lg=logarithm of density of states; emin, ebin, nbin= binning of energy
			nt = number of temperatures to output
			temp = array of temperatures
			if temp==NULL, scan all temperatures using following parameters (otherwise useless)
			tmin, tbin = binning of temperatures to output

	output: out[0][it] = average energy			(if out==NULL do not calculate)
			out[1][it] = standard deviation of energy
			out[2][it] = Cv (if KB=1 it is a pure number)
			out[3][it] = F
			prob[it][ie] = Boltzman's probability (if prob==NULL do not caluclate)
 *****************************************************************************/
void CalculateThermodynamics(double *lg, double emin, double ebin, int nbin,
			double tmin, double tbin, int nt, double **out, double kb, double *temp, double **prob, int *binok)
{
	int ie,it,chk;
	double em,em2,z,e,kt,expmax;

	// <E> = \sum_E E g(E) exp[-(E-Emin)/T] / \sum_E g(E) exp[-(E-Emin)/T]

	for (it=0;it<nt;it++)						// loop over temperatures
	{
		em = 0;
		em2 = 0;
		z = 0;
		expmax = -9E19;
		chk = 0;

		if (temp==NULL)								// which temperature
			kt = kb*(tmin + (double) it * tbin);	// scan all temperatures...
		else
			kt = kb*temp[it];						// ...or take from temp

		for (ie=0;ie<nbin;ie++)						// find maxium exponent to prevent overflows
			if ( binok[ie] == 1)
			{
				e = emin + (double) ebin * ie;		// energy from energy bin
				if (-e/kt + lg[ie] > expmax) {expmax = -e/kt + lg[ie]; chk = 1;}
			}
		if (chk==0) {
			fprintf(stderr,"expmax not found. aborting. T[%d]=%lf\n",it,kt);
			exit(1);
		}

		for (ie=0;ie<nbin;ie++)					// loop over energies
		{
			if (binok[ie] == 1)
			{
				e = emin + (double) ebin * ie;		// energy from energy bin
				if (prob != NULL) prob[it][ie] = exp(-e/kt + lg[ie] - expmax);
				if (out != NULL)
				{
					em += e * exp(-e/kt + lg[ie] - expmax);
					em2 += e * e * exp(-e/kt + lg[ie] - expmax);
				}
				z += exp(-e/kt + lg[ie] - expmax);
			}
			else if (prob != NULL) prob[it][ie]=0.;
		}

		if (out != NULL)
		{
			out[0][it] = em / z;								// average
			out[1][it] = sqrt(em2/z - out[0][it]*out[0][it]);	// standard deviation
			out[2][it] = kb * out[1][it] * out[1][it] / kt / kt;				// Cv
			out[3][it] = -kt * log(z) -kt * expmax ;
		}

		if (prob != NULL)
			for (ie=0;ie<nbin;ie++) prob[it][ie] /= z;

	}
}

void PrintThermodynamics(FILE *fout, double **in, double tmin, double tbin, int nt)
{
	int i;
	double t;

	fprintf(fout,"# TEMP\t<E>\tsigma_E\tCv\tF(T)\tS(T)\n");
	for (i=0;i<nt;i++)
	{
		t = tmin+tbin*i;
		fprintf(fout,"%lf\t\t%lf\t%lf\t%lf\t%lf\t%lf\n",t,in[0][i],in[1][i],in[2][i],in[3][i],(in[0][i]-in[3][i])/t);
	}
}

void FatalError(char *text)
{
  fprintf(stderr,"\n* Fatal Error:\n");
  fprintf(stderr,"  %s\n",text);

  exit(1);
}

double FreeEnergy(double *lg, double emin, double ebin, int nbin, double temp, double kb, int *binok)
{
	int i;
	double z,e,expmax=-9E19;

	for (i=0;i<nbin;i++)							// maximum exponent to prevent overflows
			if (binok[i] == 1)
			{
				e = emin + (double) ebin * i;		// energy from energy bin
				if ( -e/(kb*temp) + lg[i] > expmax ) expmax = -e/(kb*temp) + lg[i];
			}

	z = 0;

	for (i=0;i<nbin;i++)						// loop over energies
		if (binok[i] == 1)
		{
			e = emin + (double) ebin * i;		// energy from energy bin
			z += exp(-e/(kb*temp) + lg[i] - expmax);
		}

	return -kb * temp * log(z) - kb * temp * expmax;
}

double AverageEnergy(double *lg, double emin, double ebin, int nbin, double temp, double kb, int *binok)
{
	int i;
	double em,z,e,expmax=-9E19;

	for (i=0;i<nbin;i++)							// maximum exponent to prevent overflows
			if (binok[i] == 1)
			{
				e = emin + (double) ebin * i;		// energy from energy bin
				if ( -e/(kb*temp) + lg[i] > expmax ) expmax = -e/(kb*temp) + lg[i];
			}

	em = 0;
	z = 0;
	for (i=0;i<nbin;i++)							// loop over energies
			if (binok[i]==1)
			{
				e = emin + (double) ebin * i;		// energy from energy bin
				em += e * exp(-e/(kb*temp) + lg[i] - expmax);
				z += exp(-e/(kb*temp) + lg[i] - expmax);
			}

	return em / z;
}

double SigmaEnergy(double *lg, double emin, double ebin, int nbin, double temp, double kb, int *binok)
{
	int i;
	double em,z,e,em2,expmax=-9E19;

	for (i=0;i<nbin;i++)							// maximum exponent to prevent overflows
			if (binok[i]==1)
			{
				e = emin + (double) ebin * i;		// energy from energy bin
				if ( -e/(kb*temp) + lg[i] > expmax ) expmax = -e/(kb*temp) + lg[i];
			}

	em = 0;
	em2 = 0;
	z = 0;
	for (i=0;i<nbin;i++)							// loop over energies
			if (binok[i]==1)
			{
				e = emin + (double) ebin * i;		// energy from energy bin
				em += e * exp(-e/(kb*temp) + lg[i] - expmax);
				em2 += e * e * exp(-e/(kb*temp) + lg[i] - expmax);
				z += exp(-e/(kb*temp) + lg[i] - expmax);
			}

		return sqrt(em2/z - em*em/z/z);

}

/*****************************************************************************
 Solution algorithms of the GMH equations
 returns 1 if ok, 0 if got stucked, 1 if couln't solve
 *****************************************************************************/

// ALGORITHM WITH ANALYTICAL DERIVATIVE
#ifdef GSL_DER
int SolveEquations(double **h, double *t, double *f, double emin, double ebin, int ntemp, int nbin,
		int debug, double kb, int paranoid, int *useit, int *binok)
{
	int iter=0,status,i,j,actual_ntemp=0;
    const gsl_multiroot_fdfsolver_type *T;
    gsl_multiroot_fdfsolver *s;
    double f0;
    const size_t n = ntemp;
    gsl_vector *ff = gsl_vector_alloc (n);
	double delta = 0;


    for (i=0;i<ntemp;i++)
    	gsl_vector_set (ff, i, f[i]);

    for (i=0;i<ntemp;i++) actual_ntemp += useit[i];		// actual number of used equations (for convergence test)
    if (actual_ntemp==0) FatalError("actual_ntemp=0 in SolveEquations");

    struct rparams p = {h,t,emin,ebin,ntemp,nbin,kb,debug,paranoid,useit,binok};

    gsl_multiroot_function_fdf eqq;
    eqq.f=&GMHequation;
    eqq.df=&GMHequation_df;
    eqq.fdf=&GMHequation_fdf;
    eqq.n=n;
    eqq.params=&p;

    //T = gsl_multiroot_fdfsolver_gnewton;
    T = gsl_multiroot_fdfsolver_hybridj;
    s = gsl_multiroot_fdfsolver_alloc (T, n);
    if (s==NULL) FatalError("allocating memory for equation solver!");
    gsl_multiroot_fdfsolver_set (s, &eqq, ff);

    if (debug>2) print_state (iter, s, ntemp, useit);

    do {
			iter++;

			// solve
			status = gsl_multiroot_fdfsolver_iterate (s);

			// shift all f[i]
		/*	f0 = gsl_vector_get (s->x, 0);
			for (i=0;i<ntemp;i++)
			{
				f[i] = gsl_vector_get (s->x, i);
				f[i] -= f0;
				gsl_vector_set (s->x, i, f[i]);
			}

			if (debug>2) { fprintf(stderr,"f0=%lf\n",f0); print_state (iter, s, ntemp, useit); }
		 */
			if (debug>2) print_state (iter, s, ntemp, useit);
			if (status) break;  /* check if solver is stuck */

			// check exit condition
			status = gsl_multiroot_test_residual (s->f, MHE_RESIDUAL*actual_ntemp);
      }
    while (status == GSL_CONTINUE && iter < NITERHIS);

    if (debug>1)
    {
    	for (i=0;i<ntemp;i++)
    		if (useit[i]==1) delta += fabs(gsl_vector_get (s->f, i));
    	fprintf (stderr,"\tntemp=%d status = %s (delta=%.5e threshold=%lf)\n",ntemp, gsl_strerror (status),delta,MHE_RESIDUAL*actual_ntemp);
    }

    // set f[i] and check for NaN
    f0 = gsl_vector_get (s->x, 0);
    if (debug>2)  fprintf(stderr,"f0=%lf\n",f0);
	for (i=0;i<ntemp;i++)
	{
		f[i] = gsl_vector_get (s->x, i) - f0;					// shift to have f[0]=0

		if (isnan(f[i]))
		{
			if (debug>0)
			{
				for (j=0;j<ntemp;j++) fprintf(stderr,"f[%d]=%lf ",j,f[j]);
				fprintf(stderr,"\nWARNING: f[%d] is nan in SolveEquation\n",i);
				gsl_multiroot_fdfsolver_free (s);
				gsl_vector_free (ff);
				return -1;
			}
			if (paranoid==1) FatalError("f[i] is nan in SolveEquations");
		}
	}

	i = 1;
	if (status != GSL_SUCCESS && gsl_multiroot_test_residual(s->f, MHE_RESIDUAL2*actual_ntemp) == GSL_SUCCESS) i = 0;	// got stucked
	else if (status != GSL_SUCCESS) i = -1;			// couldn't solve

	gsl_multiroot_fdfsolver_free (s);
	gsl_vector_free (ff);

	return i;
}

void print_state (size_t iter, gsl_multiroot_fdfsolver *s, int ntemp, int *useit)
{
	int i,j;
	double a;

	fprintf (stderr,"iter = %3d\n",(int)iter);
	for (i=0;i<ntemp;i++)
		if (useit[i])
			fprintf (stderr,"f[%d]=%.3f\t",i,gsl_vector_get (s->x, i));
	fprintf(stderr,"\nJ=\n");

	for (i=0;i<ntemp;i++)
		if (useit[i])
		{
			fprintf(stderr,"%2d] ",i);
			for (j=0;j<ntemp;j++)
				if (useit[j])
				{
					 a = gsl_matrix_get(s->J,i,j);
					 fprintf(stderr,"%3.7lf ",a);
				}
			 fprintf(stderr,"\n");
		}

	for (i=0;i<ntemp;i++)
		if (useit[i])
			fprintf (stderr,"eq[%d](x)=%.3e\t",i,gsl_vector_get (s->f, i));

	fprintf(stderr,"\n\n");
}

#else
// ALGORITHM WITH NUMERICAL DERIVATIVE
int SolveEquations(double **h, double *t, double *f, double emin, double ebin, int ntemp, int nbin,
		int debug, double kb, int paranoid, int *useit, int *binok)
{
	int iter=0,status,i,j,actual_ntemp=0;
    const gsl_multiroot_fsolver_type *T;
    gsl_multiroot_fsolver *s;
    double f0;
    const size_t n = ntemp;
    gsl_vector *ff = gsl_vector_alloc (n);
	double delta = 0;


    for (i=0;i<ntemp;i++)
    	gsl_vector_set (ff, i, f[i]);

    for (i=0;i<ntemp;i++) actual_ntemp += useit[i];		// actual number of used equations (for convergence test)
    if (actual_ntemp==0) FatalError("actual_ntemp=0 in SolveEquations");

    struct rparams p = {h,t,emin,ebin,ntemp,nbin,kb,debug,paranoid,useit,binok};

    gsl_multiroot_function eqq;
    eqq.f=&GMHequation;
    eqq.n=n;
    eqq.params=&p;

    T = gsl_multiroot_fsolver_hybrids;
    s = gsl_multiroot_fsolver_alloc (T, n);
    gsl_multiroot_fsolver_set (s, &eqq, ff);

    if (debug>2) print_state (iter, s, ntemp);

    do {
			iter++;

			// solve
			status = gsl_multiroot_fsolver_iterate (s);

/*			// shift all f[i]
			f0 = gsl_vector_get (s->x, 0);
			for (i=0;i<ntemp;i++)
			{
				f[i] = gsl_vector_get (s->x, i);
				f[i] -= f0;
				gsl_vector_set (s->x, i, f[i]);
			}
*/
			if (debug>2) { fprintf(stderr,"f0=%lf\n",f0); print_state (iter, s, ntemp); }
			if (status) break;  /* check if solver is stuck */

			// check exit condition
			status = gsl_multiroot_test_residual (s->f, MHE_RESIDUAL*actual_ntemp);
      }
    while (status == GSL_CONTINUE && iter < NITERHIS);

    if (debug>1)
    {
    	for (i=0;i<ntemp;i++)
    		if (useit[i]==1) delta += fabs(gsl_vector_get (s->f, i));
	fprintf (stderr,"\tntemp=%d status = %s (delta=%.5e threshold=%lf)\n",ntemp, gsl_strerror (status),delta,MHE_RESIDUAL*actual_ntemp);
    }

      // set f[i] and check for NaN
	f0 = gsl_vector_get (s->x, 0);
	for (i=0;i<ntemp;i++)
	{
		f[i] = gsl_vector_get (s->x, i) - f0;
		if (isnan(f[i]))
		{
			if (debug>0)
			{
				for (j=0;j<ntemp;j++) fprintf(stderr,"f[%d]=%lf ",j,f[j]);
				fprintf(stderr,"\nWARNING: f[%d] is nan in SolveEquation\n",i);
				gsl_multiroot_fsolver_free (s);
				gsl_vector_free (ff);
				return -1;
			}
			if (paranoid==1) FatalError("f[i] is nan in SolveEquations");
		}
	}

	i = 1;
	if (status != GSL_SUCCESS && gsl_multiroot_test_residual(s->f, MHE_RESIDUAL2*actual_ntemp) == GSL_SUCCESS) i = 0;	// got stucked
	else if (status != GSL_SUCCESS) i = -1;			// couldn't solve

	gsl_multiroot_fsolver_free (s);
	gsl_vector_free (ff);

	return i;
}

void print_state (size_t iter, gsl_multiroot_fsolver *s, int ntemp)
{
	int i;

	fprintf (stderr,"iter = %3d\n",(int)iter);
	for (i=0;i<ntemp;i++)
		fprintf (stderr,"f[%d]=%.3f\t",i,gsl_vector_get (s->x, i));
	fprintf(stderr,"\n");

	for (i=0;i<ntemp;i++)
		fprintf (stderr,"eq[%d](x)=%.3e\t",i,gsl_vector_get (s->f, i));

	fprintf(stderr,"\n\n");
}
#endif

/*****************************************************************************
 Defines the GMH equations in GSL formalism
 ff are the input free energies, eqq are the output functions whose zeros we seek
 use only those temperatures which has useit[it]==1
 *****************************************************************************/
int GMHequation (const gsl_vector *ff, void *params, gsl_vector *eqq)
{
  int i,ie,it,nbin,ntemp,*useit,*binok;//,paranoid,debug;
  double num,den,fmax,e,emin,ebin,kb;
  double f[NTEMPMAX], eq[NTEMPMAX], *T, **h, n;

  // The GMH equations are (Eq. A.1 jesper's thesis):
  //                                                 \sum_t h_t(E)
  // \sum_E chi_T(E) -----------------------------------------------------------------------------  - 1 = 0
  //                  exp[Fmax] \sum_t chi_t(E) n_t exp[ -E/t + E/T + log(Z_T) - log(Z_t) - Fmax]

  // Take parameters
  h = ((struct rparams *) params)->h;
  T = ((struct rparams *) params)->T;
  emin = ((struct rparams *) params)->emin;
  ebin = ((struct rparams *) params)->ebin;
  ntemp = ((struct rparams *) params)->ntemp;
  nbin = ((struct rparams *) params)->nbin;
  kb = ((struct rparams *) params)->kb;
  //debug = ((struct rparams *) params)->debug;
  //paranoid = ((struct rparams *) params)->paranoid;
  useit = ((struct rparams *) params)->useit;
  binok = ((struct rparams *) params)->binok;

  // Take unknowns from GSL vector ff
  for (i=0;i<ntemp;i++)
	  f[i] = gsl_vector_get (ff, i);


  // Put in eq[i] the system of functions whose zeros must be found
  for (i=0;i<ntemp;i++)									// sum over the ntemp equations
	  if (useit[i]==1)
	  {
		  // normalization constant
		  n = 0;
		  for (ie=0;ie<nbin;ie++)
			  if (binok[ie] == 1 ) n += h[i][ie];

		  // Find maximum exponent to prevent overflows
		   fmax = -9E19;
		   for (ie=0;ie<nbin;ie++)
			   if (binok[ie] == 1)
				for (it=0;it<ntemp;it++)
					if (useit[it]==1)
						{
							e = emin + (double) ebin * ie;		// energy from energy bin
							if ( (h[it][ie]>0) && (h[i][ie]>0) && (-e/(kb*T[it])+e/(kb*T[i])-f[i]+f[it] > fmax) )
								fmax = -e/(kb*T[it])+e/(kb*T[i])-f[i]+f[it];
						}
		   if (fmax <-8E19) fmax = 0;

		  eq[i] = 0;
		  for (ie=0;ie<nbin;ie++)							// sum over energies
			if (h[i][ie]>0 && binok[ie] == 1)									// Ei must belong to the support of the histogram
			{
				e = emin + (double) ebin * ie;		// energy from energy bin
				num = 0;
				den = 0;

				for (it=0;it<ntemp;it++)					// sum over temperatures in the denominator
					if (useit[it]==1)
					{
						num += h[it][ie];						// sum numerator
						if (h[it][ie]>0)						// chi_t(e)
							den += n * exp(-e/(kb*T[it])+e/(kb*T[i])-f[i]+f[it]-fmax);
					}
				eq[i] += exp(-fmax) * num / den ;
			}
		  eq[i] -= 1.;										// the equations contain a "- 1"
	  }

 // put the functions in GSL vector eqq
  for (i=0;i<ntemp;i++)
  {
 	 if (useit[i]==1)
 		 gsl_vector_set (eqq, i, eq[i]);
 	 else
 		 gsl_vector_set (eqq, i, 0);						// if doesn't use this equation, make sure that it appears as satisfied
  }

  return GSL_SUCCESS;
}

/*****************************************************************************
 Defines the derivatives of the GMH equation
 *****************************************************************************/
int GMHequation_df (const gsl_vector *ff, void *params, gsl_matrix *J)
{
  int i,ie,iT,nbin,ntemp,*useit, *binok;//,paranoid,debug;
  double den2[NEBINMAX],fmax[NEBINMAX],e,emin,ebin,kb,sum,d,fmax2;
  double f[NTEMPMAX], deriv, *T, **h, n[NTEMPMAX],s[NEBINMAX];
  gsl_vector *eqq;

  // The GMH derivatives are (Eq. A.1.3 jesper's thesis):
  //  d F_T                                                    \exp[fmax2] \exp[-E/i-E/T-fmax2] \chi_i(E)\chi_T(E) \sum_t h_t(E)
  //  ----- = \delta(T,j) (F_t + 1) - \exp(f_T+f_i) n_i \sum_E -------------------------------------------------------------------
  //  d f_i                                                        \exp[2*fmax]( \sum_t /chi_t(E) n_t \exp[-E/t+f_t-fmax] )^2

  // Take parameters
  h = ((struct rparams *) params)->h;
  T = ((struct rparams *) params)->T;
  emin = ((struct rparams *) params)->emin;
  ebin = ((struct rparams *) params)->ebin;
  ntemp = ((struct rparams *) params)->ntemp;
  nbin = ((struct rparams *) params)->nbin;
  kb = ((struct rparams *) params)->kb;
  //debug = ((struct rparams *) params)->debug;
  //paranoid = ((struct rparams *) params)->paranoid;
  useit = ((struct rparams *) params)->useit;
  binok = ((struct rparams *) params)->binok;

  for (ie=0;ie<nbin;ie++) { s[ie] = 0; den2[ie] = 0; }

  for (i=0;i<ntemp;i++)
	  if (useit[i])
	  {
		  // Take unknowns from GSL vector ff
		  f[i] = gsl_vector_get (ff, i);

		  // normalization constants
		  n[i] = 0;
		  for (ie=0;ie<nbin;ie++)
			  if ( binok[ie] == 1 )
			  {
				  n[i] += h[i][ie];
				  s[ie] += h[i][ie];
			  }
	  }

  // denominator square
  for (ie=0;ie<nbin;ie++)
	  if ( binok[ie] == 1 )
	  {
		  fmax[ie] = -9E19;
		  for (i=0;i<ntemp;i++)					// maximum exponent to prevent overflows
			  if (h[i][ie]>0 && useit[i])
			  {
				e = emin + (double) ebin * ie;		// energy from energy bin
				if (-e/(kb*T[i]) + f[i] > fmax[ie]) fmax[ie] = -e/(kb*T[i]) + f[i];
			  }

		  d = 0;						// square denominator is den2[ie] * exp[2*fmax]
		  for (i=0;i<ntemp;i++)
			  if (h[i][ie]>0 && useit[i])
			  {
				e = emin + (double) ebin * ie;		// energy from energy bin
				d += n[i] * exp(-e/(kb*T[i]) + f[i] - fmax[ie]);
			  }
		  den2[ie] = d*d;
	  }

  // evaulate the equation itself
  eqq = gsl_vector_alloc (ntemp);
  GMHequation (ff,params,eqq);


  for (iT=0;iT<ntemp;iT++)
	  if (useit[iT])
		  for (i=0;i<ntemp;i++)
			  if (useit[i])
			  {
				  fmax2 = -9E19;							// maximum exponent to prevent overflows
				  for (ie=0;ie<nbin;ie++)
				  	if (h[iT][ie]>0 && h[i][ie]>0 && binok[ie] == 1)
				  	{
						e = emin + (double) ebin * ie;		// energy from energy bin
				  		if (-e/(kb*T[i]) - e/(kb*T[iT]) + f[i] + f[iT] - 2.*fmax[ie] > fmax2)
				  			fmax2 = -e/(kb*T[i]) - e/(kb*T[iT]) + f[i] + f[iT] - 2.*fmax[ie];
				  	}

				  if (iT == i) deriv = gsl_vector_get(eqq,iT) + 1.;
				  else deriv = 0;

				  sum = 0;									// actual sum is sum * exp[fmax2]
				  for (ie=0;ie<nbin;ie++)
					  if (h[iT][ie]>0 && h[i][ie]>0 && binok[ie] == 1 )
					  {
						e = emin + (double) ebin * ie;		// energy from energy bin
						sum += s[ie] * exp(-e/(kb*T[i]) - e/(kb*T[iT]) + f[i] + f[iT] - 2.*fmax[ie] - fmax2) / den2[ie];
					  }

				  deriv -= n[i] * sum * exp(fmax2);
				  //printf("fmax2=%lf fmax=%lf\n",fmax2,fmax);
				  gsl_matrix_set (J, iT, i, deriv);
			  }

  gsl_vector_free(eqq);

  return GSL_SUCCESS;
}

int GMHequation_fdf (const gsl_vector *ff, void *params, gsl_vector *eqq, gsl_matrix *J)
{
	GMHequation (ff, params, eqq);
	GMHequation_df (ff, params, J);

    return GSL_SUCCESS;
}
