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
 * AdjustTemperatures.c
 *
 *  Created on: Nov 24, 2009
 *      Author: guido
*/
#include "stempering.h"
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>

//#include <gsl/gsl_rng.h>	// per poter includere global.h che ha dei puntatori a tipi gsl
//#include "../global.h"	// per poter usare G_mcstep (x fini di debug)

struct gsl_param
{
	int ntemp;
	double tfirst;
	double tlast;
	double *g;
	double **boltzp;
	double ebin;
	double emin;
	int nbin;
	int debug;
	double *lg;
	int *binok;
	int gmethod;
};

/*****************************************************************************
 Changes the temperatures of tempering in the loop, adding also a new one.
 Needs InitSTempering to have been called before the loop
 Exits 1 if ok, 0 if failed and should redo, -1 if failed and should go on
 *****************************************************************************/
int AdjustSTempering(struct st_stuff *p, double step)
{
	int i,it,j;
	double tmin,tbin;
	FILE *fthe;

	// print acceptance statistics
	PrintStatistics(p, step);

	// print current dump energies
	if (strcmp(p->st_nfdumb2,""))
	{
		fthe = fopen(p->st_nfdumb2,"w");
		if (!fthe) FatalError("Cannot open file for writing dumb2 averages");
		PrintAverageEnergies(fthe,p->st_h,p->st_temp,p->st_ntemp,p->st_emin,p->st_ebin,p->st_nbin);
		fclose(fthe);
	}

	// if activated, do not remove replicas if lowest temperature is below the threshold
	if (p->st_tonlyadd>0)
		if (p->st_temp[p->st_ntemp-1]<p->st_tonlyadd)
		{
			p->st_removet = 0;
		}

	if (p->st_debug>0) fprintf(stderr,"\nADJUSTING %d TEMPERATURES (step=%.0lf, iter=%d)\n",p->st_ntemp,step,p->st_iter);
	if (p->st_debug>0 && p->st_removet>0) fprintf(stderr,"(trying to reduce the number of temperatures)\n");

	// Check if all temperatures are sampled. If not, redo simulation
	if (p->st_iter>0)
	{

		if (p->st_failure<p->st_nfail1)		// if fails for the first nfail1 times, fill with T optimized from reliable_lg
			j = CheckHistogramsAdjust(p->st_counts,p->st_temp,p->st_g,p->st_reliable_lg,&(p->st_ntemp),p->st_emin,p->st_ebin,p->st_nbin,p->st_hthresh,p->st_debug,p->st_ntempmax,p->st_h,
							p->st_prob_up,p->st_prob_down,p->st_phthresh,p->st_binok,p->st_gmethod);
		else if (p->st_failure<p->st_nfail2)	// if fails for the first nfail2 times, fill with T from extrapolated free energies
				j =CheckHistogramsAdjust3(p->st_counts,p->st_temp,p->st_g,&(p->st_ntemp),p->st_nbin,p->st_hthresh,p->st_debug,p->st_ntempmax,p->st_h,
						p->st_prob_up,p->st_prob_down,p->st_phthresh);
		else{ 		// if previous run failed as well, use reliable_t and reliable_g
				//MODIFICHE: se sono qui allora tornero allo step precedente e aggiungero una temperatura, quindi devo mettere p->st_ttarget_harvest=0
				p->st_ttarget_harvest=0;
				//FINE MODIFICHE
				j = CheckHistogramsAdjust2(p->st_counts,p->st_temp,p->st_g,&(p->st_ntemp),p->st_hthresh,p->st_debug,p->st_reliable_t,
						p->st_reliable_g,&(p->st_nreliable_t),p->st_prob_up,p->st_prob_down,p->st_phthresh);
				p->st_failure = 0;
				}

		// write on accept.dat
		//fprintf(p->st_fpacc,"%lf %d\n",step,j); fflush(NULL);

		if (j!=1) (p->st_failure) ++;						// counts how many consecutive failed runs
		if (j==-1) return -1;							// continue with same temperature, without resetting histograms
		else if (j==0)									// continue with same temperature, reset histograms
		{	
			for (it=0;it<p->st_ntemp;it++)
				for (i=0;i<p->st_nbin;i++) p->st_h[it][i]=0;

			for (it=0;it<p->st_ntemp;it++)
			{
				p->st_prob_up[it]=0;
				p->st_prob_down[it]=0;
				p->st_counts[it]=0;
			}
			return 0;
		}
	}
               
	//MODIFICHE: se la simulazione con la temperatura target inclusa e' buona, faccio simulated tempering normale
       	if(p->st_ttarget_harvest==1)       p->st_nm=1;
       	//FINE MODIFICHE
	
	// if it is here, it means that previous simulation was ok.

	// If used, add last histograms to the pile
	if (p->st_oldh != NULL)
		AddHistogramPile(p->st_ntemp,p->st_temp,p->st_h,p->st_oldh,p->st_oldt,p->st_oldt_iter,&(p->st_noldt),p->st_nbin,p->st_nkeepold,p->st_iter,p->st_debug,p->st_sum);

	// store a reliable density of states
	if (p->st_iter>0)
		for (i=0;i<p->st_nbin;i++)
			p->st_reliable_lg[i] = p->st_current_lg[i];
	else													// in the first run, there is not a current_lg
    {
		FilterHistograms(p->st_h,p->st_htmp,p->st_ntemp,p->st_nbin,p->st_binthresh);
		p->st_reliable_lg = MHistogram(p->st_htmp,p->st_temp,p->st_ntemp,p->st_nbin,p->st_emin,p->st_ebin,p->st_debug,KB,p->st_reliable_lg,
				p->st_f,p->st_paranoid,p->st_ignoreb,NULL,p->st_deltat,p->st_binok);
    }
	// store reliable temperatures and weights
	for (it=0;it<p->st_ntemp;it++)
	{
		p->st_reliable_t[it] = p->st_temp[it];
		p->st_reliable_g[it] = p->st_g[it];
	}
	p->st_nreliable_t = p->st_ntemp;

	// print restart
	if(p->st_restart>-1)
        ManageRestart(p,step);

	// advance clock of number of adjsutments
	p->st_iter ++;
	p->st_failure = 0;

	// Calculate density of states
	if (p->st_oldh == NULL)						// uses only last histograms
	{
		p->st_ntemp = OrderTemperatures(p->st_ntemp,p->st_temp,p->st_g,p->st_prob_up,p->st_prob_down,p->st_h,p->st_nbin,p->st_debug,p->st_counts);
		FilterHistograms(p->st_h,p->st_htmp,p->st_ntemp,p->st_nbin,p->st_binthresh);
		p->st_current_lg = MHistogram(p->st_htmp,p->st_temp,p->st_ntemp,p->st_nbin,p->st_emin,p->st_ebin,p->st_debug,KB,p->st_current_lg,
				p->st_f,p->st_paranoid, p->st_ignoreb,NULL,p->st_deltat,p->st_binok);
	}
	else									// uses all past histograms
	{
		p->st_ntemp = OrderTemperatures(p->st_noldt,p->st_oldt,NULL,NULL,NULL,p->st_oldh,p->st_nbin,p->st_debug,p->st_oldt_iter);
		FilterHistograms(p->st_oldh,p->st_htmp,p->st_noldt,p->st_nbin,p->st_binthresh);
		p->st_current_lg = MHistogram(p->st_htmp,p->st_oldt,p->st_noldt,p->st_nbin,p->st_emin,p->st_ebin,p->st_debug,KB,p->st_current_lg,p->st_f,
				p->st_paranoid,p->st_ignoreb,NULL,p->st_deltat,p->st_binok);
	}

	// print thermodynamics
	if (p->st_debug>1) fprintf(stderr,"\nWriting thermodynamics to file...\n");
	fthe = fopen(p->st_nfthe,"w");
	if (!fthe) FatalError("Cannot open file for writing thermodynamics");
	tmin = p->st_temp[p->st_ntemp-1]/2;
	tbin = (p->st_temp[0]*3/2 - tmin)/NTBIN;
	CalculateThermodynamics(p->st_current_lg,p->st_emin,p->st_ebin,p->st_nbin,tmin,tbin,NTBIN,p->st_out,KB,NULL,NULL,p->st_binok);
	PrintThermodynamics(fthe,p->st_out,tmin,tbin,NTBIN);
	fclose(fthe);

	// print density of states
//	if (p->st_debug>1) fprintf(stderr,"\nWriting density of states to file...\n");
//	fthe = fopen(p->st_nfdos,"w");
//	if (!fthe) FatalError("Cannot open file for writing density of states");
//	for (i=0;i<p->st_nbin;i++) fprintf(fthe,"%lf\t%lf\n",p->st_emin+p->st_ebin*i,p->st_reliable_lg[i]);
//	fclose(fthe);

	// print histograms
	if (strcmp(p->st_nfhisto,""))
	{
		if (p->st_oldh == NULL)
			PrintHistogram(p->st_nfhisto,p->st_htmp,p->st_temp,p->st_ntemp,p->st_nbin,p->st_emin,p->st_ebin);
		else
			PrintHistogram(p->st_nfhisto,p->st_htmp,p->st_oldt,p->st_noldt,p->st_nbin,p->st_emin,p->st_ebin);
	}

	// print dumb energies
	if (strcmp(p->st_nfdumb,""))
	{
		fthe = fopen(p->st_nfdumb,"w");
		if (!fthe) FatalError("Cannot open file for writing dumb averages");
		if (p->st_oldh == NULL)
			PrintAverageEnergies(fthe,p->st_h,p->st_temp,p->st_ntemp,p->st_emin,p->st_ebin,p->st_nbin);
		else
			PrintAverageEnergies(fthe,p->st_oldh,p->st_oldt,p->st_noldt,p->st_emin,p->st_ebin,p->st_nbin);
		fclose(fthe);
	}

	// check if ntmax is reached, then end
	if (p->st_ntemp == p->st_ntempmax-1)
	{
		fprintf(stderr,"ntmax reached\n");
		fprintf(stderr,"ENDING THE SIMULATION\n\n\n");
		exit(0);
	}

	OptimalWeights(p->st_ntemp,p->st_temp,p->st_g,p->st_emin,p->st_ebin,p->st_nbin,p->st_current_lg,p->st_debug,p->st_binok,p->st_gmethod);

	// adjust temperatures from 1 to ntemp-1, calculating optimal weights
	if (p->st_removet>0)
	{
		if (p->st_ntemp>3)
			p->st_ntemp = ReduceTemperatures(p->st_ntemp,p->st_temp,p->st_g,p->st_emin,p->st_ebin,p->st_nbin,p->st_debug,p->st_current_lg,
					p->st_pthresh,p->st_ntempmax,p->st_removet,p->st_boltzp,p->st_binok,p->st_gmethod);

		if (p->st_ntemp==0) return 0;
	}
	else
	{
		if (p->st_ntemp>2)
			MaximizeTempProb(p->st_ntemp,p->st_temp,p->st_g,p->st_emin,p->st_ebin,p->st_nbin,p->st_debug,p->st_current_lg,p->st_boltzp,p->st_binok,p->st_gmethod);
	}

	// add new temperature
	// MODIFICHE: se ho raggiunto la temperatura target e la simulazione e' andata male, non aggiungo altre temperature
	if(p->st_ttarget_harvest<1){
	if (p->st_k>0)								// proportionally to energy sigma
		p->st_ntemp = AddNewTemperature(p->st_ntemp,p->st_temp,p->st_g,p->st_current_lg,p->st_emin,p->st_ebin,p->st_nbin,p->st_debug,p->st_k,
				p->st_paranoid,p->st_binok,p->st_gmethod);
	else if (p->st_p_new<0)						// at wished exchange probability
		p->st_ntemp = AddNewTemperature2(p->st_ntemp,p->st_temp,p->st_g,p->st_current_lg,p->st_emin,p->st_ebin,p->st_nbin,p->st_debug,p->st_p_new,
				p->st_paranoid,p->st_boltzp,p->st_binok,p->st_gmethod);

	if (p->st_debug>0)
	{
		fprintf(stderr,"Now there are %d temperatures:\n",p->st_ntemp);
		for (j=0;j<p->st_ntemp;j++) fprintf(stderr,"%lf (g=%lf) ",p->st_temp[j],p->st_g[j]);
		fprintf(stderr,"\nBACK TO SIMULATION\n\n");
	}

	p->st_ntemp = OrderTemperatures(p->st_ntemp,p->st_temp,p->st_g,p->st_prob_up,NULL,NULL,p->st_nbin,p->st_debug,NULL);
	}
	//FINE MODIFICHE
	// reset counts
	for (it=0;it<p->st_ntemp;it++)
	{
		p->st_prob_up[it]=0;
		p->st_prob_down[it]=0;
		p->st_counts[it]=0;

		for(i=0;i<p->st_nbin;i++)
		p->st_h[it][i] = 0;
	}


	return 1;
}

/*****************************************************************************
 Adds a new temperature below the others, at some fraction of the energy stdev
 returns the number of temperatures
 *****************************************************************************/
int AddNewTemperature(int ntemp, double *temp, double *g, double *lg, double emin, double ebin, int nbin,
		int debug, double k, int paranoid, int *binok, int gmethod)
{
	int it;
	double t,em,newtemp=-1,em0,sigma0;

	// calculate average energy and stddev of the minimum temperature
	em0 = AverageEnergy(lg,emin,ebin,nbin,temp[ntemp-1],KB,binok);
	sigma0 = SigmaEnergy(lg,emin,ebin,nbin,temp[ntemp-1],KB,binok);

	// find temperature corresponding to (em0 - k*sigma0)
	for (it=NTBIN-1;it>0;it--)
	{
		t = (double) it * temp[ntemp-1]/NTBIN;				// run over T from lowest to zero
		em = AverageEnergy(lg,emin,ebin,nbin,t,KB,binok);

		if (debug>3) fprintf(stderr,"\tit=%d T=%lf em=%lf\n",it,t,em);

		if (em <= em0-k*sigma0)
		{
			newtemp = t;
			break;
		}
	}

	// if it cannot find
	if (newtemp<0)
		{
			if (paranoid==1)
			{
				fprintf(stderr,"em0=%lf sigma0=%lf em0-k*sigma0=%lf\n",em0,sigma0,em0-k*sigma0);
				FatalError("Cannot find temperature for new replica");
			}
			else
			{
				if (debug>0)
				{
					fprintf(stderr,"em0=%lf sigma0=%lf em0-k*sigma0=%lf\n",em0,sigma0,em0-k*sigma0);
					fprintf(stderr,"WARNING: cannot add new temperature (continuing)\n");
				}
				return ntemp;
			}
		}

	temp[ntemp] = newtemp;								// set new temperature

	// assign weight to new temperature
	g[ntemp] = CalculateG(lg, emin, ebin, nbin, temp[ntemp], KB, binok, gmethod) - CalculateG(lg, emin, ebin, nbin, temp[0], KB, binok, gmethod);

	ntemp ++;

	if (debug>0) {
		fprintf(stderr,"Added new temperature: T=%lf (ntemp=%d)\n",newtemp,ntemp);
		fprintf(stderr,"\t<E(%lf)>=%lf+-%lf\t<E(%lf)>=%lf\tg=%lf\n",temp[ntemp-2],em0,sigma0,temp[ntemp-1],em,g[ntemp-1]);
	}

	return ntemp;
}

/*****************************************************************************
 Calculates the weights of simulated tempering which make most uniform
 the exchange rates. It uses the density of states lg.
 *****************************************************************************/
void OptimalWeights(int ntemp, double *t, double *g, double emin, double ebin, int nbin, double *lg, int debug, int *binok, int gmethod)
{
	int it;
	double g0;

	// If g(T) = F(T)/T, the exhange probability is mostly homogeneous, i.e. p(T->T') = p(T'->T)

	if (debug>0) fprintf(stderr,"  Calculate weights...\n");

	// highest temperature
	g0 = CalculateG(lg, emin, ebin, nbin, t[0], KB, binok, gmethod);
	g[0]=0.;
	if (debug>0) fprintf(stderr,"  T[0]=%lf g[0]=%lf  ",t[0],g[0]);

	// other temperatures (to zero-th reference temperature)
	for (it=1;it<ntemp;it++)
	{
		g[it] = CalculateG(lg, emin, ebin, nbin, t[it], KB, binok, gmethod) - g0;
		if (debug>0) fprintf(stderr,"  T[%d]=%lf g[%d]=%lf  ",it,t[it],it,g[it]);
	}

}

/*****************************************************************************
 Average over p(E) of the exchange probability.
 Used by GSL optimizer.
 *****************************************************************************/
double EstimatedJumpProb (const gsl_vector *v, void *params)
{
	int i,it;
	double beta1,beta2,delta,prob,lptot=0.,energy;
	struct gsl_param *gp = (struct gsl_param *) params;
	double t[NTEMPMAX];

	// copy the set of temperatures in *t
	for(i=0;i<gp->ntemp-2;i++)
		t[i+1] = gsl_vector_get(v, i);
	t[0] = gp->tfirst;
	t[gp->ntemp-1] = gp->tlast;

	//				/
	//	p(i->j) =  / dE min(1,exp[(beta_i - beta_j) E - g_i + g_j]) * p(E,beta_i)
	//			  /E

	//
	//  L = \prod_i,j p(i->j)
	//

	// calculate optimal weights g
	OptimalWeights(gp->ntemp,t,gp->g,gp->emin,gp->ebin,gp->nbin,gp->lg,0,gp->binok,gp->gmethod);

	// Calculates p(E) at all T
	CalculateThermodynamics(gp->lg,gp->emin,gp->ebin,gp->nbin,0,0,gp->ntemp,NULL,KB,t,gp->boltzp,gp->binok);

	// probabilities of decreasing T
	if (gp->debug>3)
		{
			fprintf(stderr,"   Jump probabilities:\n   T=");
			for (it=0;it<gp->ntemp;it++) fprintf(stderr,"%lf ",t[it]);
			fprintf(stderr,"\n");
		}

	for (it=0;it<gp->ntemp-1;it++)		// taking away last
	{
		prob = 0.;					// calculates p(i,i+1)
		for (i=0;i<gp->nbin;i++)			// loop over energies
		{
			if (gp->binok[i]==1) 
			{
				beta1 = 1./(KB*t[it]);
				beta2 = 1./(KB*t[it+1]);
				energy = gp->emin + (double) gp->ebin * i;
				delta = (beta1-beta2) * energy + gp->g[it+1] - gp->g[it];
				if (delta>0) prob += gp->boltzp[it][i];
				else prob += exp(delta) * gp->boltzp[it][i];
			}
		}
		if (gp->debug>3) fprintf(stderr,"   %8lf",prob);
		lptot += log(prob);				// log of overall probability of changes down
	}
	if (gp->debug>3) fprintf(stderr,"\n           ");

	// probabilities of increasing T
	for (it=1;it<gp->ntemp;it++)			// taking away first
	{
		prob = 0.;						// calculates p(i,i+1)
		for (i=0;i<gp->nbin;i++)			// loop over energies
		{
			if (gp->binok[i]==1) 
			{
				beta1 = 1./(KB*t[it]);
				beta2 = 1./(KB*t[it-1]);
				energy = gp->emin + (double) gp->ebin * i;
				delta = (beta1-beta2) * energy + gp->g[it-1] - gp->g[it];
				if (delta>0) prob += gp->boltzp[it][i];
				else prob += exp(delta) * gp->boltzp[it][i];
			}
		}
		if (gp->debug>3) fprintf(stderr,"   %8lf",prob);
		lptot += log(prob);				// log of overall probability of changes up	
	}

	if (gp->debug>3) fprintf(stderr,"\n   Estimated log jump probability = %lf\n",lptot);

  return -lptot;									// maximize, not minimize
}

/*****************************************************************************
 Changes the temperatures in-between to optimize the exchange rate.
 Needs an estimate of the density of states lg
 Uses GSL multioptimizer (without derivative)
 returns the log of estimated probability
 *****************************************************************************/
double MaximizeTempProb(int ntemp, double *temp, double *g, double emin, double ebin, int nbin, int debug, double *lg,
		double **boltzp, int *binok, int gmethod)
{
	int i;
	size_t iter = 0;
	int status;
	double lprob = -9999.;
	double oldt[NTEMPMAX],oldg[NTEMPMAX];

	for (i=0;i<ntemp;i++)
	{
		oldt[i] = temp[i];
		oldg[i] = g[i];
	}

	gsl_multimin_function minex_func;
	const gsl_multimin_fminimizer_type *T;
	gsl_multimin_fminimizer *s = NULL;
	gsl_vector *x,*ss;

	double size=0.;

     if (debug>0) fprintf(stderr,"Rearrange %d temperatures...\n", ntemp);
     if (debug>1) fprintf(stderr,"Starting T:\n");

     // allocate and fill struct for paramaters of function to be maximized
     struct gsl_param gp = {ntemp,temp[0],temp[ntemp-1],g,boltzp,ebin,emin,nbin,debug,lg,binok,gmethod};

     minex_func.n = ntemp-2;
     minex_func.f = &EstimatedJumpProb;
     minex_func.params = &gp;

     // starting point
     OptimalWeights(ntemp,temp,g,emin,ebin,nbin,lg,0,binok,gmethod);
     x = gsl_vector_alloc (ntemp-2);
     for (i=0;i<ntemp-2;i++)
     {
    	 gsl_vector_set (x, i, temp[i+1]);
    	 if (debug>1) fprintf(stderr," T[%d]=%lf\n",i+1,temp[i+1]);
     }

     if (debug>1) fprintf(stderr," log(p)=%lf\n",-EstimatedJumpProb(x,&gp));

     // set initial step sizes to 1
     ss = gsl_vector_alloc (ntemp-2);
     gsl_vector_set_all (ss, 1.0);

     T = gsl_multimin_fminimizer_nmsimplex;
     s = gsl_multimin_fminimizer_alloc (T, ntemp-2);

     gsl_multimin_fminimizer_set (s, &minex_func, x, ss);

     // maximization loop
     do
       {
         iter++;
         status = gsl_multimin_fminimizer_iterate (s);

         if (status)
           break;

         size = gsl_multimin_fminimizer_size (s);
         status = gsl_multimin_test_size (size, CONVERGENCE);
         if (debug>3)
        	 {
				 fprintf(stderr,"\n * iter=%d\tNew temperatures:\n",(int)iter);
				 for (i=0;i<ntemp-2;i++) fprintf(stderr," T[%d]=%lf\n",i+1,gsl_vector_get (s->x, i));
				 fprintf(stderr," status = %s\n",gsl_strerror (status));
        	 }

    	 lprob = -s->fval;				// outcoming log(p)

       }
     while (status == GSL_CONTINUE && iter < NITERMAX);

     // after NITERMAX, accept the optimization if this weaker condition is met
     status = gsl_multimin_test_size (size, CONVER_WEAK);

     if (debug>0)
     {
    	 if (status == GSL_SUCCESS)
    		 fprintf (stderr,"Minimum found:\t\t(size=%e)\n",size);
    	 else
    		 fprintf(stderr,"Minimum NOT found ");
    	 fprintf (stderr,"after %5d iterations,  log(p)/n=%lf size=%e\n", (int)iter, (double)lprob/ntemp,size);
    	 fprintf (stderr,"T[0]=%lf ",temp[0]);
    	 for (i=0;i<ntemp-2;i++) fprintf (stderr,"T[%d]=%lf ",i+1,gsl_vector_get (s->x, i) );
    	 fprintf(stderr,"T[%d]=%lf\n",ntemp-1,temp[ntemp-1]);
    	 fprintf (stderr,"g[0]=%lf ",g[0]);
    	 for (i=1;i<ntemp-1;i++) fprintf (stderr,"g[%d]=%lf ",i,g[i]);
    	 fprintf(stderr,"g[%d]=%lf\n",ntemp-1,g[ntemp-1]);
     }

     // copy back new temperatures
	 if (status == GSL_SUCCESS)
	 {
		 for (i=0;i<ntemp-2;i++)
			 temp[i+1] = gsl_vector_get (s->x, i);
		 OptimalWeights(ntemp,temp,g,emin,ebin,nbin,lg,debug,binok,gmethod);
		 if (debug>0) fprintf(stderr,"Change accepted.\n");
	 }
	 else
	{
		 for (i=0;i<ntemp;i++)
		 {
			 temp[i] = oldt[i];
			 g[i] = oldg[i];
		 }
		 if (debug>0) fprintf(stderr,"Change rejected.\n");
		 lprob = -999999.;
	}

     gsl_multimin_fminimizer_free (s);
     gsl_vector_free (x);
     gsl_vector_free (ss);

     return lprob;
}

/*****************************************************************************
 Reduces the number f temperatures keeping the root of the exchange probablity
 above a threshold.
 Needs an estimate of the density of states lg
 Uses GSL multioptimizer (without derivative)
 returns the new number of temperatures
 *****************************************************************************/
int ReduceTemperatures(int ntemp, double *temp, double *g, double emin, double ebin, int nbin, int debug, double *lg,
		double lpthresh, int ntmax, int removet, double **boltzp, int *binok, int gmethod)
{
	int i,ntold;
	double lprob,oldt[NTEMPMAX],oldg[NTEMPMAX];

	if (debug>0) fprintf(stderr,"Reducing temperatures\n");

	for (i=0;i<ntemp;i++)
	{
		oldt[i] = temp[i];
		oldg[i] = g[i];
	}
	ntold = ntemp;

	// automatically set lpthresh to last value
	if (lpthresh>8)
	{
		gsl_vector *x;
		struct gsl_param gp = {ntemp,temp[0],temp[ntemp-1],g,boltzp,ebin,emin,nbin,debug,lg,binok,gmethod};
		x = gsl_vector_alloc (ntemp-2);
		for (i=0;i<ntemp-2;i++) gsl_vector_set (x, i, temp[i+1]);
		lpthresh = -EstimatedJumpProb(x,&gp)/ntemp;
		gsl_vector_free (x);

		if (debug>0) fprintf(stderr," Target lp=%lf\n",lpthresh);
	}

	// increase number of temperatures until prob>pthresh
	if (removet==99)					// reduce to any ntemp
	{
		ntemp = 2;
		temp[1] = oldt[ntold-1];
	}
	else if (ntemp-removet>2)			// reduce at most of removet
	{
		ntemp -= removet;
		temp[ntemp-1] = oldt[ntold-1];
	}
	else
	{
		ntemp = 2;
		temp[1] = oldt[ntold-1];
	}

	do 	{
			temp [ntemp] = temp[ntemp-1];						// the first and last T should not change
			temp [ntemp-1] = (temp[ntemp-2] + temp[ntemp])/2.;	// the new one is put halfway between the last two
			ntemp ++;

			if (ntemp == ntmax) FatalError("NTMAX reached in ReduceTemperatures");

			lprob = MaximizeTempProb(ntemp,temp,g,emin,ebin,nbin,debug,lg,boltzp,binok,gmethod)/ntemp;

			if (debug>1) fprintf(stderr," Trying with %d temperatures --> log(p)/n=%lf\n",ntemp,lprob);

			if (ntemp >= ntold)								// ntemp cannot be larger than before
			{
				for (i=0;i<ntold;i++)
				{
					temp[i] = oldt[i];
					g[i] = oldg[i];
				}
				ntemp = ntold;
				if (debug>0) fprintf(stderr,"...nothing done.\n");
				return ntemp;
			}

		} while (lprob < lpthresh);

	if (debug>0)
		{
			fprintf(stderr,"\nChosen %d temperatures (log(p)/n=%lf):\n",ntemp,lprob);
			for (i=0;i<ntemp;i++) fprintf(stderr,"T[%d] = %lf (g=%lf)\n",i,temp[i],g[i]);
			fprintf(stderr,"\n");
		}

	return ntemp;
}

/*****************************************************************************
 Adds a new temperature below the others, at wished estimated exchange probability
 returns the number of temperatures
 *****************************************************************************/
int AddNewTemperature2(int ntemp, double *temp, double *g, double *lg, double emin, double ebin, int nbin,
		int debug, double p_new, int paranoid, double **boltzp, int *binok, int gmethod)
{
	int it;
	double t,newtemp=-1,p1,p2;


	// loop over decreasing temperatures
	for (it=NTBIN-1;it>0;it--)
	{
		t = (double) it * temp[ntemp-1]/NTBIN;				// run over T from lowest to zero

		p1 = EstimatedJumpSingleProbability(temp[ntemp-1],t,lg,emin,ebin,nbin,debug,boltzp,binok,gmethod);
		p2 = EstimatedJumpSingleProbability(t,temp[ntemp-1],lg,emin,ebin,nbin,debug,boltzp,binok,gmethod);
		if (debug>3) fprintf(stderr,"\tit=%d Tnew=%lf log(p)=%lf+%lf\n",it,t,p1,p2);

		if (it == NTBIN-1 && p1+p2<p_new+p_new)			// if the first is already below p_new
		{
			if (paranoid == 1 || debug>0){
				fprintf(stderr,"In adding new temperature, the smallest annealing (T=%lf) already gives \n probabilities log(p)=%lf+%lf below the threshold %lf.\n",
						t,p1,p2,p_new);
						return ntemp;
						}
			if (paranoid == 1)
				FatalError("Cannot add new temperature");
		}

		if (p1+p2<p_new+p_new)						// if the probability decreases below the threshold
		{
			newtemp = (double) (it+1) * temp[ntemp-1]/NTBIN; // i.e. T at the previous step
			break;
		}
	}

	// if it cannot find
	if (newtemp<0)
		{
			newtemp = (double) 2. * temp[ntemp-1]/NTBIN;
		}
	temp[ntemp] = newtemp;								// set new temperature

	// assign weight to new temperature
	g[ntemp] = CalculateG(lg, emin, ebin, nbin, temp[ntemp], KB, binok, gmethod) - CalculateG(lg, emin, ebin, nbin, temp[0], KB, binok, gmethod);

	ntemp ++;

	if (debug>0)
	{
		fprintf(stderr,"Added new temperature: T=%lf (ntemp=%d)\n",newtemp,ntemp);
		fprintf(stderr,"\tp_exch=%lf+%lf\n",p1,p2);
	}

	return ntemp;
}

/*****************************************************************************
 Average over p(E) of the log of exchange probability, knowing t1 to find t2
 *****************************************************************************/
double EstimatedJumpSingleProbability(double t1, double t2, double *lg, double emin, double ebin, int nbin,
		int debug, double **boltzp, int *binok, int gmethod)
{
	int i;
	double beta1,beta2,delta,prob,lptot=0.,energy,g1,g2;


	//  		    /
	//	p(i->j) =  / dE min(1,exp[(beta_i - beta_j) E - g_i + g_j]) * p(E,beta_i)
	//	          /E


	// Calculates p(E) at t1
	CalculateThermodynamics(lg,emin,ebin,nbin,0,0,1,NULL,KB,&t1,boltzp,binok);

	// Calculates weights
	g1 = CalculateG(lg, emin, ebin, nbin, t1, KB, binok, gmethod);
	g2 = CalculateG(lg, emin, ebin, nbin, t2, KB, binok, gmethod);

	if (debug>3)
			fprintf(stderr,"   Jump probabilities for T=%lf -> T=%lf\n",t1,t2);

	prob = 0.;
	for (i=0;i<nbin;i++)			// loop over energies
	{
		if (binok[i] == 1){
			beta1 = 1./(KB*t1);
			beta2 = 1./(KB*t2);
			energy = emin + (double)ebin * i;
			delta = (beta1-beta2) * energy + g2 - g1;
			if (delta>0) prob += boltzp[0][i];
			else prob += exp(delta) * boltzp[0][i];
		}
	}

	lptot = log(prob);					// log of overall probability of changes down

	if (debug>3) fprintf(stderr,"\n   Estimated log jump probability = %lf) \n",lptot);

  return lptot;									// maximize, not minimize
}



/*****************************************************************************
 Check if all temperatures are sampled, without changing anything
 if not, returns 0
 if yes, returns 1
 *****************************************************************************/
int CheckHistograms(int *counts, double *temp, int ntemp, double hthresh, int debug)
{
	int it;
	double z[NTEMPMAX],ztot=0;

	if (debug>0) fprintf(stderr,"Check histogram counts:    (thresh=%lf)\n",hthresh/ntemp);

	for (it=0;it<ntemp;it++) ztot += counts[it];

	for (it=0;it<ntemp;it++)
	{
		z[it] = (double) counts[it] / ztot;
		if (debug>0) fprintf(stderr,"z(T=%lf)=%lf ",temp[it],z[it]);
	}
	if (debug>0) fprintf(stderr,"\n");

	for (it=0;it<ntemp;it++)
		if (z[it] <= hthresh/ntemp) return 0;

	return 1;
}

/*****************************************************************************
 Add histograms to the pile of old histograms
 if n=0, add all
 if n>0, keep only those belonging to the last n runs
 *****************************************************************************/
void AddHistogramPile(int ntemp, double *temp, double **h, double **oldh, double *oldt, int *oldt_iter, int *noldt, int nbin, int n, int iter, int debug, int sum)
{
	int it,it2,discard,i,did=0,ie,iold;

	if ((*noldt)+ntemp >= NHISTOMAX) FatalError("NHISTOMAX too small");

	// add all histograms to oldhisto
	for (it=0;it<ntemp;it++)
	{
		iold = it + *noldt;					// index of new histogram
		for (i=0;i<nbin;i++) oldh[iold][i] = h[it][i];			// copy histogram
		oldt[iold] = temp[it];					// copy temperature
		oldt_iter[iold] = iter;					// label each histo according to when it was collected
	}

	(*noldt) += ntemp;
	if (debug>0) fprintf(stderr,"Adding %d new histograms to the pile\n",ntemp);

	// if two histograms are at the same temperature, sum together
	if (sum==1)
	{
		did = 0;
		for (it=0;it<(*noldt);it++)
			for (it2=it+1;it2<(*noldt);it2++)
				if (oldt[it]<=oldt[it2]+EPSILON && oldt[it]>=oldt[it2]-EPSILON)		// if two temperatures are equal
				{
					for (i=0;i<nbin;i++)
						{
							oldh[it][i] += oldh[it2][i];							// add latter to former
							oldh[it2][i] = oldh[(*noldt)-1][i];						// delete latter
						}
					oldt[it2] = oldt[(*noldt)-1];
					oldt_iter[it2] = oldt_iter[(*noldt)-1];
					(*noldt) --;
					it2 --;
					did ++;
				}
		if (debug>0) fprintf(stderr,"Summing %d histogram to those with same temperature\n",did);
	}
	// if two histograms are at the same temperature, delete the older
	else if (sum==2)
	{
		did = 0;
		for (it=0;it<(*noldt);it++)
			for (it2=it+1;it2<(*noldt);it2++)
				if (oldt[it]<=oldt[it2]+EPSILON && oldt[it]>=oldt[it2]-EPSILON)		// if two temperatures are equal
				{
					if (oldt_iter[it]<oldt_iter[it2])							// if it is older than it2
					{
						for (i=0;i<nbin;i++)
						{
							oldh[it][i] = oldh[it2][i];							// overwrite it
							oldh[it2][i] = oldh[(*noldt)-1][i];					// delete last
						}
						oldt[it2] = oldt[(*noldt)-1];
						oldt_iter[it2] = oldt_iter[(*noldt)-1];
						it2 --;
					}
					else														// if it2 is older than it
					{
						for (i=0;i<nbin;i++)
						{
							oldh[it2][i] = oldh[it][i];							// overwrite it2
							oldh[it][i] = oldh[(*noldt)-1][i];					// delete last
						}
						oldt[it] = oldt[(*noldt)-1];
						oldt_iter[it] = oldt_iter[(*noldt)-1];
						it--;
					}
					(*noldt) --;
					did ++;
				}
		if (debug>0) fprintf(stderr,"Deleted %d older histograms because of redundant temperatures\n",did);
	}

	// if enabled, discard histograms older than n iterations
	did = 0;
	if (n>0)
	{
		discard = iter - n;
		for (i=0;i<*noldt;i++)
			if (oldt_iter[i] <= discard)
			{
				for (ie=0;ie<nbin;ie++) oldh[i][ie] = oldh[(*noldt)-1][ie];
				oldt[i] = oldt[(*noldt)-1];
				oldt_iter[i] = oldt_iter[(*noldt)-1];
				(*noldt) --;
				i--;
				did++;
			}
		if (debug>0) fprintf(stderr,"Removing %d old histograms\n",did);
	}

	if (debug>0) fprintf(stderr,"Now there are %d histogram in the pile\n",*noldt);

}


/*****************************************************************************
 Check if all temperatures are sampled; if not, return to last reliable
 temperatures, extrapolating the weight for the new lowest T.
 if yes, returns 1
 if there are holes among temperatures returns 0
 *****************************************************************************/
int CheckHistogramsAdjust2(int *counts, double *temp, double *g, int *ntemp, double hthresh, int debug, double *reliable_t,
		double *reliable_g, int *nreliable_t, double *prob_up, double *prob_down, double pthresh)
{
	int it,j,nbadh=0,ngoodh=0,ok;
	double z[NTEMPMAX],ztot=0,lastt;

	// obtain data from *counts
	for (it=0;it<*ntemp;it++)
	{
		z[it] = (double) counts[it];
		ztot += (double) counts[it];
	}

	// normalize counts over all temperatures
	if (debug>0) fprintf(stderr,"Check histogram counts (2): (threshold=%lf)\n",hthresh/(*ntemp));
	for (it=0;it<*ntemp;it++)
	{
		z[it] /= ztot;
		ok=1;

		// locate holes according to z and p_up/p_down
		if (it==0)
		{
			if (z[it]<=hthresh/(*ntemp) || prob_down[it]<pthresh*counts[it]) ok=0;
		}
		else if (it>0 && it<(*ntemp)-1)
		{
			if (z[it]<=hthresh/(*ntemp) || prob_up[it]<pthresh*counts[it] || prob_down[it]<counts[it]*pthresh) ok=0;
		}
		else if (it==(*ntemp)-1)
		{
			if (z[it]<=hthresh/(*ntemp) || prob_up[it]<counts[it]*pthresh ) ok=0;
		}

		// count how many holes
		if (ok==0)
		{
			nbadh ++;
		}
		else
		{
			ngoodh ++;
		}
		if (debug>0) fprintf(stderr,"z(T=%lf)=%lf\t",temp[it],z[it]);
	}

	if (debug>0) fprintf(stderr,"\nTo define a hole it uses the threshold %lf on z and %lf on prob_up/down\n",hthresh/ngoodh,pthresh);


	// if some histogram is empty
	if (nbadh>0)
	{
		lastt = temp[(*ntemp)-1];

		// reset temperatures and weights to the last reliable
		for (it=1;it<*nreliable_t;it++)
		{
			temp[it] = reliable_t[it];
			g[it] = reliable_g[it];
		}

		// get lowest temperature closer (the lowest temperature of the last (failed) run is indeed the lowest temperature)
		temp[(*nreliable_t)] = ( lastt + temp[(*nreliable_t)-1] ) / 2;

		if (*nreliable_t>1)
			g[(*nreliable_t)] = ExtrapolateLinear(temp[(*nreliable_t)],temp[(*nreliable_t)-1],g[(*nreliable_t)-1],
				temp[(*nreliable_t)-2],g[(*nreliable_t)-2]);
		else g[(*nreliable_t)] = g[(*nreliable_t)-1];

		*ntemp = *nreliable_t + 1;

			if (debug>0)
			{
				fprintf(stderr,"Resetting to last reliable temperatures and getting closer the lowest:\n");
				for (j=0;j<*ntemp;j++) fprintf(stderr,"%lf (g=%lf)\t",temp[j],g[j]);
				fprintf(stderr,"\nBACK TO SIMULATION\n\n");
			}

			return 0;
	}

	return 1;
}

double ExtrapolateLinear(double newt, double t1, double g1, double t2, double g2)
{
	double m,g0;

	m = (g1-g2) / (t1-t2);
	g0 = g2 - m*t2;

	return m * newt + g0;
}

/*****************************************************************************
 Same as CheckHistogramAdjust, but the new weights are obtained by extrapolation
 if yes, returns 1
 if there are holes among temperatures returns 0
 (not used) returns -1 to carry on the same simulation
 *****************************************************************************/
int CheckHistogramsAdjust3(int *counts, double *temp, double *g, int *ntemp, int nbin, double hthresh, int debug, int ntmax,
		double **h, double *prob_up, double *prob_down, double pthresh)
{
	int it,i,j,hole[NTEMPMAX],nbadh=0,nholes=0,ngoodh=0,ok;
	double z[NTEMPMAX],ztot=0;

	// obtain data from *counts
	for (it=0;it<*ntemp;it++)
	{
		z[it] = (double) counts[it];
		ztot += (double) counts[it];
	}

	// normalize counts over all temperatures
	if (debug>0) fprintf(stderr,"Check histogram counts: (threshold=%lf)\n",hthresh/(*ntemp));
	for (it=0;it<*ntemp;it++)
	{
		z[it] /= ztot;
		hole[it]=0;
		ok=1;

		// locate holes according to z and p_up/p_down
		if (it==0)
		{
			if (z[it]<=hthresh/(*ntemp) || prob_down[it]<counts[it]*pthresh) ok=0;
		}
		else if (it>0 && it<(*ntemp)-1)
		{
			if (z[it]<=hthresh/(*ntemp) || prob_up[it]<counts[it]*pthresh || prob_down[it]<counts[it]*pthresh) ok=0;
		}
		else if (it==(*ntemp)-1)
		{
			if (z[it]<=hthresh/(*ntemp) || prob_up[it]<counts[it]*pthresh ) ok=0;
		}

		// count how many holes
		if (ok==0)
		{
			nbadh ++;
		}
		else
		{
			ngoodh ++;
		}
		if (debug>0) fprintf(stderr,"z(T=%lf)=%lf\t",temp[it],z[it]);
	}

	if (debug>0) fprintf(stderr,"\nTo define a hole it uses the threshold %lf on z and %lf on prob_up/down\n",hthresh/ngoodh,pthresh);


	// if some histogram is empty
	if (nbadh>0)
	{
		// find holes between temperatures (calculating z only on the visited replica)
		for (it=1;it<*ntemp;it++)
			if ( (z[it] >= hthresh/ngoodh && z[it-1] < hthresh/ngoodh) ||
					(z[it-1] >= hthresh/ngoodh && z[it] < hthresh/ngoodh) ||
					prob_up[it]<counts[it]*pthresh || prob_down[it-1]<counts[it-1]*pthresh)
			{
				hole[nholes] = it;
				nholes ++;

				if (debug>0) fprintf(stderr,"There is a hole between T=%lf and T=%lf\n",temp[it-1],temp[it]);
			}


			for (i=0;i<nholes;i++)
			{
				temp[*ntemp] = (temp[hole[i]-1] + temp[hole[i]])/2.;

				g[*ntemp] = ExtrapolateLinear(temp[*ntemp],temp[hole[i]-1],g[hole[i]-1],temp[hole[i]],g[hole[i]]);

				(*ntemp) ++;
				if (*ntemp>ntmax)
				{
					fprintf(stderr,"It cannot add more temperatures because ntmax is reached (ntemp=%d)\n.ENDING SIMULATION\n",*ntemp);
					exit(0);
				}
				if (debug>0) fprintf(stderr,"Adding new temperature %lf (g=%lf)\n",temp[(*ntemp)-1],g[(*ntemp)-1]);
			}

			*ntemp = OrderTemperatures(*ntemp,temp,g,prob_up,prob_down,h,nbin,debug,counts);

			if (debug>0)
			{
				fprintf(stderr,"Retrying with new temperatures:\n");
				for (j=0;j<*ntemp;j++) fprintf(stderr,"%lf (g=%lf)\t",temp[j],g[j]);
				fprintf(stderr,"\nBACK TO SIMULATION\n\n");
			}

			return 0;
	}

	return 1;
}

/*****************************************************************************
 Check if all temperatures are sampled;
 if yes, returns 1
 if there are holes among temperatures returns 0
 (not used) returns -1 to carry on the same simulation
 *****************************************************************************/
int CheckHistogramsAdjust(int *counts, double *temp, double *g, double *reliable_lg, int *ntemp, double emin, double ebin,
		int nbin, double hthresh, int debug, int ntmax, double **h, double *prob_up, double *prob_down, double pthresh, int *binok, int gmethod)
{
	int it,i,j,hole[NTEMPMAX],nbadh=0,nholes=0,ngoodh=0,ok;
	double z[NTEMPMAX],ztot=0;

	// obtain data from *counts
	for (it=0;it<*ntemp;it++)
	{
		z[it] = (double) counts[it];
		ztot += (double) counts[it];
	}

	// normalize counts over all temperatures
	if (debug>0) fprintf(stderr,"Check histogram counts: (threshold=%lf)\n",hthresh/(*ntemp));
	for (it=0;it<*ntemp;it++)
	{
		z[it] /= ztot;
		hole[it]=0;
		ok=1;

		// locate holes according to z and p_up/p_down
		if (it==0)
		{
			if (z[it]<=hthresh/(*ntemp) || prob_down[it]<counts[it]*pthresh) ok=0;
		}
		else if (it>0 && it<(*ntemp)-1)
		{
			if (z[it]<=hthresh/(*ntemp) || prob_up[it]<counts[it]*pthresh || prob_down[it]<counts[it]*pthresh) ok=0;
		}
		else if (it==(*ntemp)-1)
		{
			if (z[it]<=hthresh/(*ntemp) || prob_up[it]<counts[it]*pthresh ) ok=0;
		}

		// count how many holes
		if (ok==0)
		{
			nbadh ++;
		}
		else
		{
			ngoodh ++;
		}
		if (debug>0) fprintf(stderr,"z(T=%lf)=%lf\t",temp[it],z[it]);
	}

	if (debug>0) fprintf(stderr,"\nTo define a hole it uses the threshold %lf on z and %lf on prob_up/down\n",hthresh/ngoodh,pthresh);


	// if some histogram is empty
	if (nbadh>0)
	{
		// find holes between temperatures (calculating z only on the visited replica)
		for (it=1;it<*ntemp;it++)
			if ( (z[it] >= hthresh/ngoodh && z[it-1] < hthresh/ngoodh) ||
					(z[it-1] >= hthresh/ngoodh && z[it] < hthresh/ngoodh) ||
					prob_up[it]<counts[it]*pthresh || prob_down[it-1]<counts[it-1]*pthresh)
			{
				hole[nholes] = it;
				nholes ++;

				if (debug>0) fprintf(stderr,"There is a hole between T=%lf and T=%lf\n",temp[it-1],temp[it]);
			}

		for (i=0;i<nholes;i++)
		{
				temp[*ntemp] = (temp[hole[i]-1] + temp[hole[i]])/2.;
				g[*ntemp] = CalculateG(reliable_lg, emin, ebin, nbin, temp[*ntemp], KB, binok, gmethod) -
						CalculateG(reliable_lg, emin, ebin, nbin, temp[0], KB, binok, gmethod);
				(*ntemp) ++;

				if (*ntemp>ntmax)
				{
					fprintf(stderr,"It cannot add more temperatures because ntmax is reached (ntemp=%d)\n.ENDING SIMULATION\n",*ntemp);
					exit(0);
				}
				if (debug>0) fprintf(stderr,"Adding new temperature %lf (g=%lf)\n",temp[(*ntemp)-1],g[(*ntemp)-1]);
		}

		*ntemp = OrderTemperatures(*ntemp,temp,g,prob_up,prob_down,h,nbin,debug,counts);


		if (debug>0)
		{
				fprintf(stderr,"Retrying with new temperatures:\n");
				for (j=0;j<*ntemp;j++) fprintf(stderr,"%lf (g=%lf)\t",temp[j],g[j]);
				fprintf(stderr,"\nBACK TO SIMULATION\n\n");
		}

			return 0;
	}

	return 1;
}

void ManageRestart(struct st_stuff *p, double step)
{
	int i,j;
	FILE *fr;

	// shift the stack by one, loosing the last
	for (i=NRESTARTS-1;i>0;i--)
		if ( ((p->st_st_restart)+i-1)->step >-1 )
		{
			((p->st_st_restart)+i)->step = ((p->st_st_restart)+i-1)->step;
			((p->st_st_restart)+i)->ntemp = ((p->st_st_restart)+i-1)->ntemp;
			for (j=0;j<NTEMPMAX;j++)
			{
				(((p->st_st_restart)+i)->temp)[j] = (((p->st_st_restart)+i-1)->temp)[j];
				(((p->st_st_restart)+i)->g)[j] = (((p->st_st_restart)+i-1)->g)[j];
			}
			for (j=0;j<p->st_nbin;j++)
				(((p->st_st_restart)+i)->lg)[j] = (((p->st_st_restart)+i-1)->lg)[j];
		}

	// write the current data on the 0-th element of the stack
	((p->st_st_restart)+0)->step = step;
	((p->st_st_restart)+0)->ntemp = p->st_nreliable_t;
	for (j=0;j<NTEMPMAX;j++)
	{
		(((p->st_st_restart)+0)->temp)[j] = p->st_reliable_t[j];
		(((p->st_st_restart)+0)->g)[j] = p->st_reliable_g[j];
	}
	for (j=0;j<p->st_nbin;j++)
		(((p->st_st_restart)+0)->lg)[j] = p->st_reliable_lg[j];

	// print restart
	fr = fopen("RESTART_ST","w");
	if (!fr) FatalError("Cannot open restart file for writing");

	for (i=0;i<NRESTARTS;i++)
		if ( ((p->st_st_restart)+i)->step >-1 )
		{
			fprintf(fr,"restart %d\n",i);
			fprintf(fr,"%lf\t%d\n",((p->st_st_restart)+i)->step,((p->st_st_restart)+i)->ntemp);
			for (j=0;j<((p->st_st_restart)+i)->ntemp;j++) fprintf(fr,"%lf\t%lf\n",(((p->st_st_restart)+i)->temp)[j],(((p->st_st_restart)+i)->g)[j]);
			fprintf(fr,"%lf %lf %d\n",p->st_emin,p->st_ebin,p->st_nbin);
			for (j=0;j<p->st_nbin;j++) fprintf(fr,"%lf\n",(((p->st_st_restart)+i)->lg)[j]);
		}

	fclose(fr);
}

void FilterHistograms(double **h, double **hout, int ntemp, int nbin, double binthresh)
{
	int i,it;

	for (it=0;it<ntemp;it++)
		for (i=0;i<nbin;i++)
		{
			if (h[it][i]>binthresh) hout[it][i] = h[it][i];
			else hout[it][i] = 0;
		}
}

double CalculateG(double *lg, double emin, double ebin, int nbin, double temp, double kb, int *binok, int method)
{
	if (method==0) return FreeEnergy(lg, emin, ebin, nbin, temp, KB, binok) / (KB * temp);
	else if (method==1) return AverageEnergy(lg, emin, ebin, nbin, temp, KB, binok) / (KB * temp);
	else if (method==2) return 0;
	FatalError("I don't know how to calculate g");
	return -1;
}
