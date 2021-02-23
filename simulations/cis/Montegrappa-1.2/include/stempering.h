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
 * adjust_st.h
 *
 *  Created on: Nov 24, 2009
 *      Author: guido
 */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_multimin.h>


#define KB              1.0


#define NTEMPMAX	1000
#define NTBIN           200             // # of bins to invert <E(T)> in AddNewTemperature
#define CONVERGENCE     1E-9    // convergence criterion for probability maximization
#define CONVER_WEAK 1E-7        // after NITERMAX, accept if this weaker condition is met
#define NITERMAX        5000    // maximum number of iterations       "
#define NHISTOMAX       150             // maximum number of histograms to remember
#define NEBINMAX        10001   // maximum number of energy bins
#define NRESTARTS       10              // number of last snapshots to restart

#include "do_mhistogram.h"	//must be after #include "define.h" (if it is included)
//#define NEEDMEMORY

        struct st_restart
        {
        double step;
        int ntemp;
        double *temp;
        double *g;
        double *lg;
        };


        struct st_stuff
        {
        char st_method[15];
        int st_nm;					// shortening for method (1=st, 2=adaptive)
        int st_ntempmax;
        int st_ntemp;
        int st_nstep;			// attempts to change temperature every nstep
        int st_nstep_adj;			// adjust temperatures and weights every nstep_adj
        int st_nprint;
        int st_npre;				// equilibration state before everything
        double st_emin;
        double st_emax;
        double st_ebin;
        int st_nbin;
        char st_nftout[50];
        double st_k;				// new temperature put at k-times the standard deviation of the lowest
        double st_pthresh;			// probability threshold used to decrease the number of temperatures
        int st_debug;
        int st_keepall;			// keep all history of histograms
        char st_nfthe[50];			// output file for thermodynamics
        char st_nfdos[50];			// output file for density of states
        char st_nfdumb[50];		// output file for reliable dumb averages
        char st_nfdumb2[50];		// output file for current dumb averages
        char st_nfhisto[50];		// histogram file
        int st_nprintt;			// print thermodynamics every nprintt steps
        int st_paranoid;			// 0=give only warnings when something goes wrong, 1=stop
        double st_hthresh;			// threshold on fraction of time in a single temperature
        char st_anneal[15];		// how to set lower temperature (sigma/prob)
        double st_p_new;			// wished exchange probability
        int st_restart;			// if to restart
        double st_binthresh;		// threshold on (normalized) content of single bin
        double st_tonlyadd;		// below this temperature, only add temperatures
        int st_nkeepold;			// keep only last histograms
        int st_sum;				// 1=in keepall, sum together histograms with same temperature
        int st_removet;			// 0=do not remove temperatures, 99=any removal, ?=maximum decrease of ntemp in a shot
        double st_tstop;			// stop the simulation when this temperature is reached
        double st_tnorm;			// when this temperature is reached do nrmal stempering
        int st_nfail1;				// number of failures to fill holes with data from lg
        int st_nfail2;				// number of failures to fill holes with extrapolated free energies


         double st_phthresh;                // threshold on prob_up and prob_down to accept sampling
        int st_ignoreb;                    // in mhistogram ignore temperatures which prevent equation solving
        double st_deltat;                  // in mhistogram discard histograms which are closer than deltat in temperature
        int st_gmethod;                    // how to calculate optimal weights g   

        double *st_temp;                   // temperatures
        double *st_g;                      // weigths
        int st_itemp;                      // actual temperature
        double *st_prob_up;
        double *st_prob_down;
        int *st_counts;
        FILE *st_ftout;
        double **st_h;                     // the histogram h[temp][e]
        double **st_htmp;                  // temporary histogram h[temp][e] for normalization purposes
        double *st_f;                      // log Z, only needed to mhistogram
        double *st_current_lg;             // current log of density of states
        double *st_reliable_lg;            // last reliable log of density of states
        double *st_reliable_t;             //      "       temperatures
        double *st_reliable_g;             //      "       weights
        int st_nreliable_t;
        int st_count_st;                   // counter to attempt temperature change
        int st_count_adj;                  // counter to adjust temperatures/weights
        int st_count_print;                // counter to print temperatures
        int st_count_printt;               // counter to print thermodynamics
        int st_iter;                       // current iteration of adjust_st
        int *st_oldt_iter;                 // label each old histo with when it was collected
        int st_failure;                    // how many consecutive failed runs
        double st_energy;                  // potential energy of the system
        int *st_binok;                     // if to use a given energy bin (set in MHistogram)

        double **st_oldh;                  // past histograms
        double *st_oldt;                   // related temperatures
        int st_noldt;                              // how many old histograms
        int st_noadjust;                   // if 1, do not adjust temperatures
        double **st_out;                   // thermodynamics output
        double **st_boltzp;                // boltzmann distribution
        FILE *st_fpacc;
        double st_ttarget;                 // below this temperature, set to this and make last run
        int st_ttarget_harvest;            // if 1, print trajectory
        int st_printpdb;
        FILE *st_pdbf;
        char st_pdbnf[500];
        struct st_restart *st_st_restart;
        };





// stempering.c
void FreeStempering(struct st_stuff *x);
int STempering(double energy, double step, struct st_stuff *p);
int DoSTempering(double energy, struct st_stuff *p);
int OrderTemperatures(int ntemp, double *temp, double *g, double *prob_up,
		double *prob_down, double **h, int nbin, int debug, int *counts);
void PrintAverageEnergies(FILE *fout, double **h, double *temp, int ntemp, double emin, double ebin, int nbin);
void PrintHistogram(char *filename, double **h, double *t, int ntemp, int nbin, double emin, double ebin);
void Restart(struct st_stuff *p);
int FindClosestT(double t, int ntemp, double *temp);
void PrintStatistics(struct st_stuff *p, double step);

// adjust_st.c
void ManageRestart(struct st_stuff *p, double step);
int AdjustSTempering(struct st_stuff *p, double step);
int AddNewTemperature(int ntemp, double *temp, double *g, double *lg, double emin, double ebin,
		int nbin, int debug, double k, int paranoid, int *binok, int gmethod);
void OptimalWeights(int ntemp, double *t, double *g, double emin, double ebin, int nbin, double *lg, int debug, int *binok, int gmethod);
double EstimatedJumpProb (const gsl_vector *v, void *params);
double MaximizeTempProb(int ntemp, double *temp, double *g, double emin, double ebin, int nbin, int debug, double *lg,
		double **boltzp, int *binok, int gmethod);
int ReduceTemperatures(int ntemp, double *temp, double *g, double emin, double ebin, int nbin, int debug, double *lg,
		double pthresh, int ntmax, int removet, double **boltzp, int *binok, int gmethod);
int AddNewTemperature2(int ntemp, double *temp, double *g, double *lg, double emin, double ebin, int nbin,
		int debug, double p_new, int paranoid, double **boltzp, int *binok, int gmethod);
double EstimatedJumpSingleProbability(double t1, double t2, double *lg, double emin, double ebin, int nbin,
		int debug, double **boltzp, int *binok, int gmethod);
int CheckHistogramsAdjust(int *counts, double *temp, double *g, double *reliable_lg, int *ntemp, double emin, double ebin,
		int nbin, double hthresh, int debug, int ntmax, double **h, double *prob_up, double *prob_down, double pthresh, int *binok,
		int gmethod);
void NormalizeHistograms(double **h, double **htmp, int ntemp, int nbin, double binthresh, int gmethod);
int CheckHistograms(int *counts, double *temp, int ntemp, double htresh, int debug);
void AddHistogramPile(int ntemp, double *temp, double **h, double **oldh, double *oldt, int *oldt_iter, int *noldt, int nbin, int n, int iter, int debug, int sum);
int CheckHistogramsAdjust2(int *counts, double *temp, double *g, int *ntemp, double hthresh, int debug,
		double *reliable_t, double *reliable_g, int *nreliable_t, double *prob_up, double *prob_down, double phthresh);
double ExtrapolateLinear(double newt, double t1, double g1, double t2, double g2);
int CheckHistogramsAdjust3(int *counts, double *temp, double *g, int *ntemp, int nbin, double hthresh, int debug, int ntmax,
		double **h, double *prob_up, double *prob_down, double phthresh);
void FilterHistograms(double **h, double **hout, int ntemp, int nbin, double binthresh);
void CutHistograms(double **h, int ntemp, double emin, double ebin, int nbin, int *nbinnew, double *eminnew, int debug);
double CalculateG(double *lg, double emin, double ebin, int nbin, double temp, double kb, int *binok, int method);

// do_mhistogram.c
double *MHistogram(double **h, double *t, int ntemp, int nbin, double emin, double ebin, int debug, double kb, double *,
		double *, int paranoid, int ignoreb, FILE *ff, double deltaT, int *binok);

// memory1.c
double **AlloDoubleMat(int l, int m);
#ifdef NEEDMEMORY
double *AlloDouble(int l);
int *AlloInt(int l);
#endif
struct st_restart *AlloRestart(int ntemp, int nbin, int nres);
void FreeRestart(struct st_restart *x,int nres);
void FreeDoubleMat(double **x, int l);
