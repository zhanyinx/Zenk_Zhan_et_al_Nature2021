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
 * do_mhistogram.h
 *
 *  Created on: Nov 18, 2009
 *      Author: guido
 */

#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>

  #ifndef NTEMPMAX
  #define NTEMPMAX 		1000	// Maximum number of histograms to be read
  #endif
#ifndef NEBINMAX
#define NEBINMAX 		10001
#endif
#ifndef	EPSILON
#define EPSILON			0.000000001	// a small number
#endif

#define MHE_RESIDUAL	0.00001	// Residual required for the solution of MHE equations
#define MHE_RESIDUAL2	0.001	// Residual required for the solution of MHE equations for savings
#define NITERHIS		5000

#define GSL_DER			// use explicit-derivative algorithm to solve MHE equations

struct rparams
{
	double **h;
	double *T;
	double emin;
	double ebin;
	int ntemp;
	int nbin;
	double kb;
	int debug;
	int paranoid;
	int *useit;
	int *binok;
};


double *MHistogram(double **h, double *t, int ntemp, int nbin, double emin, double ebin,
		int debug, double kb,double *, double *, int paranoid, int ignoreb, FILE *ff, double deltaT, int *binok);
void DensityOfStates(double **h, double *t, double *f, double emin, double ebin, int ntemp,
		int nbin, double *g, int debug, double kb, int *useit, int *binok);
int SolveEquations(double **h, double *t, double *f, double emin, double ebin, int ntemp, int nbin,
		int debug, double kb, int paranoid, int *useit, int *binok);
void PrintThermodynamics(FILE *fout, double **in, double tmin, double tbin, int nt);
void CalculateThermodynamics(double *g, double emin, double ebin, int nbin,
								double tmin, double tbin, int nt, double **out, double kb, double *temp, double **prob, int *binok);
int GMHequation(const gsl_vector *, void *, gsl_vector *);
#ifdef GSL_DER
void print_state (size_t iter, gsl_multiroot_fdfsolver * s, int ntemp, int *useit);
#else
void print_state (size_t iter, gsl_multiroot_fsolver * s, int ntemp);
#endif
void FatalError(char *text);
double FreeEnergy(double *lg, double emin, double ebin, int nbin, double temp, double kb, int *binok);
double SigmaEnergy(double *lg, double emin, double ebin, int nbin, double temp, double kb, int *binok);
double AverageEnergy(double *lg, double emin, double ebin, int nbin, double temp, double kb, int *binok);
int GMHequation_df (const gsl_vector *ff, void *params, gsl_matrix *J);
int GMHequation_fdf (const gsl_vector *ff, void *params, gsl_vector *eqq, gsl_matrix *J);


// memory.c
int **AlloIntMat(int l, int m);
double **AlloDoubleMat(int l, int m);
double *AlloDouble(int l);
int *AlloInt(int l);
