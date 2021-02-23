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
 * optimizepot.h
 *
 *  Created on: Jan 31, 2011
 *      Author: guido
 */



#define OP_NCONTMAX	500	// maximum number of contacts
#define OP_NRESMAX	100	// maximum number of restrains
#define REALLOCATE		// if the number of contacts exceed OP_NCONTMAX, reallocate the structure
//#define OP_DEBUG

struct s_optimizepot *InitializeOptimizePot(struct s_mc_parms *parms, int ntypes, struct s_potential *u, FILE *fproc, int iproc, FILE *fp, int *nrestr);
struct s_optimizepot *AlloOptimizePot(struct s_mc_parms *parms, int ntypes, int nrestr, FILE *fproc);
struct s_optimizepot_input *ReadOPRestrains(struct s_mc_parms *parms);
void OP_GetRestrain(int it, struct s_polymer *p, int ipol, struct s_optimizepot_input *op_input);
void OP_AddEnergy(struct s_polymer *p, int a1, int a2, double mul);
double Chi2(double *x, double *xexp, double *sigma, int n);
double OP_function(double **e, struct s_optimizepot *x, struct s_optimizepot_input *in, struct s_mc_parms *parms);
void OP_SamplePotential(struct s_polymer *p, struct s_mc_parms *parms, int ntypes, double **emat, int irun, double **ematold);
double OP_functionDiff(int iw, double deltae, struct s_optimizepot *x, struct s_optimizepot_input *in, struct s_mc_parms *parms);
void FreeOpt(struct s_optimizepot *x,struct s_mc_parms *parms,int nrestr);
void FreeOptInput(struct s_optimizepot_input *x);
