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



struct s_optimizepot
{
	int *it1;		// this itype...
	int *it2;		// ...interacts with itype
	double **mul;		// the interaction is multiplied by this geometric factor in each snapshot
	int ncontacts;		// this is the dimension of the arrays above
	int nallocont;		// size of it1, it2 and mul in terms of ncontacts, in terms of OP_NCONTMAX
	double *eold;		// the original energy of each recorded structure
	double *efix;		// the non-optimizible part of the energy (dihedrals, angles, etc.)

	double **restrain;	// the raw value of the restrains in each recorded conformation

				// the following are for internal calculations:
	double *rest_av;	// the calculated thermal average of restrains
	double *enew;		// new energies
	double *t;		// the temperature of the ith snapshot

	int nframes;		// number of recorded conformations
	int nframesmax;		// size of the allocated structure
	int icount;

	int record;		// 1= record the energy when calling EnergyPairs
};

struct s_optimizepot_input
{
	int ndata;		// number of restrain data
	double *expdata;
	double *sigma;		// experimental error (i.e. weight of the data)
	int *datatype;		// the kind of data (0=contact, 1=distance, 2=1/distance^6, 3=block contact)
	int *i1;
	int *i2;
	int *i3;
	int *i4;
};


