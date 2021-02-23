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

#include <mpi.h>
//MPI Types

void Create_rot_datatype(MPI_Datatype *Rottype);
void Create_side_datatype(MPI_Datatype *Sidetype);
void Create_back_datatype(MPI_Datatype *Backtype);
void Create_pot_datatype(MPI_Datatype *Pottype);
void Create_parms_datatype(MPI_Datatype *Parmstype);

/*
struct s_mpi_parms
{
	int my_rank;
	int nprocs;
	MPI_Datatype Rot_mpi;
	MPI_Datatype Parms_mpi;
	MPI_Datatype Back_mpi;
	MPI_Datatype Side_mpi;
	MPI_Datatype Pot_mpi;
	MPI_Status astatus;

};
*/

//Send-Receive
struct s_mc_parms *send_parms(int iproc, int nprocs, MPI_Datatype Parmstype, struct s_mc_parms *parms, MPI_Status astatus);//,double *tmin);
void send_struct(int *nback, int iproc, int nprocs, int *nat, int *ntypes, MPI_Status astatus);
struct s_polymer *send_pol(int iproc, int nprocs, int nback, MPI_Datatype Backtype, MPI_Datatype Sidetype,  MPI_Datatype Rottype, struct s_polymer *startp, MPI_Status astatus, int npol, int shell, int nosidechains);
struct s_potential *send_pot(int nat, int ntypes, int noangpot, int nodihpot, int hb, int iproc, int nprocs, MPI_Datatype Pottype, struct s_potential *u, MPI_Status astatus);
void send_double_matrix(int length1, int length2, int iproc, double **m, int source);
void send_int_matrix(int length1, int length2, int iproc, int **m, int source);

//exchange
int ExchangePol(struct s_polymer *polymer, struct s_polymer *replica, struct s_polymer *oldp, struct s_mc_parms *parms, struct s_potential *u, int iproc, int ntemp, int even, int *ex_count, int *ex_acc, MPI_Datatype Backtype, MPI_Datatype Sidetype, MPI_Datatype Rottype, MPI_Status astatus,unsigned long long istep);
//Potential-related
#ifdef OPTIMIZEPOT
struct s_optimizepot_input *Allo_op_input(int nrestr);
struct s_optimizepot_input *send_op_input(int nrestr,struct s_optimizepot_input *in);
#endif
