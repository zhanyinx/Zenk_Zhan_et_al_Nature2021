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
 * montegrappa.h
 *
 *  Created on: Sep 15, 2010
 *      Author: guido
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <signal.h>


        #define NVER		1
        #define NSUBVER	 	2	

        #define NRESMAX		1000		// maximum number of backbone atoms
        #define NATOMMAX	60000		// maximum number of atoms
        #define NSIDEMAX	15		// maximum number of sidechain atoms for a given backbone atom
        #define NROTMAX		40		// maximum number of rotamers allowed
        #define NCONTMAX	200		// maximum number of other atoms in contact with a given one (in polymer structure)
        #define NCONTMAX2	10000		// twice the number of other atoms in contact with a given one (in the moving structure, i.e. >>NCONTMAX)
        #define NSHELLMAX	1000		// maximum number of atoms in the shell of a given one
        #define NMOVES		11		// number of types of allowed moves
        #define NDIHFMAX	361		// maximum number of bins in tabled dihedral potential

        #define NREPMAX		40
        #define NAAMAX		1000
        #define NCHAINMAX	25


	//#define DEBUG_SHELL

        // tables to calculate fast trigonometrics, exp, etc.
        // #define TABLED_F			// pay attention when activate it !!!!
        					// square root (bin=0.0025)
        #define FSQRT_L		80000		// # of elements in the table
        #define FSQRT_BINRC	800.		// recirpocal of binning of the table
        #define FSQRT_MAX	100		// = L / BIN
        #define FSQRT_HBIN	0.000625	// half bin
        #define FTRIG_L		72000		// # of elements in the table
        #define FTRIG_BINRC	200.		// reciprocal of binning of the table
        #define FTRIG_HBIN	0.0025		// half bin
        #define FEXP_L		10000		// # of elements in the table
        #define FEXP_BINRC	1000.		// reciprocal of binning of the table
        #define FEXP_MAX	10		// = L / BIN
        #define FEXP_HBIN	0.0005		// half bin
        #define FACOS_L		40000		// # of elements in the table
        #define FACOS_BINRC	20000.		// reciprocal of binning of the table
        #define FACOS_HBIN	0.000025	// half bin

        //LM_DELTAA for BGS steps
        #define LM_DELTAA	0.001
	#define DELTAOMEGA	20.		//omega dihedral tolerance



        #define NSTART 4
        #define NCHECK 1

        //#define DEBUG

        #define PI 3.14159265
        #define EPSILON 0.0001
        #define EPSILON2 0.0001
        #define LARGE 9E9

#define OPTIMIZEPOT

#ifdef STEMPERING
#include "stempering.h"
#endif

#include "struct.h"
#include "struct_op.h"



#ifdef OPTIMIZEPOT
#include "optimizepot.h"
#endif



#ifdef ACTIVE_MPI

#include <mpi.h>
#include "MPIfunc.h"
#define buffer_max 5000000
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


#else
struct s_mpi_parms
{

};
#endif


 

void Welcome(FILE *fp);

// io.c
void Error(char *text);
void ReadPolymer(char *fname, struct s_polymer *p, FILE *flog, int npol, int debug, int *iamax, int *itypemax);
#ifdef STEMPERING
void ResetStStuff(struct s_mc_parms *x,char *fname);
#endif
int FindKeyword(char *string, char *keyword);
void PrintPolymer(char *fname, struct s_polymer *p, int nchains);
void PrintPDBStream(struct s_polymer *p, int npol, FILE *fp);
struct s_mc_parms *ReadMcParms(char *fname);
int ReadPotential(char *fname, struct s_potential *u, struct s_mc_parms *parms, int na, int ntype);
void SetLookbackTables(struct s_polymer *p, int nc);
void PrintVector(FILE *fp, char c[], struct vector x);
void PrintAngles(FILE *fp, char c[], struct angles x);
void PrintStructure(struct s_polymer *p, int npol, FILE *fp, int shell);
void ReadParD(char *s, char key[20], int *par);
void ReadParL(char *s, char key[20], long *par);
void ReadParF(char *s, char key[20], double *par);
void ReadParS(char *s, char key[20], char *par);
void ReadParN(char *s, char key[20], int *par);
void ReadParLLU(char *s, char key[20], unsigned long long *par);
void PrintPotential(struct s_potential *u, char *eoutfile, int nat, int ntypes, int noangpot, int nodihpot, int hb);
void ComputeIAPDB(struct s_polymer *p,struct s_mc_parms *mc_parms);


// memory.c
struct s_polymer *AlloPolymer(int npol, int n, int nside, int nrot, int natoms, int shell, int noside, FILE *flog);
struct s_tables *InitTables(FILE *fp);
struct s_potential *AlloPotential(int natoms, int ntypes, int noangpot, int nodihpot, int hb);
int **AlloIntMatrix(int l, int m);
int *AlloInt(int l);
double **AlloDoubleMatrix(int l, int m);
double *AlloDouble(int l);
 void FreePolymer(struct s_polymer *p,int npol, int nback, int nside, int shell, int noside);

 void FreeDoubleMatrix(double **x,int l);

void FreePotential(struct s_potential *x, int ntypes, int noangpot, int nodihpot, int hb);

void FreeTables(struct s_tables *t);


// mc.c
void Do_MC(struct s_polymer *p, struct s_polymer *fragment, struct s_polymer *replica, struct s_polymer *native, struct s_potential *pot, struct s_mc_parms *parms, FILE *ftrj, FILE *fe, struct s_polymer *oldp, FILE *fproc, int irun,struct s_mpi_parms *mpiparms);
void CopyResiduePositions(struct s_back *from, struct s_back *to);
void CopyResiduePositions_NOCONT(struct s_back *from, struct s_back *to);
int MoveBackboneFlip(struct s_polymer *p,struct s_polymer *oldp, struct s_mc_parms *mc_parms, struct s_potential *u,unsigned long long step, int debug, double t);
int MoveBackbonePivot(struct s_polymer *p,struct s_polymer *oldp, struct s_potential *pot, struct s_mc_parms *mc_parms, double t);
int MoveMultiplePivot(struct s_polymer *p, struct s_polymer *oldp, struct s_potential *pot, int nmul, struct s_mc_parms *mc_parms, double t);
int Metropolis(double deltaE, double T, struct s_tables *t);
void UpdateMonomer(struct s_polymer *from, struct s_polymer *to, int w, int n, int shell);
int MoveSidechain(struct s_polymer *p,struct s_polymer *oldp, struct s_mc_parms *mc_parms, struct s_potential *pot,
		unsigned long long istep, int debug, double t);
void UpdateMonomerRange(struct s_polymer *from, struct s_polymer *to, int wfrom, int wto, int p, int shell);
int MoveLoosePivot(struct s_polymer *p, struct s_polymer *oldp, struct s_potential *pot, int nmul, struct s_mc_parms *parms, double t);
void SoftExit();
int MoveMultipleFlip(struct s_polymer *p,struct s_polymer *oldp, struct s_mc_parms *parms, struct s_potential *pot,
		unsigned long long istep, int debug, double t);
int MoveCoM(struct s_polymer *p,struct s_polymer *oldp, struct s_mc_parms *parms, struct s_potential *pot,
		unsigned long long istep, int debug, double t);
int MoveRotation(struct s_polymer *p,struct s_polymer *oldp, struct s_mc_parms *parms, struct s_potential *pot, unsigned long long istep, int debug, double t);

int MoveClusterCoM(struct s_polymer *p,struct s_polymer *oldp, struct s_mc_parms *parms, struct s_potential *pot,int istep, int debug, double t, FILE *fproc, int iproc);
int MoveClusterRot(struct s_polymer *p,struct s_polymer *oldp, struct s_mc_parms *parms, struct s_potential *pot,int istep,int npol, int debug, double t, FILE *fproc, int iproc);



void Anneal(struct s_mc_parms *p, double *t, int *counter, int *status, int *ok, int *ishell, int *mcount);

//bias_mp.c
void InvertTriang( double **Inv, double **Mat, int n );
double Squared_n_Norma( double *vect, int dim );
void TransposedMatOnVect( double **Mat, double *vect, double *risu, int n);
int print_square_matrix(double **m,int n);
int print_vector(double *v,int n);
int Gaussian_Angles(double *angles,int n);
void MatA( double **A, double **G, int dim,double a,double b);
void Cholesky_2( double **L, double **A, int dim);


void CopyResidueCoordinates(struct s_back *from,struct s_back *to);
int ComputeG(double **g,struct s_polymer *fragment,struct s_polymer *p,int ip,int k,int n,int npol,struct s_mc_parms *parms);
int MoveBiasedGaussian(struct s_polymer *p,struct s_polymer *oldp,struct s_polymer *fragment,struct s_potential *pot,int nmul,struct s_mc_parms *mc_parms,double t);
int MoveReallyLocal(struct s_polymer *p,struct s_polymer *oldp,struct s_potential *pot,int nmul,struct s_mc_parms *mc_parms,double t);
int CopyFragment(struct s_polymer *p,struct s_polymer *f,int iw,int nmul,int natom_fragment,int ip);
int B_Metropolis(double deltaE,double T,double WN,double WD,struct s_tables *t);
//local_move.c
struct s_polymer *Allo_Fragment(int npol, int nback,int nang, FILE *flog);
int LocalMove(struct s_polymer *p, struct s_polymer *oldp,struct s_polymer *fragment,struct s_potential *pot,int nmul,struct s_mc_parms *parms, double t);
int Compute_G(struct s_polymer *fragment,struct s_polymer *p,int ip,int istart,int natom_fragment,int nang,struct s_mc_parms *parms);
void FreeFragment(struct s_polymer *f,struct s_mc_parms *parms);

//misc.c
int irand(int r);
double frand(void);
long Randomize(int n);
double FastSqrt(double x, struct s_tables *t);
double FastSin(double x, struct s_tables *t);
double FastCos(double x, struct s_tables *t);
double FastExp(double x, struct s_tables *t);
double FastAcos(double x, struct s_tables *t);
double Norm2(struct vector a);
void MatrixVectorProduct(double T[][3], struct vector in, struct vector *out);
double DAbs(double x);
int Abs(int x);
void CopyDoubleMatrix(double **from, double **to, int n, int m);
double gauss(double av, double sig);
void CopyVector(struct vector *from, struct vector *to);

// geometry.c
int Flip(struct s_polymer *p, int iw, double dw);
int MoveHead(struct s_polymer *p, struct s_mc_parms *mc_parms);
int MoveTail(struct s_polymer *p, struct s_mc_parms *mc_parms);
struct angles Cartesian2Spherical( struct vector A, struct vector B,
                                       struct vector C, struct vector D, struct s_tables *tables, int *out);
struct vector Spherical2Cartesian(struct vector A, struct vector B,
                                  struct vector C, struct angles W, struct s_tables *tables, int *out);
struct vector RotateVector(struct vector v, double theta, int w, struct s_tables *tables);
int PivotForward(struct s_polymer *p, int iw, double dw, int n, struct s_mc_parms *parms);
int PivotBackward(struct s_polymer *p, int iw, double dw, int n, struct s_mc_parms *parms);
int Pivot(struct s_polymer *p, int iw, double dw, int ip, struct s_mc_parms *parms);
double Dist2(struct vector a, struct vector b);
double Dist(struct vector a, struct vector b);
double Angle( struct vector B, struct vector C, struct vector D, struct s_tables *tables, int *out);
double Dihedral( struct vector A, struct vector B,
                                       struct vector C, struct vector D, struct s_tables *tables, int *out);
int AddSidechain(struct s_polymer *p, int istart, int istop, int ipol);
int FlipFragment(struct s_polymer *p, int ifrom, int ito, double dw);
void CopyPolymer(struct s_polymer *from, struct s_polymer *to, int cfrom, int cto, int noside, int norot);
void CopyAllPolymers(struct s_polymer *from, struct s_polymer *to, int n, int noside, int norot);
void DisplaceCoM(struct s_polymer *p, int ip, double dx, double dy, double dz);
double CosAngle( struct vector B, struct vector C, struct vector D, struct s_tables *tables, int *out);
void RotationX(struct s_polymer *p,int ip,double dtheta);
void RotationY(struct s_polymer *p,int ip,double dtheta);
void RotationZ(struct s_polymer *p,int ip,double dtheta);
void RotationClusterX(struct s_polymer *p,int ip,double dtheta,int icluster,int cluster[NCHAINMAX][NCHAINMAX],int npol_cluster);
void RotationClusterY(struct s_polymer *p,int ip,double dtheta,int icluster,int cluster[NCHAINMAX][NCHAINMAX],int npol_cluster);
void RotationClusterZ(struct s_polymer *p,int ip,double dtheta,int icluster,int cluster[NCHAINMAX][NCHAINMAX],int npol_cluster);

// potential.c
double TotalEnergy(struct s_polymer *p, struct s_potential *u, struct s_mc_parms *parms, int npol, int update, int sidechains, int debug, int iproc);
double EnergyMonomer(struct s_polymer *p, struct s_potential *u, int i, int ci, int npol, int update, int shell, int sidechains, int disentangle, int hb);
void AddContact(struct s_polymer *p, int i, int j, int ci, int cj, double e);
void AddShell(struct s_polymer *p, int i, int j, int ci, int cj);
void ResetContactsMonomer(struct s_polymer *p, int i, int ci);
double GetEnergyMonomer(struct s_polymer *p, int ip, int iw);
void PrintContacts(FILE *fp, struct s_polymer *p,int ip, unsigned long long step);
void CountContacts(FILE *fp,struct s_polymer *polymer,struct s_mc_parms *parms,unsigned long long step);
double EnergyMonomerShell(struct s_polymer *p, struct s_potential *u, int i,int ip, int update);
void UpdateShell(struct s_polymer *p, struct s_mc_parms *parms);
void CopyShell(struct s_polymer *from,struct s_polymer *to, struct s_mc_parms *parms);
double GetEnergyMonomerRange(struct s_polymer *p, int from, int to, int ip);
double EnergyMonomerRange(struct s_polymer *p, struct s_potential *u, int from, int to, int ip, int npol, int shell, int update, int nosidechains, int disentangle, int hb);
double EnergyPair(struct s_polymer *p, struct s_potential *u, int i, int j, int ci, int cj, int update, int sidechains, int disentangle, int tooclose, int hb);
double EnergyAngles(struct s_polymer *p, struct s_potential *u, int iw, int ic, int update);
double EnergyDihedrals(struct s_polymer *p, struct s_potential *u, int iw, int ic, int update);
void PrintEnergies(FILE *fp, int nc, unsigned long long step, struct s_polymer *p);
void PrintEnergies_Parallel(FILE *fp, int nc, unsigned long long step, struct s_polymer *p,int my_rank);
void CopyPotential(struct s_potential *from, struct s_potential *to, int nat, int ntypes);
int CheckOverlaps(struct s_polymer *p, struct s_potential *u, struct s_mc_parms *parms, int npol, int nosidechains, int pr,FILE *fproc);
void CompareStructures(struct s_polymer *a, struct s_polymer *b, int nc, int natoms);
int EnergyBox(struct s_polymer *p, struct s_potential *u, int iw, int ic);
int EnergyBoxPolymer(struct s_polymer *p, struct s_potential *u, int ic);


// constants for gaussian random generator
float ran1(long *idum);
double gasdev(long *idum);

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)
