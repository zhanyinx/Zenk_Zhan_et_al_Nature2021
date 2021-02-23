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
 * struct.h
 *
 *  Created on: Sep 15, 2010
 *      Author: guido
 */
//#include "stempering.h"




struct vector { double x,y,z; };

struct angles { double ang;
                double dih;
                double r;
              };


struct s_polymer
{
	int nback;			// number of backbone atoms
	char title[80];
	struct s_back *back;
	struct vector **vback;		// lookback table (given an atom, gives the address of the cartesian position; set by ReadPolymer)
	struct s_tables *tables;	// tables of sin,cos,etc. ( needs InitTables before use)
	double etot;			// total energy
	double t;			// its temperature (necessary for op)
      double **A,**G,**L,**Y;      //local move
      double *g_ang,*d_ang;        //loval move



	#ifdef OPTIMIZEPOT
	struct s_optimizepot *op;
	#endif
};

struct s_back
{
	int ia;
	int iapdb;			//->WARNING: see ComputeIAPDB and BGS
	int itype;
	char aa[5];
	char type[5];
	int iaa;
	struct vector pos;		// cartesian position
	struct angles sph;		// position in spherical coordinates based on the 3 preceding backbone atoms
	int nside;			// number of atoms in its sidechain (0=none)
	int move;			// 0 means that it should't be moved
	int irot;			// current rotamer id
	int nrot;			// total number of rotamers

	struct s_side *side;		// sidechain atoms

	int ncontacts;			// number of backbone atoms whose atom (including sidechain) are in contact with this
	int *contacts;			// list of atoms in contact
	int *contacts_p;		// ... which chain it belongs to
	double *e;			// interaction energy of this with the others residues
	double e_ang;			// angular energy
	double e_dih;			// dihedral energy
	int nshell;			// number of backbone atoms in the shell
	int *shell;			// atoms in the shell
	int *shell_p;			// ... which chain it belongs to

	double d2_next;			// square distance to the next backbone
	double a_next;			// angle
	double d_next;			// dihedral
};

struct s_side
{
	int ia;
	int iapdb;
	int itype;
	char type[5];
	struct vector pos;		// cartesian position

	struct s_rotamers *rot;		// rotamers
};

struct s_rotamers
{
	int b1;						// the atoms which act as reference state
	int b2;
	int b3;
	struct angles ang;				// the spherical coordinate of the sidechain atom
};

struct s_potential
{
	double **e;					// energy of the square well
	double **r_2;					// square width of the square well
	double **r0_2;					// square width of hard core repulsion

	double *e_ang;					// energy constant of angular potential
	double *ang0;					// equilibrium angle

	int dih_periodic;				// 1=activate periodic dihedrals
	double *e_dih1;					// energy constant of dihedral potential with multeplicity 1
	double *dih01;					// phase shift of diherdrals with multeplicity 1
	double *e_dih3;					// energy constant of dihedral potential with multeplicity 3
	double *dih03;					// phase shift of diherdrals with multeplicity 3

	int dih_tabled;					// 1=activate tabled dihedrals
	int *dih_which;					// 0=phi, 1=psi
	double *dih_pa;					// weight of alpha-helix in tabled dihedral
	double *dih_pb;					// weight of beta-sheet in tabled dihedrals
	double *dih_f_psi_a;				// terms of the potential
	double *dih_f_phi_a;
	double *dih_f_psi_b;
	double *dih_f_phi_b;

	double g_r0hard;				// global hardcore repulsion
	double g_ehomo;					// global homopolymer e
	double g_rhomo;
	double g_anglek;
	double g_angle0;
	double g_dihk1;
	double g_dih01;
	double g_dihk3;
	double g_dih03;
	int g_imin;					// imin
	double g_dihbin;				// bin of tabled dihedral potential
	double g_dihe;					// energy for tabled dihedral potential
	char boxtype;					// n=none, c=cubic, s=spherical
	double boxsize;

	double lbox;

	double kr2_splice;				// splice energy well into two parts, the distance being kr * r
	double ke_splice;				// and the depth is ke * e
	int splice;						// =1 to activate

	int *hc_type;					// hardcore parameters
	double *hc_r0;
	int hc_number;

	int *hb;					// hydrogen bond: 0=none, 1=donor, 2=acceptor
	int *hb_iam;					// ia of the atom preceding donor/acceptor

	double **sigma;				//energy dihedrals with minimum in Ramachandran dihedrals
	int **dih0;
	int dih_ram;
	double e_dihram;
	double **ab_propensity;
};

////////struct st_restart
////////{
////////double step;
////////int ntemp;
////////double *temp;
////////double *g;
////////double *lg;
////////};


////////struct st_stuff
////////{
////////char st_method[15];
////////int st_nm;					// shortening for method (1=st, 2=adaptive)
////////int st_ntempmax;
////////int st_ntemp;
////////int st_nstep;			// attempts to change temperature every nstep
////////int st_nstep_adj;			// adjust temperatures and weights every nstep_adj
////////int st_nprint;
////////int st_npre;				// equilibration state before everything
////////double st_emin;
////////double st_emax;
////////double st_ebin;
////////int st_nbin;
////////char st_nftout[50];
////////double st_k;				// new temperature put at k-times the standard deviation of the lowest
////////double st_pthresh;			// probability threshold used to decrease the number of temperatures
////////int st_debug;
////////int st_keepall;			// keep all history of histograms
////////char st_nfthe[50];			// output file for thermodynamics
////////char st_nfdos[50];			// output file for density of states
////////char st_nfdumb[50];		// output file for reliable dumb averages
////////char st_nfdumb2[50];		// output file for current dumb averages
////////char st_nfhisto[50];		// histogram file
////////int st_nprintt;			// print thermodynamics every nprintt steps
////////int st_paranoid;			// 0=give only warnings when something goes wrong, 1=stop
////////double st_hthresh;			// threshold on fraction of time in a single temperature
////////char st_anneal[15];		// how to set lower temperature (sigma/prob)
////////double st_p_new;			// wished exchange probability
////////int st_restart;			// if to restart
////////double st_binthresh;		// threshold on (normalized) content of single bin
////////double st_tonlyadd;		// below this temperature, only add temperatures
////////int st_nkeepold;			// keep only last histograms
////////int st_sum;				// 1=in keepall, sum together histograms with same temperature
////////int st_removet;			// 0=do not remove temperatures, 99=any removal, ?=maximum decrease of ntemp in a shot
////////double st_tstop;			// stop the simulation when this temperature is reached
////////double st_tnorm;			// when this temperature is reached do nrmal stempering
////////int st_nfail1;				// number of failures to fill holes with data from lg
////////int st_nfail2;				// number of failures to fill holes with extrapolated free energies
////////	 double st_phthresh;                // threshold on prob_up and prob_down to accept sampling
////////int st_ignoreb;                    // in mhistogram ignore temperatures which prevent equation solving
////////double st_deltat;                  // in mhistogram discard histograms which are closer than deltat in temperature
////////int st_gmethod;                    // how to calculate optimal weights g   

////////double *st_temp;                   // temperatures
////////double *st_g;                      // weigths
////////int st_itemp;                      // actual temperature
////////double *st_prob_up;
////////double *st_prob_down;
////////int *st_counts;
////////FILE *st_ftout;
////////double **st_h;                     // the histogram h[temp][e]
////////double **st_htmp;                  // temporary histogram h[temp][e] for normalization purposes
////////double *st_f;                      // log Z, only needed to mhistogram
////////double *st_current_lg;             // current log of density of states
////////double *st_reliable_lg;            // last reliable log of density of states
////////double *st_reliable_t;             //      "       temperatures
////////double *st_reliable_g;             //      "       weights
////////int st_nreliable_t;
////////int st_count_st;                   // counter to attempt temperature change
////////int st_count_adj;                  // counter to adjust temperatures/weights
////////int st_count_print;                // counter to print temperatures
////////int st_count_printt;               // counter to print thermodynamics
////////int st_iter;                       // current iteration of adjust_st
////////int *st_oldt_iter;                 // label each old histo with when it was collected
////////int st_failure;                    // how many consecutive failed runs
////////double st_energy;                  // potential energy of the system
////////int *st_binok;                     // if to use a given energy bin (set in MHistogram)

////////double **st_oldh;                  // past histograms
////////double *st_oldt;                   // related temperatures
////////int st_noldt;                              // how many old histograms
////////int st_noadjust;                   // if 1, do not adjust temperatures
////////double **st_out;                   // thermodynamics output
////////double **st_boltzp;                // boltzmann distribution
////////FILE *st_fpacc;
////////double st_ttarget;                 // below this temperature, set to this and make last run
////////int st_ttarget_harvest;            // if 1, print trajectory
////////int st_printpdb;
////////FILE *st_pdbf;
////////char st_pdbnf[500];
////////struct st_restart *st_st_restart;
////////};









struct s_mc_parms
{
	int npol;		// number of chains
	unsigned long long nstep;
	long seed;
	double dw_flip;		// maximum angle allowed in a flip
	double dw_pivot;	// maximum angle allowed in a flip
	double dw_mpivot;	// maximum angle allowed in a flip
	double dw_lpivot;	// maximum angle allowed in a flip
	double dw_mflip;	// maximum angle allowed in a multiple flip
	char fntrj[50];		// name of trajectory file
	char fne[50];		// name of energy file
	char flastp[50];	// name of last conformation;
	char fnproc[50];	// name of the log file
	int nprinttrj;		// when to print trajectory file
	int nprintlog;		// when to print log
	int nprinte;		// when to print energy
	FILE *flog;
	int shell;		// 1=activate shell
	int nshell;		// update shell every %d step
	double r2shell;		// square distance to be in the shell of a backbone atom
	int ntemp;		// number of temperatures (replicas)	
	#ifdef ACTIVE_MPI
	double T[NREPMAX];	// temperatures
	#else
	double T;
	#endif
	int randdw;		// 1=dw from flat distribution, 2=dw from gaussian distribution
	int debug;
	int movetype[NMOVES];	//activate(0/1) move: 0=flip, 1=pivot, 2=multipivot, 3=sidechain, 4=loose multipivot
	int nmul_mpivot;	// how many backbone atoms to move in multiple pivot move
	int nmul_lpivot;	// how many backbone atoms to move in loose pivot move
	int nmul_mflip;		// how many backbone atoms to move at most in multiple flip move
	int nosidechains;	//0=there are, 1=no sidechains
	int noangpot;		//0=there is a potential on angles, 1; there is not
	int nodihpot;		//0=there is a potential on dihedrals, 1; there is not
	int nrun;		// number of repetitions of the MC
	int always_restart;	//1=in different irun, start always from input structure
	int record_native; //1=records the input (native( structure as first snapshot
	int acc;		// # of accepted moves
	int mov;		// # of attempted moves
	int disentangle;	//0=if two atoms are closer than r_hardcore, then do not calculate the energy
	int stempering;		//1=do stempering
	double dx_com;		//displacemente oc center of mass
	double dx_clm;		//displacement of cluster
	double dtheta;
	double r_cloose;	// constrains in the bond distance,
	double a_cloose;	// angles
	double d_cloose;	// and dihedral, relative to the initial position (-1 to disable)
	int hb;			// activate hydrogen bonds
	#ifdef OPTIMIZEPOT
	 struct s_optimizepot_input *op_input;
	 char fnop[50];		// name of file of experimental restrains
	 char op_minim[50];	// optimize procedure: none, sample, steep
	 int op_itermax;	// # of iteration steps
	 double op_step;
	 double op_T;		// temperature corresponding to the experimental data (can be different from T)
	 int op_deltat;		// = nstep/op_nframes
	 double op_stop;	// stop condition
	 int op_print;		// print log every * step
	 double op_emin;	// matrix elements cannot go below this limit
	 double op_emax;	// matrix elements cannot go below this limit
	 int op_wait;		// discard the first steps
	 double op_r;		// default well width;
	 double op_r0;		// default well hardcore;
	#endif
	int anneal;		// 1=anneal,0=not 
	int anneal_often;	// every * istep 
	int anneal_step;	// carry out * step at higher temperature (* decreases to zero)
	double anneal_t;	// higher temperature
	int anneal_recov;	// do nothing for * steps after returning to actual temperature
        int nconf;              //number of conformation that want to save
	int nstep_exchange;	
	int nmul_local;

	double bgs_a;
        double bgs_b;


	int chi2start;
	
	int ishell;

	#ifdef STEMPERING
	struct st_stuff	*p;
	#endif

	int iT_bias; 	//if replica index > iT_bias, no local move -> no shells


};

struct s_tables
{
	double *fast_sqrt;
	double *fast_sin;
	double *fast_cos;
	double *fast_expp;		// exp of positive numbers
	double *fast_expm;		// exp of negative numbers
	double *fast_acos;		// arccos

};
