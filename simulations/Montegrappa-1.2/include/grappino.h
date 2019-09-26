/*
 * grappino.h
 *
 *  Created on: Sep 24, 2010
 *      Author: guido
 */

#define ROTMAX          40
#define ROT_ATMAX       20
#define AAKINDMAX		30
#define ROT_CB_ANG 		120.
#define ROT_CB_DIH 		180.
#define ROT_CB_RADIUS 	1.51
#define ROT_CB_ANG_PRO 		120.
#define ROT_CB_DIH_PRO 		180.
#define ROT_CB_RADIUS_PRO 	1.51
#define NPOLMAX				150

#define NTYPEMAX			200

struct s_parms
{
	char pdbfile[200];
	char poutfile[200];
	char eoutfile[200];
	char cntfile[200];
	char ab_propensityfile[200];

	char potential[50];
	char atomtypes[50];
	char model[50];
	int hydrogens;
	double rnat;
	double rhard;
	double kr_splice;			// double <1 to separate the box-potential in 2 (width)
	double ke_splice;			// double <1 to separate the box-potential in 2 (height)
	int splice;				// =1 to activate
	int use_nativedist;
	int nosidechain;			
	double k_native_r;
	double k_native_hc;
	char back_a[100][100];
    	int move_a[100];
	int n_back_a;
	double tthresh;				// threshold to define a bonded interaction
	int debug;
	int imin;				// minimum distance along the chain
	double goe;				// go energy
	double e_homo;				// energy of global homopolymeric interaction
	double r_homo;				// its radius
	int rotamers;				// how to deal with rotamers (0=off)
	int sidescratch;				// build sidechain from rotamers instead than reading from pdb
	char rotfile[200];			// file with rotamer library
	int cb_pdb;				// wether to keep CB from PDB
	int pdb_rot;				// whether to add the pdb rotamers or not
	double e_dih1;				// energy dihedrals
	double e_dih3;
	double e_ang;				// energy angles
	int go_dih;				// activate dihedrals go-style
	int go_ang;				// activate ang go-style
	double cys;				// makee al cysteines interact with some energy
	double cysr;				// .. its radius
	int go_noise;				// add noise to the go model (1=only to native contacts, 2=to any contact)
	double go_noise_sigma;			// stdev of the noise distribution

	char op_file[100];			//file to eventually write optimize potential datas for montegrappa
	char op_kind[100];
	int hydro_at[NTYPEMAX];
	double hydro_e;
	double hydro_r;
	char typesfile[100];
    
    int dih_ram;				//energy dihedrals with minimum in Ramachandran dihedrals
	double e_dihram;
	double sig_a_phi;
	double sig_b_phi;
	double sig_a_psi;
	double sig_b_psi;
	int phi_0_a;
	int phi_0_b;
	int psi_0_a;
	int psi_0_b;


};


struct atom_s
{
	char atom[5];
	char aa[5];
	int chain;
	int iaa;
	struct vector pos;
};

struct rot_atom_s {
        char     A1[5];
        char     A2[5];
        char     A3[5];
        char     kind[5];
        double  ang;
        double  dih;
        double  r;
};

struct rot_input_s {
        int     nrot;
        int     nat;
        char	aa[5];
        struct  rot_atom_s **atom;
        double 	weight[ROTMAX];
};

struct aa_presence
{
	int aa_isthere[50];
	char aa[50][50];
	double **aa_en;
};

struct atom_gromacs
{
	char aa[10];
	int isthere;
	char pdb_name[50][50];
	char gro_name[50][50];
	double sigma[50];
	double epsilon[50];
	double charge[50];
	int itype[50];
	int natom;
	double mass[50];
};

void Parse(FILE *fp, struct s_parms *p);
void AppendPotentialComments(struct s_parms *p, char *eoutfile);
struct atom_s *AlloAtoms(int n);
void Help(FILE *fp);
void GrappinoWelcome(FILE *fp);

double **ReadEMat(char *fname, struct s_parms *p);
double Dist(struct vector a, struct vector b);
double **AlloDoubleMatrix(int l, int m);
int **AlloIntMatrix(int l, int m);

// pdb.c
int PDB2CACB(struct atom_s *pdb, struct atom_s *ca, int n);
int PDB2CA(struct atom_s *pdb, struct atom_s *ca, int n);
int PDB2NCAC(struct atom_s *pdb, struct atom_s *ca, int n);
int ReadPDB(struct atom_s *x, char *pdbfile, int hydrogens, int *nchain, int *nbackmax, struct s_parms *p);
int IsBackbone(char *atom, struct s_parms *p);
void CreateTopology(struct atom_s *a, int n, int **top, double thresh, int debug);
void FindSidechain(struct atom_s *a, int natoms, int ipol, int iback, int w, int **top, struct s_polymer *p, struct s_parms *parms,
		int b1, int b2, int b3);
double Dist(struct vector a, struct vector b);
int CreateSidechainFromPDB(struct atom_s *a, int nchains, int na, struct s_polymer *polymer, struct s_parms *parms, int iamax);
int CreateBackboneFromPDB(struct atom_s *a, int nchains, int na, struct s_polymer *polymer, struct s_parms *parms, int *ib);
int Pdb2Polymer(struct atom_s *a, int nchains, int na, struct s_polymer *polymer, struct s_parms *parms, struct rot_input_s *rotamers, int nrot_kinds);
double DumbRMSD2(struct vector *a, struct vector *b, int n);
void SetRotamersSimilarToPDB(int nc, struct s_polymer *p, struct atom_s *a, int napdb);
void CopyPDB(struct atom_s *from, struct atom_s *to, int n);
int SimplifyPDB(struct atom_s *x, int n, char *model);
void ReadPropensity(char *fname, struct s_potential *u);


// energy.c
double **ContactMap(struct s_parms *parms, struct s_polymer *p, int nchains, int nat, int debug);
void Go_Pairs(struct s_parms *parms, struct s_polymer *p, double **e, double **r, double **r0, int nchains, int natoms, double **cm);
void Go_Dihedrals(struct s_parms *parms, struct s_polymer *p, int nc, double *dih01, double *dih03, struct s_potential *u);
void Go_Angles(struct s_parms *parms, struct s_polymer *p, int nc, double *ang, struct s_potential *u);
void Ram_Dihedrals(struct s_parms *p, struct s_potential *u);
int SetGoTypes(struct s_polymer *p, int nchains, int nat);
void DisulphideBonds(struct s_parms *parms, struct s_polymer *p, double **e, double **r2, double **r02, int nchains, int nat);
void OP_AddEnergy(struct s_polymer *p, int a1, int a2, double mul);
int ReadTypes(struct s_polymer *p, int nchains, int nat, char *nfile);
void PrintOpGoFile(char *nfile, double **cm, struct s_polymer *p, int nchains, char *kind, int imin);

// rotamers.c
struct rot_input_s *AlloRotamer(int naa, FILE *flog);
int ReadRotamers(char *fname, struct rot_input_s *x, int debug);
int Rot2PolymerScratch(int nrot_kinds, int nchains, struct s_polymer *p, struct rot_input_s *r, struct s_parms *parms, int iamax);
void SetDefaultCB(struct rot_input_s *CB, int nrot);
void SetPdbCB(struct s_back *b,struct rot_input_s *CB);
int AddPdbRot(struct s_back *b,int nrot);
int AddDefaultCB(struct s_polymer *p, int nc, int *N, int *CA, int *C, int naa, int na, int debug);
void SubstituteDefaultCB(struct s_polymer *p, int nc);
int Rot2Polymer(int nrot_kinds, int nchains, struct s_polymer *p, struct rot_input_s *r, struct s_parms *parms);

