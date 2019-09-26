 /*
 * Grappino.c
 *
 *  Created on: Sep 24, 2010
 *      Author: guido
 */

#include "montegrappa.h"
#include "grappino.h"

int main(int argc, char *argv[])
{
	struct s_parms *p;                      // parameters structure
	struct atom_s *pdb;                     // atom structure
	struct s_polymer *polymer;              // polymer structure
	int npdb, nchain, i, nrot_kinds=0, ntypes=0, nat, nbackmax;
	struct rot_input_s *rotamers;           // rotamers structure
	struct s_potential *u;                  // potential structure
	double **cm;                            // ???
	FILE *fin;

    GrappinoWelcome(stderr);
    
	if (argc != 2) {
        fprintf(stderr,"usage: grappino file.in\n");
        exit(1);
    }
	
	// open input file
	fin = fopen(argv[1],"r");
	if (!fin) Error("Cannot open inputfile");

	// allocate and read parameters from input file
	p = (struct s_parms *) calloc(1,sizeof(struct s_parms));
	if (!p) Error("Cannot allocate parameter structure");
	Parse(fin,p);
	fclose(fin);

	// get pdb and save the number of atoms (npdb)
	fprintf(stderr,"Getting pdb...\n");
	pdb = AlloAtoms(NATOMMAX);
	npdb = ReadPDB(pdb,p->pdbfile,p->hydrogens,&nchain,&nbackmax,p);

    // if needed, simplify the pdb (e.g. CA, CACB, NCAC)
	if (!strcmp(p->model,"CA") || !strcmp(p->model,"CACB") || !strcmp(p->model,"NCAC"))
            npdb = SimplifyPDB(pdb,npdb,p->model);
    
	// load rotamers
    
	if (p->rotamers)
	{
		rotamers = AlloRotamer(AAKINDMAX,stderr);
		nrot_kinds = ReadRotamers(p->rotfile,rotamers,p->debug);
	}

	

	///////////////
	// POLYMER FILE
	///////////////

	// allocate polymer and tables
	fprintf(stderr,"Allocate polymer\n");
	polymer = AlloPolymer(nchain,nbackmax,NSIDEMAX,NROTMAX,npdb,0,0,stderr);
	polymer->tables = InitTables(stderr);
	for (i=1;i<nchain;i++) (polymer+i)->tables = polymer->tables;			// tables of all polymers point to same address

	// transfer pdb to polymer structure
	nat = Pdb2Polymer(pdb,nchain,npdb,polymer,p,rotamers,nrot_kinds);		// nat is the total number of atoms

	/////////////////
	// POTENTIAL FILE
	/////////////////

	// set atomtypes
	if (!strcmp(p->atomtypes,"go")) ntypes = SetGoTypes(polymer,nchain,npdb);
	else ntypes = ReadTypes(polymer,nchain,npdb,p->atomtypes);
	
    
	cm = ContactMap(p,polymer,nchain,ntypes,p->debug);

	u = AlloPotential(npdb,ntypes,0,0,0);
	u->splice=0;
	u->g_imin = p->imin;
	fprintf(stderr,"imin = %d\n",u->g_imin);

	// define go potential
	if (!strcmp(p->potential,"go")) {
		Go_Pairs(p,polymer,u->e,u->r_2,u->r0_2,nchain,ntypes,cm);
		u->g_r0hard = p->rhard;
	}

	if (p->splice) {
		u->splice=1;
		u->ke_splice = p->ke_splice;
		u->kr2_splice = p->kr_splice;	// warning: here kr2_splice is not squared, it's gonna be squared in montegrappa!
	}

	// cysteine bridge
	if (p->cys<-EPSILON || p->cys>EPSILON)
		DisulphideBonds(p,polymer,u->e,u->r_2,u->r0_2,nchain,ntypes);


	// global energy terms
	if (p->e_homo<-EPSILON || p->e_homo> EPSILON) {
        u->g_ehomo = p->e_homo;
        u->g_rhomo = p->r_homo;
    }
    
    // go dihedrals and go angles
	if (p->go_dih)
        Go_Dihedrals(p,polymer,nchain,u->dih01,u->dih03,u);
	if (p->go_ang)
        Go_Angles(p,polymer,nchain,u->ang0,u);
    
    // ramachandran dihedrals
    if (p->dih_ram)
	{
		Ram_Dihedrals(p,u);
        fprintf(stderr,"Opening propensity file %s...\n",p->ab_propensityfile);
        fflush(stderr);
		ReadPropensity(p->ab_propensityfile,u);
	}



	// print output
	PrintPolymer(p->poutfile,polymer,nchain);
	PrintPotential(u,p->eoutfile,npdb,ntypes,0,0,0);
	AppendPotentialComments(p, p->eoutfile);

	//if you want a file.op to optimize the potential 
	if(strcmp(p->op_file,""))
        PrintOpGoFile(p->op_file,cm,polymer,nchain,p->op_kind,p->imin/3); // imin in aa

	exit(0);
}


void Parse(FILE *fp, struct s_parms *p)
{
	char aux[500];
	int i,j,sum = 0;

	// Set default values
	p->hydrogens = 0;
	p->rnat = 3.5;
	p->rhard = 2.;
	p->tthresh = 2.;
	p->imin = 9;
	strcpy(p->potential,"go");
	p->goe = -1.;
	p->e_homo = 0.;
	p->r_homo = 0.;
	strcpy(p->poutfile,"polymer.pol");
	strcpy(p->eoutfile,"potential.pot");
	strcpy(p->cntfile,"");
	p->n_back_a = 3;
	strcpy(p->back_a[0],"N");
    strcpy(p->back_a[1],"CA");
    strcpy(p->back_a[2],"C");
	for(i=0;i<100;i++) p->move_a[i]=1;
	p->debug=0;
	p->rotamers=0;
	p->cb_pdb = 0;
	p->pdb_rot = 0;
	p->go_dih = 0;
	p->go_ang = 0;
	p->k_native_r=1.;
	p->k_native_hc=0.5;
	p->go_noise = 0;
	p->go_noise_sigma = 1.;
	p->nosidechain=0;
	p->sidescratch = 0;
	strcpy(p->op_file,"");
	strcpy(p->op_kind,"GO_DIST_CA");
    
	p->dih_ram = 0;
	p->e_dihram =1.00;
	p->sig_a_phi = 25.;
	p->sig_b_phi = 30.;
	p->sig_a_psi = 30.;
	p->sig_b_psi = 35.;
	p->phi_0_a = -57;
	p->phi_0_b = -129;
	p->psi_0_a = -47;
	p->psi_0_b = 124;
    
    // read parameters file
	while(fgets(aux,500,fp)!=NULL)
    {
		ReadParS(aux,"pdbfile",p->pdbfile);
		ReadParS(aux,"polfile",p->poutfile);
		ReadParS(aux,"potfile",p->eoutfile);
		ReadParS(aux,"contactfile",p->cntfile);
        ReadParS(aux,"propensityfile",p->ab_propensityfile);
		ReadParS(aux,"model",p->model);
		ReadParN(aux,"hydrogens",&(p->hydrogens));
		ReadParF(aux,"r_hardcore",&(p->rhard));
		ReadParF(aux,"r_native",&(p->rnat));
		ReadParN(aux,"use_nativedist",&(p->use_nativedist));
		ReadParF(aux,"k_native_r",&(p->k_native_r));
		ReadParF(aux,"k_native_hc",&(p->k_native_hc));
		ReadParF(aux,"r_homo",&(p->r_homo));
		ReadParF(aux,"e_homo",&(p->e_homo));
		ReadParF(aux,"r_bonded",&(p->tthresh));
		ReadParD(aux,"debug",&(p->debug));
		ReadParD(aux,"imin",&(p->imin));
		ReadParF(aux,"go-energy",&(p->goe));
		ReadParS(aux,"potential",p->potential);
		ReadParS(aux,"atomtypes",p->atomtypes);
		ReadParN(aux,"rotamers",&(p->rotamers));
		ReadParS(aux,"rotfile",p->rotfile);
		ReadParN(aux,"cb_pdb",&(p->cb_pdb));
		ReadParN(aux,"pdb_rot",&(p->pdb_rot));
		ReadParN(aux,"go_dihedrals",&(p->go_dih));
		ReadParN(aux,"go_angles",&(p->go_ang));
		ReadParN(aux,"nosidechain",&(p->nosidechain));
		ReadParN(aux,"sidescratch",&(p->sidescratch));
		ReadParF(aux,"e_dih1",&(p->e_dih1));
		ReadParF(aux,"e_dih3",&(p->e_dih3));
		ReadParF(aux,"e_ang",&(p->e_ang));
		ReadParF(aux,"cys_e",&(p->cys));
		ReadParF(aux,"cys_r",&(p->cysr));
		ReadParD(aux,"go_noise",&(p->go_noise));
		ReadParF(aux,"go_nsigma",&(p->go_noise_sigma));

		ReadParS(aux,"op_file",p->op_file);
		ReadParS(aux,"op_kind",p->op_kind);
        
        ReadParN(aux,"dih_ram",&(p->dih_ram));
		ReadParF(aux,"e_dihram",&(p->e_dihram));
		ReadParF(aux,"sig_a_phi",&(p->sig_a_phi));
		ReadParF(aux,"sig_b_phi",&(p->sig_b_phi));
		ReadParF(aux,"sig_a_psi",&(p->sig_a_psi));
		ReadParF(aux,"sig_b_psi",&(p->sig_b_psi));
		ReadParD(aux,"phi_0_a",&(p->phi_0_a));
		ReadParD(aux,"phi_0_b",&(p->phi_0_b));
		ReadParD(aux,"psi_0_a",&(p->psi_0_a));
		ReadParD(aux,"psi_0_b",&(p->psi_0_b));

        // read explicitly declared backbone atoms
		if (!strncmp(aux,"backbone_atoms",14)) {
            (p->n_back_a)=0;
            ReadParD(aux,"backbone_atoms",&(p->n_back_a));
            if ( sscanf(aux,"backbone_atoms %*d %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s",
					  p->back_a[0],p->back_a[1],p->back_a[2],p->back_a[3],p->back_a[4],p->back_a[5],
                      p->back_a[6],p->back_a[7],p->back_a[8],p->back_a[9],p->back_a[10],p->back_a[11],
                      p->back_a[12],p->back_a[13],p->back_a[14],p->back_a[15],p->back_a[16],p->back_a[17],
                      p->back_a[18],p->back_a[19])<1 ) Error("Cannot read backbone_atoms in inputfile");
            // backbone atoms manually entered are unlocked by default,
            // one must specify if they need to be locked
            for(i=0;i<p->n_back_a;i++) p->move_a[i]=1;
        }


        //locked atoms: only backbone atoms can be locked, allowed moves of side chain atoms are in rotamers libraries
        if (!strncmp(aux,"locked_atoms",12)) {
            char locked_a[100][100];
            if ( sscanf(aux,"locked_atoms %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s",
                locked_a[0],locked_a[1],locked_a[2],locked_a[3],locked_a[4],locked_a[5],
                locked_a[6],locked_a[7], locked_a[8],locked_a[9],locked_a[10],locked_a[11],
                locked_a[12],locked_a[13],locked_a[14],locked_a[15],locked_a[16],locked_a[17],
                locked_a[18],locked_a[19])>=1 ) {
                    for (j=0;j<(p->n_back_a);j++)
                        for (i=0;i<(p->n_back_a);i++)
                            if (!strcmp(locked_a[j],p->back_a[i])) p->move_a[i]=0;
            }
            else Error("Cannot read locked_atoms in inputfile");
        }

        //splice option in potential
		if (!strncmp(aux,"splice",6)) {
            p->splice=1;
            if ( sscanf(aux,"splice %lf %lf",&(p->kr_splice), &(p->ke_splice)) < 1 )
                Error("Cannot read splice option in inputfile (splice <double> <double>)");
        }

    }
    
    // print backbone atoms
	fprintf(stderr, "%d backbone atoms found. They are: ",p->n_back_a);
	for (i=0;i<(p->n_back_a);i++)
        fprintf(stderr, "%s",p->back_a[i]);
    // print locked atoms
    fprintf(stderr, "\n these backbone atoms are locked: ");
	for (i=0;i<(p->n_back_a);i++)
        if (p->move_a[i]==0) fprintf(stderr, "	%s",p->back_a[i]);
    
	for (j=0;j<(p->n_back_a);j++)
        sum += p->move_a[j];
    if (sum==p->n_back_a)
        fprintf(stderr,"\n\n** WARNING: all backbone atoms are UNLOCKED ==> EVERYTHING can be MOVED\n"
                "** remember that the keyword locked_atoms must be written below backbone_atoms not above!!\n");
    if (sum==0)
        fprintf(stderr,"\n\n** WARNING: all backbone atoms are LOCKED ==> NOTHING can be MOVED\n");

    if (p->splice) {
        fprintf(stderr,"\n splice option enabled with kr_splice = %lf\t ke_splice = %lf\n",p->kr_splice,p->ke_splice);
        if (p->kr_splice>1)
        Error("in spice option kr_splice can't be bigger than 1\n");
		if (p->ke_splice>1) fprintf(stderr,"WARNING: in splice option ke_splice > 1 sounds strange, it makes long range interactions stronger than short range ones!!\n");
		if (p->use_nativedist && p->kr_splice*p->k_native_r*p->rnat< p->k_native_hc*p->rnat) fprintf(stderr,"\n ** WARNING **: kr_splice*k_native_r*r_native < k_native_hc*r_native ==> spliced distance smaller than hardcore repulsion!!\n");
            fprintf(stderr,"# splice option enabled. Size of the potential well:\n\n");
			fprintf(stderr,"# E \t |\n");
			fprintf(stderr,"# | \t |\t\t__________ 0\n");
			fprintf(stderr,"# | \t |\t\t|\n");
			fprintf(stderr,"# | \t |\t________| %4.2lf *E_native\n",p->ke_splice);
			fprintf(stderr,"# | \t |\t|\n");
			fprintf(stderr,"# | \t |______| E_native\n");
			fprintf(stderr,"# |______________________________> r\n");
			fprintf(stderr,"# 0 \t%4.2lf\t%4.2lf\t%4.2lf\n\n",p->k_native_hc*p->rnat,p->kr_splice*p->k_native_r*p->rnat,p->k_native_r*p->rnat);
	 }

	 fprintf(stderr, "\n \n");
}

/****************************************************************************
 Append potential parameters
 *****************************************************************************/
void AppendPotentialComments(struct s_parms *p, char *eoutfile)
{
	FILE *fout;

		// open file
		if (strcmp(eoutfile,"stdout"))
		{
			fout = fopen(eoutfile,"a");
			if (!fout) Error("Cannot open file to write potential comments");
		}
		else fout = stdout;

		fprintf(fout,"[ comments ]\n");

		fprintf(fout,"# potential type %s\n",p->potential);
		if(p->go_dih) fprintf(fout,"# energy of dihedral e_dih1 = %lf\te_dih3 = %lf\n",p->e_dih1,p->e_dih3);
        if(p->dih_ram) fprintf(fout,"# energy of ramachandran dihedrals e_dihram = %lf\n",p->e_dihram);
		if(p->go_ang) fprintf(fout,"# energy of angles e_ang = %lf\n",p->e_ang);
		fprintf(fout,"# energy of global hohomopolymeric interaction e_homo = %lf\t and its radius r_homo = %lf\n",p->e_homo,p->r_homo);
		fprintf(fout,"# threshold to define a bonded interaction = %lf\n",p->tthresh);
		fprintf(fout,"# native distance = %lf\n",p->rnat);
		fprintf(fout,"# hard core distance = %lf\n",p->rhard);
		if(p->use_nativedist)
		{
			fprintf(fout,"# use_nativedistance is on with parameters k_native_r = %lf\t and k_native_hc = %lf\n",p->k_native_r,p->k_native_hc);
		if(p->splice)
			{
			 fprintf(fout,"# splice option enabled. Size of the potential well:\n\n");
			 fprintf(fout,"# E \t |\n");
			 fprintf(fout,"# | \t |\t\t__________ 0\n");
			 fprintf(fout,"# | \t |\t\t|\n");
			 fprintf(fout,"# | \t |\t________| %4.2lf *E_native\n",p->ke_splice);
			 fprintf(fout,"# | \t |\t|\n");
			 fprintf(fout,"# | \t |______| E_native\n");
			 fprintf(fout,"# |______________________________> r\n");
			 fprintf(fout,"# 0 \t%4.2lf\t%4.2lf\t%4.2lf\n\n",p->k_native_hc*p->rnat,p->kr_splice*p->k_native_r*p->rnat,p->k_native_r*p->rnat);
			}
		}

		if (strcmp(eoutfile,"stdout")) fclose(fout);
}



/****************************************************************************
 Allocate atom structure
 *****************************************************************************/
struct atom_s *AlloAtoms(int n)
{
	struct atom_s *x;

	x = (struct atom_s *) calloc(n,sizeof(struct atom_s));
	if (!x) Error("Cannot allocate atom structure");

	return x;
}


void GrappinoWelcome(FILE *fp) {
    
    fprintf(fp,"\n\n");
    fprintf(fp,"       [=]\n");
    fprintf(fp,"       | |\t*Grappino*\n");
    fprintf(fp,"       }@{\t\tv1.0\n");
    fprintf(fp,"      /   \\ \n");
    fprintf(fp,"     /     \\ \n");
    fprintf(fp,"    /       \\ \n");
    fprintf(fp,"    :_______;\n");
    fprintf(fp,"    | MONTE |   \\~~~/\n");
    fprintf(fp,"    |GRAPPA |    \\_/\n");
    fprintf(fp,"    |-------|     |\n");
    fprintf(fp,"    '_______'   __|__\n");
    fprintf(fp,"\n\n");
    fprintf(fp,"\nG. Tiana, 2010\n");
    
}
