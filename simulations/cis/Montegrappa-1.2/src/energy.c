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
 * energy.c
 *
 *  Created on: Oct 11, 2010
 *      Author: guido
 */

#include "montegrappa.h"
#include "grappino.h"

/*****************************************************************************
 Go two-body potential (in this case ntypes=natoms)
 *****************************************************************************/
void Go_Pairs(struct s_parms *parms, struct s_polymer *p, double **e, double **r2, double **r02, int nchains, int nat, double **cm)
{
	int i,j,k=0;
	double r;
	FILE *fout;

	if (parms->go_noise) 
	{
		fprintf(stderr,"Adding noise (kind %d, sigma=%lf) to the Go potential\n",parms->go_noise,parms->go_noise_sigma);
		srand(time(0));	
	}
	if (parms->debug>1)  fprintf(stderr,"\n\nGo contacts:\n");

	for (i=0;i<nat;i++)
		for (j=i;j<nat;j++)
		{
			if (cm[i][j]>0)									// if it is a native contact
			{
				e[i][j] = parms->goe;
				if (parms->go_noise) e[i][j] += gauss(0.,parms->go_noise_sigma);	// if you want to add noise to the go potential
				e[j][i] = e[i][j];							// make the matrix symmetric
				if (parms->use_nativedist) r = cm[i][j] * cm[i][j] * parms->k_native_r * parms->k_native_r;
				else r = parms->rnat * parms->rnat;

				if (r2[i][j]>r) r = r2[i][j];					// if already set, use the maximum between the two
			
				r2[i][j] = r;
				r2[j][i] = r2[i][j];

				if (parms->use_nativedist) r = cm[i][j] * cm[i][j] * parms->k_native_hc * parms->k_native_hc;
				else r = parms->rhard * parms->rhard;

				if (r02[i][j]>0 && r02[i][j]<r) r=r02[i][j];
				r02[i][j] = r;
				r02[j][i] = r02[i][j];
				if (parms->debug>1)  fprintf(stderr,"%d-%d\t",i,j);
				k++;
			}

			if (parms->go_noise==2)								// if you want to add non-go noise
			{
				e[i][j] = gauss(0.,parms->go_noise_sigma);
				e[j][i] = e[i][j];
				r2[i][j] =  parms->rnat * parms->rnat;
				r2[j][i] = r2[i][j];
				r02[i][j] = parms->rhard * parms->rhard;
				r02[j][i] = r02[i][j];
			}
		}
	if (parms->debug>1)  fprintf(stderr,"\n");


	// if required, print native contacts to contactfile
	if (strcmp(parms->cntfile,""))
	{
		fout = fopen(parms->cntfile,"w");
		if (!fout) Error("Cannot open file to write native contacts");
		for (i=0;i<nat;i++)
			for (j=0;j<nat;j++)
				if (cm[i][j]>0) fprintf(fout,"%4d %4d\n",i,j);


		fclose (fout);
	}
	if (parms->debug>0) fprintf(stderr,"Found %d Go contacts\n",k);

}

void Go_Dihedrals(struct s_parms *parms, struct s_polymer *p, int nc, double *dih01, double *dih03, struct s_potential *u)
{
	int i,j,ia,out;

	for (i=0;i<nc;i++)
		for (j=2;j<(p+i)->nback-1;j++)
		{
			ia =  (((p+i)->back)+j)->ia;
			u->dih01[ia] = Dihedral( (((p+i)->back)+j-2)->pos, (((p+i)->back)+j-1)->pos, 
						(((p+i)->back)+j)->pos, (((p+i)->back)+j+1)->pos, p->tables,&out);
			u->e_dih1[ia] = parms->e_dih1;
			u->dih03[ia] = u->dih01[ia];
			u->e_dih3[ia] = parms->e_dih3;
		}
}

void Go_Angles(struct s_parms *parms, struct s_polymer *p, int nc, double *ang, struct s_potential *u)
{
	int i,j,ia,out;

	for (i=0;i<nc;i++)
		for (j=1;j<(p+i)->nback-1;j++)
		{
			ia =  (((p+i)->back)+j)->ia;
			u->ang0[ia] = Angle( (((p+i)->back)+j-1)->pos, (((p+i)->back)+j)->pos, (((p+i)->back)+j+1)->pos, p->tables,&out );
			u->e_ang[ia] = parms->e_ang;
		}
}

void Ram_Dihedrals(struct s_parms *p, struct s_potential *u)
{
	if(p->dih_ram!=0)
	{
		u->dih_ram = p->dih_ram;
		u->e_dihram = p->e_dihram;
		u->sigma[0][0] = p->sig_a_phi;
		u->sigma[1][0] = p->sig_b_phi;
		u->sigma[0][1] = p->sig_a_psi;
		u->sigma[1][1] = p->sig_b_psi;
		u->dih0[0][0] = p->phi_0_a;
		u->dih0[1][0] = p->phi_0_b;
		u->dih0[0][1] = p->psi_0_a;
		u->dih0[1][1] = p->psi_0_b;
	}
}

/*****************************************************************************
 Contact map between atoms sticking out of different backbone atoms
 (nat x nat)
 *****************************************************************************/
double **ContactMap(struct s_parms *parms, struct s_polymer *p, int nchains, int nat,int debug)
{
	double **cm;
	int ib1,ib2,is1,is2,c1,c2;
	double d;

	cm = AlloDoubleMatrix(nat,nat);
	if (debug>2) fprintf(stderr,"Contact map:\n");

	// loop on backbone atoms
	for (c1=0;c1<nchains;c1++)
		for (c2=c1;c2<nchains;c2++)
			for (ib1=0;ib1<(p+c1)->nback;ib1++)
				for (ib2=0;ib2<(p+c2)->nback;ib2++)
					if (c1!=c2 || ib1<=ib2-parms->imin)
					{
						// backbone-backbone interactions
						d = Dist( (((p+c1)->back)+ib1)->pos, (((p+c2)->back)+ib2)->pos );
						if ( d < parms->rnat )
						{
							cm[(((p+c1)->back)+ib1)->itype][(((p+c2)->back)+ib2)->itype] = d;
							cm[(((p+c2)->back)+ib2)->itype][(((p+c1)->back)+ib1)->itype] = d;
							if (debug>2) fprintf(stderr,"\t%d%s (%d%s) - %d%s (%d%s)\td=%lf\n",(((p+c1)->back)+ib1)->ia,(((p+c1)->back)+ib1)->type,
									(((p+c1)->back)+ib1)->iaa,(((p+c1)->back)+ib1)->aa,(((p+c2)->back)+ib2)->ia,(((p+c2)->back)+ib2)->type,
									(((p+c2)->back)+ib2)->iaa,(((p+c2)->back)+ib2)->aa,d);
						}

						// sidechain-backbone interactions
						for (is1=0;is1<(((p+c1)->back)+ib1)->nside;is1++)
						{
							d = Dist( (((((p+c1)->back)+ib1)->side)+is1)->pos, (((p+c2)->back)+ib2)->pos );
							if ( d < parms->rnat )
							{
								cm[(((((p+c1)->back)+ib1)->side)+is1)->itype][(((p+c2)->back)+ib2)->itype] = d;
								cm[(((p+c2)->back)+ib2)->itype][(((((p+c1)->back)+ib1)->side)+is1)->itype] = d;
								if (debug>2) fprintf(stderr,"\t%d%s (%d%s) - %d%s (%d%s)\td=%lf\n",(((((p+c1)->back)+ib1)->side)+is1)->ia,(((((p+c1)->back)+ib1)->side)+is1)->type,
									(((p+c1)->back)+ib1)->iaa,(((p+c1)->back)+ib1)->aa,(((p+c2)->back)+ib2)->ia,(((p+c2)->back)+ib2)->type,
									(((p+c2)->back)+ib2)->iaa,(((p+c2)->back)+ib2)->aa,d);

							}
						}
						for (is2=0;is2<(((p+c2)->back)+ib2)->nside;is2++)
						{
							d = Dist( (((((p+c2)->back)+ib2)->side)+is2)->pos, (((p+c1)->back)+ib1)->pos );
							if ( d < parms->rnat )
							{
								cm[(((((p+c2)->back)+ib2)->side)+is2)->itype][(((p+c1)->back)+ib1)->itype] = d;
								cm[(((p+c1)->back)+ib1)->itype][(((((p+c2)->back)+ib2)->side)+is2)->itype] = d;
								if (debug>2) fprintf(stderr,"\t%d%s (%d%s) - %d%s (%d%s)\td=%lf\n",(((p+c1)->back)+ib1)->ia,(((p+c1)->back)+ib1)->type,
									(((p+c1)->back)+ib1)->iaa,(((p+c1)->back)+ib1)->aa,(((((p+c2)->back)+ib2)->side)+is2)->ia,(((((p+c2)->back)+ib2)->side)+is2)->type,
									(((p+c2)->back)+ib2)->iaa,(((p+c2)->back)+ib2)->aa,d);

							}
						}
						// sidechain-sidechain interactions
						for (is1=0;is1<(((p+c1)->back)+ib1)->nside;is1++)
							for (is2=0;is2<(((p+c2)->back)+ib2)->nside;is2++)
							{
								d = Dist( (((((p+c1)->back)+ib1)->side)+is1)->pos, (((((p+c2)->back)+ib2)->side)+is2)->pos );
								if ( d < parms->rnat )
								{
									cm[(((((p+c1)->back)+ib1)->side)+is1)->itype][(((((p+c2)->back)+ib2)->side)+is2)->itype] = d;
									cm[(((((p+c2)->back)+ib2)->side)+is2)->itype][(((((p+c1)->back)+ib1)->side)+is1)->itype] = d;
									if (debug>2) fprintf(stderr,"\t%d%s (%d%s) - %d%s (%d%s)\td=%lf\n",(((((p+c1)->back)+ib1)->side)+is1)->ia,(((((p+c1)->back)+ib1)->side)+is1)->type,
										(((p+c1)->back)+ib1)->iaa,(((p+c1)->back)+ib1)->aa,(((((p+c2)->back)+ib2)->side)+is2)->ia,(((((p+c2)->back)+ib2)->side)+is2)->type,
										(((p+c2)->back)+ib2)->iaa,(((p+c2)->back)+ib2)->aa,d);

								}
							}
					}

			return cm;
}

/*****************************************************************************
 Assign as atom types the atom id (and ntypes=natoms), as it is in the Go model
 *****************************************************************************/
int SetGoTypes(struct s_polymer *p, int nchains, int nat)
{
	int ic,i,j;

	fprintf(stderr,"Set Go types\n");

	for (ic=0;ic<nchains;ic++)
		for (i=0;i<(p+ic)->nback;i++)
		{
			(((p+ic)->back)+i)->itype = (((p+ic)->back)+i)->ia;
			for (j=0;j<(((p+ic)->back)+i)->nside;j++) (((((p+ic)->back)+i)->side)+j)->itype = (((((p+ic)->back)+i)->side)+j)->ia;
		}

	return nat;
}

/*****************************************************************************
 Read atom types from a library file
 *****************************************************************************/
int ReadTypes(struct s_polymer *p, int nchains, int nat, char *nfile)
{
	int i,ic,j,k,it=0,type[NTYPEMAX],nt=0;
	char aux[500],atom[NTYPEMAX][5],aa[NTYPEMAX][5];
	FILE *fp;

	// read from file
	fprintf(stderr,"Reading atom types from file %s\n",nfile);

	fp = fopen(nfile,"r");
	if (!fp) Error("Cannot open file for reading atom types, check atomtypes directive");

	while(fgets(aux,500,fp)!=NULL)
	{
		if ( sscanf(aux,"%d %s %s",&(type[it]),atom[it],aa[it]) == 3) it++;
		if (it>=NTYPEMAX) Error("NTYPEMAX too small in ReadTypes");
		if (it>nt) nt=it;
	}	

	fprintf(stderr,"Read %d atomtypes\n",it);
	fclose(fp);

	// assign types
	for (ic=0;ic<nchains;ic++)
		for (i=0;i<(p+ic)->nback;i++)
		{
			// backbone
			(((p+ic)->back)+i)->itype = -1;
			for (k=0;k<it;k++)
				if (  ( !strcmp((((p+ic)->back)+i)->type, atom[k]) || !strcmp( "*", atom[k]) )
					&&  ( !strcmp( (((p+ic)->back)+i)->aa, aa[k]) || !strcmp( "*", aa[k])) )
						(((p+ic)->back)+i)->itype = type[k];
			if ( (((p+ic)->back)+i)->itype == -1)
			{
				fprintf(stderr,"Cannot assign atomtype to chain %d, back=%d (aa=%s atom=%s)\n",ic,i,
					(((p+ic)->back)+i)->aa,(((p+ic)->back)+i)->type);
				Error("Cannot assign atom type");
			}

			// sidechain
			for (j=0;j<(((p+ic)->back)+i)->nside;j++)
			{
				(((((p+ic)->back)+i)->side)+j)->itype = -1;
				for (k=0;k<it;k++)
					if ( ( !strcmp((((((p+ic)->back)+i)->side)+j)->type, atom[k]) || !strcmp( "*", atom[k]) )
						&& ( !strcmp( (((p+ic)->back)+i)->aa, aa[k]) || !strcmp( "*", aa[k])) )
							(((((p+ic)->back)+i)->side)+j)->itype = type[k];


				if ( (((((p+ic)->back)+i)->side)+j)->itype == -1)
				{
					fprintf(stderr,"Cannot assign atomtype to chain %d, back=%d (aa=%s atom=%s)\n",ic,i,
						(((p+ic)->back)+i)->aa,(((((p+ic)->back)+i)->side)+j)->type);
					Error("Cannot assign atom type");
				}
			}
		}

	return nt;
}


/*****************************************************************************
 Add a cys-cys interaction
 *****************************************************************************/
void DisulphideBonds(struct s_parms *parms, struct s_polymer *p, double **e, double **r2, double **r02, int nchains, int nat)
{
	int i,j,ic,jc,ki,kj,it1,it2,cnt=0;

	fprintf(stderr,"Adding disulphide bonds...\n");

	for (ic=0;ic<nchains;ic++)
		for (jc=0;jc<nchains;jc++)
			for (i=0;i<(p+ic)->nback;i++)
				for (j=0;j<(p+jc)->nback;j++)
					if ( !strcmp( (((p+ic)->back)+i)->aa , "CYS" ) )
						if ( !strcmp( (((p+jc)->back)+j)->aa , "CYS" ) )
							for (ki=0;ki<(((p+ic)->back)+i)->nside;ki++)
								for (kj=0;kj<(((p+jc)->back)+j)->nside;kj++)
									if ( !strncmp( (((((p+ic)->back)+i)->side)+ki)->type , "S" ,1) )
										if ( !strncmp( (((((p+jc)->back)+j)->side)+kj)->type , "S" ,1) )
											if (ic!=jc || i!=j)
												{
													it1 = (((((p+ic)->back)+i)->side)+ki)->itype;
													it2 = (((((p+jc)->back)+j)->side)+kj)->itype;
													e[it1][it2] = parms->cys;
													e[it2][it1] = e[it1][it2];
													r2[it1][it2] = parms->cysr * parms->cysr;
													r2[it2][it1] = r2[it1][it2];
													if (r02[it1][it2]<EPSILON)
													{
														r02[it1][it2] = parms->rhard;
														r02[it2][it1] = r02[it1][it2];
													}
													cnt++;
												}
	fprintf(stderr,"Added %d disulphide bond interactions.\n",cnt);
}

void HydrophobicInteraction(struct s_parms *parms, struct s_polymer *p, double **e, double **r2, double **r02, int natypes)
{
	int i,j,cnt=0;

	fprintf(stderr,"Adding hydrophobic interaction...\n");

	for (i=0;i<natypes;i++)
		for (j=0;j<natypes;j++)
			if (parms->hydro_at[i]==1 && parms->hydro_at[j]==1)
			{
				e[i][j] = parms->hydro_e;
				e[j][i] = e[i][j];
				r2[i][j] = parms->hydro_r * parms->hydro_r;
				r2[j][i] = r2[i][j];
				if (r02[i][j]<EPSILON)
				{
					r02[i][j] = parms->rhard;
					r02[j][i] = r02[i][j];
				}
				cnt++;
			}
}

/**********************************************
 Create OP file
 **********************************************/
void PrintOpGoFile(char *nfile, double **cm, struct s_polymer *p, int nchains, char *kind, int imin)
{
	int ic,jc,i,j,nCA;
	FILE *fp;
    
	fprintf(stderr,"Write OP file %s for potential optimization (kind=%s)\n",nfile,kind);
    
	fp = fopen(nfile,"w");
	if (!fp) Error("Cannot open op file for writing");
    
	// distances between native CA (numbers are ia)
	if (!strcmp(kind,"GO_DIST_CA"))
	{
        int nrest = 0;
        for (ic=0;ic<nchains;ic++)
			for (i=0;i<(p+ic)->nback;i++)
                for (jc=ic;jc<nchains;jc++) {
                    // same chain
                    if (ic == jc) {
                        nCA = 0;
                        for (j=i+1;j<(p+jc)->nback;j++) {
                            // counts CA atoms to check imin condition
                            if  (!strcmp((((p+jc)->back)+j)->type,"CA") )
                                nCA++;
                            if ( !strcmp((((p+ic)->back)+i)->type,"CA") && !strcmp((((p+jc)->back)+j)->type,"CA") && ( nCA >= imin ) )
                                nrest++;
                        }
                    }
                    // different chains
                    else {
                        for (j=0;j<(p+jc)->nback;j++) {
                            if ( !strcmp((((p+ic)->back)+i)->type,"CA") && !strcmp((((p+jc)->back)+j)->type,"CA") )
                                nrest++;
                        }
                    }
                }
        
		fprintf(fp,"ndata %d\n",nrest);
        
        for (ic=0;ic<nchains;ic++)
			for (i=0;i<(p+ic)->nback;i++)
                for (jc=ic;jc<nchains;jc++) {
                    if (ic == jc) {
                        nCA = 0;
                        for (j=i+1;j<(p+jc)->nback;j++) {
                            if  (!strcmp((((p+jc)->back)+j)->type,"CA") )
                                nCA++;
                            if ( !strcmp((((p+ic)->back)+i)->type,"CA") && !strcmp((((p+jc)->back)+j)->type,"CA") && ( nCA >= imin ) )
                                fprintf(fp,"%d\t%d\t\t1\t%lf\t0.5\n",(((p+ic)->back)+i)->ia,(((p+jc)->back)+j)->ia,
                                        Dist( (((p+ic)->back)+i)->pos, (((p+jc)->back)+j)->pos ) );
                        }
                    }
                    else {
                        for (j=0;j<(p+jc)->nback;j++) {
                            if ( !strcmp((((p+ic)->back)+i)->type,"CA") && !strcmp((((p+jc)->back)+j)->type,"CA") )
                                fprintf(fp,"%d\t%d\t\t1\t%lf\t0.5\n",(((p+ic)->back)+i)->ia,(((p+jc)->back)+j)->ia,
                                        Dist( (((p+ic)->back)+i)->pos, (((p+jc)->back)+j)->pos ) );
                        }
                    }
                }
	}
    //
	if(!strcmp(kind,"GO_DIST_ALLATOM"))
	{
        int nrest = 0;
        int iaa,jaa,iside,jside;
        for (ic=0;ic<nchains;ic++)
            for (i=0;i<(p+ic)->nback;i++) {
                iaa = i/3 + 1;  //define aminoacid id
                for (jc=ic;jc<nchains;jc++) {
                    if (ic == jc) {    // same chain
                        nCA = 0;
                        for (j=i;j<(p+jc)->nback;j++) {
                            jaa = j/3 + 1;
                            if ((jaa-iaa)%2 == 0) { // only counts restraints between even aa and between even aa
                                if  (!strcmp((((p+jc)->back)+j)->type,"CA") ) // counts CA atoms to check imin condition
                                    nCA++;
                                if ( !strcmp((((p+ic)->back)+i)->type,"CA") && !strcmp((((p+jc)->back)+j)->type,"CA") && ( nCA >= imin ) ){
                                    if( strcmp(((((p+ic)->back)+i)->aa),"GLY") && strcmp(((((p+jc)->back)+j)->aa),"GLY") )
                                        nrest+=4;   // one for CA-CA, one for sidechain/sidechain and two for CA/sidechain
                                    if( strcmp(((((p+ic)->back)+i)->aa),"GLY") && !strcmp(((((p+jc)->back)+j)->aa),"GLY") )
                                        nrest+=2;
                                    if( !strcmp(((((p+ic)->back)+i)->aa),"GLY") && strcmp(((((p+jc)->back)+j)->aa),"GLY") )
                                        nrest+=2;
                                    if( !strcmp(((((p+ic)->back)+i)->aa),"GLY") && !strcmp(((((p+jc)->back)+j)->aa),"GLY") )
                                        nrest+=1;
                                }
                            }
                        }
                    }
                    // different chains
                    else {
                        for (j=0;j<(p+jc)->nback;j++) {
                            jaa = j/3 + 1;
                            if ((jaa-iaa)%2 == 0) { // only counts restraints between even aa and between even aa
                                if ( !strcmp((((p+ic)->back)+i)->type,"CA") && !strcmp((((p+jc)->back)+j)->type,"CA") ) {
                                    if( strcmp(((((p+ic)->back)+i)->aa),"GLY") && strcmp(((((p+jc)->back)+j)->aa),"GLY") )
                                        nrest+=4;   // one for CA-CA, one for sidechain/sidechain and two for CA/sidechain
                                    if( strcmp(((((p+ic)->back)+i)->aa),"GLY") && !strcmp(((((p+jc)->back)+j)->aa),"GLY") )
                                        nrest+=2;
                                    if( !strcmp(((((p+ic)->back)+i)->aa),"GLY") && strcmp(((((p+jc)->back)+j)->aa),"GLY") )
                                        nrest+=2;
                                    if( !strcmp(((((p+ic)->back)+i)->aa),"GLY") && !strcmp(((((p+jc)->back)+j)->aa),"GLY") )
                                        nrest+=1;
                                }
                            }
                        }
                    }
                }
            }
        
        fprintf(fp,"ndata %d\n",nrest);
        
        for (ic=0;ic<nchains;ic++)
            for (i=0;i<(p+ic)->nback;i++) {
                iaa = i/3 + 1;  //define aminoacid id
                for (jc=ic;jc<nchains;jc++) {
                    if (ic == jc) {    // same chain
                        nCA = 0;
                        for (j=i;j<(p+jc)->nback;j++) {
                            jaa = j/3 + 1;
                            if (abs((iaa-jaa))%2 == 0) { // only counts restraints between even aa and between even aa
                                if  (!strcmp((((p+jc)->back)+j)->type,"CA") ) // counts CA atoms to check imin condition
                                    nCA++;
                                if ( !strcmp((((p+ic)->back)+i)->type,"CA") && !strcmp((((p+jc)->back)+j)->type,"CA") && ( nCA >= imin ) ) {
                                    fprintf(fp,"%d\t%d\t\t1\t%lf\t0.5\n",(((p+ic)->back)+i)->ia,(((p+jc)->back)+j)->ia,
                                            Dist( (((p+ic)->back)+i)->pos, (((p+jc)->back)+j)->pos ) );   // Backbone-Backbone
                                    iside = ((((p+ic)->back)+i)->nside)-1;
                                    jside = ((((p+jc)->back)+j)->nside)-1;
                                    if( !strcmp(((((p+ic)->back)+i)->aa),"GLY") && strcmp(((((p+jc)->back)+j)->aa),"GLY") )
                                        fprintf(fp,"%d\t%d\t\t1\t%lf\t0.5\n",(((p+ic)->back)+i)->ia,(((((p+jc)->back)+j)->side)+jside)->ia,
                                                Dist( (((p+ic)->back)+i)->pos, (((((p+jc)->back)+j)->side)+jside)->pos ) );    // Backbone-Sidechain
                                    if( strcmp(((((p+ic)->back)+i)->aa),"GLY") && !strcmp(((((p+jc)->back)+j)->aa),"GLY") )
                                        fprintf(fp,"%d\t%d\t\t1\t%lf\t0.5\n",(((((p+ic)->back)+i)->side)+iside)->ia,(((p+jc)->back)+j)->ia,
                                                Dist( (((((p+ic)->back)+i)->side)+iside)->pos, (((p+jc)->back)+j)->pos ) );   // Sidechain-Backbone
                                    if( strcmp(((((p+ic)->back)+i)->aa),"GLY") && strcmp(((((p+jc)->back)+j)->aa),"GLY") ) {
                                        fprintf(fp,"%d\t%d\t\t1\t%lf\t0.5\n",(((p+ic)->back)+i)->ia,(((((p+jc)->back)+j)->side)+jside)->ia,
                                                Dist( (((p+ic)->back)+i)->pos, (((((p+jc)->back)+j)->side)+jside)->pos ) );    // Backbone-Sidechain
                                        fprintf(fp,"%d\t%d\t\t1\t%lf\t0.5\n",(((((p+ic)->back)+i)->side)+iside)->ia,(((p+jc)->back)+j)->ia,
                                                Dist( (((((p+ic)->back)+i)->side)+iside)->pos, (((p+jc)->back)+j)->pos ) );    // Sidechain-Backbone
                                        fprintf(fp,"%d\t%d\t\t1\t%lf\t0.5\n",(((((p+ic)->back)+i)->side)+iside)->ia,(((((p+jc)->back)+j)->side)+jside)->ia,
                                                Dist( (((((p+ic)->back)+i)->side)+iside)->pos, (((((p+jc)->back)+j)->side)+jside)->pos ) );    // Sidechain-Sidechain
                                    }
                                }
                            }
                        }
                    }
                    // different chains
                    else {
                        for (j=0;j<(p+jc)->nback;j++) {
                            jaa = j/3 + 1;
                            if (abs((iaa-jaa))%2 == 0) { // only counts restraints between even aa and between odd aa
                                if ( !strcmp((((p+ic)->back)+i)->type,"CA") && !strcmp((((p+jc)->back)+j)->type,"CA") ) {
                                    fprintf(fp,"%d\t%d\t\t1\t%lf\t0.5\n",(((p+ic)->back)+i)->ia,(((p+jc)->back)+j)->ia,
                                            Dist( (((p+ic)->back)+i)->pos, (((p+jc)->back)+j)->pos ) );   // Backbone-Backbone
                                    iside = ((((p+ic)->back)+i)->nside)-1;
                                    jside = ((((p+jc)->back)+j)->nside)-1;
                                    if( !strcmp(((((p+ic)->back)+i)->aa),"GLY") && strcmp(((((p+jc)->back)+j)->aa),"GLY") )
                                        fprintf(fp,"%d\t%d\t\t1\t%lf\t0.5\n",(((p+ic)->back)+i)->ia,(((((p+jc)->back)+j)->side)+jside)->ia,
                                                Dist( (((p+ic)->back)+i)->pos, (((((p+jc)->back)+j)->side)+jside)->pos ) );    // Backbone-Sidechain
                                    if( strcmp(((((p+ic)->back)+i)->aa),"GLY") && !strcmp(((((p+jc)->back)+j)->aa),"GLY") )
                                        fprintf(fp,"%d\t%d\t\t1\t%lf\t0.5\n",(((((p+ic)->back)+i)->side)+iside)->ia,(((p+jc)->back)+j)->ia,
                                                Dist( (((((p+ic)->back)+i)->side)+iside)->pos, (((p+jc)->back)+j)->pos ) );   // Sidechain-Backbone
                                    if( strcmp(((((p+ic)->back)+i)->aa),"GLY") && strcmp(((((p+jc)->back)+j)->aa),"GLY") ) {
                                        fprintf(fp,"%d\t%d\t\t1\t%lf\t0.5\n",(((p+ic)->back)+i)->ia,(((((p+jc)->back)+j)->side)+jside)->ia,
                                                Dist( (((p+ic)->back)+i)->pos, (((((p+jc)->back)+j)->side)+jside)->pos ) );    // Backbone-Sidechain
                                        fprintf(fp,"%d\t%d\t\t1\t%lf\t0.5\n",(((((p+ic)->back)+i)->side)+iside)->ia,(((p+jc)->back)+j)->ia,
                                                Dist( (((((p+ic)->back)+i)->side)+iside)->pos, (((p+jc)->back)+j)->pos ) );    // Sidechain-Backbone
                                        fprintf(fp,"%d\t%d\t\t1\t%lf\t0.5\n",(((((p+ic)->back)+i)->side)+iside)->ia,(((((p+jc)->back)+j)->side)+jside)->ia,
                                                Dist( (((((p+ic)->back)+i)->side)+iside)->pos, (((((p+jc)->back)+j)->side)+jside)->pos ) );    // Sidechain-Sidechain
                                    }
                                }
                            }
                        }
                    }
                }
            }
        
    }
    
	
    
	fclose(fp);
	
    
}
