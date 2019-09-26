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
 * pdb.c
 *
 *  Created on: Oct 11, 2010
 *      Author: guido
 */

#include "montegrappa.h"
#include "grappino.h"

/*****************************************************************************
 Reads a pdbfile and puts the information in the atom structure
 returns the number of atoms
	 if hydrogens==0, do not read them
 *****************************************************************************/
int ReadPDB(struct atom_s *x, char *pdbfile, int hydrogens, int *nchain, int *nbackmax, struct s_parms *parms)
{
  int i,k,ok;
  //int firstaa=-1;
  char aux[1000],aux2[15];
  char atom[5],aa[4],chainname,oldchainname;
  int iatom=0,iaa,iaa0=-1;
  double px,py,pz;
  FILE *fp;

  *nchain=0;
  *nbackmax=0;
  chainname=' ';
  oldchainname=' ';

  fp = fopen(pdbfile,"r");
  if (!fp) {
      fprintf(stderr,"Problem with file %s.\n",pdbfile);
      Error("Cannot find PDB file ");
  }

  while(fgets(aux,500,fp)!=NULL)
    //if the first four letters of aux are not ATOM
    if (!strncmp(aux,"ATOM",4)){
      ok=1;
      // Copy atom name, positioned between 12 and 16 on the pdb lines
      k=0;
      for(i=12;i<=15;i++) {
	    if (aux[i]!=' ') {
          atom[k]=aux[i];
	      k++;
	    }
      }
      atom[k]='\0';

	  if (aux[16]!=' ') Error("Grappino can't handle alternate positions of atoms ");
      // remove hydrogens if requested
	  if (hydrogens==0 && atom[0]=='H') ok=0;
        for(i=0;i<parms->n_back_a;i++)
		  if (!strcmp(atom,parms->back_a[i])) (*nbackmax) ++;
        
        // Copy amino acid name
      k=0;
	  for(i=17;i<=19;i++) {
	    if (aux[i]!=' ') {
          aa[k]=aux[i];
	      k++;
        }
      }
      aa[k]='\0';

	  // Copy chain
      if (aux[21]!=' ') chainname=aux[21];

	  // Copy amino acid number
	  k=0;
      for(i=22;i<=25;i++) {
        if (aux[i]!=' ') {
	    		   aux2[k]=aux[i];
	    		   k++;
	    }
      }
      aux2[k]='\0';
	  iaa = atoi(aux2);
	  if (iaa0==-1)
        iaa0=iaa;

      // Copy x coordinate
	  k=0;
      for(i=30;i<=37;i++) {
        if (aux[i]!=' ') {
          aux2[k]=aux[i];
	      k++;
	    }
      }
      aux2[k]='\0';
	  px = atof(aux2);

	  // Copy y coordinate
	  k=0;
	  for(i=38;i<=45;i++)
        if (aux[i]!=' ') {
	      aux2[k]=aux[i];
	      k++;
	    }
      aux2[k]='\0';
	  py = atof(aux2);

	  // Copy z coordinate
	  k=0;
	  for(i=46;i<=53;i++)
        if (aux[i]!=' ') {
          aux2[k]=aux[i];
	      k++;
	    }
      aux2[k]='\0';
	  pz = atof(aux2);

	  // copy everything to structure
	  strcpy( (x+iatom)->atom,atom);
	  strcpy( (x+iatom)->aa,aa);
	  (x+iatom)->iaa = iaa;
	  ((x+iatom)->pos).x = px;
	  ((x+iatom)->pos).y = py;
	  ((x+iatom)->pos).z = pz;
			 
	  if (iatom>0)
	    if (chainname != oldchainname ) {
          (*nchain) ++;
		  iaa0=-1;
	    }
      (x+iatom)->chain = *nchain;
	  oldchainname = chainname;

      // accept atom
      if (ok==1) iatom ++;
      if (iatom>NATOMMAX)
        Error("NATOMMAX too small");
      }


   (*nchain)++;

   fprintf(stderr,"Read %d atoms in %d chains (nbackmax=%d)\n",iatom,*nchain,*nbackmax);

   return iatom;
}

/*****************************************************************************
 Turn a full pdb into a CA model
 *****************************************************************************/
int PDB2CA(struct atom_s *pdb, struct atom_s *ca, int n)
{
   int i,k=0;

   for(i=0;i<n;i++)
	 if ( !strcmp( (pdb+i)->atom, "CA" ) )
	 {
	    strcpy((ca+k)->aa,(pdb+i)->aa);
	    strcpy((ca+k)->atom,(pdb+i)->atom);
	    (ca+k)->chain = (pdb+i)->chain;
	    (ca+k)->iaa = (pdb+i)->iaa;
	    (ca+k)->pos = (pdb+i)->pos;
        k++;
	 }

   // check if every residue has a CA
   for(i=0;i<k;i++)
	 if ( (ca+i)->iaa != i+1) fprintf(stderr,"WARNING: missing CA for residue %d\n",(ca+i)->iaa);

	return k;
}

/*****************************************************************************
 Turn a full pdb into a CA-CB model
 *****************************************************************************/
int PDB2CACB(struct atom_s *pdb, struct atom_s *ca, int n)
{
   int i,k=0,j,nacb;
   struct vector cm;

   for(i=0;i<n;i++)
	 if ( !strcmp( (pdb+i)->atom, "CA" ) )
	 {
		// add CA
	    strcpy((ca+k)->aa,(pdb+i)->aa);
	    strcpy((ca+k)->atom,(pdb+i)->atom);
	    (ca+k)->chain = (pdb+i)->chain;
	    (ca+k)->iaa = (pdb+i)->iaa;
	    (ca+k)->pos = (pdb+i)->pos;
	    k++;

	    // add CB
        // If the residue is a GLY, use H instead
	    if (strcmp((ca+k-1)->aa,"GLY"))
	    {
	       nacb=0;
	       cm.x=0; cm.y=0; cm.z=0;
	       for(j=0;j<n;j++)
	    	 if ( (pdb+j)->iaa == (pdb+i)->iaa )
	    	   if ( strcmp( (pdb+j)->atom,"C") && strcmp( (pdb+j)->atom,"O") && strcmp( (pdb+j)->atom,"N") && strcmp( (pdb+j)->atom,"H") && strcmp( (pdb+j)->atom,"HA") )
	    	   {
	    		   cm.x += ((pdb+j)->pos).x;
	    		   cm.y += ((pdb+j)->pos).y;
	    		   cm.z += ((pdb+j)->pos).z;
	    		   nacb ++;
	    	   }
		    strcpy((ca+k)->aa,(pdb+i)->aa);
		    strcpy((ca+k)->atom,"CB");
		    (ca+k)->chain = (pdb+i)->chain;
		    (ca+k)->iaa = (pdb+i)->iaa;
		    cm.x /= nacb;
		    cm.y /= nacb;
		    cm.z /= nacb;
		    (ca+k)->pos = cm;
		    k++;
	    }

	 }

	return k;
}

/*****************************************************************************
 Turn a full pdb into a NCAC model
 *****************************************************************************/
int PDB2NCAC(struct atom_s *pdb, struct atom_s *ca, int n)
{
    int i,k=0;
    
    for(i=0;i<n;i++)
        if ( !strcmp( (pdb+i)->atom, "N" ) || !strcmp( (pdb+i)->atom, "CA" ) || \
            !strcmp( (pdb+i)->atom, "C" ))
        {
            strcpy((ca+k)->aa,(pdb+i)->aa);
            strcpy((ca+k)->atom,(pdb+i)->atom);
            (ca+k)->chain = (pdb+i)->chain;
            (ca+k)->iaa = (pdb+i)->iaa;
            (ca+k)->pos = (pdb+i)->pos;
            k++;
        }
    
	return k;
}

/****************************************************************************
 Fill an empty polymer file with atom positions and return the total atom number 
 *****************************************************************************/
int Pdb2Polymer(struct atom_s *a, int nchains, int na, struct s_polymer *polymer, struct s_parms *parms, struct rot_input_s *rotamers, int nrot_kinds)
{
	int i,ib[NPOLMAX],iamax=0;

	// generate the backbone, putting 
 	iamax = CreateBackboneFromPDB(a,nchains,na,polymer,parms,ib);

	fprintf(stderr,"\nPdb2Polymer (nchains=%d natoms=%d):\n",nchains,na);
                                   
	if (!parms->nosidechain)
	{
		// or just using rotamer library
		if (parms->sidescratch)
		{
			iamax = Rot2PolymerScratch(nrot_kinds,nchains,polymer,rotamers,parms,iamax);
			SetRotamersSimilarToPDB(nchains,polymer,a,na);
		}
		// create sidechain from pdb
		else 	
		{
			iamax = CreateSidechainFromPDB(a,nchains,na,polymer,parms,iamax);
			if (parms->rotamers) Rot2Polymer(nrot_kinds,nchains,polymer,rotamers,parms);
			if (!parms->pdb_rot) SetRotamersSimilarToPDB(nchains,polymer,a,na);
		}
	}

	fprintf(stderr,"Backbone atoms found for the different chains: ");
	for(i=0;i<nchains;i++) fprintf(stderr,"%d ",(polymer+i)->nback);
	fprintf(stderr,"\n");

	fprintf(stderr,"Transferred in polymer %d atoms\n",iamax);


	return iamax;
}


/****************************************************************************
 Generate backbone from pdb
 *****************************************************************************/
int CreateBackboneFromPDB(struct atom_s *a, int nchains, int na, struct s_polymer *polymer, struct s_parms *parms, int *ib)
{
	int i,ic,j;//iamax=0;

	for(i=0;i<nchains;i++) ib[i]=0;

	// add backbone
	for(i=0;i<na;i++)
	{
		if (IsBackbone((a+i)->atom,parms))											// if i is backbone atom
		{
			ic = (a+i)->chain;
			for(j=0;j<parms->n_back_a;j++)											// check if that atom can be moved
			if (!strcmp((a+i)->atom,parms->back_a[j])) (((polymer+ic)->back)+ib[ic])->move=parms->move_a[j];

			(((polymer+ic)->back)+ib[ic])->iaa = (a+i)->iaa;
			strcpy((((polymer+ic)->back)+ib[ic])->aa,(a+i)->aa);
			strcpy((((polymer+ic)->back)+ib[ic])->type,(a+i)->atom);
			(((polymer+ic)->back)+ib[ic])->ia = i;									// the atom id is taken from the pdb
			((((polymer+ic)->back)+ib[ic])->pos).x = ((a+i)->pos).x;
			((((polymer+ic)->back)+ib[ic])->pos).y = ((a+i)->pos).y;
			((((polymer+ic)->back)+ib[ic])->pos).z = ((a+i)->pos).z;

			*(((polymer+ic)->vback)+i) = &((((polymer+ic)->back)+ib[ic])->pos);

			if (parms->debug>1)
				fprintf(stderr,"\tbackbone=%3d\tatom=%3d %s\tch=%d\taa=%3d%s\tcor=%lf %lf %lf\n",ib[ic],i,(a+i)->atom,
						(a+i)->chain,(a+i)->iaa,(a+i)->aa,((a+i)->pos).x,((a+i)->pos).y,((a+i)->pos).z);

			(((polymer+ic)->back)+ib[ic])->nside = 0;
			(((polymer+ic)->back)+ib[ic])->nrot = 0;
			(((polymer+ic)->back)+ib[ic])->irot = 0;

			ib[ic]++;
		}
	}
	
	for(ic=0;ic<nchains;ic++)
		(polymer+ic)->nback = ib[ic];

	return i;
}

/****************************************************************************
 Generate sidechain from the atoms of the PDB (fills rotamer 0)
 *****************************************************************************/
int CreateSidechainFromPDB(struct atom_s *a, int nchains, int na, struct s_polymer *polymer, struct s_parms *parms, int iamax)
{
	int **top,i,ic,ib;

	// get topology of bonded interactions
	top = AlloIntMatrix(na,na);
	CreateTopology(a,na,top,parms->tthresh,parms->debug);		// matrix indexed on structure atom_s

	// add sidechains
	for(ic=0;ic<nchains;ic++)
	{
		ib = (polymer+ic)->nback;
		for(i=0;i<(polymer+ic)->nback;i++)
		{
			if (parms->debug>1) fprintf(stderr,"\tside of chain=%d back=%d atom=%d\n",ic,i,(((polymer+ic)->back)+i)->ia);

			if (i>0 && i<ib-1)  	// the reference system for the first sidechain (usually CB) is (i+1),(i-1),i
				FindSidechain(a,na,ic,i,(((polymer+ic)->back)+i)->ia,top,polymer,parms,(((polymer+ic)->back)+i+1)->ia,
						(((polymer+ic)->back)+i-1)->ia,(((polymer+ic)->back)+i)->ia);
			else if (i<ib-1)	 	// for the first is (i+2),(i+1),i
				FindSidechain(a,na,ic,i,(((polymer+ic)->back)+i)->ia,top,polymer,parms,(((polymer+ic)->back)+i+2)->ia,
						(((polymer+ic)->back)+i+1)->ia,(((polymer+ic)->back)+i)->ia);
			else if (i==ib-1)	// for the last is (nback-3),(nback-2),(nback-1)
				FindSidechain(a,na,ic,i,(((polymer+ic)->back)+i)->ia,top,polymer,parms,(((polymer+ic)->back)+ ib-3 )->ia,
						(((polymer+ic)->back)+ ib-2 )->ia, (((polymer+ic)->back)+ ib-1 )->ia);

		}
	}

	for (i=0;i<na;i++) free(*(top+i));

	free(top);

	return na;
}

int IsBackbone(char *atom, struct s_parms *p)
{
	int i;

	for(i=0;i<p->n_back_a;i++)
		if (!strcmp(atom,p->back_a[i])) return 1;

	return 0;
}

void CreateTopology(struct atom_s *a, int n, int **top, double thresh, int debug)
{
	int i,j,c=0;

	fprintf(stderr,"Create topology of bonded interactions\n");
	if (debug>2) fprintf(stderr,"Bonded topology:\n");

	for(i=0;i<n;i++)
		for(j=i+1;j<n;j++)
			if ( (a+i)->chain == (a+j)->chain )					// must be the same chain
			 if (Dist((a+i)->pos,(a+j)->pos)<thresh)			// bonded if dist<thresh
				 if ( !(!strcmp((a+i)->aa,"PRO") && !strcmp((a+i)->atom,"N") && !strcmp((a+j)->atom,"CD")) &&
						  !(!strcmp((a+i)->aa,"PRO") && !strcmp((a+j)->atom,"N") && !strcmp((a+i)->atom,"CD")) &&	// if PRO, assign backbone to CA even if it comes later
						  !(!strcmp((a+i)->aa,"TRP") && !strcmp((a+i)->atom,"CD2") && !strcmp((a+j)->atom,"CE3")) &&
						  !(!strcmp((a+i)->aa,"TRP") && !strcmp((a+j)->atom,"CD2") && !strcmp((a+i)->atom,"CE3"))  ) 	// if TRP, turn from the side of CE3 (like rotamers.lib)
				 {
					top[i][j]=1;
					top[j][i]=1;
					c++;

					if (debug>2) fprintf(stderr,"\t%d%s(%d%s) - %d%s(%d%s)\n",i,(a+i)->atom,(a+i)->iaa,(a+i)->aa,j,(a+j)->atom,(a+j)->iaa,(a+j)->aa);
				 }

	fprintf(stderr,"Found %d covalent bond to define the molecule topology\n",c);
}

/****************************************************************************
 Given the backbone atom w from a pdb, build recursively the rotamer
 structure in the s_polymer, according to the topology top.
 ib1,ib2,ib3 are the basis set for the next sidechain atom to set.
 *****************************************************************************/
void FindSidechain(struct atom_s *a, int natoms, int ipol, int iback, int w, int **top, struct s_polymer *p, struct s_parms *parms,
		int ib1, int ib2, int ib3)
{
	int i,j,ok,out;
	struct vector b1,b2,b3;

	for(i=0;i<natoms;i++)
		if (top[w][i]==1)
			if ( (a+i)->iaa == (a+w)->iaa && (a+i)->chain == (a+w)->chain )
				if ( IsBackbone((a+i)->atom,parms)==0 )
				{
					// check if already added
					ok=1;
					for(j=0;j<(((p+ipol)->back)+iback)->nside;j++)
						if ( (((((p+ipol)->back)+iback)->side)+j)->ia == i ) ok = 0;

					if (ok==1)
					{
						// fill s_polymer with it
						(((((p+ipol)->back)+iback)->side)+ (((p+ipol)->back)+iback)->nside )->ia = i;
						((((((p+ipol)->back)+iback)->side)+ (((p+ipol)->back)+iback)->nside )->pos).x = ((a+i)->pos).x;
						((((((p+ipol)->back)+iback)->side)+ (((p+ipol)->back)+iback)->nside )->pos).y = ((a+i)->pos).y;
						((((((p+ipol)->back)+iback)->side)+ (((p+ipol)->back)+iback)->nside )->pos).z = ((a+i)->pos).z;
						strcpy( (((((p+ipol)->back)+iback)->side)+ (((p+ipol)->back)+iback)->nside )->type, (a+i)->atom );
						*(((p+ipol)->vback)+i) = &((((((p+ipol)->back)+iback)->side)+ (((p+ipol)->back)+iback)->nside )->pos);

						// set in 0th rotamer the three preceding atoms
						(((((((p+ipol)->back)+iback)->side)+ (((p+ipol)->back)+iback)->nside )->rot)+0)->b1 = ib1;
						(((((((p+ipol)->back)+iback)->side)+ (((p+ipol)->back)+iback)->nside )->rot)+0)->b2 = ib2;
						(((((((p+ipol)->back)+iback)->side)+ (((p+ipol)->back)+iback)->nside )->rot)+0)->b3 = ib3;

						if (parms->debug>1)
							fprintf(stderr,"\t\tside atom=%d %s\tbasis=%d %d %d\n",i,(a+i)->atom,ib1,ib2,ib3);

						// set the spherical coordinates
						b1 = (a+ib1)->pos;
						b2 = (a+ib2)->pos;
						b3 = (a+ib3)->pos;
						if (parms->debug>1)
							fprintf(stderr,"\t\t\tx=(%lf %lf %lf)\n\t\t\tb1=(%lf %lf %lf)\n\t\t\tb2=(%lf %lf %lf\n\t\t\tb3=(%lf %lf %lf)\n",
									((a+i)->pos).x,((a+i)->pos).y,((a+i)->pos).z,b1.x,b1.y,b1.z,b2.x,b2.y,b2.z,b3.x,b3.y,b3.z);

						(((((((p+ipol)->back)+iback)->side)+ (((p+ipol)->back)+iback)->nside )->rot)+0)->ang = Cartesian2Spherical(b1,b2,
                                													b3, (a+i)->pos, (p+ipol)->tables, &out);
						if (out==0) Error("Cannot get spherical coordinates in Findsidechain");

						if (parms->debug>1)
							fprintf(stderr,"\t\t\tang=(%lf %lf %lf)\n",((((((((p+ipol)->back)+iback)->side)+ 
										(((p+ipol)->back)+iback)->nside )->rot)+0)->ang).ang,
										((((((((p+ipol)->back)+iback)->side)+ (((p+ipol)->back)+iback)->nside )->rot)+0)->ang).dih,
										((((((((p+ipol)->back)+iback)->side)+ (((p+ipol)->back)+iback)->nside )->rot)+0)->ang).r);

						(((p+ipol)->back)+iback)->nside ++;
						if ( (((p+ipol)->back)+iback)->nside >= NSIDEMAX ) Error("NSIDEMAX too small");
						(((p+ipol)->back)+iback)->nrot = 1;

						// find its neighbours
						FindSidechain(a,natoms,ipol,iback,i,top,p,parms,ib2,ib3,i);
					}
				}

}

/****************************************************************************
 Distance (SLOW!!! DO not use in dynamics)
 *****************************************************************************/
double Dist(struct vector a, struct vector b)
{
	return sqrt((a.x-b.x)*(a.x-b.x) + (a.y-b.y)*(a.y-b.y) + (a.z-b.z)*(a.z-b.z));
}


/****************************************************************************
 Assign rotamer to be close to that in PDB, and calcualte cartesian coordinates 
 *****************************************************************************/
void SetRotamersSimilarToPDB(int nc, struct s_polymer *p, struct atom_s *a, int napdb)
{
	int ic,i,ir=0,irotmin,iaa,is,ipdb,found;	
	struct vector pdb[ROT_ATMAX],cpol[ROT_ATMAX];
	double rmsd2min=99999.,x;

	fprintf(stderr,"Set rotamers close to those of PDB\n");

	for(ic=0;ic<nc;ic++)
		for(i=0;i<(p+ic)->nback;i++)
			if (!strcmp( (((p+ic)->back)+i)->type , "CA" ) && (((p+ic)->back)+i)->nside>0)		// loop on all CA atoms
			{
				iaa = (((p+ic)->back)+i)->iaa;

				// find positions in the pdb
				for(is=1;is< (((p+ic)->back)+i)->nside; is++)				// loop on sides of polymer, except CB
				{
					found=0;
					for(ipdb=0;ipdb<napdb;ipdb++)
					{
						if ( (a+ipdb)->chain == ic && (a+ipdb)->iaa == iaa && 
							!strcmp((a+ipdb)->atom,((((((p+ic)->back)+i)->side)+is)->type) ) )
						{ CopyVector(&((a+ipdb)->pos),pdb+is-1); found=1; }		// this contains the corresponding positions in the pdb
					}
					if (found==0) Error("Cannot find atom in SetRotamersSimilarToPDB");
				}
				// find the best rotamer
				rmsd2min=99999.;
				irotmin=0;
				for(ir=0;ir< (((p+ic)->back)+i)->nrot;ir++)
				{
					(((p+ic)->back)+i)->irot = ir;

					AddSidechain(p,i,i,ic);
					for(is=1;is<(((p+ic)->back)+i)->nside; is++) CopyVector( &((((((p+ic)->back)+i)->side)+is)->pos),cpol+is-1 );
					x = DumbRMSD2( pdb, cpol, (((p+ic)->back)+i)->nside-1 );
					if (x<rmsd2min) { rmsd2min=x; irotmin=ir; };

				}	
				// set rotamer and cartesian coordinates
				(((p+ic)->back)+i)->irot = irotmin;	
				fprintf(stderr,"%d %d irotmin=%d/%d\n",ic,i,irotmin,(((p+ic)->back)+i)->nrot);
				AddSidechain(p,i,i,ic);
			}

}

/****************************************************************************
 Calculates RMSD between two sets of vectors without roto-translation 
 *****************************************************************************/
double DumbRMSD2(struct vector *a, struct vector *b, int n)
{
	int i;
	double rmsd2=0;

	for(i=0;i<n;i++)
	{
		rmsd2 += Dist2(a[i],b[i]);
	}

	return rmsd2/n;
}

/****************************************************************************
 Copy PDB structure
 *****************************************************************************/
void CopyPDB(struct atom_s *from, struct atom_s *to, int n)
{
	int i;

	for(i=0;i<n;i++)
	{
		strcpy( (to+i)->atom, (from+i)->atom );
		strcpy( (to+i)->aa, (from+i)->aa );
		(to+i)->chain = (from+i)->chain;
		(to+i)->iaa = (from+i)->iaa;
		CopyVector( &((from+i)->pos), &((to+i)->pos) );	
	}
}

/****************************************************************************
 Generate a CA, CACB, NCAC model for a pdb
 *****************************************************************************/
int SimplifyPDB(struct atom_s *x, int n, char *model)
{
	struct atom_s *y;
	int m;

	y = AlloAtoms(n);

	if (!strcmp(model,"CA")) m = PDB2CA(x,y,n);
	else if (!strcmp(model,"CACB")) m = PDB2CACB(x,y,n);
    else if (!strcmp(model,"NCAC")) m = PDB2NCAC(x,y,n);
	else Error("Model not defined in SimplifyPDB");

	CopyPDB(y,x,m);

	free(y);
	return m;
}

/****************************************************************************
 Read the file of the propensity alfa/beta
 *****************************************************************************/

void ReadPropensity(char *fname, struct s_potential *u)
{
    FILE *fp;
    fp = fopen(fname,"r");
    if (!fp) Error("Cannot open propensity file");
    fprintf(stderr,"Reading Propensity File....\n");
    
    int naa=0, i;
    char n[5], t[5], aux[50];
    float coil, alfa, beta;
    
    while(fgets(aux,500,fp)!=NULL)
    {
        if(fscanf(fp,"%d %s %s %f %f %f",&i,n,t,&coil,&alfa,&beta)==6)
    	{
            u->ab_propensity[0][naa]=alfa;
            u->ab_propensity[1][naa]=beta;
            naa++;
    	}
    }
    fclose(fp);
    
    return;
}
