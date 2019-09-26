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


#include "montegrappa.h"
#include "grappino.h"

/***********************************************************************
 * Allocate structure for temporary reading rotamers from file
 ***********************************************************************/
struct rot_input_s *AlloRotamer(int naa, FILE *flog)
{
  int i,j,size=0;
  struct rot_input_s *x;

  x = (struct rot_input_s *) calloc(naa,sizeof(struct rot_input_s));
  size += sizeof(struct rot_input_s);
  if (!x) Error("Cannot allocate rotamers");

  for (i=0;i<naa;i++)
    {
       ((x+i)->atom) = (struct rot_atom_s **) calloc(ROTMAX,sizeof(struct rot_atom_s *));
       size += sizeof(struct rot_atom_s *);
       if (!(((x+i)->atom))) Error("Cannot allocate rotamers (2)");
       for (j=0;j<ROTMAX;j++)
           {
              *(((x+i)->atom)+j) = (struct rot_atom_s *) calloc(ROT_ATMAX,sizeof(struct rot_atom_s));
              size += sizeof(struct rot_atom_s);
              if (!(*((x+i)->atom)+j)) Error("Cannot allocate rotamers (3)");
           }
    }

  fprintf(flog,"Allocate rotamer structure of size  %lf MB\n",(double)size/1024/1024);

  return x;
}

/***********************************************************************
 * Read rotamers from file, putting in rot_input_s and returns
 * the number of amino acid types
 ***********************************************************************/
int ReadRotamers(char *fname, struct rot_input_s *x, int debug)
{
  int i,j,nrot,nat,naa=0;
  double ang,dih,r;
  char aa[10],at[5],a1[5],a2[5],a3[5];
  FILE *fp;

  fp = fopen(fname,"r");
  if (!fp) Error("Cannot open rotamer file");

  if ( debug > 1)
    {
       fprintf(stderr,"Rotamers:\n");
       fprintf(stderr,"a   rot atm 	   A1  A2   A3	atm		ang	dih		r\n");
    }

  while (fscanf(fp,"%s %d %d",aa,&nrot,&nat)==3)
    {
       (x+naa)->nrot = nrot;
       (x+naa)->nat = nat;
       strcpy( (x+naa)->aa,aa);
       for (i=0;i<nrot;i++)
        for (j=0;j<nat;j++)
          {
             if (fscanf(fp,"%s %s %s %s %lf %lf %lf",a1,a2,a3,at,&ang,&dih,&r)!=7) Error("Error in reading rotamer file");
             else
             {
				 strcpy( (*(((x+naa)->atom)+i)+j)->A1 ,a1);
				 strcpy( (*(((x+naa)->atom)+i)+j)->A2 ,a2);
				 strcpy( (*(((x+naa)->atom)+i)+j)->A3 ,a3);
				 strcpy( (*(((x+naa)->atom)+i)+j)->kind, at);
				 (*(((x+naa)->atom)+i)+j)->ang = ang;
				 (*(((x+naa)->atom)+i)+j)->dih = dih;
				 (*(((x+naa)->atom)+i)+j)->r = r;
             }

          }

       if ( debug > 1)
         for (i=0;i<nrot;i++)
           for (j=0;j<nat;j++)
                fprintf(stderr,"%3s %2d %2d\t%4s %4s %4s\t%4s     %11lf %11lf %11lf\n",(x+naa)->aa,i,j,(*(((x+naa)->atom)+i)+j)->A1,(*(((x+naa)->atom)+i)+j)->A2,(*(((x+naa)->atom)+i)+j)->A3,(*(((x+naa)->atom)+i)+j)->kind,(*(((x+naa)->atom)+i)+j)->ang,(*(((x+naa)->atom)+i)+j)->dih,(*(((x+naa)->atom)+i)+j)->r);

       naa++;
       if (naa>=AAKINDMAX) Error("AAKINDMAX too small");
    }

  return naa;
}


/***********************************************************************
 * Add sidechains from scratch
 ***********************************************************************/
int Rot2PolymerScratch(int nrot_kinds, int nchains, struct s_polymer *p, struct rot_input_s *r, struct s_parms *parms, int nat)
{
	int i,is,ic,N[NATOMMAX],C[NATOMMAX],CA[NATOMMAX],natomback=0,aamax=-1;
	int k,rot_kind,irot,iaa=0,b1,b2,b3,m,iaa0;

	fprintf(stderr,"\nAdd sidechains from rotamers library (neglecting pdb)...\n");
	
	for (i=0;i<NRESMAX;i++) CA[i] = -1;
	for (ic=0;ic<nchains;ic++) natomback += (p+ic)->nback;

	// find main backbone atoms
	for (ic=0;ic<nchains;ic++)
	{
		iaa0 = (((p+ic)->back)+0)->iaa;			// id of first aa of this chain
		for (i=0;i<NRESMAX;i++) { CA[i] = -1; C[i] = -1; N[i] = -1; } 

		// find the iback corresponding to N, CA, C
		for (i=0;i<(p+ic)->nback;i++)
		{
			if ( (((p+ic)->back)+i)->iaa -iaa0 > NRESMAX ) Error("NRESMAX too small in Rot2PolymerScratch");
			if (!strcmp( (((p+ic)->back)+i)->type,"N") )  N[ (((p+ic)->back)+i)->iaa -iaa0] =  (((p+ic)->back)+i)->ia;
			if (!strcmp( (((p+ic)->back)+i)->type,"C") )  C[ (((p+ic)->back)+i)->iaa -iaa0] =  (((p+ic)->back)+i)->ia;
			if (!strcmp( (((p+ic)->back)+i)->type,"CA") )  CA[ (((p+ic)->back)+i)->iaa -iaa0] =  (((p+ic)->back)+i)->ia;
			aamax = (((p+ic)->back)+i)->iaa - iaa0;
		}
		if (parms->debug>2) fprintf(stderr,"ic=%d\tiaa0=%d naa=%d\n",ic,iaa0,aamax+1);

		// check if all aa have N, CA, C
		for (i=0;i<=aamax;i++)
		{
			if ( N[ (((p+ic)->back)+i)->iaa -iaa0 ] == -1 ) { fprintf(stderr,"chain=%d amino acid=%d\n",ic,(((p+ic)->back)+i)->iaa); 
										Error("Cannot find N atom in Rot2PolymerScratch"); }
			if ( C[ (((p+ic)->back)+i)->iaa -iaa0] == -1 ) { fprintf(stderr,"chain=%d amino acid=%d\n",ic,(((p+ic)->back)+i)->iaa); 
										Error("Cannot find C atom in Rot2PolymerScratch"); }
			if ( CA[ (((p+ic)->back)+i)->iaa-iaa0] == -1 ) { fprintf(stderr,"chain=%d amino acid=%d\n",ic,(((p+ic)->back)+i)->iaa); 
										Error("Cannot find CA atom in Rot2PolymerScratch"); }
		}

		for (i=0;i<(p+ic)->nback;i++)
		{
			(((p+ic)->back)+i)->nside = 0;
		}

		// add the CB
		nat = AddDefaultCB(p,ic,N,CA,C,aamax+1,nat,parms->debug);

		// add sidechain atoms other than CB
		fprintf(stderr,"Adding the rotamers of the atoms other than CB\n");
		for (i=0;i<(p+ic)->nback;i++)
			if (!strcmp( (((p+ic)->back)+i)->type , "CA" ))										// loop on all CA atoms
			{
				if (parms->debug>1) fprintf(stderr,"Back=%3d\tchain=%d\taa=%s\t",i,ic,(((p+ic)->back)+i)->aa);

				// find which amino acid (i.e. which element of rotamer structure)
				iaa = (((p+ic)->back)+i)->iaa -iaa0;

				if (iaa>NRESMAX) Error("NRESMAX too small in Rot2PolymerScratch");
				rot_kind=-1;
				for (k=0;k<nrot_kinds;k++)
					if (!strcmp((((p+ic)->back)+i)->aa,(r+k)->aa)) rot_kind = k;
				if (rot_kind==-1 && parms->debug >1) fprintf(stderr,"  no rotamer");
				if (parms->debug>1) fprintf(stderr,"\tRotamer id is %d\n",rot_kind);


				//propagate CB, which is defined only for rotamer 0, to all rotamers
				if (strcmp( (((p+ic)->back)+i)->aa , "GLY" ))
					for (irot=1;irot< (r+rot_kind)->nrot ;irot++)
					{
						(((((((p+ic)->back)+i)->side)+0)->rot)+irot)->b1 = (((((((p+ic)->back)+i)->side)+0)->rot)+0)->b1;
						(((((((p+ic)->back)+i)->side)+0)->rot)+irot)->b2 = (((((((p+ic)->back)+i)->side)+0)->rot)+0)->b2;
						(((((((p+ic)->back)+i)->side)+0)->rot)+irot)->b3 = (((((((p+ic)->back)+i)->side)+0)->rot)+0)->b3;
						((((((((p+ic)->back)+i)->side)+0)->rot)+irot)->ang).ang = ((((((((p+ic)->back)+i)->side)+0)->rot)+0)->ang).ang;
						((((((((p+ic)->back)+i)->side)+0)->rot)+irot)->ang).dih = ((((((((p+ic)->back)+i)->side)+0)->rot)+0)->ang).dih;
						((((((((p+ic)->back)+i)->side)+0)->rot)+irot)->ang).r = ((((((((p+ic)->back)+i)->side)+0)->rot)+0)->ang).r;
					}


				// add sidechain atoms
				if (rot_kind>-1)
				for (is=0;is<(r+rot_kind)->nat;is++)												// over atoms to be put
				{
					// append new sidechain	
					(((((p+ic)->back)+i)->side)+is+1)->ia = nat;
					*(((p+ic)->vback)+nat) = &((((((p+ic)->back)+i)->side)+is+1)->pos);				// is+1 because is=0 is CB

					strcpy( (((((p+ic)->back)+i)->side)+is+1)->type , (*(((r+rot_kind)->atom)+0)+is)->kind );		// assumes that all rotamers have the same atom (the 0th)
					if (parms->debug>1) fprintf(stderr,"\tAtom=%d (%d of %d)\tnside=%d",nat,is,(r+rot_kind)->nat,(((p+ic)->back)+i)->nside);

					for (irot=0;irot<(r+rot_kind)->nrot;irot++)									// over rotamers
					{ 
						((((((((p+ic)->back)+i)->side)+is+1)->rot)+irot)->ang).ang = (*(((r+rot_kind)->atom)+irot)+is)->ang;
						((((((((p+ic)->back)+i)->side)+is+1)->rot)+irot)->ang).dih = (*(((r+rot_kind)->atom)+irot)+is)->dih;
						((((((((p+ic)->back)+i)->side)+is+1)->rot)+irot)->ang).r = (*(((r+rot_kind)->atom)+irot)+is)->r;
						if (parms->debug>1) fprintf(stderr,"\t\tirot=%d\tang=%lf dih=%lf r=%lf\n",irot,(*(((r+rot_kind)->atom)+irot)+is)->ang,
								(*(((r+rot_kind)->atom)+irot)+is)->dih,(*(((r+rot_kind)->atom)+irot)+is)->r);

						b1 = -1;					// find b1
						if ( !strcmp((*(((r+rot_kind)->atom)+irot)+is)->A1 , "N") ) b1 = N[iaa];
						else if ( !strcmp((*(((r+rot_kind)->atom)+irot)+is)->A1 , "CA") ) b1 = CA[iaa];
						else if ( !strcmp((*(((r+rot_kind)->atom)+irot)+is)->A1 , "C") ) b1 = C[iaa];
						else
							for (m=0;m<(((p+ic)->back)+i)->nside;m++)							// over atoms of the sidechain
								if (!strcmp( (*(((r+rot_kind)->atom)+irot)+is)->A1, (((((p+ic)->back)+i)->side)+m)->type ))
									b1 = (((((p+ic)->back)+i)->side)+m)->ia;					// b1 is the atom id of the basis set
						if (b1 == -1) { fprintf(stderr,"ic=%d i=%d is=%d irot=%d A1=%s\n",ic,i,is,irot,(*(((r+rot_kind)->atom)+irot)+is)->A1); Error("Cannot find atom to build rotamer (1)"); }
						if (parms->debug>1) fprintf(stderr,"\t\tA1=%s b1=%d\n",(*(((r+rot_kind)->atom)+irot)+is)->A1,b1);

						b2 = -1;					// find b2
						if ( !strcmp((*(((r+rot_kind)->atom)+irot)+is)->A2 , "N") ) b2 = N[iaa];
						else if ( !strcmp((*(((r+rot_kind)->atom)+irot)+is)->A2 , "CA") ) b2 = CA[iaa];
						else if ( !strcmp((*(((r+rot_kind)->atom)+irot)+is)->A2 , "C") ) b2 = C[iaa];
						else
							for (m=0;m<(((p+ic)->back)+i)->nside;m++)							// over atoms of the sidechain
								if (!strcmp( (*(((r+rot_kind)->atom)+irot)+is)->A2, (((((p+ic)->back)+i)->side)+m)->type ))
									b2 = (((((p+ic)->back)+i)->side)+m)->ia;					// b2 is the atom id of the basis set
						if (b2 == -1) { fprintf(stderr,"ic=%d i=%d is=%d irot=%d A2=%s\n",ic,i,is,irot,(*(((r+rot_kind)->atom)+irot)+is)->A2); Error("Cannot find atom to build rotamer (2)"); }
						if (parms->debug>1) fprintf(stderr,"\t\tA2=%s b2=%d\n",(*(((r+rot_kind)->atom)+irot)+is)->A2,b2);

						b3 = -1;					// find b3
						if ( !strcmp((*(((r+rot_kind)->atom)+irot)+is)->A3 , "N") ) b3 = N[iaa];
						else if ( !strcmp((*(((r+rot_kind)->atom)+irot)+is)->A3 , "CA") ) b3 =CA[iaa];
						else if ( !strcmp((*(((r+rot_kind)->atom)+irot)+is)->A3 , "C") ) b3 = C[iaa];
						else
							for (m=0;m<(((p+ic)->back)+i)->nside;m++)							// over atoms of the sidechain
								if (!strcmp( (*(((r+rot_kind)->atom)+irot)+is)->A3, (((((p+ic)->back)+i)->side)+m)->type ))
										b3 = (((((p+ic)->back)+i)->side)+m)->ia;				// b3 is the atom id of the basis set
						if (b3 == -1) { fprintf(stderr,"ic=%d i=%d is=%d irot=%d A3=%s\n",ic,i,is,irot,(*(((r+rot_kind)->atom)+irot)+is)->A3); Error("Cannot find atom to build rotamer (3)"); }
						if (parms->debug>1) fprintf(stderr,"\t\tA3=%s b3=%d\n",(*(((r+rot_kind)->atom)+irot)+is)->A3,b3);
						if (parms->debug>1) fprintf(stderr,"\t\ttype=%s\n",(((((p+ic)->back)+i)->side)+is+1)->type);

						(((((((p+ic)->back)+i)->side)+is+1)->rot)+irot)->b1 = b1;
						(((((((p+ic)->back)+i)->side)+is+1)->rot)+irot)->b2 = b2;
						(((((((p+ic)->back)+i)->side)+is+1)->rot)+irot)->b3 = b3;
					}
					(((p+ic)->back)+i)->nside ++;
					nat ++;
				}

				(((p+ic)->back)+i)->nrot = (r+rot_kind)->nrot;
				if (parms->debug>1) fprintf(stderr,"\tnrot=%d\tnside=%d (rot id=%d)\n",(r+rot_kind)->nrot,(((p+ic)->back)+i)->nside,rot_kind );
			}
	}
	// careful : still needs to be converted in cartesian coordinates

	return nat;
}

/***********************************************************************
 * Add CB to CA according to default polar coordinates (only in rotamer 0)
 ***********************************************************************/
int AddDefaultCB(struct s_polymer *p, int ic, int *N, int *CA, int *C, int naa, int nat, int debug)
{
	int i,out;

	fprintf(stderr,"Building default CB for chain %d\n",ic);

	for (i=0;i<naa;i++)
		if (CA[i] != -1 && strcmp((((p+ic)->back)+CA[i])->aa,"GLY") )
		{
			(((((p+ic)->back)+CA[i])->side)+ 0)->ia = nat;
			strcpy((((((p+ic)->back)+CA[i])->side)+ 0)->type,"CB");
			if (!strcmp((((p+ic)->back)+CA[i])->aa,"PRO"))
			{
				((((((((p+ic)->back)+CA[i])->side)+ 0)->rot)+0)->ang).ang = ROT_CB_ANG_PRO;
				((((((((p+ic)->back)+CA[i])->side)+ 0)->rot)+0)->ang).dih = ROT_CB_DIH_PRO;
				((((((((p+ic)->back)+CA[i])->side)+ 0)->rot)+0)->ang).r = ROT_CB_RADIUS_PRO;
			}
			else
			{
				((((((((p+ic)->back)+CA[i])->side)+ 0)->rot)+0)->ang).ang = ROT_CB_ANG;
				((((((((p+ic)->back)+CA[i])->side)+ 0)->rot)+0)->ang).dih = ROT_CB_DIH;
				((((((((p+ic)->back)+CA[i])->side)+ 0)->rot)+0)->ang).r = ROT_CB_RADIUS;
			}
			(((p+ic)->back)+CA[i])->nside = 1;
			(((p+ic)->back)+CA[i])->nrot = 1;
			(((((((p+ic)->back)+CA[i])->side)+ 0)->rot)+0)->b1 = (((p+ic)->back)+C[i])->ia;
			(((((((p+ic)->back)+CA[i])->side)+ 0)->rot)+0)->b2 = (((p+ic)->back)+N[i])->ia;
			(((((((p+ic)->back)+CA[i])->side)+ 0)->rot)+0)->b3 = (((p+ic)->back)+CA[i])->ia;
                                                                                                                                                   
			// set cartesian coordinates for CB
			(((((p+ic)->back)+CA[i])->side)+ 0)->pos = Spherical2Cartesian( (((p+ic)->back)+C[i])->pos ,
					(((p+ic)->back)+N[i])->pos , (((p+ic)->back)+CA[i])->pos , (((((((p+ic)->back)+CA[i])->side)+ 0)->rot)+0)->ang, p->tables,&out);
			if (out==0) Error("Cannot find cartesian coordinates of CB");

			// set position lookback table for CB
			*(((p+ic)->vback)+nat) = &( (((((p+ic)->back)+CA[i])->side)+0)->pos );


			if (debug>1) { fprintf(stderr,"aa=%d\tCB  b1=%3d b2=%3d b3=%3d   ia=%d  ",i,(((((((p+ic)->back)+CA[i])->side)+ 0)->rot)+0)->b1,
				(((((((p+ic)->back)+CA[i])->side)+ 0)->rot)+0)->b2,(((((((p+ic)->back)+CA[i])->side)+ 0)->rot)+0)->b3,(((((p+ic)->back)+CA[i])->side)+ 0)->ia);
				PrintVector(stderr,"CB", (((((p+ic)->back)+CA[i])->side)+0)->pos); }

			nat ++;
		}

	return nat;
}

/***********************************************************************
 * Substitute default CB ine xisting PDB (only in rotamer 0)
 ***********************************************************************/
void SubstituteDefaultCB(struct s_polymer *p, int nc)
{
	int i,ic,is,out;

	fprintf(stderr,"Substituting default CB\n");

	for (ic=0;ic<nc;ic++)
		for (i=0;i<(p+ic)->nback;i++)
			for (is=0;is<(((p+ic)->back)+i)->nside;is++)
				if (!strcmp((((((p+ic)->back)+i)->side)+is)->type,"CB"))
				{
					if (!strcmp((((p+ic)->back)+i)->aa,"PRO"))
					{
						((((((((p+ic)->back)+i)->side)+ is)->rot)+0)->ang).ang = ROT_CB_ANG_PRO;
						((((((((p+ic)->back)+i)->side)+ is)->rot)+0)->ang).dih = ROT_CB_DIH_PRO;
						((((((((p+ic)->back)+i)->side)+ is)->rot)+0)->ang).r = ROT_CB_RADIUS_PRO;
					}
					else
					{
						((((((((p+ic)->back)+i)->side)+ is)->rot)+0)->ang).ang = ROT_CB_ANG;
						((((((((p+ic)->back)+i)->side)+ is)->rot)+0)->ang).dih = ROT_CB_DIH;
						((((((((p+ic)->back)+i)->side)+ is)->rot)+0)->ang).r = ROT_CB_RADIUS;
					}

					(((((p+ic)->back)+i)->side)+is)->pos = Spherical2Cartesian( (((p+ic)->back)+i+1)->pos ,
							(((p+ic)->back)+i-1)->pos , (((p+ic)->back)+i)->pos , (((((((p+ic)->back)+i)->side)+is)->rot)+0)->ang, p->tables,&out);
					if (out==0) Error("Cannot find cartesian coordinates of CB");
				}

	return ;
}

/***********************************************************************
 *
 ***********************************************************************/
int Rot2Polymer(int nrot_kinds, int nchains, struct s_polymer *p, struct rot_input_s *r, struct s_parms *parms)
{
	int i,j,ic,N[NATOMMAX],C[NATOMMAX],CA[NATOMMAX],nat=0;
	int k,rot_kind,irot,iaa,b1,b2,b3,m,iar,keeppdb=0;

	fprintf(stderr,"\nAdd rotamers from rotamers library...\n");

	for (i=0;i<NRESMAX;i++) CA[i] = -1;

	// find main backbone atoms
	for (ic=0;ic<nchains;ic++)
		for (i=0;i<(p+ic)->nback;i++)
		{
			if (!strcmp( (((p+ic)->back)+i)->type,"N") ) { N[ (((p+ic)->back)+i)->iaa ] = (((p+ic)->back)+i)->ia; }
			if (!strcmp( (((p+ic)->back)+i)->type,"C") ) { C[ (((p+ic)->back)+i)->iaa ] = (((p+ic)->back)+i)->ia; }
			if (!strcmp( (((p+ic)->back)+i)->type,"CA") ) { CA[ (((p+ic)->back)+i)->iaa ] = (((p+ic)->back)+i)->ia; }
		}

	// if required, change CB
	if (!parms->cb_pdb) SubstituteDefaultCB(p,nchains);

	// if required, keep the sidechain of the pdb as rotamer 0
	if (parms->pdb_rot) keeppdb=1;

	// add sidechain atoms other than CB
	fprintf(stderr,"Adding the rotamers of the atoms other than CB\n");
	for (ic=0;ic<nchains;ic++)
		for (i=0;i<(p+ic)->nback;i++)
			if (!strcmp( (((p+ic)->back)+i)->type , "CA" ))										// loop on all CA atoms
			{
				if (parms->debug>1) fprintf(stderr,"Back=%3d\taa=%s\t",i,(((p+ic)->back)+i)->aa);
				// find which amino acid (i.e. which element of rotamer structure)
				iaa = (((p+ic)->back)+i)->iaa;
				rot_kind=-1;
				for (k=0;k<nrot_kinds;k++)
					if (!strcmp((((p+ic)->back)+i)->aa,(r+k)->aa)) rot_kind = k;
				if (rot_kind==-1 && parms->debug >1) fprintf(stderr,"  no rotamer");
				if (parms->debug>1 && (r+rot_kind)->nat >0) fprintf(stderr,"\tRotamer id is %d\n",rot_kind);

				//propagate CB, which in the pdb are defined only for rotamer 0, to all  rotamers
				if ( rot_kind>-1)
					for (irot=1;irot< (r+rot_kind)->nrot + keeppdb;irot++)
					{
						(((((((p+ic)->back)+i)->side)+0)->rot)+irot)->b1 = (((((((p+ic)->back)+i)->side)+0)->rot)+0)->b1;
						(((((((p+ic)->back)+i)->side)+0)->rot)+irot)->b2 = (((((((p+ic)->back)+i)->side)+0)->rot)+0)->b2;
						(((((((p+ic)->back)+i)->side)+0)->rot)+irot)->b3 = (((((((p+ic)->back)+i)->side)+0)->rot)+0)->b3;
						((((((((p+ic)->back)+i)->side)+0)->rot)+irot)->ang).ang = ((((((((p+ic)->back)+i)->side)+0)->rot)+0)->ang).ang;
						((((((((p+ic)->back)+i)->side)+0)->rot)+irot)->ang).dih = ((((((((p+ic)->back)+i)->side)+0)->rot)+0)->ang).dih;
						((((((((p+ic)->back)+i)->side)+0)->rot)+irot)->ang).r = ((((((((p+ic)->back)+i)->side)+0)->rot)+0)->ang).r;
					}

				(((p+ic)->back)+i)->nrot = (r+rot_kind)->nrot + keeppdb;

				// add new spherical coordinates to rotamers
				for (j=1;j<(((p+ic)->back)+i)->nside;j++)						// over atoms to be put except CB (j=0)
					if (strcmp((((((p+ic)->back)+i)->side)+j)->type,"CB") && rot_kind > -1)
					{
						//if (parms->debug>1) fprintf(stdout,"\tSide atom=%d (%d of %d)\tnside=%d",nat,j,(r+rot_kind)->nat,(((p+ic)->back)+i)->nside);

						for (irot=0;irot<(r+rot_kind)->nrot;irot++)									// over rotamers
						{
							// find atom within rotamer
							iar=-1;
							for (m=0;m<(r+rot_kind)->nat;m++)										// over atoms of the rotamer
								if (!strcmp( (*(((r+rot_kind)->atom)+irot)+m)->kind, (((((p+ic)->back)+i)->side)+j)->type ))
									iar = m;														// iar is the correct element of rotamer array
							if (iar == -1)
							{
								fprintf(stderr,"WARNING: cannot find atom %s of amino acid %s in rotamer library\n",(((((p+ic)->back)+i)->side)+j)->type,(((p+ic)->back)+i)->aa);
								exit(1);
							}

							((((((((p+ic)->back)+i)->side)+j)->rot)+irot+keeppdb)->ang).ang = (*(((r+rot_kind)->atom)+irot)+iar)->ang;
							((((((((p+ic)->back)+i)->side)+j)->rot)+irot+keeppdb)->ang).dih = (*(((r+rot_kind)->atom)+irot)+iar)->dih;
							((((((((p+ic)->back)+i)->side)+j)->rot)+irot+keeppdb)->ang).r = (*(((r+rot_kind)->atom)+irot)+iar)->r;

							if (parms->debug>1) fprintf(stderr,"\t\tirot=%d\tang=%lf dih=%lf r=%lf\n",irot+keeppdb,(*(((r+rot_kind)->atom)+irot+keeppdb)+iar)->ang,
									(*(((r+rot_kind)->atom)+irot+keeppdb)+iar)->dih,(*(((r+rot_kind)->atom)+irot+keeppdb)+iar)->r);

							b1 = -1;					// find b1
							if ( !strcmp((*(((r+rot_kind)->atom)+irot)+iar)->A1 , "N") ) b1 = N[iaa];
							else if ( !strcmp((*(((r+rot_kind)->atom)+irot)+iar)->A1 , "CA") ) b1 =CA[iaa];
							else if ( !strcmp((*(((r+rot_kind)->atom)+irot)+iar)->A1 , "C") ) b1 =C[iaa];
							else
								for (m=0;m<(((p+ic)->back)+i)->nside;m++)							// over atoms of the sidechain
									if (!strcmp( (*(((r+rot_kind)->atom)+irot)+iar)->A1, (((((p+ic)->back)+i)->side)+m)->type ))
										b1 = (((((p+ic)->back)+i)->side)+m)->ia;					// b1 is the atom id of the basis set
							if (b1 == -1) { fprintf(stderr,"ic=%d i=%d iar=%d irot=%d A1=%s\n",ic,i,iar,irot,(*(((r+rot_kind)->atom)+irot)+iar)->A1); Error("Cannot find atom to build rotamer (1)"); }
							if (parms->debug>1) fprintf(stderr,"\t\tA1=%s b1=%d\n",(*(((r+rot_kind)->atom)+irot)+iar)->A1,b1);

							b2 = -1;					// find b2
							if ( !strcmp((*(((r+rot_kind)->atom)+irot)+iar)->A2 , "N") ) b2 = N[iaa];
							else if ( !strcmp((*(((r+rot_kind)->atom)+irot)+iar)->A2 , "CA") ) b2 = CA[iaa];
							else if ( !strcmp((*(((r+rot_kind)->atom)+irot)+iar)->A2 , "C") ) b2 = C[iaa];
							else
								for (m=0;m<(((p+ic)->back)+i)->nside;m++)							// over atoms of the sidechain
									if (!strcmp( (*(((r+rot_kind)->atom)+irot)+iar)->A2, (((((p+ic)->back)+i)->side)+m)->type ))
										b2 = (((((p+ic)->back)+i)->side)+m)->ia;					// b1 is the atom id of the basis set
							if (b2 == -1) { fprintf(stderr,"ic=%d i=%d iar=%d irot=%d A2=%s\n",ic,i,iar,irot,(*(((r+rot_kind)->atom)+irot)+iar)->A2); Error("Cannot find atom to build rotamer (1)"); }
							if (parms->debug>1) fprintf(stderr,"\t\tA2=%s b2=%d\n",(*(((r+rot_kind)->atom)+irot)+iar)->A2,b2);


							b3 = -1;					// find b3
							if ( !strcmp((*(((r+rot_kind)->atom)+irot)+iar)->A3 , "N") ) b3 = N[iaa];
							else if ( !strcmp((*(((r+rot_kind)->atom)+irot)+iar)->A3 , "CA") ) b3 = CA[iaa];
							else if ( !strcmp((*(((r+rot_kind)->atom)+irot)+iar)->A3 , "C") ) b3 = C[iaa];
							else
								for (m=0;m<(((p+ic)->back)+i)->nside;m++)							// over atoms of the sidechain
									if (!strcmp( (*(((r+rot_kind)->atom)+irot)+iar)->A3, (((((p+ic)->back)+i)->side)+m)->type ))
										b3 = (((((p+ic)->back)+i)->side)+m)->ia;					// b1 is the atom id of the basis set
							if (b3 == -1) { fprintf(stderr,"ic=%d i=%d iar=%d irot=%d A3=%s\n",ic,i,iar,irot,(*(((r+rot_kind)->atom)+irot)+iar)->A3); Error("Cannot find atom to build rotamer (1)"); }
							if (parms->debug>1) fprintf(stderr,"\t\tA3=%s b3=%d\n",(*(((r+rot_kind)->atom)+irot)+iar)->A3,b3);

							(((((((p+ic)->back)+i)->side)+j)->rot)+irot+keeppdb)->b1 = b1;
							(((((((p+ic)->back)+i)->side)+j)->rot)+irot+keeppdb)->b2 = b2;
							(((((((p+ic)->back)+i)->side)+j)->rot)+irot+keeppdb)->b3 = b3;

					}
				}

				if (parms->debug>1) fprintf(stderr,"\tnrot=%d\tnside=%d (rot id=%d)\n",(((p+ic)->back)+i)->nrot,(((p+ic)->back)+i)->nside,rot_kind );
			}

	return nat;
}
