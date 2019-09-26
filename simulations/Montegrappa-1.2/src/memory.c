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
 * memory.c
 *
 *  Created on: Sep 15, 2010
 *      Author: guido
 */

#include "montegrappa.h"

/*****************************************************************************
 Allocate npol polymer structures
 *****************************************************************************/
struct s_polymer *AlloPolymer(int npol, int nback, int nside, int nrot, int natoms, int shell, int noside, FILE *flog)
{
	int ires,iside,ipol,size;
	struct s_polymer *p;

	p = (struct s_polymer *) calloc(npol,sizeof(struct s_polymer));
	size = npol*sizeof(struct s_polymer);
	if (!p) Error("Cannot allocate polymer");

	for (ipol=0;ipol<npol;ipol++)
	{
		(p+ipol)->back = (struct s_back *) calloc(nback,sizeof(struct s_back));
		if (!(p+ipol)->back) Error("Cannot allocate backbone");
		size += nback * sizeof(struct s_back);

		(p+ipol)->vback = (struct vector **) calloc(natoms,sizeof(struct vector *));
		if (!(p+ipol)->vback) Error("Cannot allocate vback");
		size += nback * sizeof(struct vector *);

		// allocate sidechains
		for (ires=0;ires<nback;ires++)
		{
			(((p+ipol)->back)+ires)->nside = 0;
			(((p+ipol)->back)+ires)->nrot = 0;
			if (!noside)
			{
				(((p+ipol)->back)+ires)->side = (struct s_side *) calloc(nside,sizeof(struct s_side));
				size += nside*sizeof(struct s_side);
				if (!(((p+ipol)->back)+ires)->side) Error("Cannot alocate  sidechain");
			}
			(((p+ipol)->back)+ires)->contacts = (int *) calloc(NCONTMAX,sizeof(int));
			size += NCONTMAX*sizeof(int);
			(((p+ipol)->back)+ires)->contacts_p = (int *) calloc(NCONTMAX,sizeof(int));
			size += NCONTMAX*sizeof(int);
			(((p+ipol)->back)+ires)->e = (double *) calloc(NCONTMAX,sizeof(double));
			size += NCONTMAX*sizeof(double);
			if (shell)
			{
				(((p+ipol)->back)+ires)->shell = (int *) calloc(NSHELLMAX,sizeof(int));
				size += NCONTMAX*sizeof(int);
				(((p+ipol)->back)+ires)->shell_p = (int *) calloc(NSHELLMAX,sizeof(int));
				size += NCONTMAX*sizeof(int);
			}
			(((p+ipol)->back)+ires)->ncontacts = 0;
			(((p+ipol)->back)+ires)->nshell = 0;
			(((p+ipol)->back)+ires)->move = 1;


			if (!noside)
			{
				for (iside=0;iside<nside;iside++)
				{
					(((((p+ipol)->back)+ires)->side)+iside)->rot = (struct s_rotamers *) calloc(nrot,sizeof(struct s_rotamers));
					size += nrot*sizeof(struct s_rotamers);
				}

				for (iside=0;iside<natoms;iside++) ((p+ipol)->vback)[iside] = NULL;
			}
		}

	}

	if (noside) { nside=0; nrot=0; }
//	fprintf(stderr,"Allocate polymer structure of size  %lf MB (nchains=%d, nback=%d, nside=%d, nrot=%d)\n",(double)size/1024/1024,npol,nback,nside,nrot);


	return p;
}

/*****************************************************************************
 Initialize tables to calculate sqrt, sin, cos, exp
 *****************************************************************************/
struct s_tables *InitTables(FILE *fp)
{
      struct s_tables *t;
      int i,size=0;

      // Allocate overall table structure
      t = (struct s_tables *) calloc(1,sizeof(struct s_tables));
      if (!t) Error("cannot allocate table");

      // Sqrt
      t->fast_sqrt = (double *) calloc(FSQRT_L,sizeof(double));
      if (!t->fast_sqrt) Error("cannot allocate sqrt in table");
      size += FSQRT_L*sizeof(double);

      for (i=0;i<FSQRT_L;i++)  *((t->fast_sqrt)+i) = sqrt( (double)i / FSQRT_BINRC );

      // Sin and Cos
      t->fast_sin = (double *) calloc(FTRIG_L,sizeof(double));
      if (!t->fast_sin) Error("cannot allocate sin in table");
      size += FTRIG_L*sizeof(double);
      t->fast_cos = (double *) calloc(FTRIG_L,sizeof(double));
      if (!t->fast_cos) Error("cannot allocate cos in table");
      size += FTRIG_L*sizeof(double);

      for (i=0;i<FTRIG_L;i++)  *((t->fast_sin)+i) = sin( (double)i / FTRIG_BINRC / 180. * PI);
      for (i=0;i<FTRIG_L;i++)  *((t->fast_cos)+i) = cos( (double)i / FTRIG_BINRC / 180. * PI);

      // Exp
      t->fast_expp = (double *) calloc(FEXP_L,sizeof(double));
      if (!t->fast_expp) Error("cannot allocate expp in table");
      size += FEXP_L*sizeof(double);
      t->fast_expm = (double *) calloc(FEXP_L,sizeof(double));
      if (!t->fast_expm) Error("cannot allocate expm in table");
      for (i=0;i<FEXP_L;i++)  *((t->fast_expm)+i) = exp( -(double)i / FEXP_BINRC );
      size += FEXP_L*sizeof(double);

      // Acos
      t->fast_acos = (double *) calloc(FACOS_L,sizeof(double));
      if (!t->fast_sin) Error("cannot allocate sin in table");
      for (i=0;i<FACOS_L;i++)  *((t->fast_acos)+i) = 180. * acos( (double) i / FACOS_BINRC - 1.) / PI;
      size += FACOS_L*sizeof(double);

      fprintf(fp,"Allocate tables of size %lf MB\n",(double)size/1024/1024);

      return t;
}


//ASTEMPERING
void FreePolymer(struct s_polymer *p,int npol, int nback, int nside, int shell, int noside){
 
int ires,iside,ipol;



for (ipol=0;ipol<npol;ipol++)
	{
		for (ires=0;ires<nback;ires++)
		{
		if (!noside)
			{
			for (iside=0;iside<nside;iside++)
			{
				free((((((p+ipol)->back)+ires)->side)+iside)->rot);

			}


		}

if (shell)
{
free((((p+ipol)->back)+ires)->shell);
free((((p+ipol)->back)+ires)->shell_p);
}


if (!noside)
free((((p+ipol)->back)+ires)->side);

free((((p+ipol)->back)+ires)->contacts);
free((((p+ipol)->back)+ires)->contacts_p);

free((((p+ipol)->back)+ires)->e) ;



}

free((p+ipol)->vback);
free((p+ipol)->back);

	}



free(p);





}


//ASTEMPERING
void FreeTables(struct s_tables *t){

free(t->fast_acos);
free(t->fast_expm);
free(t->fast_expp);
free(t->fast_cos);
free(t->fast_sin);
free(t->fast_sqrt);
free(t);

}



struct s_potential *AlloPotential(int natoms, int ntypes, int noangpot, int nodihpot, int hb)
{
	int i,k;
	struct s_potential *x;

	x = (struct s_potential *) calloc(1,sizeof(struct s_potential));
	if (!x) Error("Cannot allocate potential structure");

	// allocate pair potentials
	x->e = (double **) calloc(ntypes,sizeof(double *));
	if (!(x->e)) Error("Cannot allocate e in potential structure");
	x->r_2 = (double **) calloc(ntypes,sizeof(double *));
	if (!(x->r_2)) Error("Cannot allocate r in potential structure");
	x->r0_2 = (double **) calloc(ntypes,sizeof(double *));
	if (!(x->r0_2)) Error("Cannot allocate r0 in potential structure");

	for (i=0;i<ntypes;i++)
	{
		*((x->e)+i) = (double *) calloc(ntypes,sizeof(double));
		if (!(*((x->e)+i))) Error("Cannot allocate e in potential structure");
		*((x->r_2)+i) = (double *) calloc(ntypes,sizeof(double));
		if (!(*((x->r_2)+i))) Error("Cannot allocate r in potential structure");
		*((x->r0_2)+i) = (double *) calloc(ntypes,sizeof(double));
		if (!(*((x->r0_2)+i))) Error("Cannot allocate r0 in potential structure");

		for (k=0;k<ntypes;k++)
		{
			(x->e)[i][k] = 0;
			(x->r_2)[i][k] = 0;
			(x->r0_2)[i][k] = 0;

		}
	}

	// allocate angular and dihedrals
	if (!noangpot)
	{
		x->e_ang = (double *) calloc(natoms,sizeof(double));
		if (!(x->e_ang)) Error("Cannot allocate e_ang in potential structure");
		x->ang0 = (double *) calloc(natoms,sizeof(double));
		if (!(x->ang0)) Error("Cannot allocate ang0 in potential structure");
	}

	if (!nodihpot)
      {
            x->e_dih1 = (double *) calloc(natoms,sizeof(double));
            if (!(x->e_dih1)) Error("Cannot allocate e_dih1 in potential structure");
            x->dih01 = (double *) calloc(natoms,sizeof(double));
            if (!(x->dih01)) Error("Cannot allocate dih01 in potential structure");
            x->e_dih3 = (double *) calloc(natoms,sizeof(double));
            if (!(x->e_dih3)) Error("Cannot allocate e_dih3 in potential structure");
            x->dih03 = (double *) calloc(natoms,sizeof(double));
            if (!(x->dih03)) Error("Cannot allocate dih03 in potential structure");
            x->dih_pa = (double *) calloc(natoms,sizeof(double));
            if (!(x->dih_pa)) Error("Cannot allocate dihpha in potential structure");
            x->dih_pb = (double *) calloc(natoms,sizeof(double));
            if (!(x->dih_pb)) Error("Cannot allocate dihphb in potential structure");
            x->dih_which = (int *) calloc(natoms,sizeof(int));
            if (!(x->dih_which)) Error("Cannot allocate dig_which in potential structure");
            x->dih_f_psi_a =  (double *) calloc(NDIHFMAX,sizeof(double));
            if (!(x->dih_f_psi_a)) Error("Cannot allocate dih_f_psi_a in potential structure");
            x->dih_f_psi_b =  (double *) calloc(NDIHFMAX,sizeof(double));
            if (!(x->dih_f_psi_b)) Error("Cannot allocate dih_f_psi_b in potential structure");
            x->dih_f_phi_a =  (double *) calloc(NDIHFMAX,sizeof(double));
            if (!(x->dih_f_phi_a)) Error("Cannot allocate dih_f_phi_a in potential structure");
            x->dih_f_phi_b =  (double *) calloc(NDIHFMAX,sizeof(double));
            if (!(x->dih_f_phi_b)) Error("Cannot allocate dih_f_phi_b in potential structure");
            x->sigma = (double **) calloc(2,sizeof(double *));
            if (!(x->sigma)) Error("Cannot allocate sigma in potential structure");
            x->dih0 = (int **) calloc(2,sizeof(int *));
            if (!(x->dih0)) Error("Cannot allocate dih0 in potential structure");
            x->ab_propensity = (double **) calloc(2,sizeof(double *));
            if (!(x->ab_propensity)) Error("Cannot allocate ab_propensity in potential structure");

            for (i=0;i<2;i++)
            {
                  *((x->sigma)+i) = (double *) calloc(2,sizeof(double));
                  if (!(*((x->sigma)+i))) Error("Cannot allocate sigma in potential structure");
                  *((x->dih0)+i) = (int *) calloc(2,sizeof(int));
                  if (!(*((x->dih0)+i))) Error("Cannot allocate dih0 in potential structure");
                  *((x->ab_propensity)+i) = (double *) calloc(NAAMAX,sizeof(double));
                  if (!(*((x->sigma)+i))) Error("Cannot allocate ab-propensity in potential structure");
                  for (k=0;k<2;k++)
                  {
                        (x->sigma)[i][k] = 0;
                        (x->dih0)[i][k] = 0;
                  }
                  for (k=0;k<NAAMAX;k++)
                        (x->ab_propensity)[i][k] = LARGE;
            }
      }


	x->hc_type = (int *) calloc(ntypes,sizeof(int));
	if (!(x->hc_type)) Error("Cannot allocate hc_type in potential structure");
	x->hc_r0 = (double *) calloc(ntypes,sizeof(double));
	if (!(x->hc_r0)) Error("Cannot allocate hc_r0 in potential structure");
	x->hc_number=0;
	x->g_imin=0;
	x->boxtype = 'n';

	if (hb)
	{
		x->hb =  (int *) calloc(natoms,sizeof(int));
		if (!(x->hb)) Error("Cannot allocate hb in potential structure");
		x->hb_iam =  (int *) calloc(natoms,sizeof(int));
		if (!(x->hb_iam)) Error("Cannot allocate hb_iam in potential structure");
		for (i=0;i<natoms;i++) x->hb[i]=0;
	}

	return x;
}

//ASTEMPERING
void FreePotential(struct s_potential *x, int ntypes, int noangpot, int nodihpot, int hb){

int i;

if (hb)
         {
                free(x->hb);
free(x->hb_iam);
}
free(x->hc_type);
free(x->hc_r0);



if (!nodihpot)
      {
            free(x->e_dih1);
            free(x->dih01);
            free(x->e_dih3);
            free(x->dih03);
            free(x->dih_pa);
            free(x->dih_pb);
            free(x->dih_which);
            free(x->dih_f_psi_a);
            free(x->dih_f_psi_b);
            free(x->dih_f_phi_a);
            free(x->dih_f_phi_b);

            for (i=0;i<2;i++)
            {
                  free(*((x->sigma)+i));	
                  free(*((x->dih0)+i));		
                  free(*((x->ab_propensity)+i));
            }
            free(x->sigma);
            free(x->dih0);
            free(x->ab_propensity);
      }

if (!noangpot)
{
free(x->e_ang);
free(x->ang0);
}




	for (i=0;i<ntypes;i++)
{
free(*((x->e)+i));
free(*((x->r_2)+i));
free(*((x->r0_2)+i));

}

free(x->e);
free(x->r_2);
free(x->r0_2);

free(x);



}






/*****************************************************************************
 Allocate int matrix in the form x[l][m]
 *****************************************************************************/
int **AlloIntMatrix(int l, int m)
{
  int **x;
  int i,j;

  x = (int **) calloc(l,sizeof(int *));
  if (!x) Error("Cannot allocate int matrix");

  for (i=0;i<l;i++)
    {
      *(x+i) = calloc(m,sizeof(int));
      if (!(*(x+i))) Error("Cannot allocate int matrix");
      for (j=0;j<m;j++) *(*(x+i)+j) = 0;
    }

  return x;
}

/*****************************************************************************
 Allocate double matrix in the form x[l][m]
 *****************************************************************************/
double **AlloDoubleMatrix(int l, int m)
{
  double **x;
  int i,j;

  x = (double **) calloc(l,sizeof(double *));
  if (!x) Error("Cannot allocate double matrix");

  for (i=0;i<l;i++)
    {
      *(x+i) = calloc(m,sizeof(double));
      if (!(*(x+i))) Error("Cannot allocate double matrix");
      for (j=0;j<m;j++) *(*(x+i)+j) = 0.;
    }

  return x;
}

//ASTEMPERING
void FreeDoubleMatrix(double **x,int l){
int i;
for (i=0;i<l;i++) free(*(x+i));

free(x);

}



int *AlloInt(int n)
{
	int *x,i;
	
	x = (int *) calloc(n,sizeof(int));
	if (!x) Error("Cannot allocate int vector");

	for (i=0;i<n;i++) x[i]=0;

	return x;
}

double *AlloDouble(int n)
{
	double *x;
	int i;
	
	x = (double *) calloc(n,sizeof(double));
	if (!x) Error("Cannot allocate double vector");

	for (i=0;i<n;i++) x[i]=0.;

	return x;
}




