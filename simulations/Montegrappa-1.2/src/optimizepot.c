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
 * optimizepot.h
 *
 *  Created on: Sep 19, 2011
 *      Author: guido
 */

#include "montegrappa.h"


/***************************************************
 Allocate op structure and read constrain from file
 ***************************************************/
struct s_optimizepot *InitializeOptimizePot(struct s_mc_parms *parms, int ntypes, struct s_potential *u, FILE *fproc, int iproc, FILE *fp, int *nrestr)
{
	int i,j;
	struct s_optimizepot *x;

	// if it does need to optimize potential, allocate the structure and set record=0
	if (!strcmp(parms->op_minim,"none"))
	{
		x = (struct s_optimizepot *) calloc(1,sizeof(struct s_optimizepot));
		if (x==NULL) Error("Cannot allocate s_optimizepot");
		x->record= 0;
		return x;
	}

	fprintf(fproc,"\nActivating module OPTIMIZE POTENTIAL\n\nInitialize structures\n");
	
	// read restrains from file (and allocates input structure)
	if(iproc==0)
	{
		parms->op_input = ReadOPRestrains(parms);
		fprintf(stderr," ntypes=%d\n OP_NCONTMAX=%d\n nrestr=%d\n nframes=%llu\n",ntypes,OP_NCONTMAX,parms->op_input->ndata,parms->nstep/parms->op_deltat+1);
	
	*nrestr=parms->op_input->ndata;
	}	


	

	// allocate structure for recording conformations
	#ifdef ACTIVE_MPI
	MPI_Bcast(nrestr,1,MPI_INT,0,MPI_COMM_WORLD);
	if(iproc!=0){ 
	parms->op_input = Allo_op_input(*nrestr);
	}
	parms->op_input = send_op_input((*nrestr),parms->op_input);
	fprintf(fproc," ntypes=%d\n OP_NCONTMAX=%d\n nrestr=%d\n nframes=%llu\n",ntypes,OP_NCONTMAX,parms->op_input->ndata,parms->nstep/parms->op_deltat+1);
	
	#endif
	x = AlloOptimizePot(parms,ntypes,parms->op_input->ndata,fproc);

	#ifdef OP_DEBUG
	 fprintf(fp,"LOG of OPTIMIZEPOTENTIAL\n");
	#endif

	// set default hardcores and width
	if (parms->op_r>0)
		for (i=0;i<ntypes;i++)
			for (j=i+u->g_imin;j<ntypes;j++)
				if ( u->r_2[i][j] <0 )
	  			{
					 u->r_2[i][j] = parms->op_r * parms->op_r;
					 u->r_2[j][i] = parms->op_r * parms->op_r;
				  	 u->r0_2[i][j] = parms->op_r0 * parms->op_r0;
					 u->r0_2[j][i] = parms->op_r0 * parms->op_r0;
				}

	return x;
}

/***************************************************
 Allocate op structure
 NCONTMAX can be small; if reached, the array is reallocated
 ***************************************************/
struct s_optimizepot *AlloOptimizePot(struct s_mc_parms *parms, int ntypes, int nrestr, FILE *fproc)
{
	struct s_optimizepot *x;
	int nframes;

	nframes = (parms->nstep / parms->op_deltat) + 1;

	x = (struct s_optimizepot *) calloc(1,sizeof(struct s_optimizepot));
	if (x==NULL) Error("Cannot allocate s_optimizepot");

	x->it1 = AlloInt(OP_NCONTMAX);
	x->it2 = AlloInt(OP_NCONTMAX);
	x->mul = AlloDoubleMatrix(nframes,OP_NCONTMAX);
	x->eold = AlloDouble(nframes);
	x->enew = AlloDouble(nframes);	
	x->efix = AlloDouble(nframes);
	x->t = AlloDouble(nframes);
	x->restrain = AlloDoubleMatrix(nframes,nrestr);
	x->rest_av = AlloDouble(nrestr);

	fprintf(fproc,"Allocate OP structures of size %lf MB\n",((double)ntypes*ntypes*sizeof(double)+2*OP_NCONTMAX*sizeof(int)+nframes*OP_NCONTMAX*sizeof(double)+
					3*nframes*sizeof(double)+nframes*nrestr*sizeof(double)+nrestr*sizeof(double))/1024/1024);
	
	#ifdef ACTIVE_MPI
	fprintf(stderr,"Allocate OP structures of size %lf MB\n",((double)ntypes*ntypes*sizeof(double)+2*OP_NCONTMAX*sizeof(int)+nframes*OP_NCONTMAX*sizeof(double)+
															 3*nframes*sizeof(double)+nframes*nrestr*sizeof(double)+nrestr*sizeof(double))/1024/1024);
	#endif
	
	x->nframesmax = nframes;
	x->nframes = 0;
	x->icount = 0;
	x->nallocont = 1;		// i.e., how many multiples of OP_NCONTMAX
	x->record= 0;
	x->ncontacts = 0;

	return x;
}


void FreeOpt(struct s_optimizepot *x,struct s_mc_parms *parms,int nrestr){

    int nframes;
 
    nframes = (parms->nstep / parms->op_deltat) + 1;


    free(x->it1);
    free(x->it2);
    
        free(x->eold);
    free(x->enew);
    free(x->efix);
    free(x->t);
    free(x->rest_av);
	
	FreeDoubleMatrix((x->mul),nframes);
	FreeDoubleMatrix((x->restrain),nframes);


free(x);




}

/***************************************************
 Read list of restrains from file.
 Format:  i1  i2  i3  i4  type  value  sigma 
 ***************************************************/

struct s_optimizepot_input *ReadOPRestrains(struct s_mc_parms *parms)
{
	struct s_optimizepot_input *x;
	int i;
	char aux[500];
	FILE *fp;
	
	// allocate input structure
	x = (struct s_optimizepot_input *) calloc(1,sizeof(struct s_optimizepot_input));
	if (x==NULL) Error("Cannot allocate s_optimizepot_input");

	// open restrain file
	fp = fopen(parms->fnop,"r");
	if (!fp) Error("Cannot open restrain file");	
	
	//read first lines
	fgets(aux,500,fp);
	if (sscanf(aux,"ndata %d",&(x->ndata)) != 1) Error("First line of restrain data should be: ndata <int>");

	// allocate arrays
	x->expdata = AlloDouble(x->ndata);
	x->sigma = AlloDouble(x->ndata);
	x->datatype = AlloInt(x->ndata);
	x->i1 = AlloInt(x->ndata);
	x->i2 = AlloInt(x->ndata);
	x->i3 = AlloInt(x->ndata);
	x->i4 = AlloInt(x->ndata);

	// read file
	for (i=0;i<x->ndata;i++)
	{
		fgets(aux,500,fp);

		if ( sscanf(aux,"%d %d %d %d %d %lf %lf",&(x->i1[i]),&(x->i2[i]),&(x->i3[i]),&(x->i4[i]),&(x->datatype[i]),&(x->expdata[i]),&(x->sigma[i])) != 7 )
		{
			if ( sscanf(aux,"%d %d %d %lf %lf",&(x->i1[i]),&(x->i2[i]),&(x->datatype[i]),&(x->expdata[i]),&(x->sigma[i])) != 5 )
			{
				fprintf(stderr,"line: %d string: %s\n",i,aux); 
				Error("Error in reading restrain file"); 
			}
		}
		else if ( x->datatype[i] != 3) Error("If it reads i3 and i4 in restrain file, datatype must be 3");
	}

	fprintf(stderr," %d restrains read\n",x->ndata);
	fclose(fp);

	return x;
}

void FreeOptInput(struct s_optimizepot_input *x)
{
	free(x->expdata);
	free(x->sigma);
	free(x->datatype);
	free(x->i1);
	free(x->i2);
	free(x->i3);
	free(x->i4);
	free(x);

}
/***************************************************
 Get the actual values of the resrtain from chain ipol
 of polymer p
 ***************************************************/
void OP_GetRestrain(int it, struct s_polymer *p, int ipol, struct s_optimizepot_input *opi)
{
	int i,k,l,m,i1,i2,i3,i4;
	double xraw,d;

	// loop on all restrains 
	for (i=0;i<opi->ndata;i++)
	{
		
		i1 = opi->i1[i];
		i2 = opi->i2[i];
		xraw = 0.;
		
		// contact between amino acids
		if ( opi->datatype[i] == 0 )
		{
			for (k=0;k<(((p+ipol)->back)+i1)->ncontacts;k++)
				if ( *(((((p+ipol)->back)+i1)->contacts)+k) == i2 ) xraw = 1.;
		}
		// distance between atoms
		else if ( opi->datatype[i] == 1 )
		{
				
		xraw = FastSqrt( Dist2( *(*(((p+ipol)->vback)+i1)) , *(*(((p+ipol)->vback)+i2)) ), p->tables );
					
		}
		// 1/r^6 between atoms
	 	else if ( opi->datatype[i] == 2 )
		{
			d = Dist2( *(*(((p+ipol)->vback)+i1)) , *(*(((p+ipol)->vback)+i2)) );
			xraw = 1. / (d*d*d);
		}
		// contact between blocks of amino acids
		else if ( opi->datatype[i] == 3 )
		{
			i3 = opi->i3[i];		
			i4 = opi->i4[i];		
			for (k=i1;k<=i2;k++)
				for (l=i3;l<=i4;l++)
					for (m=0;m<(((p+ipol)->back)+k)->ncontacts;m++)
						if ( *(((((p+ipol)->back)+k)->contacts)+m) == l ) xraw = 1.;
		}
		
		p->op->restrain[it][i] = xraw;
	}

	
	#ifdef OP_DEBUG
	 FILE *fp;
	 fp = fopen("op_log","a");
	 for (i=0;i<opi->ndata;i++)
	 	fprintf(fp,"restrains: snap=%5d restr=%4d  value=%lf  T\n",it,i,p->op->restrain[it][i]);
	 fclose(fp);
	#endif

}

/***************************************************
 Add energy to energy table 
 ***************************************************/
void OP_AddEnergy(struct s_polymer *p, int a1, int a2, double mul)
{
	int iframe,icon,i;

	iframe = p->op->nframes;
	icon = p->op->ncontacts;
                                                                          
	// if there is already a contact between types a1 and a2, sum mul
	for (i=0;i<icon;i++)
		if ( (a1 == p->op->it1[i] && a2 == p->op->it2[i]) ||
			(a1 == p->op->it2[i] && a2 == p->op->it1[i]) ) 
		{
			p->op->mul[iframe][i] += mul;
			return;
		}

	// otherwise add to the list
	p->op->it1[icon] = a1;
	p->op->it2[icon] = a2;
	p->op->mul[iframe][icon] = mul;
	p->op->ncontacts ++;


	// if number of contacts exceed allocated size, reallocate
	#ifdef REALLOCATE
	if ( p->op->ncontacts >= OP_NCONTMAX * p->op->nallocont )
	{
		p->op->nallocont ++;	

		fprintf(stderr,"Reallocate OP_contact structure of size %ld MB (for %d contacts)\n",(2*OP_NCONTMAX * p->op->nallocont * sizeof(int) +
				OP_NCONTMAX * p->op->nallocont * p->op->nframesmax*sizeof(double))/1024/1024,OP_NCONTMAX * p->op->nallocont); 
		p->op->it1 = (int *) realloc(p->op->it1,( OP_NCONTMAX * p->op->nallocont * sizeof(int)));
		if (!p->op->it1) Error("Cannot reallocate it1");
		p->op->it2 = (int *) realloc(p->op->it2,( OP_NCONTMAX * p->op->nallocont * sizeof(int)));
		if (!p->op->it2) Error("Cannot reallocate it2");
		for (iframe=0;iframe<p->op->nframesmax;iframe++)
		{
			p->op->mul[iframe] = (double *) realloc( (p->op->mul[iframe]) , ( OP_NCONTMAX * p->op->nallocont * sizeof(double) ) );
			if (!p->op->mul[iframe]) Error("Cannot reallocate mul");
		}
		for (iframe=0;iframe<p->op->nframesmax;iframe++)
			for (i=p->op->ncontacts;i<OP_NCONTMAX * p->op->nallocont;i++) p->op->mul[iframe][i] = 0;
	}
	#else
	if ( p->op->ncontacts >= OP_NCONTMAX) Error("OP_NCONTMAX too small");
	#endif
			
}

/***************************************************
 Random search in energy space
 ***************************************************/
void OP_SamplePotential(struct s_polymer *p, struct s_mc_parms *parms, int ntypes, double **u, int irun, double **ematold)
{
	int istep,iacc=0,a1,a2,iw,ifr,i,j,n=0,nstep=0;
	double chi2,chi2old,deltae,eav=0,eav2=0,de2=0;
	char aux[500];
	FILE *fp, *fc;
	

	if ( p->op->nframes == 0 ) Error("Cannot optimize potential: no frame recorded for OP_SamplePotential");

	//keep old matrix for comparison
	CopyDoubleMatrix(u,ematold,ntypes,ntypes);


	// initial chi2
	chi2old = OP_function(u,p->op,parms->op_input,parms);
	
//	fprintf(stderr,"Potential Optimization \n> Initial value of Chi2 = %lf | Reduced = %lf\n",chi2old,chi2old/(parms->op_input->ndata));
	fprintf(stderr,"Potential Optimization \n> Initial value of Chi2 %lf\n",chi2old/(parms->op_input->ndata));

	fprintf(stderr,"> History\t\t\tChi2");
	//print chi2
	if(parms->nrun>1)
	{
		sprintf(aux,"chi2.dat");
		fc = fopen(aux,"a");
		fprintf(fc,"%d\t%lf\n",irun, chi2old/(parms->op_input->ndata));
		fclose(fc);
	}
	
	// minimize
	for (istep=0;istep<parms->op_itermax;istep++)
	{
		iw = irand(p->op->ncontacts);
		deltae = (frand()-0.5) * parms->op_step;

		do { deltae = (frand()-0.5) * parms->op_step; }					// generate a random deltae
		while ( u[p->op->it1[iw]][p->op->it2[iw]]+deltae < parms->op_emin ||
			 u[p->op->it1[iw]][p->op->it2[iw]]+deltae > parms->op_emax );		// check if matrix element within boundaries

		chi2 = OP_functionDiff(iw,deltae,p->op,parms->op_input,parms);

		// accept
		if (chi2 < chi2old)
		{
			a1 = p->op->it1[iw];
			a2 = p->op->it2[iw];
			u[a1][a2] += deltae;
			u[a2][a1] += deltae;
			chi2old = chi2;
			iacc ++;
		}
		// reject
		else
		{	
			for (ifr=0;ifr<p->op->nframes;ifr++)			// return to previous value
				p->op->enew[ifr] -= deltae * p->op->mul[ifr][iw];	
		}

		nstep ++;

		// print something
	//	if (!(istep%parms->op_print)) fprintf(stderr,"   %5d\t\t%lf\t%lf\n",istep,chi2old, chi2old/(parms->op_input->ndata));
		if (!(istep%parms->op_print)) fprintf(stderr,"   %5d\t\t%lf\n",istep,chi2old/(parms->op_input->ndata));
		// exit condition
		if (chi2<parms->op_stop) break;
	}

	// calculate the properties of the new matrix
	for (i=0;i<ntypes;i++)
		for (j=i+1;j<ntypes;j++) 
		{ 
			eav += u[i][j]; 
			eav2 += u[i][j]*u[i][j]; 
			de2 += (u[i][j]-ematold[i][j])*(u[i][j]-ematold[i][j]);
			n++; 
		}


	fprintf(stderr,"> Result of the sampling: Chi2=%lf\n",chi2old/parms->op_input->ndata);
	fprintf(stderr,"> accepted changes = %lf\n",(double)iacc/nstep);
	fprintf(stderr,"> average E = %lf +/- %lf\n",eav/n,sqrt(eav2/n-eav*eav/n/n));
	fprintf(stderr,"> rms difference = %lf\n\n",sqrt(de2/n));
	

	// print to file
	sprintf(aux,"restraints_%d.dat",irun);
	fp = fopen(aux,"w");
	for (i=0;i<parms->op_input->ndata;i++)
		fprintf(fp,"%d\t%lf\t%lf\t%lf\n",i,p->op->rest_av[i],parms->op_input->expdata[i],parms->op_input->sigma[i]);
	fclose(fp);

	// reset structures
	for (i=0;i<p->op->ncontacts;i++)
	{
		p->op->it1[i] = 0;
		p->op->it2[i] = 0;
		for (ifr=0;ifr<p->op->nframes;ifr++)  p->op->mul[ifr][i] = 0;
	}
	for (ifr=0;ifr<p->op->nframes;ifr++)
	{
		p->op->eold[ifr] = 0;
		p->op->efix[ifr] = 0;
		p->op->enew[ifr] = 0;
		p->op->t[ifr] = -1;
		for (i=0;i<parms->op_input->ndata;i++)  p->op->restrain[ifr][i] = 0;
	}
	for (i=0;i<parms->op_input->ndata;i++) p->op->rest_av[i] = 0;
	p->op->ncontacts = 0;
	p->op->nframes = 0;
	p->op->icount = 0;

	// reset structures						
	for (i=0;i<p->op->ncontacts;i++)
	{
		p->op->it1[i] = 0;
		p->op->it2[i] = 0;
		for (ifr=0;ifr<p->op->nframes;ifr++)  p->op->mul[ifr][i] = 0;
	}
	for (ifr=0;ifr<p->op->nframes;ifr++)
	{
		p->op->eold[ifr] = 0;
		p->op->efix[ifr] = 0;
		p->op->enew[ifr] = 0;
		p->op->t[ifr] = -1;
		for (i=0;i<parms->op_input->ndata;i++)  p->op->restrain[ifr][i] = 0;
	}
	for (i=0;i<parms->op_input->ndata;i++) p->op->rest_av[i] = 0;
	p->op->ncontacts = 0;
	p->op->nframes = 0;
	p->op->icount = 0;	
	
	
}


/***************************************************
 Calculate the chi2 with experimental data, updating
 the vector rest_av and e_new
 ***************************************************/
double OP_function(double **e, struct s_optimizepot *x, struct s_optimizepot_input *in, struct s_mc_parms *parms)
{
	double chi2,enew,eold,mul,emax=-9E19,z=0.,boltz;
	int ifr,icont,it1,it2,ir;

	for (ir=0;ir<in->ndata;ir++) x->rest_av[ir]=0.;
                             
	// calculate energies
	for (ifr=0;ifr<x->nframes;ifr++)
	{
		x->enew[ifr] = x->efix[ifr];			// one-body energy
		for (icont=0;icont<x->ncontacts;icont++)	// calculate new two-body energy
			if (x->mul[ifr][icont]>0)		// if in frame ifr there is the icont contact
			{
				it1 = x->it1[icont];
				it2 = x->it2[icont];
				mul = x->mul[ifr][icont];
				x->enew[ifr] += e[it1][it2] * mul;
			}
                                              
		// calculate emax to avoid overflow in the exp
		if (-x->enew[ifr]/parms->op_T + x->eold[ifr]/x->t[ifr] > emax) 
						emax = -x->enew[ifr]/parms->op_T + x->eold[ifr]/x->t[ifr];
	}
                           
	// calculate thermal averages
	for (ifr=0;ifr<x->nframes;ifr++)
	{
		eold = x->eold[ifr];						// total energy with the original matrix
		enew = x->enew[ifr];
		boltz = exp(-enew/parms->op_T + eold/x->t[ifr] - emax);

		// make the average
		for (ir=0;ir<in->ndata;ir++) 
			x->rest_av[ir] +=  x->restrain[ifr][ir] * boltz;
                                                                                                                            
		z += boltz;
	}


	for (ir=0;ir<in->ndata;ir++) x->rest_av[ir] /= z;
	
	chi2 = Chi2(x->rest_av,in->expdata,in->sigma,in->ndata);

	return chi2;
}

/***************************************************
 Calculate chi2 between two vectors
 ***************************************************/
double Chi2(double *x, double *xexp, double *sigma, int n)
{
	int i;
	double chi2=0.;

	for (i=0;i<n;i++)
		chi2 += (x[i] - xexp[i])*(x[i] - xexp[i]) / ( sigma[i] * sigma[i] );

	return chi2;
}

/***************************************************
 Calculate the chi2 with experimental data, given
 the energy difference in the contact iw
 it needs x->enew and changes x->enew[ifr]
 ***************************************************/
double OP_functionDiff(int iw, double deltae, struct s_optimizepot *x, struct s_optimizepot_input *in, struct s_mc_parms *parms)
{
	double chi2,enew,eold,emax=-9E19,z=0.,boltz;
	int ifr,ir;

	for (ir=0;ir<in->ndata;ir++) x->rest_av[ir]=0.;


	// calculate energies
	for (ifr=0;ifr<x->nframes;ifr++)
	{
		// if the iw contact is present in the ifr frame, then its energy is affected
		if ( x->mul[ifr][iw] > 0 )
			x->enew[ifr] += deltae * x->mul[ifr][iw];

		// calculate emax to avoid overflow in the exp
		if (-x->enew[ifr]/parms->op_T + x->eold[ifr]/x->t[ifr] > emax) 
						emax = -x->enew[ifr]/parms->op_T + x->eold[ifr]/x->t[ifr];
	}

	// calculate thermal averages
	for (ifr=0;ifr<x->nframes;ifr++)
	{
		eold = x->eold[ifr];						// total energy with the original matrix
		enew = x->enew[ifr];
		boltz = exp(-enew/parms->op_T + eold/x->t[ifr] - emax);

		// make the average
		for (ir=0;ir<in->ndata;ir++)
			x->rest_av[ir] += x->restrain[ifr][ir] * boltz;
		z += boltz;
	}

	for (ir=0;ir<in->ndata;ir++) x->rest_av[ir] /= z;
	
	chi2 = Chi2(x->rest_av,in->expdata,in->sigma,in->ndata);

	#ifdef OP_DEBUG
	 FILE *fp;
	 fp = fopen("op_log","a");
	 fprintf(fp,"minim:  iw=%d deltae=%lf\n",iw,deltae);
 	 for (ifr=0;ifr<x->nframes;ifr++)
	 	fprintf(fp,"\teold[%d]=%lf enew[%d]=%lf efix[%d]=%lf mul[%d]=%lf",ifr,x->eold[ifr],ifr,x->enew[ifr],ifr,x->efix[ifr],ifr,x->mul[ifr][iw]);
         fprintf(fp,"\nemax=%lf\n",emax);
         for (ir=0;ir<in->ndata;ir++)
		fprintf(fp,"\t <x>[%d]=%lf\t",ir,x->rest_av[ir]);
         fprintf(fp,"\n");
         fprintf(fp,"\n");
	 fclose(fp);
	#endif

	return chi2;
}

