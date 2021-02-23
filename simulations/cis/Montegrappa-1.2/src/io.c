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
 * io.c
 *
 *  Created on: Sep 15, 2010
 *      Author: guido
 */
#include "montegrappa.h"
	
/*****************************************************************************
 Standard error function which stops the execution
 *****************************************************************************/
void Error(char *text)
{
  time_t rawtime;

  time(&rawtime);
  fprintf(stderr,"---------------------------------------------\n");
  fprintf(stderr,"Error:\n");
  fprintf(stderr,"* %s\n",text);
  fprintf(stdout,"after %3.1f s. at %s\n",(double) ( clock() / CLOCKS_PER_SEC ),
                                    asctime(localtime(&rawtime)));
  fprintf(stderr,"---------------------------------------------\n\n");

  exit(1);
}

/********************************************************************
 Read polymer structure from file
 	 returns the number of atoms read
 ********************************************************************/
void ReadPolymer(char *fname, struct s_polymer *p, FILE *flog, int npol, int debug, int *iamax, int *itypemax)
{
	int iback=0,iab,irot,iar,b1,b2,b3,iat,iaa,ia,ch,ichmax=-1,move,rr=0,itype;
	int readback=0,readside=0,readrot=0,nat_back=0;//wback[NATOMMAX];//cback[NATOMMAX];
	double ang,dih,r,x,y,z;
	char aux[500],keyword[100],at[5],aa[5],type[5];
	FILE *fp;


	fprintf(stderr,"read polymer from %s\n",fname);
	#ifdef DEBUG
	 fprintf(flog,"read polymer from %s\n",fname); fflush(flog);
	#endif

	fp = fopen(fname,"r");
	if (!fp) Error("Cannot open polymer file");

	*iamax = -1;
	*itypemax = -1;
	for (ch=0;ch<npol;ch++) (p+ch)->nback = 0;

	while(fgets(aux,500,fp)!=NULL)
    	{
	     if ( FindKeyword(aux,keyword)==1 && !strcmp(keyword,"backbone") ) { readback=1; readside=0; readrot=0;}
	     if ( FindKeyword(aux,keyword)==1 && !strcmp(keyword,"sidechains") ) { readback=0; readside=1; readrot=0;}
	     if ( FindKeyword(aux,keyword)==1 && !strcmp(keyword,"rotamers") ) { readback=0; readside=0; readrot=1;}

	     // Read backbone information
	     if ( readback==1 )
	     {
	    	 if ( sscanf( aux, "%d %d %s %d %s %d %d %lf %lf %lf %d", &iback, &ia,type,&itype,aa,&iaa,&ch,&x,&y,&z,&move) == 11 )
	    	 {
	    		 if (debug>3) fprintf(flog,"iback=%d ia=%3d type=%3s itype=%3d aa=%3s iaa=%3d ch=%2d x=%lf y=%lf z=%lf move=%d\n",iback,ia,type,itype,aa,iaa,ch,x,y,z,move);

	    		 if (iback >= NRESMAX) Error("NRESMAX too small");
	    		 (((p+ch)->back)+iback)->ia = ia;
	    		 strcpy( (((p+ch)->back)+iback)->type,type );
	    		 (((p+ch)->back)+iback)->itype = itype;
	    		 strcpy( (((p+ch)->back)+iback)->aa,aa );
	    		 (((p+ch)->back)+iback)->iaa = iaa;
	    		 ((((p+ch)->back)+iback)->pos).x = x;
	    		 ((((p+ch)->back)+iback)->pos).y = y;
	    		 ((((p+ch)->back)+iback)->pos).z = z;
	    		 (((p+ch)->back)+iback)->move = move;

//			fprintf(stderr, "%d\t %d\t%s\t%d\t%s\t%d\t%d\t%lf\t%lf\t%lf\t%d\n",iback,(((p+ch)->back)+iback)->ia,(((p+ch)->back)+iback)->type,(((p+ch)->back)+iback)->itype,(((p+ch)->back)+iback)->aa,
  //                                      (((p+ch)->back)+iback)->iaa,ch,((((p+ch)->back)+iback)->pos).x,((((p+ch)->back)+iback)->pos).y,((((p+ch)->back)+iback)->pos).z, (((p+ch)->back)+iback)->move );


	    		 //wback[ (((p+ch)->back)+iback)->ia ] = iback;				// record backbone id in the lookback table
	    		 //cback[ (((p+ch)->back)+iback)->ia ] = ch;					// record chain in the lookback table
	    		 if (ch > ichmax) ichmax = ch;								// maximum chain-id read
	    		 if (ichmax+1 > npol) Error("NPOLMAX too small in reading polymers from file");

	    		 ((p+ch)->vback)[ (((p+ch)->back)+iback)->ia ] = &(((((p+ch)->back)+iback)->pos));		// record address of cartesian coordinates
	    		 	 	 	 	 	 	 	 	 	 	 	 	 	 	 	 	 	 	 	 	 	 	// of each atom
	    		 if ( ((p->back)+iback)->ia > *iamax ) *iamax = ((p->back)+iback)->ia;					// maximum atom-id overall
	    		 if ( ((p->back)+iback)->itype > *itypemax ) *itypemax = ((p->back)+iback)->itype;		// maximum atom-type overall

	    		 (p+ch)->nback =iback+1;
	    		 nat_back ++;
	    	 }
	     }

	     // Read rotamer information
	     if ( readrot==1 )
	     {
	    	 if ( sscanf(aux, "%d %d %d %d %d %d %d %d %s %d %lf %lf %lf",&iab,&ch,&irot,&iar,&b1,&b2,&b3,&iat,at,&itype,&ang,&dih,&r) == 13 )
	    	 {
	    		 if (debug>3) fprintf(flog,"iback=%3d ch=%d irot=%2d iside=%2d b1=%3d b2=%3d b3=%3d ia=%3d type=%3s itype=%3d ang=%lf %lf %lf (nback=%d)\n",iab,ch, irot,iar,b1,b2,b3,iat,at,itype,ang,dih,r,(p+ch)->nback);
	    		 if (ch>ichmax+1) Error("Reading sidechain of non-existing polymer");
	    		 if (iab>(p+ch)->nback) Error("Reading sidechain of non-existing backbone");
	    		 (((((p+ch)->back)+iab)->side)+iar)->ia = iat;
	    		 (((((p+ch)->back)+iab)->side)+iar)->itype = itype;
	    		 if ( itype > *itypemax ) *itypemax = itype;
	    		 if (iat > *iamax ) *iamax = iat;
	    		 strcpy( (((((p+ch)->back)+iab)->side)+iar)->type , at );
	    		 (((((((p+ch)->back)+iab)->side)+iar)->rot)+irot)->b1 = b1;
	    		 (((((((p+ch)->back)+iab)->side)+iar)->rot)+irot)->b2 = b2;
	    		 (((((((p+ch)->back)+iab)->side)+iar)->rot)+irot)->b3 = b3;
	    		 ((((((((p+ch)->back)+iab)->side)+iar)->rot)+irot)->ang).ang = ang;
	    		 ((((((((p+ch)->back)+iab)->side)+iar)->rot)+irot)->ang).dih = dih;
	    		 ((((((((p+ch)->back)+iab)->side)+iar)->rot)+irot)->ang).r = r;
	    		 if ( irot+1 > (((p+ch)->back)+iab)->nrot ) (((p+ch)->back)+iab)->nrot = irot+1; 		// set # of rotamers
	    		 if ( iar+1 > (((p+ch)->back)+iab)->nside ) (((p+ch)->back)+iab)->nside = iar+1;		// set number of atoms in the sidechain
	    		 rr++;
	    	 }
	     }

	     // Read sidechain information based on atom number (by default rotamer is in rotameric state 0)
	     if ( readside == 1 )
	     {
	    	 if ( sscanf(aux, "%d %d %d",&iab,&ch,&irot) == 3)
	    	 {
	    		 if (debug>3) fprintf(flog,"iback=%d ch=%d irot=%d\n",iab,ch,irot);
	    		 if (iab>(p+ch)->nback) Error("Reading rotameric state for non-existing atom");
	    		 (((p+ch)->back)+iab)->irot = irot;
	    	 }
	     }
    }

	// check number of chains
	if (ichmax+1 != npol) Error("Number of chains contained in polymer file is different from that indicated in parameter file");

	// assign addresses of sidechain position to lookback table
	for (ch=0;ch<=ichmax;ch++)
		for (iab=0;iab<(p+ch)->nback;iab++)
			for (iat=0;iat<(((p+ch)->back)+iab)->nside;iat++)
				((p+ch)->vback)[ (((((p+ch)->back)+iab)->side)+iat)->ia ] = &((((((p+ch)->back)+iab)->side)+iat)->pos);

	(*iamax) ++;		// number of atoms and of atomtypes
	(*itypemax)++;

	for (ch=0;ch<ichmax;ch++)

	if (p->nback==0) Error("No backbone atoms read");
	fprintf(flog,"Read %d backbone atoms in %s\n",nat_back,fname);
	fprintf(flog,"Read %d rotamers for sidechains\n",rr);
	fprintf(flog,"There are %d disjoint chains in the system.\n",ichmax+1);
	fprintf(flog,"There are %d atoms\n",*iamax);
	fprintf(flog,"There are %d atom types\n",*itypemax);

	return;
}

void SetLookbackTables(struct s_polymer *p, int nc)
{
	int iab,iat;

	// assign addresses of backbone position
	for (iab=0;iab<(p+nc)->nback;iab++)
	{
		// given the atom, find the pointer to the position vector
		((p+nc)->vback)[ (((p+nc)->back)+iab)->ia ] = &(((((p+nc)->back)+iab)->pos));
		for (iat=0;iat<(((p+nc)->back)+iab)->nside;iat++)
			((p+nc)->vback)[ (((((p+nc)->back)+iab)->side)+iat)->ia ] = &((((((p+nc)->back)+iab)->side)+iat)->pos);
	}

}

/********************************************************************
 Given a string, returns the content of square brackets, excluding
 spaces and returns 1,
 if there are no square brackets returns 0
 ********************************************************************/
int FindKeyword(char *string, char *keyword)
{
   int i=0,key=0,k=0;
   char c;

   do {
		   c = string[i];
		   if (c=='[') key=1;
		   if (c==']') key=0;
		   if (key==1 && c!=' ' && c!='[' && c!='\t')
		   {
			   keyword[k] = c;
			   k++;
		   }
		   i++;

	   } while ( c != '\0' && i<500 );
	keyword[k] = '\0';

	if (!strcmp(keyword,"")) return 0;
	else return 1;
}

/********************************************************************
 Print polymer structure to file (or stdout or stderr)
 ********************************************************************/
void PrintPolymer(char *fname, struct s_polymer *p, int nchains)
{
	int i,cl=0,j,k,ic;
	FILE *fout;
	
	// decide where to write output
	if (!strcmp(fname,"stdout")) fout = stdout;
	else if (!strcmp(fname,"stderr")) fout = stderr;
	else
	{
		fout = fopen(fname,"w");
		if (!fout) Error("Cannot open polymer file for reading");
		cl=1;
	}

	// write it
	fprintf(fout,"[ backbone ]\n");
	fprintf(fout,"back\tia\ttype\titype\taa\tiaa\tch\tx\t\ty\t\tz\t\ttomove\n");
	
	for (ic=0;ic<nchains;ic++)
		for (i=0;i<(p+ic)->nback;i++)
		{
			fprintf(fout, "%d\t %d\t%s\t%d\t%s\t%d\t%d\t%lf\t%lf\t%lf\t%d\n",i,(((p+ic)->back)+i)->ia,(((p+ic)->back)+i)->type,(((p+ic)->back)+i)->itype,(((p+ic)->back)+i)->aa,
					(((p+ic)->back)+i)->iaa,ic,((((p+ic)->back)+i)->pos).x,((((p+ic)->back)+i)->pos).y,((((p+ic)->back)+i)->pos).z, (((p+ic)->back)+i)->move );
//			fprintf(stderr, "%d\t %d\t%s\t%d\t%s\t%d\t%d\t%lf\t%lf\t%lf\t%d\n",i,(((p+ic)->back)+i)->ia,(((p+ic)->back)+i)->type,(((p+ic)->back)+i)->itype,(((p+ic)->back)+i)->aa,
  //                                      (((p+ic)->back)+i)->iaa,ic,((((p+ic)->back)+i)->pos).x,((((p+ic)->back)+i)->pos).y,((((p+ic)->back)+i)->pos).z, (((p+ic)->back)+i)->move );

//			#ifdef	DEBUG_SHELL
			//fprintf(stderr,"ic = %d  \t i = %d \t nshell = %d\t",ic,i,(((p+ic)->back)+i)->nshell);	
//			ntshell+=(((p+ic)->back)+i)->nshell;
//			#endif
		}

//	fprintf(stderr,"NTSHELL = %d \n",ntshell);

	fprintf(fout,"\n[ rotamers ]\n");
	fprintf(fout,"back ch rot at\t b1  b2   b3\tia type itype\tang\t    dih\t      r\n");
	for (ic=0;ic<nchains;ic++)
		for (i=0;i<(p+ic)->nback;i++)
			for (j=0;j<(((p+ic)->back)+i)->nside;j++)
				for (k=0;k<(((p+ic)->back)+i)->nrot;k++)
					fprintf(fout, "%3d %3d %2d %2d\t%3d %3d %3d\t%3d %3s %3d\t%lf %lf %lf\n",i,ic,k,j,(((((((p+ic)->back)+i)->side)+j)->rot)+k)->b1,
							(((((((p+ic)->back)+i)->side)+j)->rot)+k)->b2,(((((((p+ic)->back)+i)->side)+j)->rot)+k)->b3,
							(((((p+ic)->back)+i)->side)+j)->ia,(((((p+ic)->back)+i)->side)+j)->type,(((((p+ic)->back)+i)->side)+j)->itype,
							((((((((p+ic)->back)+i)->side)+j)->rot)+k)->ang).ang,((((((((p+ic)->back)+i)->side)+j)->rot)+k)->ang).dih,
							((((((((p+ic)->back)+i)->side)+j)->rot)+k)->ang).r );

	fprintf(fout,"\n[ sidechains ]\n");
	fprintf(fout,"back\tch\tirot\n");
	for (ic=0;ic<nchains;ic++)
		for (i=0;i<(p+ic)->nback;i++)
			if ((((p+ic)->back)+i)->nrot > 0)
				fprintf(fout,"%d\t%d\t%d\n",i,ic,(((p+ic)->back)+i)->irot);


	if (cl==1) fclose(fout);
}

/********************************************************************
 Print pdb to stream
 ********************************************************************/
void PrintPDBStream(struct s_polymer *p, int npol, FILE *fp)
{
	int i,j,ch;
	char cc;

	fprintf(fp,"TITLE %s\n",p->title);

	for (ch=0;ch<npol;ch++)
	{
		for (i=0;i<(p+ch)->nback;i++)
	      {
			  if (npol<24) cc = ch+65;
			  else cc=' ';

	          fprintf(fp,"ATOM  %5d %-4s%3s %c%4d    %8.3lf%8.3lf%8.3lf\n",(((p+ch)->back)+i)->ia,(((p+ch)->back)+i)->type,
	        		  (((p+ch)->back)+i)->aa,cc,(((p+ch)->back)+i)->iaa,((((p+ch)->back)+i)->pos).x,((((p+ch)->back)+i)->pos).y,
	        		  ((((p+ch)->back)+i)->pos).z );

	          for (j=0;j<(((p+ch)->back)+i)->nside;j++)
	        	  fprintf(fp,"ATOM  %5d %-4s%3s %c%4d    %8.3lf%8.3lf%8.3lf\n", (((((p+ch)->back)+i)->side)+j)->ia,
	        			  (((((p+ch)->back)+i)->side)+j)->type, (((p+ch)->back)+i)->aa,cc,(((p+ch)->back)+i)->iaa,
	        			  ((((((p+ch)->back)+i)->side)+j)->pos).x,((((((p+ch)->back)+i)->side)+j)->pos).y,
	        			  ((((((p+ch)->back)+i)->side)+j)->pos).z);
	      }
        if (npol>=24) fprintf(fp,"ENDMDL\n");
	}
	fprintf(fp,"ENDMDL\n");

	#ifdef DEBUG
	 fflush(fp);
	#endif
}

/********************************************************************
 Allocate structure and read MC parameters
 ********************************************************************/
struct s_mc_parms *ReadMcParms(char *fname)
{
	int i,chk;
	struct s_mc_parms *x;
	char aux[1000],nlog[1000];
	FILE *fp;
	
	x = calloc(1,sizeof(struct s_mc_parms));
	if (!x) Error("Cannot allocate mc_parms");

	#ifdef STEMPERING
	int r;
	double input_temp[NTEMPMAX],input_g[NTEMPMAX],dumb;
	char aux_gm[1000];
	x->p = (struct st_stuff *) calloc(1,sizeof(struct st_stuff));

	strcpy((x->p)->st_nftout,"temperatures.dat");
        strcpy((x->p)->st_nfdos,"dos.dat");
        strcpy((x->p)->st_nfthe,"");
        strcpy((x->p)->st_nfdumb,"");
        strcpy((x->p)->st_nfdumb2,"");
        strcpy((x->p)->st_nfhisto,"");
        strcpy((x->p)->st_anneal,"");
        strcpy((x->p)->st_nfthe,"thermdodyn.dat");
        strcpy((x->p)->st_nftout,"temperatures.dat");
        strcpy((x->p)->st_nfdumb,"dumb.dat");
        strcpy((x->p)->st_pdbnf,"harvest.pdb");
	
	(x->p)->st_nprint = 1;
        (x->p)->st_k = 1.;
        (x->p)->st_p_new = 1.;          // >0 means not using it
        (x->p)->st_ntempmax = 30;
        (x->p)->st_ntemp = 1;
        (x->p)->st_debug = 0;
        (x->p)->st_pthresh = 0;
        (x->p)->st_nstep_adj = 100000;
        (x->p)->st_nprintt = 500000;
        (x->p)->st_keepall = 0;
        (x->p)->st_paranoid = 0;
        (x->p)->st_oldh = NULL;
        (x->p)->st_npre = 0;
        (x->p)->st_hthresh = 0.70;
        (x->p)->st_restart = -1;
        (x->p)->st_binthresh = 0;
        (x->p)->st_tonlyadd = 0;
        (x->p)->st_nkeepold = 0;
        (x->p)->st_sum = 0;
        (x->p)->st_removet = 99;
        (x->p)->st_tstop=0;
        (x->p)->st_tnorm=0;
        (x->p)->st_noadjust = 0;
        (x->p)->st_nfail1=0;
        (x->p)->st_nfail2=20;
        (x->p)->st_phthresh = 0.1;
        (x->p)->st_ignoreb=0;
        (x->p)->st_deltat=0;
        (x->p)->st_gmethod=0;
        (x->p)->st_ttarget = 0;
        (x->p)->st_ttarget_harvest = 0;
        (x->p)->st_printpdb=1000;

	#endif

	fp = fopen(fname,"r");
	if (!fp) Error("Cannot open parameter file");

	// defaults

	x->npol = -1;
	x->nstep = 100000;
	x->seed = -1;
	x->dw_flip = 30;
	x->dw_pivot = 10.;
	x->dw_mpivot = 1.;
	x->dw_lpivot = 1.;
	x->dw_mflip = 30.;
	x->dx_com = 1.;
	x->dx_clm = 1.;
	x->dtheta=1.;
	x->randdw = 1;
	x->nprinttrj = 1000;
	x->nprintlog = 1000;
	strcpy(x->fntrj,"traj");
	x->r2shell = 6.;
	x->ntemp = 0;
	x->debug = 0;
	x->shell = 0;
	x->nshell = 10;
	x->nmul_mpivot = 3;
	x->nmul_lpivot = 3;
	x->nmul_mflip=100;
	for (i=0;i<NMOVES;i++) x->movetype[i] = -1;
	x->nosidechains=0;
	strcpy(nlog,"montegrappa.log");
	x->noangpot=0;
	x->nodihpot=0;
	strcpy(x->fne,"energy");
	strcpy(x->flastp,"last");
	strcpy(x->fnproc,"proc");	
	x->nprinte = 1000;
	x->nrun = 1;
	x->always_restart = 0;
	x->record_native = 0;
	x->disentangle = 0;
	x->r_cloose = 0.5;
	x->a_cloose = -1;
	x->d_cloose = -1;
	x->chi2start =0;
	x->ishell =0;
	x->bgs_a=2;
	x->bgs_b=1;
	#ifdef OPTIMIZEPOT
	strcpy(x->fnop,"");
	strcpy(x->op_minim,"none");
	x->op_deltat = 10000;
	x->op_itermax = 1000;
	x->op_step=1.;
	x->op_stop=0.1;
	x->op_T = 1;
	x->op_print = 100;
	x->op_emin = -9E19;
	x->op_emax = 9E19;
	x->op_wait = 0;
	x->op_r = 3.;
	x->op_r0 = 1.;
	#endif
	x->hb=0;
	#ifdef ACTIVE_MPI
	x->nstep_exchange = 10000;
	#endif

//ASTEMPERING	
	x->nconf=-1;



	while(fgets(aux,1000,fp)!=NULL)
	{
                ReadParD(aux,"nconf",&(x->nconf));

		ReadParD(aux,"nchains",&(x->npol));
		ReadParLLU(aux,"nstep",&(x->nstep));
		ReadParL(aux,"seed",&(x->seed));
		ReadParF(aux,"dw_flip",&(x->dw_flip));
		ReadParF(aux,"dw_pivot",&(x->dw_pivot));
		ReadParF(aux,"dw_mpivot",&(x->dw_mpivot));
		ReadParF(aux,"dw_lpivot",&(x->dw_lpivot));
		ReadParF(aux,"dw_mflip",&(x->dw_mflip));
		ReadParF(aux,"dx_com",&(x->dx_com));
		ReadParF(aux,"dtheta",&(x->dtheta));
		ReadParF(aux,"dx_clm",&(x->dx_com));
		ReadParD(aux,"randdw",&(x->randdw));
		ReadParD(aux,"nprinttrj",&(x->nprinttrj));
		ReadParD(aux,"nprintlog",&(x->nprintlog));
		ReadParD(aux,"nprinte",&(x->nprinte));
		ReadParD(aux,"chi2start",&(x->chi2start));
		ReadParS(aux,"traj",x->fntrj);
		ReadParS(aux,"logfile",nlog);
		ReadParS(aux,"lastp",x->flastp);
		ReadParS(aux,"efile",x->fne);
		ReadParS(aux,"procfile",x->fnproc);
		ReadParF(aux,"r2shell",&(x->r2shell));
		ReadParD(aux,"ntemp",&(x->ntemp));
		ReadParD(aux,"debug",&(x->debug));
		ReadParN(aux,"shell",&(x->shell));
		ReadParD(aux,"nshell",&(x->nshell));
		ReadParD(aux,"nmul_mpivot",&(x->nmul_mpivot));
		ReadParD(aux,"nmul_lpivot",&(x->nmul_lpivot));
		ReadParD(aux,"nmul_mflip",&(x->nmul_mflip));
		ReadParD(aux,"flip",&(x->movetype[0]));
		ReadParD(aux,"pivot",&(x->movetype[1]));
		ReadParD(aux,"mpivot",&(x->movetype[2]));
		ReadParD(aux,"sidechain",&(x->movetype[3]));
		ReadParD(aux,"lpivot",&(x->movetype[4]));
		ReadParD(aux,"mflip",&(x->movetype[5]));
		ReadParD(aux,"movecom",&(x->movetype[6]));
		ReadParD(aux,"movebias",&(x->movetype[7]));
		ReadParD(aux,"moverot",&(x->movetype[8]));
		ReadParD(aux,"comcluster",&(x->movetype[9]));
		ReadParD(aux,"rotcluster",&(x->movetype[10]));

		ReadParF(aux,"r_cloose",&(x->r_cloose));
		ReadParF(aux,"a_cloose",&(x->a_cloose));
		ReadParF(aux,"d_cloose",&(x->d_cloose));
		ReadParN(aux,"nosidechain",&(x->nosidechains));
		ReadParN(aux,"noangpot",&(x->noangpot));
		ReadParN(aux,"stempering",&(x->stempering));
		ReadParN(aux,"nodihpot",&(x->nodihpot));
		ReadParN(aux,"disentangle",&(x->disentangle));
		ReadParD(aux,"nrun",&(x->nrun));
		ReadParN(aux,"always_restart",&(x->always_restart));
		ReadParN(aux,"record_native",&(x->record_native));
		ReadParN(aux,"hb",&(x->hb));
		ReadParN(aux,"do_anneal",&(x->anneal));
		ReadParD(aux,"anneal_often",&(x->anneal_often));
		ReadParD(aux,"anneal_step",&(x->anneal_step));
		ReadParD(aux,"anneal_recov",&(x->anneal_recov));
		ReadParF(aux,"anneal_t",&(x->anneal_t));
		 ReadParD(aux,"iT_bias",&(x->iT_bias));
	#ifdef STEMPERING

	if (!strncmp(aux,"st_method",9))
        if ( !sscanf(aux,"st_method %s",(x->p)->st_method) ) FatalError("Cannot read method in inputfile");

       if (!strncmp(aux,"st_anneal",9))
         if ( !sscanf(aux,"st_anneal %s",(x->p)->st_anneal) ) FatalError("Cannot read anneal in inputfile");

       if (!strncmp(aux,"st_ntmax",8))
        if ( !sscanf(aux,"st_ntmax %d",&(x->p)->st_ntempmax) ) FatalError("Cannot read ntmax in inputfile");

       if (!strncmp(aux,"st_nstep",8))
         if ( !sscanf(aux,"st_nstep %d",&(x->p)->st_nstep) ) FatalError("Cannot read nstep in inputfile");

       if (!strncmp(aux,"st_preamble",11))
          if ( !sscanf(aux,"st_preamble %d",&(x->p)->st_npre) ) FatalError("Cannot read nstep in inputfile");

       if (!strncmp(aux,"st_nsadj",8))
         if ( !sscanf(aux,"st_nsadj %d",&(x->p)->st_nstep_adj) ) FatalError("Cannot read nsadj in inputfile");

       if (!strncmp(aux,"st_nprint",9))
         if ( !sscanf(aux,"st_nprint %d",&(x->p)->st_nprint) ) FatalError("Cannot read nprint in inputfile");

        if (!strncmp(aux,"st_printpdb",11))
         if ( !sscanf(aux,"st_printpdb %d",&(x->p)->st_printpdb) ) FatalError("Cannot read printpdb in inputfile");

       if (!strncmp(aux,"st_nprntt",9))
          if ( !sscanf(aux,"st_nprntt %d",&(x->p)->st_nprintt) ) FatalError("Cannot read nprntt in inputfile");

       if (!strncmp(aux,"st_ntemp",8))
         if ( !sscanf(aux,"st_ntemp %d",&(x->p)->st_ntemp) ) FatalError("Cannot read ntemp in inputfile");

       if (!strncmp(aux,"st_debug",8))
          if ( !sscanf(aux,"st_debug %d",&(x->p)->st_debug) ) FatalError("Cannot read debug in inputfile");

	 if (!strncmp(aux,"st_temperatures",15))
       {
           if ((x->p)->st_ntemp == 0) FatalError("You must specify ntemp before listing the temperatures");
           for (i=0;i<(x->p)->st_ntemp;i++)
           {
                   r = fscanf(fp,"%lf %lf",&input_temp[i],&dumb);
                   if (r==2) input_g[i] = dumb;
                   else if (r==1) input_g[i] = 0.;
                   else FatalError("Cannot read temperature in inputfile");
           }
       }

       if (!strncmp(aux,"st_emin",7))
         if ( !sscanf(aux,"st_emin %lf",&(x->p)->st_emin) ) FatalError("Cannot read emin in inputfile");

       if (!strncmp(aux,"st_emax",7))
          if ( !sscanf(aux,"st_emax %lf",&(x->p)->st_emax) ) FatalError("Cannot read emax in inputfile");

       if (!strncmp(aux,"st_ebin",7))
          if ( !sscanf(aux,"st_ebin %lf",&(x->p)->st_ebin) ) FatalError("Cannot read bin in inputfile");

       if (!strncmp(aux,"st_tfile",8))
          if ( !sscanf(aux,"st_tfile %s",(x->p)->st_nftout) ) FatalError("Cannot read tfile in inputfile");

       if (!strncmp(aux,"st_thefile",10))
          if ( !sscanf(aux,"st_thefile %s",(x->p)->st_nfthe) ) FatalError("Cannot read thefile in inputfile");

       if (!strncmp(aux,"st_dosfile",10))
          if ( !sscanf(aux,"st_dosfile %s",(x->p)->st_nfdos) ) FatalError("Cannot read dosfile in inputfile");

       if (!strncmp(aux,"st_histofile",12))
           if ( !sscanf(aux,"st_histofile %s",(x->p)->st_nfhisto) ) FatalError("Cannot read histofile in inputfile");

       if (!strncmp(aux,"st_dumbfile",11))
          if ( !sscanf(aux,"st_dumbfile %s",(x->p)->st_nfdumb) ) FatalError("Cannot read dumbfile in inputfile");

       if (!strncmp(aux,"st_currentdumbfile",18))
          if ( !sscanf(aux,"st_currentdumbfile %s",(x->p)->st_nfdumb2) ) FatalError("Cannot read currentdumbfile in inputfile");

       if (!strncmp(aux,"st_k_new",8))
           if ( !sscanf(aux,"st_k_new %lf",&(x->p)->st_k) ) FatalError("Cannot read k_new in inputfile");

       if (!strncmp(aux,"st_lpthresh",11))
           if ( !sscanf(aux,"st_lpthresh %lf",&(x->p)->st_pthresh) ) FatalError("Cannot read lpthresh in inputfile");

 if (!strncmp(aux,"st_hthresh",10))
                 if ( !sscanf(aux,"st_hthresh %lf",&(x->p)->st_hthresh) ) FatalError("Cannot read hthresh in inputfile");

       if (!strncmp(aux,"st_keepall",10)) (x->p)->st_keepall = 1;

       if (!strncmp(aux,"st_paranoid",11)) (x->p)->st_paranoid = 1;

       if (!strncmp(aux,"st_restart",10))
           if ( !sscanf(aux,"st_restart %d",&(x->p)->st_restart) ) FatalError("Cannot read restart in inputfile");

       if (!strncmp(aux,"st_lp_new",9))
           if ( !sscanf(aux,"st_lp_new %lf",&(x->p)->st_p_new) ) FatalError("Cannot read lp_new in inputfile");

       if (!strncmp(aux,"st_t_only_add",13))
                 if ( !sscanf(aux,"st_t_only_add %lf",&(x->p)->st_tonlyadd) ) FatalError("Cannot read t_only_add in inputfile");

       if (!strncmp(aux,"st_binthresh",12))
                        if ( !sscanf(aux,"st_binthresh %lf",&(x->p)->st_binthresh) ) FatalError("Cannot read binthresh in inputfile");

       if (!strncmp(aux,"st_nkeepold",11))
           if ( !sscanf(aux,"st_nkeepold %d",&(x->p)->st_nkeepold) ) FatalError("Cannot read nkeepold in inputfile");

       if (!strncmp(aux,"st_sumoldhisto",14))
           if ( !sscanf(aux,"st_sumoldhisto %d",&(x->p)->st_sum) ) FatalError("Cannot read sumoldhisto in inputfile");

       if (!strncmp(aux,"st_removet",10))
           if ( !sscanf(aux,"st_removet %d",&(x->p)->st_removet) ) FatalError("Cannot read removet in inputfile");

       if (!strncmp(aux,"st_tstop",8))
           if ( !sscanf(aux,"st_tstop %lf",&(x->p)->st_tstop) ) FatalError("Cannot read tstop in inputfile");


 if (!strncmp(aux,"st_ttarget",10))
           if ( !sscanf(aux,"st_ttarget %lf",&(x->p)->st_ttarget) ) FatalError("Cannot read ttarget in inputfile");

       if (!strncmp(aux,"st_tnorm",8))
           if ( !sscanf(aux,"st_tnorm %lf",&(x->p)->st_tnorm) ) FatalError("Cannot read tnorm in inputfile");

        if (!strncmp(aux,"st_nfail1",9))
           if ( !sscanf(aux,"st_nfail1 %d",&(x->p)->st_nfail1) ) FatalError("Cannot read nfail1 in inputfile");

        if (!strncmp(aux,"st_nfail2",9))
           if ( !sscanf(aux,"st_nfail2 %d",&(x->p)->st_nfail2) ) FatalError("Cannot read nfail2 in inputfile");

        if (!strncmp(aux,"st_phthresh",11))
               if ( !sscanf(aux,"st_phthresh %lf",&(x->p)->st_phthresh) ) FatalError("Cannot read phthresh in inputfile");

        if (!strncmp(aux,"st_ignoreb",10))
                if ( !sscanf(aux,"st_ignoreb %d",&(x->p)->st_ignoreb) ) FatalError("Cannot read ignoreb in inputfile");

        if (!strncmp(aux,"st_deltat",9))
                if ( !sscanf(aux,"st_deltat %lf",&(x->p)->st_deltat) ) FatalError("Cannot read deltat in inputfile");
        if (!strncmp(aux,"st_gmethod",10))
                if ( !sscanf(aux,"st_gmethod %s",aux_gm) ) FatalError("Cannot read gmethod in inputfile");
        if (!strncmp(aux,"st_pdbnf",8))
                if ( !sscanf(aux,"st_pdbnf %s",(x->p)->st_pdbnf) ) FatalError("Cannot read pdbnf in inputfile");


    

        if (!strcmp((x->p)->st_anneal,"sigma")) (x->p)->st_p_new = 0;
        else if (!strcmp((x->p)->st_anneal,"prob"))     (x->p)->st_k = 0;




		#endif


		#ifdef OPTIMIZEPOT
		 ReadParS(aux,"op_minim",(x->op_minim));
		 ReadParS(aux,"op_file",x->fnop);
		 ReadParD(aux,"op_deltat",&(x->op_deltat));
		 ReadParD(aux,"op_itermax",&(x->op_itermax));
		 ReadParD(aux,"op_print",&(x->op_print));
		 ReadParD(aux,"op_wait",&(x->op_wait));
		 ReadParF(aux,"op_step",&(x->op_step));
		 ReadParF(aux,"op_stop",&(x->op_stop));
		 ReadParF(aux,"op_T",&(x->op_T));
		 ReadParF(aux,"op_emin",&(x->op_emin));
		 ReadParF(aux,"op_emax",&(x->op_emax));
		 ReadParF(aux,"op_rw",&(x->op_r));
		 ReadParF(aux,"op_r0",&(x->op_r0));
		 ReadParD(aux,"nmul_local",&(x->nmul_local));
		 ReadParF(aux,"bgs_a",&(x->bgs_a));
		 ReadParF(aux,"bgs_b",&(x->bgs_b));
		#endif
		#ifdef ACTIVE_MPI
		 ReadParD(aux,"step_exchange",&(x->nstep_exchange));
		 
		#endif
		#ifndef ACTIVE_MPI
		ReadParF(aux,"Temp",&(x->T));
		#endif

		#ifdef ACTIVE_MPI
		if (!strncmp(aux,"temperatures",12))
       		{
    	   		if (x->ntemp==0) Error("You must specify ntemp before listing the temperatures");
			if(x->ntemp>NREPMAX) Error("NREPMAX too small");
			int r;
			fprintf(stderr, "- PT Temperatures:\n");
			for(i=0;i<x->ntemp;i++)
			{
				r = fscanf(fp,"%lf",&(x->T[i]));
				if(r!=1) Error("Cannot read temperature in inputfile");
				fprintf(stderr,"  rank %d\t%f\t",i,x->T[i]);
    			}
			fprintf(stderr,"\n");
      		}	
		#endif
	}
//	#ifndef ACTIVE_MPI
//	if((x->ntemp)!=1) Error("ntemp must be 1");
//	#endif

	if(x->nrun>1 && !strcmp(x->op_minim,"none")) fprintf(stderr,"WARNING: nrun > 1 without optimization parameters\n");

	#ifndef STEMPERING
	if(x->stempering) Error("Please compile with version=STEMPERING if you want to enable simulated tempering!!!\n");	
	#endif
	
	// check compulsory data
	if (x->npol<1) Error("no chain defined in parameter file");
	if (x->npol>NCHAINMAX) Error("Too many chains defined in parameter file");
	chk=-1;
	for (i=0;i<NMOVES;i++) if ( x->movetype[i] != -1 ) chk=1;
	if (chk==-1) Error("No move type defined");


	x->flog=fopen(nlog,"w");
	if (!x->flog) Error("Cannot open log file");
	


	#ifdef STEMPERING
	if (x->stempering){
	//restart
	if ((x->p)->st_restart > -1)
                {
                        Restart(x->p);
                        if ((x->p)->st_debug>0) fprintf(stderr,"Start the simulation...\n\n");
                   //     return x->p;
                }

        if ((x->p)->st_debug>0)
        {
                fprintf(stderr,"\n************************************\n");
                fprintf(stderr,"* SIMULATED TEMPERING              *\n");
                fprintf(stderr,"************************************\n\n");
                fflush(stderr);

                if (!strcmp((x->p)->st_method,"stempering") || !strcmp((x->p)->st_method,"adaptive"))
                        fprintf(stderr,"method =\t%s\n",(x->p)->st_method);
                else FatalError("Method not known");
                fprintf(stderr,"ntempmax = \t%d\n",(x->p)->st_ntempmax);
                fprintf(stderr,"ntemp = \t%d\n",(x->p)->st_ntemp);
                fprintf(stderr,"temperatures:\n");
                for (i=0;i<(x->p)->st_ntemp;i++)
                        fprintf(stderr,"T[%d] = %lf\n",i,input_temp[i]);
                fprintf(stderr,"weights:\n");
                for (i=0;i<(x->p)->st_ntemp;i++)
                        fprintf(stderr,"g[%d] = %lf\n",i,input_g[i]);
                fprintf(stderr,"nstep = \t%d\n",(x->p)->st_nstep);
                if (!strcmp((x->p)->st_method,"adaptive"))
                                fprintf(stderr,"nsadj = \t%d\n",(x->p)->st_nstep_adj);
                if ((x->p)->st_npre>0) fprintf(stderr,"preamble = \t%d\n",(x->p)->st_npre);
                fprintf(stderr,"nprint = \t%d\n",(x->p)->st_nprint);
                fprintf(stderr,"nprintt = \t%d\n",(x->p)->st_nprintt);
                fprintf(stderr,"debug = \t%d\n",(x->p)->st_debug);
                fprintf(stderr,"emin = \t\t%lf\n",(x->p)->st_emin);
                fprintf(stderr,"emax = \t\t%lf\n",(x->p)->st_emax);
                fprintf(stderr,"ebin = \t\t%lf\n",(x->p)->st_ebin);
                fprintf(stderr,"nbin = \t\t%d\n",(int)(((x->p)->st_emax-(x->p)->st_emin)/(x->p)->st_ebin));
                fprintf(stderr,"binthresh = \t%lf\n",(x->p)->st_binthresh);
                fprintf(stderr,"phthresh = \t%lf\n",(x->p)->st_phthresh);
                if ((x->p)->st_tstop>0) fprintf(stderr,"tstop = \t\t%lf\n",(x->p)->st_tstop);
                if ((x->p)->st_tnorm>0) fprintf(stderr,"tnorm = \t\t%lf\n",(x->p)->st_tnorm);

                if (!strcmp((x->p)->st_anneal,"sigma"))
                {
                        fprintf(stderr,"Anneal = \tproportionally to energy stdev\n");
                        fprintf(stderr,"k_new = \t%lf\n",(x->p)->st_k);
                }


	else if (!strcmp((x->p)->st_anneal,"prob"))
                {
                        fprintf(stderr,"Anneal = \tsetting wished exchange probability\n");
                        fprintf(stderr,"p_new = \t%lf\n",(x->p)->st_p_new);
                }
                else
                        fprintf(stderr,"Anneal = \tnone (fixed range of temperatures\n");

                if ((x->p)->st_removet>0)
                {
                        if ((x->p)->st_pthresh<0)
                                fprintf(stderr,"lpthresh = \t%lf (to reduce the number of temperatures)\n",(x->p)->st_pthresh);
                        else if ((x->p)->st_pthresh>8)
                                fprintf(stderr,"lpthresh = automatically calculated\n");
                        fprintf(stderr,"removet = \t%d\n",(x->p)->st_removet);
                }
                fprintf(stderr,"hthresh = \t%lf (to keep a histogram)\n",(x->p)->st_hthresh);

                if ((x->p)->st_keepall==1 && (x->p)->st_nkeepold==0) fprintf(stderr,"Calculate g(E) from all past history\n");
                else if ((x->p)->st_keepall==1) fprintf(stderr,"Calculate g(E) from last %d histograms\n",(x->p)->st_nkeepold);
                if ((x->p)->st_keepall==1 && (x->p)->st_sum==1) fprintf(stderr,"Sum together histograms with same temperature\n");
                if ((x->p)->st_keepall==1 && (x->p)->st_sum==2) fprintf(stderr,"Delete histograms with same temperature than newer\n");

                if ((x->p)->st_paranoid==1) fprintf(stderr,"Be paranoid\n");
                if ((x->p)->st_tonlyadd>0) fprintf(stderr,"Below T=%lf, only add temperatures\n",(x->p)->st_tonlyadd);
                if (strcmp((x->p)->st_nfhisto,"")) fprintf(stderr,"Write histograms to %s\n",(x->p)->st_nfhisto);
        }

        // Check the allocation dimensions

        if (NHISTOMAX>NTEMPMAX) FatalError("NHISTOMAX should be smaller than NTEMPMAX");
        if ((x->p)->st_ntempmax>NTEMPMAX) FatalError("ntmax should be smaller than NTEMPMAX");
        if ((x->p)->st_ntemp>(x->p)->st_ntempmax) FatalError("ntemp should be smaller than ntmax");

        // Initialize common stuff
        (x->p)->st_ftout = fopen((x->p)->st_nftout,"w");
        if ((x->p)->st_ftout==NULL) FatalError("Error opening 'tfile' for writing");
        fprintf((x->p)->st_ftout,"# STEP\tTEMP\tNTEMP\n");
        (x->p)->st_f = NULL;
        (x->p)->st_nbin = (int) (((x->p)->st_emax-(x->p)->st_emin)/(x->p)->st_ebin) + 1;
        
	(x->p)->st_h = AlloDoubleMat((x->p)->st_ntempmax,(x->p)->st_nbin);
        (x->p)->st_binok = AlloInt((x->p)->st_nbin);
        (x->p)->st_htmp = AlloDoubleMat(NHISTOMAX,(x->p)->st_nbin);
    fprintf(stderr,"allocate htmp %dx%d\n",NHISTOMAX,(x->p)->st_nbin);

        (x->p)->st_out = AlloDoubleMat(4,NTBIN);
        (x->p)->st_boltzp = AlloDoubleMat(NHISTOMAX,(x->p)->st_nbin);

        // Initialize standard simulated tempering
        if (!strcmp((x->p)->st_method,"stempering"))
        {
                if ((x->p)->st_ntemp<2) FatalError("You need more than two tempering to make a stempering");
                (x->p)->st_nm = 1;
                (x->p)->st_temp = AlloDouble((x->p)->st_ntemp);
                (x->p)->st_g = AlloDouble((x->p)->st_ntemp);
                (x->p)->st_prob_down = AlloDouble((x->p)->st_ntemp);
                (x->p)->st_prob_up = AlloDouble((x->p)->st_ntemp);
                (x->p)->st_counts = AlloInt((x->p)->st_ntemp);

                for (i=0;i<(x->p)->st_ntemp;i++)
                {
                        (x->p)->st_temp[i] = input_temp[i];
                        (x->p)->st_g[i] = input_g[i];
                }
        }

        // Initialize adaptive simulated tempering
         else if (!strcmp((x->p)->st_method,"adaptive"))
        {
                (x->p)->st_nm = 2;
                (x->p)->st_temp = AlloDouble((x->p)->st_ntempmax);
                (x->p)->st_g = AlloDouble((x->p)->st_ntempmax);
                (x->p)->st_prob_down = AlloDouble((x->p)->st_ntempmax);
                (x->p)->st_prob_up = AlloDouble((x->p)->st_ntempmax);
                (x->p)->st_counts = AlloInt((x->p)->st_ntempmax);
                (x->p)->st_f = AlloDouble((x->p)->st_ntempmax);
                (x->p)->st_reliable_lg = AlloDouble((x->p)->st_nbin);
                (x->p)->st_current_lg = AlloDouble((x->p)->st_nbin);
                (x->p)->st_reliable_t = AlloDouble((x->p)->st_ntempmax);
                (x->p)->st_reliable_g = AlloDouble((x->p)->st_ntempmax);

                for (i=0;i<(x->p)->st_ntemp;i++)
                        (x->p)->st_g[i] = input_g[i];

                if ((x->p)->st_keepall==1)
                {
                        (x->p)->st_oldh = AlloDoubleMat(NHISTOMAX,(x->p)->st_nbin);
                        (x->p)->st_oldt = AlloDouble(NHISTOMAX);
                        (x->p)->st_oldt_iter = AlloInt(NHISTOMAX);
                        (x->p)->st_noldt = 0;

                        if ((x->p)->st_debug>0) fprintf(stderr,"Allocate memory to keep all histograms up to %d.\n",NHISTOMAX);
                }

                //(x->p)->st_fpacc = fopen("accept.dat","w");
                //if (!(x->p)->st_fpacc) FatalError("cannot open accetp.dat");

        }

        for (i=0;i<(x->p)->st_ntemp;i++)
                        (x->p)->st_temp[i] = input_temp[i];

	// structure for restart
	if((x->p)->st_restart>-1)
	(x->p)->st_st_restart = AlloRestart(NTEMPMAX,(x->p)->st_nbin,NRESTARTS);

        // Start from highest temperature
         (x->p)->st_itemp = 0;

        // Reset counters
	 (x->p)->st_count_st = 0;
        (x->p)->st_count_adj = 0;
        (x->p)->st_count_print = 0;
        (x->p)->st_count_printt = 0;
        (x->p)->st_iter = 0;
        (x->p)->st_failure = 0;



}


	#endif	




	return x;
}


#ifdef STEMPERING
void ResetStStuff(struct s_mc_parms *x,char *fname){
	
	int i,r;
        double input_temp[NTEMPMAX],input_g[NTEMPMAX],dumb;
        char aux[500],aux_gm[500];
	FILE *fp;
        x->p = (struct st_stuff *) calloc(1,sizeof(struct st_stuff));
	

	strcpy((x->p)->st_nftout,"temperatures.dat");
        strcpy((x->p)->st_nfdos,"dos.dat");
        strcpy((x->p)->st_nfthe,"");
        strcpy((x->p)->st_nfdumb,"");
        strcpy((x->p)->st_nfdumb2,"");
        strcpy((x->p)->st_nfhisto,"");
        strcpy((x->p)->st_anneal,"");
        strcpy((x->p)->st_nfthe,"thermdodyn.dat");
        strcpy((x->p)->st_nftout,"temperatures.dat");
        strcpy((x->p)->st_nfdumb,"dumb.dat");
        strcpy((x->p)->st_pdbnf,"harvest.pdb");
	
	(x->p)->st_nprint = 1;
        (x->p)->st_k = 1.;
        (x->p)->st_p_new = 1.;          // >0 means not using it
        (x->p)->st_ntempmax = 30;
        (x->p)->st_ntemp = 1;
        (x->p)->st_debug = 0;
        (x->p)->st_pthresh = 0;
        (x->p)->st_nstep_adj = 100000;
        (x->p)->st_nprintt = 500000;
        (x->p)->st_keepall = 0;
        (x->p)->st_paranoid = 0;
        (x->p)->st_oldh = NULL;
        (x->p)->st_npre = 0;
        (x->p)->st_hthresh = 0.70;
        (x->p)->st_restart = -1;
        (x->p)->st_binthresh = 0;
        (x->p)->st_tonlyadd = 0;
        (x->p)->st_nkeepold = 0;
        (x->p)->st_sum = 0;
        (x->p)->st_removet = 99;
        (x->p)->st_tstop=0;
        (x->p)->st_tnorm=0;
        (x->p)->st_noadjust = 0;
        (x->p)->st_nfail1=0;
        (x->p)->st_nfail2=20;
        (x->p)->st_phthresh = 0.1;
        (x->p)->st_ignoreb=0;
        (x->p)->st_deltat=0;
        (x->p)->st_gmethod=0;
        (x->p)->st_ttarget = 0;
        (x->p)->st_ttarget_harvest = 0;
        (x->p)->st_printpdb=1000;

	fp = fopen(fname,"r");
	if (!fp) Error("Cannot open parameter file");

	while(fgets(aux,500,fp)!=NULL)
	{

	if (!strncmp(aux,"st_method",9))
        if ( !sscanf(aux,"st_method %s",(x->p)->st_method) ) FatalError("Cannot read method in inputfile");

       if (!strncmp(aux,"st_anneal",9))
         if ( !sscanf(aux,"st_anneal %s",(x->p)->st_anneal) ) FatalError("Cannot read anneal in inputfile");

       if (!strncmp(aux,"st_ntmax",8))
        if ( !sscanf(aux,"st_ntmax %d",&(x->p)->st_ntempmax) ) FatalError("Cannot read ntmax in inputfile");

       if (!strncmp(aux,"st_nstep",8))
         if ( !sscanf(aux,"st_nstep %d",&(x->p)->st_nstep) ) FatalError("Cannot read nstep in inputfile");

       if (!strncmp(aux,"st_preamble",11))
          if ( !sscanf(aux,"st_preamble %d",&(x->p)->st_npre) ) FatalError("Cannot read nstep in inputfile");

       if (!strncmp(aux,"st_nsadj",8))
         if ( !sscanf(aux,"st_nsadj %d",&(x->p)->st_nstep_adj) ) FatalError("Cannot read nsadj in inputfile");

       if (!strncmp(aux,"st_nprint",9))
         if ( !sscanf(aux,"st_nprint %d",&(x->p)->st_nprint) ) FatalError("Cannot read nprint in inputfile");

        if (!strncmp(aux,"st_printpdb",11))
         if ( !sscanf(aux,"st_printpdb %d",&(x->p)->st_printpdb) ) FatalError("Cannot read printpdb in inputfile");

       if (!strncmp(aux,"st_nprntt",9))
          if ( !sscanf(aux,"st_nprntt %d",&(x->p)->st_nprintt) ) FatalError("Cannot read nprntt in inputfile");

       if (!strncmp(aux,"st_ntemp",8))
         if ( !sscanf(aux,"st_ntemp %d",&(x->p)->st_ntemp) ) FatalError("Cannot read ntemp in inputfile");

       if (!strncmp(aux,"st_debug",8))
          if ( !sscanf(aux,"st_debug %d",&(x->p)->st_debug) ) FatalError("Cannot read debug in inputfile");

	 if (!strncmp(aux,"st_temperatures",15))
       {
           if ((x->p)->st_ntemp == 0) FatalError("You must specify ntemp before listing the temperatures");
           for (i=0;i<(x->p)->st_ntemp;i++)
           {
                   r = fscanf(fp,"%lf %lf",&input_temp[i],&dumb);
                   if (r==2) input_g[i] = dumb;
                   else if (r==1) input_g[i] = 0.;
                   else FatalError("Cannot read temperature in inputfile");
           }
       }

       if (!strncmp(aux,"st_emin",7))
         if ( !sscanf(aux,"st_emin %lf",&(x->p)->st_emin) ) FatalError("Cannot read emin in inputfile");

       if (!strncmp(aux,"st_emax",7))
          if ( !sscanf(aux,"st_emax %lf",&(x->p)->st_emax) ) FatalError("Cannot read emax in inputfile");

       if (!strncmp(aux,"st_ebin",7))
          if ( !sscanf(aux,"st_ebin %lf",&(x->p)->st_ebin) ) FatalError("Cannot read bin in inputfile");

       if (!strncmp(aux,"st_tfile",8))
          if ( !sscanf(aux,"st_tfile %s",(x->p)->st_nftout) ) FatalError("Cannot read tfile in inputfile");

       if (!strncmp(aux,"st_thefile",10))
          if ( !sscanf(aux,"st_thefile %s",(x->p)->st_nfthe) ) FatalError("Cannot read thefile in inputfile");

       if (!strncmp(aux,"st_dosfile",10))
          if ( !sscanf(aux,"st_dosfile %s",(x->p)->st_nfdos) ) FatalError("Cannot read dosfile in inputfile");

       if (!strncmp(aux,"st_histofile",12))
           if ( !sscanf(aux,"st_histofile %s",(x->p)->st_nfhisto) ) FatalError("Cannot read histofile in inputfile");

       if (!strncmp(aux,"st_dumbfile",11))
          if ( !sscanf(aux,"st_dumbfile %s",(x->p)->st_nfdumb) ) FatalError("Cannot read dumbfile in inputfile");

       if (!strncmp(aux,"st_currentdumbfile",18))
          if ( !sscanf(aux,"st_currentdumbfile %s",(x->p)->st_nfdumb2) ) FatalError("Cannot read currentdumbfile in inputfile");

       if (!strncmp(aux,"st_k_new",8))
           if ( !sscanf(aux,"st_k_new %lf",&(x->p)->st_k) ) FatalError("Cannot read k_new in inputfile");

       if (!strncmp(aux,"st_lpthresh",11))
           if ( !sscanf(aux,"st_lpthresh %lf",&(x->p)->st_pthresh) ) FatalError("Cannot read lpthresh in inputfile");

 if (!strncmp(aux,"st_hthresh",10))
                 if ( !sscanf(aux,"st_hthresh %lf",&(x->p)->st_hthresh) ) FatalError("Cannot read hthresh in inputfile");

       if (!strncmp(aux,"st_keepall",10)) (x->p)->st_keepall = 1;

       if (!strncmp(aux,"st_paranoid",11)) (x->p)->st_paranoid = 1;

       if (!strncmp(aux,"st_restart",10))
           if ( !sscanf(aux,"st_restart %d",&(x->p)->st_restart) ) FatalError("Cannot read restart in inputfile");

       if (!strncmp(aux,"st_lp_new",9))
           if ( !sscanf(aux,"st_lp_new %lf",&(x->p)->st_p_new) ) FatalError("Cannot read lp_new in inputfile");

       if (!strncmp(aux,"st_t_only_add",13))
                 if ( !sscanf(aux,"st_t_only_add %lf",&(x->p)->st_tonlyadd) ) FatalError("Cannot read t_only_add in inputfile");

       if (!strncmp(aux,"st_binthresh",12))
                        if ( !sscanf(aux,"st_binthresh %lf",&(x->p)->st_binthresh) ) FatalError("Cannot read binthresh in inputfile");

       if (!strncmp(aux,"st_nkeepold",11))
           if ( !sscanf(aux,"st_nkeepold %d",&(x->p)->st_nkeepold) ) FatalError("Cannot read nkeepold in inputfile");

       if (!strncmp(aux,"st_sumoldhisto",14))
           if ( !sscanf(aux,"st_sumoldhisto %d",&(x->p)->st_sum) ) FatalError("Cannot read sumoldhisto in inputfile");

       if (!strncmp(aux,"st_removet",10))
           if ( !sscanf(aux,"st_removet %d",&(x->p)->st_removet) ) FatalError("Cannot read removet in inputfile");

       if (!strncmp(aux,"st_tstop",8))
           if ( !sscanf(aux,"st_tstop %lf",&(x->p)->st_tstop) ) FatalError("Cannot read tstop in inputfile");


 if (!strncmp(aux,"st_ttarget",10))
           if ( !sscanf(aux,"st_ttarget %lf",&(x->p)->st_ttarget) ) FatalError("Cannot read ttarget in inputfile");

       if (!strncmp(aux,"st_tnorm",8))
           if ( !sscanf(aux,"st_tnorm %lf",&(x->p)->st_tnorm) ) FatalError("Cannot read tnorm in inputfile");

        if (!strncmp(aux,"st_nfail1",9))
           if ( !sscanf(aux,"st_nfail1 %d",&(x->p)->st_nfail1) ) FatalError("Cannot read nfail1 in inputfile");

        if (!strncmp(aux,"st_nfail2",9))
           if ( !sscanf(aux,"st_nfail2 %d",&(x->p)->st_nfail2) ) FatalError("Cannot read nfail2 in inputfile");

        if (!strncmp(aux,"st_phthresh",11))
               if ( !sscanf(aux,"st_phthresh %lf",&(x->p)->st_phthresh) ) FatalError("Cannot read phthresh in inputfile");

        if (!strncmp(aux,"st_ignoreb",10))
                if ( !sscanf(aux,"st_ignoreb %d",&(x->p)->st_ignoreb) ) FatalError("Cannot read ignoreb in inputfile");

        if (!strncmp(aux,"st_deltat",9))
                if ( !sscanf(aux,"st_deltat %lf",&(x->p)->st_deltat) ) FatalError("Cannot read deltat in inputfile");
        if (!strncmp(aux,"st_gmethod",10))
                if ( !sscanf(aux,"st_gmethod %s",aux_gm) ) FatalError("Cannot read gmethod in inputfile");
        if (!strncmp(aux,"st_pdbnf",8))
                if ( !sscanf(aux,"st_pdbnf %s",(x->p)->st_pdbnf) ) FatalError("Cannot read pdbnf in inputfile");


    

        if (!strcmp((x->p)->st_anneal,"sigma")) (x->p)->st_p_new = 0;
        else if (!strcmp((x->p)->st_anneal,"prob"))     (x->p)->st_k = 0;
	}

	if (x->stempering){
	//restart
	if ((x->p)->st_restart > -1)
                {
                        Restart(x->p);
                        if ((x->p)->st_debug>0) fprintf(stderr,"Start the simulation...\n\n");
                   //     return x->p;
                }

        if ((x->p)->st_debug>0)
        {
                fprintf(stderr,"\n************************************\n");
                fprintf(stderr,"* SIMULATED TEMPERING              *\n");
                fprintf(stderr,"************************************\n\n");
                fflush(stderr);

                if (!strcmp((x->p)->st_method,"stempering") || !strcmp((x->p)->st_method,"adaptive"))
                        fprintf(stderr,"method =\t%s\n",(x->p)->st_method);
                else FatalError("Method not known");
                fprintf(stderr,"ntempmax = \t%d\n",(x->p)->st_ntempmax);
                fprintf(stderr,"ntemp = \t%d\n",(x->p)->st_ntemp);
                fprintf(stderr,"temperatures:\n");
                for (i=0;i<(x->p)->st_ntemp;i++)
                        fprintf(stderr,"T[%d] = %lf\n",i,input_temp[i]);
                fprintf(stderr,"weights:\n");
                for (i=0;i<(x->p)->st_ntemp;i++)
                        fprintf(stderr,"g[%d] = %lf\n",i,input_g[i]);
                fprintf(stderr,"nstep = \t%d\n",(x->p)->st_nstep);
                if (!strcmp((x->p)->st_method,"adaptive"))
                                fprintf(stderr,"nsadj = \t%d\n",(x->p)->st_nstep_adj);
                if ((x->p)->st_npre>0) fprintf(stderr,"preamble = \t%d\n",(x->p)->st_npre);
                fprintf(stderr,"nprint = \t%d\n",(x->p)->st_nprint);
                fprintf(stderr,"nprintt = \t%d\n",(x->p)->st_nprintt);
                fprintf(stderr,"debug = \t%d\n",(x->p)->st_debug);
                fprintf(stderr,"emin = \t\t%lf\n",(x->p)->st_emin);
                fprintf(stderr,"emax = \t\t%lf\n",(x->p)->st_emax);
                fprintf(stderr,"ebin = \t\t%lf\n",(x->p)->st_ebin);
                fprintf(stderr,"nbin = \t\t%d\n",(int)(((x->p)->st_emax-(x->p)->st_emin)/(x->p)->st_ebin));
                fprintf(stderr,"binthresh = \t%lf\n",(x->p)->st_binthresh);
                fprintf(stderr,"phthresh = \t%lf\n",(x->p)->st_phthresh);
                if ((x->p)->st_tstop>0) fprintf(stderr,"tstop = \t\t%lf\n",(x->p)->st_tstop);
                if ((x->p)->st_tnorm>0) fprintf(stderr,"tnorm = \t\t%lf\n",(x->p)->st_tnorm);

                if (!strcmp((x->p)->st_anneal,"sigma"))
                {
                        fprintf(stderr,"Anneal = \tproportionally to energy stdev\n");
                        fprintf(stderr,"k_new = \t%lf\n",(x->p)->st_k);
                }


	else if (!strcmp((x->p)->st_anneal,"prob"))
                {
                        fprintf(stderr,"Anneal = \tsetting wished exchange probability\n");
                        fprintf(stderr,"p_new = \t%lf\n",(x->p)->st_p_new);
                }
                else
                        fprintf(stderr,"Anneal = \tnone (fixed range of temperatures\n");

                if ((x->p)->st_removet>0)
                {
                        if ((x->p)->st_pthresh<0)
                                fprintf(stderr,"lpthresh = \t%lf (to reduce the number of temperatures)\n",(x->p)->st_pthresh);
                        else if ((x->p)->st_pthresh>8)
                                fprintf(stderr,"lpthresh = automatically calculated\n");
                        fprintf(stderr,"removet = \t%d\n",(x->p)->st_removet);
                }
                fprintf(stderr,"hthresh = \t%lf (to keep a histogram)\n",(x->p)->st_hthresh);

                if ((x->p)->st_keepall==1 && (x->p)->st_nkeepold==0) fprintf(stderr,"Calculate g(E) from all past history\n");
                else if ((x->p)->st_keepall==1) fprintf(stderr,"Calculate g(E) from last %d histograms\n",(x->p)->st_nkeepold);
                if ((x->p)->st_keepall==1 && (x->p)->st_sum==1) fprintf(stderr,"Sum together histograms with same temperature\n");
                if ((x->p)->st_keepall==1 && (x->p)->st_sum==2) fprintf(stderr,"Delete histograms with same temperature than newer\n");

                if ((x->p)->st_paranoid==1) fprintf(stderr,"Be paranoid\n");
                if ((x->p)->st_tonlyadd>0) fprintf(stderr,"Below T=%lf, only add temperatures\n",(x->p)->st_tonlyadd);
                if (strcmp((x->p)->st_nfhisto,"")) fprintf(stderr,"Write histograms to %s\n",(x->p)->st_nfhisto);
        }

        // Check the allocation dimensions

        if (NHISTOMAX>NTEMPMAX) FatalError("NHISTOMAX should be smaller than NTEMPMAX");
        if ((x->p)->st_ntempmax>NTEMPMAX) FatalError("ntmax should be smaller than NTEMPMAX");
        if ((x->p)->st_ntemp>(x->p)->st_ntempmax) FatalError("ntemp should be smaller than ntmax");

        // Initialize common stuff
        (x->p)->st_ftout = fopen((x->p)->st_nftout,"w");
        if ((x->p)->st_ftout==NULL) FatalError("Error opening 'tfile' for writing");
        fprintf((x->p)->st_ftout,"# STEP\tTEMP\tNTEMP\n");
        (x->p)->st_f = NULL;
        (x->p)->st_nbin = (int) (((x->p)->st_emax-(x->p)->st_emin)/(x->p)->st_ebin) + 1;
        
	(x->p)->st_h = AlloDoubleMat((x->p)->st_ntempmax,(x->p)->st_nbin);
        (x->p)->st_binok = AlloInt((x->p)->st_nbin);
        (x->p)->st_htmp = AlloDoubleMat(NHISTOMAX,(x->p)->st_nbin);
    fprintf(stderr,"allocate htmp %dx%d\n",NHISTOMAX,(x->p)->st_nbin);

        (x->p)->st_out = AlloDoubleMat(4,NTBIN);
        (x->p)->st_boltzp = AlloDoubleMat(NHISTOMAX,(x->p)->st_nbin);

        // Initialize standard simulated tempering
        if (!strcmp((x->p)->st_method,"stempering"))
        {
                if ((x->p)->st_ntemp<2) FatalError("You need more than two tempering to make a stempering");
                (x->p)->st_nm = 1;
                (x->p)->st_temp = AlloDouble((x->p)->st_ntemp);
                (x->p)->st_g = AlloDouble((x->p)->st_ntemp);
                (x->p)->st_prob_down = AlloDouble((x->p)->st_ntemp);
                (x->p)->st_prob_up = AlloDouble((x->p)->st_ntemp);
                (x->p)->st_counts = AlloInt((x->p)->st_ntemp);

                for (i=0;i<(x->p)->st_ntemp;i++)
                {
                        (x->p)->st_temp[i] = input_temp[i];
                        (x->p)->st_g[i] = input_g[i];
                }
        }

        // Initialize adaptive simulated tempering
         else if (!strcmp((x->p)->st_method,"adaptive"))
        {
                (x->p)->st_nm = 2;
                (x->p)->st_temp = AlloDouble((x->p)->st_ntempmax);
                (x->p)->st_g = AlloDouble((x->p)->st_ntempmax);
                (x->p)->st_prob_down = AlloDouble((x->p)->st_ntempmax);
                (x->p)->st_prob_up = AlloDouble((x->p)->st_ntempmax);
                (x->p)->st_counts = AlloInt((x->p)->st_ntempmax);
                (x->p)->st_f = AlloDouble((x->p)->st_ntempmax);
                (x->p)->st_reliable_lg = AlloDouble((x->p)->st_nbin);
                (x->p)->st_current_lg = AlloDouble((x->p)->st_nbin);
                (x->p)->st_reliable_t = AlloDouble((x->p)->st_ntempmax);
                (x->p)->st_reliable_g = AlloDouble((x->p)->st_ntempmax);

                for (i=0;i<(x->p)->st_ntemp;i++)
                        (x->p)->st_g[i] = input_g[i];

                if ((x->p)->st_keepall==1)
                {
                        (x->p)->st_oldh = AlloDoubleMat(NHISTOMAX,(x->p)->st_nbin);
                        (x->p)->st_oldt = AlloDouble(NHISTOMAX);
                        (x->p)->st_oldt_iter = AlloInt(NHISTOMAX);
                        (x->p)->st_noldt = 0;

                        if ((x->p)->st_debug>0) fprintf(stderr,"Allocate memory to keep all histograms up to %d.\n",NHISTOMAX);
                }

                //(x->p)->st_fpacc = fopen("accept.dat","w");
                //if (!(x->p)->st_fpacc) FatalError("cannot open accetp.dat");

        }

        for (i=0;i<(x->p)->st_ntemp;i++)
                        (x->p)->st_temp[i] = input_temp[i];

	// structure for restart
	if((x->p)->st_restart>-1)
	(x->p)->st_st_restart = AlloRestart(NTEMPMAX,(x->p)->st_nbin,NRESTARTS);

        // Start from highest temperature
         (x->p)->st_itemp = 0;

        // Reset counters
	 (x->p)->st_count_st = 0;
        (x->p)->st_count_adj = 0;
        (x->p)->st_count_print = 0;
        (x->p)->st_count_printt = 0;
        (x->p)->st_iter = 0;
        (x->p)->st_failure = 0;


	}
}
#endif
/********************************************************************
 Parse input file to obtain parameters (D=decimal, F=double, S=string, N=standalone)
 ********************************************************************/
void ReadParD(char *s, char key[20], int *par)
{
	int l,r;

	l = strlen(key);
	if (!strncmp(s,key,l))
	{
		r = sscanf(s,"%*s %d",par);
		if (r<1) { fprintf(stderr,"ERROR: Cannot read %s in parameter file (r=%d)\n",key,r); exit(1); }
		fprintf(stderr,"- %s  \t%d\n",key,*par);
	}

}

void ReadParLLU(char *s, char key[20], unsigned long long *par)
{
	int l,r;

	l = strlen(key);
	if (!strncmp(s,key,l))
	{
		r = sscanf(s,"%*s %llu",par);
		if (r<1) { fprintf(stderr,"ERROR: Cannot read %s in parameter file (r=%d)\n",key,r); exit(1); }
		fprintf(stderr,"- %s  \t%llu\n",key,*par);
	}

}


void ReadParL(char *s, char key[20], long *par)
{
	int l,r;

	l = strlen(key);
	if (!strncmp(s,key,l))
	{
		r = sscanf(s,"%*s %ld",par);
		if (r<1) { fprintf(stderr,"ERROR: Cannot read %s in parameter file (r=%d)\n",key,r); exit(1); }
		fprintf(stderr,"- %s  \t%ld\n",key,*par);
	}

}

void ReadParF(char *s, char key[20], double *par)
{
	int l,r;

	l = strlen(key);
	if (!strncmp(s,key,l))
	{
		r = sscanf(s,"%*s %lf",par);
		if (r<1) { fprintf(stderr,"ERROR: Cannot read %s in parameter file\n",key); exit(1); }
		fprintf(stderr,"- %s  \t%lf\n",key,*par);
	}
}

void ReadParS(char *s, char key[20], char *par)
{
	int l,r;

	l = strlen(key);
	if (!strncmp(s,key,l))
	{
		r = sscanf(s,"%*s %s",par);
		if (r<1) { fprintf(stderr,"ERROR: Cannot read %s in parameter file\n",key); exit(1); }
		fprintf(stderr,"- %s  \t%s\n",key,par);
	}
}

void ReadParN(char *s, char key[20], int *par)
{
	int l;

	l = strlen(key);
	if (!strncmp(s,key,l))
	{
		*par=1;
		fprintf(stderr,"- %s  \n",key);
	}
}

/********************************************************************
 Read data for potential from file
	 warning: if pairs is read before global, it overrides it
 ********************************************************************/
int ReadPotential(char *fname, struct s_potential *u, struct s_mc_parms *parms, int na, int ntype)
{
	int i,j,k=0,chapt=0,ab,dihtype,iaa,dih0;
	double e,r,r0,k0,d1,d3,d01,d03,kr_splice,pa,pb,sigma;//r0vec[ntype];
	char aux[500],keyword[100],ccc,c_ab,c_dihtype;
	FILE *fp;
	int c_d=0,c_a=0;
	u->splice = 0;
	u->g_r0hard = 0;
	u->g_dihe=-1.;
	u->dih_periodic = 0;
	u->dih_tabled = 0;
	u->splice=0;
	u->dih_ram=0;
	u->e_dihram=1;

	//for (j=0;j<ntype;j++) r0vec[j]=-1;

	fprintf(parms->flog,"Read potential from %s\n",fname);

	fp = fopen(fname,"r");
	if (!fp) Error("Cannot open potential file");

	while(fgets(aux,500,fp)!=NULL)
	{
	    if ( FindKeyword(aux,keyword)==1 && !strcmp(keyword,"pairs") ) chapt = 1;
	    if ( FindKeyword(aux,keyword)==1 && !strcmp(keyword,"global") ) chapt = 4;
	    if ( FindKeyword(aux,keyword)==1 && !strcmp(keyword,"angles") )
	    	{ 
			chapt = 2;
			if (parms->noangpot) 
			{ 
				fprintf(stderr,"WARNING: you defined noangpot in par file, but pot file contains angular potential (ignoring latter)\n"); 
				chapt=0; 
			}
		}
	
	    if ( FindKeyword(aux,keyword)==1 && !strcmp(keyword,"dihedrals") )
	    	{ 
			chapt = 3;
			if (parms->nodihpot) 
			{ 
				fprintf(stderr,"WARNING: you defined nodihpot in par file, but pot file contains dihedral potential (ignoring latter)\n"); 
				chapt=0; 
			}
		}
	    if ( FindKeyword(aux,keyword)==1 && !strcmp(keyword,"comments") ) chapt=5;
	    if ( FindKeyword(aux,keyword)==1 && !strcmp(keyword,"dihedral_table") )
		{
			chapt = 7;
			if (parms->nodihpot) 
			{ 
				fprintf(stderr,"WARNING: you defined nodihpot in par file, but pot file contains dihedral potential (ignoring latter)\n"); 
				chapt=0; 
			}
		}
	    if ( FindKeyword(aux,keyword)==1 && !strcmp(keyword,"dihedral_weight") )
	    	{ 
			chapt = 6;
			if (parms->nodihpot) 
			{ 
				fprintf(stderr,"WARNING: you defined nodihpot in par file, but pot file contains dihedral potential (ignoring latter)\n"); 
				chapt=0; 
			}
		}
	    if ( FindKeyword(aux,keyword)==1 && !strcmp(keyword,"hydrogen_bonds") ) chapt = 8;
 	    if ( FindKeyword(aux,keyword)==1 && !strcmp(keyword,"Ramachandran_Dihedrals") )
	    	{ 
			chapt = 9;
			if (parms->nodihpot) 
			{ 
				fprintf(stderr,"WARNING: you defined nodihpot in par file, but pot file contains dihedral potential (ignoring latter)\n"); 
				chapt=0; 
			}
		}
	    if ( FindKeyword(aux,keyword)==1 && !strcmp(keyword,"Alfa/Beta_propensity") )
	    	{ 
			chapt = 10;
			if (parms->nodihpot) 
			{ 
				fprintf(stderr,"WARNING: you defined nodihpot in par file, but pot file contains dihedral potential (ignoring latter)\n"); 
				chapt=0; 
			}
		}
	 
	    /////////////////////////
	    // read pairs interaction
	    /////////////////////////
	    if (chapt==1)
	    {
	    	if (sscanf(aux,"%d %d %lf %lf %lf",&i,&j,&e,&r,&r0)==5)
	    	{
				if (i>ntype-1) { fprintf(stderr,"* ERROR: found atomtype %d in ReadPotential where ntype=%d\n",i,ntype); exit(1); }
                                if (j>ntype-1) { fprintf(stderr,"* ERROR: found atomtype %d in ReadPotential where ntype=%d\n",j,ntype); exit(1); }
				(u->e)[i][j] = e;
				(u->r_2)[i][j] = r*r;
				(u->r0_2)[i][j] = r0*r0;
				(u->e)[j][i] = e;				// if the whole matrix is completely defined, these are useless
				(u->r_2)[j][i] = r*r;
				(u->r0_2)[j][i] = r0*r0;
				k++;
		}
	    }

	    /////////////////////////
	    // read angular potential
	    /////////////////////////
	    if (chapt==2)
	    {
		
	    	if (sscanf(aux,"%d %lf %lf",&i,&k0,&r0)==3)
	    	{

			if (i>na-1) { fprintf(stderr,"* ERROR: found atom %d in ReadPotential where natom=%d\n",i,na); exit(1); }
	    		(u->e_ang)[i] = k0;
	    		(u->ang0)[i] = r0;
	    	}
	    }

	    //////////////////////////
	    // read dihedral potential
	    //////////////////////////
	    if (chapt==3)
	    {
	    	if (sscanf(aux,"%d %lf %lf %lf %lf",&i,&d1,&d01,&d3,&d03)==5)
	    	{
			
				if (i>na-1) { fprintf(stderr,"* ERROR: found atom %d in ReadPotential where natom=%d\n",i,na); exit(1); }
	    		(u->e_dih1)[i] = d1;
	    		(u->dih01)[i] = d01;
	    		(u->e_dih3)[i] = d3;
	    		(u->dih03)[i] = d03;
	    	}
	    	u->dih_periodic = 1;
	    }

	    ///////////////////////////
	    // read global interactions
	    ///////////////////////////
	    if (chapt==4)
	    {
	    	// hardcore
	    	if (sscanf(aux,"hardcore %lf",&(u->g_r0hard))==1)
	    	{
	    		for (i=0;i<ntype;i++)
	    			for(j=0;j<ntype;j++)
	    				/*if (i != j )*/ (u->r0_2)[j][i] = u->g_r0hard * u->g_r0hard;
	    		fprintf(parms->flog,"Global hardcore = %lf\n",u->g_r0hard);
	    	}

	    	// imin
	    	if (sscanf(aux,"imin %d",&i)==1)
	    	{
	    		u->g_imin = i;
	    		fprintf(parms->flog,"Global imin = %d\n",u->g_imin);
	    	}

	    	// non go dihedrals (tabled dih)
	    	if (sscanf(aux,"e_dih %lf",&e)==1)		//strength of dih potential
	    	{
	    		u->g_dihe = e;
	    		fprintf(parms->flog,"Energy tabled dih = %lf\n",u->g_dihe);
	    	}

	    	if (sscanf(aux,"g_dihbin %lf",&e)==1)		//bin that determines usage of dih table
	    	{
	    		u->g_dihbin = e;
	    		fprintf(parms->flog,"Bin tabled dih = %lf\n",u->g_dihbin);
			if(360./u->g_dihbin > (double)NDIHFMAX)Error("NDIHFMAX too small ");
	    	}

	    	// attractive well
	    	if (sscanf(aux,"homopolymeric %lf %lf",&(u->g_ehomo),&(u->g_rhomo))==1)
	    	{
	    		for (i=0;i<ntype;i++)
	    			for(j=0;j<ntype;j++)
	    				if (i != j )
	    				{
								(u->r_2)[j][i] = u->g_rhomo*u->g_rhomo;
								(u->e)[j][i] = u->g_ehomo;
	    				}
	    		fprintf(parms->flog,"Global homopolymeric interaction = %lf (r=%lf)\n",u->g_ehomo,u->g_rhomo);
	    	}

	    	// angular potential
	    	if (sscanf(aux,"angle %lf %lf",&(u->g_anglek),&(u->g_angle0))==2)
	    	{
	    	    for (i=0;i<na;i++)
	    	    {
	    	    	(u->e_ang)[i] = u->g_anglek;
	    	    	(u->ang0)[i] = u->g_angle0;
	    	    }
	    	    fprintf(parms->flog,"Global angular potential: k=%lf ang0=%lf\n",u->g_anglek,u->g_angle0);
	    	}

	    	// dihedral potential
	    	if (sscanf(aux,"dihedral1 %lf %lf",&(u->g_dihk1),&(u->g_dih01))==2)
	    	{
	    		for (i=0;i<na;i++)
	    		{
	    			(u->e_dih1)[i] = u->g_dihk1;
	    	    	(u->dih01)[i] = u->g_dih01;
	    		}
	    		fprintf(parms->flog,"Global dihedral potential: k=%lf dih0=%lf n=1\n",u->g_dihk1,u->g_dih01);
	    	}
	    	if (sscanf(aux,"dihedral3 %lf %lf",&(u->g_dihk3),&(u->g_dih03))==2)
	    	{
	    		for (i=0;i<na;i++)
	    		{
	    			(u->e_dih3)[i] = u->g_dihk3;
	    			(u->dih03)[i] = u->g_dih03;
	    		}
	    		fprintf(parms->flog,"Global dihedral potential: k=%lf dih0=%lf n=3\n",u->g_dihk3,u->g_dih03);
	    	}

			// ramachandran dihedral
			if (sscanf(aux,"e_dihram %lf",&e)==1)
	    	{
				u->e_dihram = e;	    	
					fprintf(parms->flog,"Global energy ramachandran dihedral: e_dihram=%lf\n",u->e_dihram);
	    	}

	    	// splice well
	    	if (sscanf(aux,"splice %lf %lf",&kr_splice,&(u->ke_splice))==2)
	    		{
	    			u->splice = 1;
		    		u->kr2_splice = kr_splice * kr_splice;
		    		fprintf(parms->flog,"Splice energy well into two parts (kr=%lf ke=%lf)\n",kr_splice,u->ke_splice);
		    		fprintf(stderr,"Splice energy well into two parts (kr=%lf ke=%lf)\n",kr_splice,u->ke_splice);
		    		if (kr_splice>=1) Error("kr_splice > 1");
		    		if (u->ke_splice>=1) fprintf(stderr,"WARNING: ke_splice > 1 sounds strange....\n");
	    		}

	    	// box
	    	if (sscanf(aux,"boxtype %c",&(u->boxtype))) fprintf(stderr,"Molecules in a box of type %c\n",u->boxtype);
	    	if (sscanf(aux,"boxsize %lf",&(u->boxsize)))
		{
			 fprintf(stderr,"Size of the box is %lf",u->boxsize);
			//set boxsize=boxsize/2 (see EnergyBox in Potential.c)
	 		u->boxsize=(u->boxsize/2.);
	
		}
	    }


	    ///////////////
	    //read comments
	    ///////////////
	    if(chapt==5)
	    {
	    	fprintf(stderr,"%s",aux);
	    }

	    ///////////////////////
	    // read dihedral tables
	    ///////////////////////
	    if (chapt==6)
	    {
	    	if (sscanf(aux,"%d %lf %lf %d",&i,&d1,&d3,&j)==4)   //last int is 0 if phi (ia CA) and 1 if psi (ia C)
	    	{
	    				if (i>na-1) { fprintf(stderr,"* ERROR: found atom %d in ReadPotential where natom=%d\n",i,na); exit(1); }
	    	    		(u->dih_pa)[i] = d1;
	    	    		(u->dih_pb)[i] = d3;
	    	    		(u->dih_which)[i] = j;
	    	}
	    	u->dih_tabled = 1;
	    }

	    if (chapt==7)
	    {
	    	if (sscanf(aux,"%d %lf %lf %lf %lf",&i,&d1,&d3,&d01,&d03)==5)		//phi alpha..
	    	{
				if (i>=NDIHFMAX) Error("NDIHFMAX too small when reading dihedral tables ");
	    	    		(u->dih_f_phi_a)[i] = d1;
	    	    		(u->dih_f_phi_b)[i] = d3;
	    	    		(u->dih_f_psi_a)[i] = d01;
	    	    		(u->dih_f_psi_b)[i] = d03;

	    	}
	    }

	    ///////////////////////
	    // hydrogen bonds
	    ///////////////////////
	    if (chapt==8)
	    {
	    	if (sscanf(aux,"%d %c %d",&i,&ccc,&j)==3)
	    	{
	    		if (ccc=='d') { (u->hb)[i]=1; c_d++; }
	    		if (ccc=='a') { (u->hb)[i]=2; c_a++; }
	    		u->hb_iam[i] = j;
	    	}
	    }
	    
	    /////////////////////////////
	    // read ramachandran dihedral
	    /////////////////////////////
	    if (chapt==9)
	    {
            u->dih_ram = 1;
	    	if (sscanf(aux,"%c %c %lf %d",&c_ab,&c_dihtype,&sigma,&dih0)==4)
	    	{
                //int is 0  if ALPHA and 1 if BETA
                if (c_ab == 'a' ) ab = 0;
                else {
                    if (c_ab == 'b' ) ab = 1;
                    else Error("Wrong name for secondary structure in dihedral potential");
                }
                
                if (c_dihtype == 'f' ) dihtype = 0;
                else {
                    if (c_dihtype == 'p' ) dihtype = 1;
                    else Error("Wrong name for dihedral name in dihedral potential");
                }
                (u->sigma)[ab][dihtype] = sigma;
                (u->dih0)[ab][dihtype] = dih0;
	       	}
	    }
        
	    if (chapt==10)
	    {
            u->dih_ram = 1;
	    	if (sscanf(aux,"%d %lf %lf",&iaa,&pa,&pb)==3)
	    	{
                (u->ab_propensity)[0][iaa] = pa;
				(u->ab_propensity)[1][iaa] = pb;
	    	}
	    }
        
	}
	fclose(fp);


	if (k==0) fprintf(parms->flog,"WARNING: no pair interaction read in potential file.\n");
	else fprintf(parms->flog,"Read %d pair interaction in potential file\n",k);

	return k;
}

/********************************************************************
 Print single elements of the structures
 ********************************************************************/

void PrintVector(FILE *fp, char c[], struct vector x)
{
	fprintf(fp,"%s = ( %lf %lf %lf )\n",c,x.x,x.y,x.z);
}

void PrintAngles(FILE *fp, char c[], struct angles x)
{
	fprintf(fp,"%s = ( %lf %lf %lf )\n",c,x.ang,x.dih,x.r);
}

/********************************************************************
 Print structures (for debugging purposes)
 ********************************************************************/
void PrintStructure(struct s_polymer *p, int npol, FILE *fp, int shell)
{
	int i,j,k;//chk=0;

	for (i=0;i<npol;i++)
	{
		fprintf(fp,"Polymer %d\n",i);
		for (j=0;j<(p+i)->nback;j++)
		{
			fprintf(fp,"\tch=%d back=%d%s at=%d%s\n",i,j,(((p+i)->back)+j)->aa,(((p+i)->back)+j)->ia,(((p+i)->back)+j)->type);
			fprintf(fp,"\t  Contacts: (%d)\t",(((p+i)->back)+j)->ncontacts);
			for (k=0;k<(((p+i)->back)+j)->ncontacts;k++)
			{
				//chk=0;
				fprintf(fp,"%d(ch=%d e=%lf) ",*(((((p+i)->back)+j)->contacts)+k),*(((((p+i)->back)+j)->contacts_p)+k),*(((((p+i)->back)+j)->e)+k) );
				//for (l=0;l<(((p+ *(((((p+i)->back)+j)->contacts_p)+k) )->back)+ *(((((p+i)->back)+j)->contacts)+k) )->ncontacts;l++)
				// if ( *(((((p+ *(((((p+i)->back)+j)->contacts_p)+k) )->back)+ *(((((p+i)->back)+j)->contacts)+k) )->contacts)+l) == j) chk=1;
				//if (chk==0) fprintf(stderr,"\nWARNING: asymmetry between residues %d and %d\n",j,*(((((p+i)->back)+j)->contacts)+k));
			}
			fprintf(fp,"\n");
			if (shell>0)
			{
				fprintf(fp,"\t  Shell: ");
				for (k=0;k<(((p+i)->back)+j)->nshell;k++) fprintf(fp,"%d(ch=%d) ",*(((((p+i)->back)+j)->shell)+k),
									*(((((p+i)->back)+j)->shell_p)+k) );
				fprintf(fp,"\n");
			}
			if ((((p+i)->back)+j)->nside>0)
			{
				fprintf(fp,"\t  Sidechains (rot=%d): ",(((p+i)->back)+j)->irot);
				for (k=0;k<(((p+i)->back)+j)->nside;k++) fprintf(fp,"%d%s ",(((((p+i)->back)+j)->side)+k)->ia,
						(((((p+i)->back)+j)->side)+k)->type );
				fprintf(fp,"\n");
			}


		}
		fprintf(fp,"\tEtot=%lf\n",(p+i)->etot);
	}
}


/********************************************************************
 Print potential file
 ********************************************************************/
void PrintPotential(struct s_potential *u, char *eoutfile, int nat, int ntypes, int noangpot, int nodihpot, int hb)
{
	int i,j;
	char ccc=' ';
	FILE *fout;

	// open file
	if (strcmp(eoutfile,"stdout"))
	{
		fout = fopen(eoutfile,"w");
		if (!fout) Error("Cannot open file to write potential");
	}
	else fout = stdout;

	// print globals
	if (u->g_r0hard>0 || u->g_rhomo>0 || u->g_anglek || u->g_dihk1 || u->g_dihk3 || (u->splice>0) || u->g_imin || u->dih_ram)
	{
		fprintf(fout,"[ global ]\n");
		if (u->g_r0hard>0) fprintf(fout,"hardcore %lf\n",u->g_r0hard);
		if (u->g_rhomo>0) fprintf(fout,"homopolymeric %lf %lf\n",u->g_ehomo,u->g_rhomo);
		if (u->g_anglek>0) fprintf(fout,"angle %lf %lf\n",u->g_anglek,u->g_angle0);
		if (u->g_dihk1>0) fprintf(fout,"dihedral1 %lf %lf\n",u->g_dihk1,u->g_dih01);
		if (u->g_dihk3>0) fprintf(fout,"dihedral3 %lf %lf\n",u->g_dihk3,u->g_dih03);
		if (u->splice>0) fprintf(fout,"splice %lf %lf \n",sqrt(u->kr2_splice),u->ke_splice);
		if (u->g_imin) fprintf(fout,"imin %d\n",u->g_imin);
		if (u->dih_ram) fprintf(fout,"e_dihram %lf\n",u->e_dihram);
		// non go dihedrals (tabled dih)
		if(u->dih_tabled)fprintf(fout,"e_dih %lf\n",u->g_dihe);		//strength of dih potential
		if(u->dih_tabled)fprintf(fout,"g_dihbin %lf\n",u->g_dihbin);//bin that determines usage of dih table
		if (u->boxtype !='n')
		{
			fprintf(fout,"boxtype %c\n",u->boxtype);
			fprintf(fout,"boxsize %lf\n",u->boxsize);
		}


		fprintf(fout,"\n");
	}

	// print (in case) hardcores
	if (u->hc_number)
	{
		fprintf(fout,"[ hardcores ]\n");

		for (i=0;i<u->hc_number;i++)
			fprintf(fout,"%d\t%lf\n",u->hc_type[i],u->hc_r0[i]);
	}

	// print pairs
	fprintf(fout,"[ pairs ]\n");

	for (i=0;i<ntypes;i++)
		for (j=i;j<ntypes;j++)
			if (u->e[i][j]<-EPSILON || u->e[i][j]>EPSILON)
			{
				fprintf(fout,"%3d\t%3d\t\t%lf\t%lf\t%lf\n",i,j,u->e[i][j],sqrt(u->r_2[i][j]),sqrt(u->r0_2[i][j]));
			}

	// angles
	if (!noangpot)
	{
		fprintf(fout,"[ angles ]\n");

		for (i=1;i<nat;i++)
			if (u->e_ang[i]<-EPSILON || u->e_ang[i]>EPSILON)
				fprintf(fout,"%d %lf %lf\n",i,u->e_ang[i],u->ang0	[i]);
	}

	if (!nodihpot)
	{
		fprintf(fout,"[ dihedrals ]\n");

		for (i=1;i<nat;i++)
			if (u->e_dih1[i]<-EPSILON || u->e_dih1[i]>EPSILON || u->e_dih3[i]<-EPSILON || u->e_dih3[i]>EPSILON)
				fprintf(fout,"%d %lf %lf %lf %lf\n",i,u->e_dih1[i],u->dih01[i],u->e_dih3[i],u->dih03[i]);
	
		if(u->dih_ram!=0)
		{
			fprintf(fout,"[ Ramachandran_Dihedrals ]\n");
            fprintf(fout,"a/b\tpsi/phi\tstdev\tpsi0/phi0\n");
			//stampo: 1. alfa=0/beta=1
			//	  2. phi=0/psi=1
			//	  3. sigma
			//	  4. dih0
			for(i=0; i<2; i++) {
				for(j=0;j<2;j++) {
                    if (i == 0) {
                        if (j == 0) fprintf(fout,"a\tf\t%f\t%d\n",u->sigma[i][j],u->dih0[i][j]);
                        else fprintf(fout,"a\tp\t%f\t%d\n",u->sigma[i][j],u->dih0[i][j]);
                    }
                    else{
                        if (j == 0) fprintf(fout,"b\tf\t%f\t%d\n",u->sigma[i][j],u->dih0[i][j]);
                        else fprintf(fout,"b\tp\t%f\t%d\n",u->sigma[i][j],u->dih0[i][j]);
                    }
                }
            }
			fprintf(fout,"[ Alfa/Beta_propensity ]\n");
            fprintf(fout,"ia\tp_a\tp_b\n");
			//stampo:	iaa	prop_a[iaa]	prop_b[iaa]
			for(i=0;i<NAAMAX;i++)
				if(u->ab_propensity[0][i]<2)
					fprintf(fout,"%d\t%lf\t%lf\n",i,u->ab_propensity[0][i],u->ab_propensity[1][i]);
		}
        
	}


	if (hb)
	{
	fprintf(fout,"[ hydrogen_bonds ]\n");

		for (i=0;i<nat;i++)
		{
			if ((u->hb)[i]==1) ccc='d';
			if ((u->hb)[i]==2) ccc='a';
			if ((u->hb)[i]!=0)  fprintf(fout,"%d \t%c\t%d\n",i,ccc,u->hb_iam[i]);
		}
	}


	if (strcmp(eoutfile,"stdout")) fclose(fout);


}







void	ComputeIAPDB(struct s_polymer *p,struct s_mc_parms *parms){
	int ipol;
	int i;
	for(ipol=0;ipol<parms->npol;ipol++)
	{
	//	fprintf(stderr,"IADB: setting up polymer n. %d\n",ipol);
		for(i=0;i<(p+ipol)->nback;i++)
		{	
			if( !strcmp( (((p+ipol)->back)+i)->type,"N") )
                        (((p+ipol)->back)+i)->iapdb=0;
			if( !strcmp( (((p+ipol)->back)+i)->type,"CA") )
			(((p+ipol)->back)+i)->iapdb=1;
			if( !strcmp( (((p+ipol)->back)+i)->type,"C") )
                        (((p+ipol)->back)+i)->iapdb=2;
		}
		for(i=0;i<(p+ipol)->nback;i++)
                {
          //              fprintf(stderr,"%d\t%d\n",i,(((p+ipol)->back)+i)->iapdb);
                }

	}










}

