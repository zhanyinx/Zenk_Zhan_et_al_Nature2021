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
 * mc.c
 *
 *  Created on: Sep 20, 2010
 *      Author: guido
 */



#include "montegrappa.h"
/*****************************************************************************
 Monte Carlo loop
 *****************************************************************************/
int softexit;

void Do_MC(struct s_polymer *p, struct s_polymer *fragment, struct s_polymer *replica, struct s_polymer *native, struct s_potential *pot, struct s_mc_parms *parms, FILE *ftrj, FILE *fe, struct s_polymer *oldp, FILE *fproc,int irun, struct s_mpi_parms *mpiparms)
{
	int i,ok=0,iprinttrj=0,iprintlog=0,iprinte=0,ntm=0,ci,mcount[NMOVES],macc[NMOVES],mdone[NMOVES];
//	int anneal_count=0, anneal_status=0;

	double t;
	unsigned long long istep=0;               
		
	if(parms->shell==1)
		parms->ishell=0;

 
	#ifdef OPTIMIZEPOT
	 struct s_optimizepot *op = p->op;
	 op->icount=0;
	 op->nframes = 0;
	 op->ncontacts = 0;
	 op->record = 0;
	#endif


//ASTEMPERING
	#ifdef STEMPERING
//	struct st_stuff *st_p;
	int st_iprint=0;
	int st_o=0;
//	if (parms->stempering) st_p = parms->p;
	int nconf;
	#endif

	

        

        #ifdef ACTIVE_MPI
	int my_rank=mpiparms->my_rank;
	MPI_Datatype Backtype=mpiparms->Back_mpi;
        MPI_Datatype Sidetype=mpiparms->Side_mpi;
        MPI_Datatype Rottype=mpiparms->Rot_mpi;
        MPI_Status astatus=mpiparms->astatus;
	
	int ptempering_count=0;
	int ex_count[(parms->ntemp)-1],ex_acc[(parms->ntemp)-1];
        t=parms->T[my_rank];
        #else
	int my_rank=0;
	t=parms->T;
	#endif
	softexit =0;

/*	#ifdef ACTIVE_MPI
	if(my_rank==0)
	for(i=0;i<4;i++)
	{
		fprintf(stderr,"rank %d\t temperature = %lf\n",i,parms->T[i]);

	}
	fprintf(stderr,"rank %d\t my temp is  = %lf\n",my_rank,t);
	#endif
*/
	// print debug info
	if (parms->debug>0)
		fprintf(stderr,"Starting step loop....\n");
	if (parms->debug>1)
	{
		fprintf(stderr,"Moves = %d %d %d %d\n",parms->movetype[0],parms->movetype[1],
				parms->movetype[2],parms->movetype[3]);
		for (ci=0;ci<parms->npol;ci++)
			for (i=0;i<(p+ci)->nback;i++) ntm += (((p+ci)->back)+i)->move;
		fprintf(stderr,"Backbone atoms which can be moved = %d\n",ntm);
	}
	
	CopyAllPolymers(p,oldp,parms->npol,parms->nosidechains,parms->nosidechains);
	parms->acc = 0;
	parms->mov = 0;
	
	for(i=0;i<NMOVES;i++)
	{
		mcount[i]=0;
		macc[i]=0;
		mdone[i]=0;
	}

	     #ifdef ACTIVE_MPI
      //for (i=0;i<(parms->ntemp)-1;i++) { ex_count[i]=0; ex_acc[i]=0;}
      #endif
  
#ifdef OPTIMIZEPOT
               if(my_rank==0)
                {
                        if (strcmp(parms->op_minim,"none")  && ok>-1){

			#ifdef STEMPERING
                        if (parms->stempering){
                        if ((parms->p)->st_ttarget_harvest==1 && st_iprint>=(parms->p)->st_printpdb && t <= (parms->p)->st_ttarget + EPSILON && (parms->p)->st_nm==1)
                           {

                                        
                                       // #ifdef ACTIVE_MPI

                                       // fprintf(stderr,"> recording polymer %d/%llu\tE=%lf\tT=%lf\n",op->nframes,
                                       //        (parms->nstep-parms->op_wait)/parms->op_deltat,p->etot,parms->T[my_rank]);
                                       // #else
                                       // fprintf(stderr,"> recording polymer %d/%llu\tE=%lf\tT=%lf\n",op->nframes,
                                       //       (parms->nstep-parms->op_wait)/parms->op_deltat,p->etot,parms->T);
                                       // #endif
                                        

                                        op->record = 1;         // fill it,it2,mul
                                        TotalEnergy(p,pot,parms,parms->npol,0,parms->nosidechains,parms->debug,my_rank);      // record energy
                                        #ifdef OP_DEBUG
                                        FILE *fp;
                                        fp = fopen("op_log","a");
                                        for (i=0;i<p->op->ncontacts;i++)
                                                fprintf(fp,"contacts: snap=%5d  %d-%d\tm=%lf\n",op->nframes,
                                        p->op->it1[i],p->op->it2[i],p->op->mul[op->nframes][i]);
                                        fclose(fp);
                                        #endif
                                        op->record = 0;
                                        OP_GetRestrain(op->nframes,p,0,parms->op_input);            // record restrains
                                        #ifdef ACTIVE_MPI       
                                        op->t[op->nframes] = parms->T[my_rank];
                                        #else
                                        op->t[op->nframes] = parms->T;
                                        #endif

                         // record temperature
                            		op->nframes ++;

					 st_iprint = 0;

                                       
                          }                               

					if ((parms->p)->st_ttarget_harvest>0 && ok>-1 && t <= (parms->p)->st_ttarget + EPSILON ) st_iprint ++;

                                         if (st_o == 1) t = (parms->p)->st_temp[(parms->p)->st_itemp];

                                        }

                                        
                                else if (op->icount == parms->op_deltat)
                                
                                #else
				if (op->icount == parms->op_deltat)
                                
        
				#endif
				{
                                        #ifdef ACTIVE_MPI
                                        
                                        fprintf(stderr,"> recording polymer %d/%llu\tE=%lf\tT=%lf\n",op->nframes,
                                               (parms->nstep-parms->op_wait)/parms->op_deltat,p->etot,parms->T[my_rank]);
                                        #else
                                        fprintf(stderr,"> recording polymer %d/%llu\tE=%lf\tT=%lf\n",op->nframes,
                                               (parms->nstep-parms->op_wait)/parms->op_deltat,p->etot,parms->T);
                                        #endif
                                        

                                        op->record = 1;         // fill it,it2,mul
                                        TotalEnergy(p,pot,parms,parms->npol,0,parms->nosidechains,parms->debug,my_rank);      // record energy
                                        #ifdef OP_DEBUG
                                        FILE *fp;
                                        fp = fopen("op_log","a");
                                        for (i=0;i<p->op->ncontacts;i++)
                                                fprintf(fp,"contacts: snap=%5d  %d-%d\tm=%lf\n",op->nframes,
                                        p->op->it1[i],p->op->it2[i],p->op->mul[op->nframes][i]);
                                        fclose(fp);
                                        #endif
                                        op->record = 0;
                                        OP_GetRestrain(op->nframes,p,0,parms->op_input);            // record restrains
                                        #ifdef ACTIVE_MPI       
                                        op->t[op->nframes] = parms->T[my_rank];
                                        #else
                                        op->t[op->nframes] = parms->T;
                                        #endif

                         // record temperature
                            op->nframes ++;
                                        op->icount=0;
                                }
                                if (istep >= parms->op_wait)
                                        op->icount ++;

                }
                }
                #endif


	////////////////////////////////////////
	// MC LOOP                            //
	/////////////////////////////////////////
	do
	{
	

		#ifdef DEBUG
		if (parms->debug>1)
			fprintf(stderr,"-----------------\nStep %llu\n",istep);
		#endif



		// make a move (one per allowed type)
		if (mcount[0] == parms->movetype[0])						// flip of backbone
		{
			ok = MoveBackboneFlip(p,oldp,parms,pot,istep,parms->debug,t);
			if (ok>-1) mcount[0] = 0;
			if (ok==1) macc[0]++;
			mdone[0]++;
		}
		
		if (mcount[1] == parms->movetype[1])						// pivot of backbone
		{

			ok = MoveBackbonePivot(p,oldp,pot,parms,t);
			if (ok>-1) mcount[1] = 0;
			if (ok==1) macc[1]++;
			mdone[1]++;
		}
	
		if (mcount[2] == parms->movetype[2])						// multiple pivot of backbone
		{
			ok = MoveMultiplePivot(p,oldp,pot,parms->nmul_mpivot,parms,t);
			if (ok>-1) mcount[2] = 0;
			if (ok==1) macc[2]++;
			mdone[2]++;
		}
		if (mcount[3] == parms->movetype[3])						// move of sidechain
		{

			ok = MoveSidechain(p,oldp,parms,pot,istep,parms->debug,t);
			if (ok>-1) mcount[3] = 0;
			if (ok==1) macc[3]++;
			mdone[3]++;
		}
	
		if (mcount[4] == parms->movetype[4])						// loose pivot
		{
			ok = MoveLoosePivot(p,oldp,pot,parms->nmul_lpivot,parms,t);
			if (ok>-1) mcount[4] = 0;
			if (ok==1) macc[4]++;
			mdone[4]++;
		}
	
		if (mcount[5] == parms->movetype[5])						// multiple flip
		{
			ok = MoveMultipleFlip(p,oldp,parms,pot,istep,parms->debug,t);
			if (ok>-1) mcount[5] = 0;
			if (ok==1) macc[5]++;
			mdone[5]++;
		}
	
		if (mcount[6] == parms->movetype[6])						// move of center of mass
		{
			ok = MoveCoM(p,oldp,parms,pot,istep,parms->debug,t);
			if (ok>-1) mcount[6] = 0;
			if (ok==1) macc[6]++;
			mdone[6]++;
		}

		if (mcount[7] == parms->movetype[7])
                {


                        ok=LocalMove(p,oldp,fragment,pot,parms->nmul_local,parms,t);
			if(ok>-1) mcount[7]=0;
                        if(ok==1) macc[7]++;
                        mdone[7]++;
		}

		if (mcount[8] == parms->movetype[8])
                {
                        ok=MoveRotation(p,oldp,parms,pot,istep,parms->debug,t);
                        if(ok>-1) mcount[8]=0;
                        if(ok==1) macc[8]++;
                        mdone[8]++;
                }
				
		// if you want to anneal
		//if (parms->anneal) Anneal(parms,&t,&anneal_count,&anneal_status,&ok,&ishell,mcount);


		if (mcount[9] == parms->movetype[9])
                {
///                        ok=MoveRotation(p,oldp,parms,pot,istep,parms->debug,t);
		       ok=MoveClusterCoM(p,oldp,parms,pot,istep,parms->debug,t,fproc,my_rank); 

                       if(ok>-1) mcount[9]=0;
                        if(ok==1) macc[9]++;
                        mdone[9]++;
                }

		if (mcount[10] == parms->movetype[10])
                {
//                        ok=MoveRotation(p,oldp,parms,pot,istep,parms->debug,t);
			ok=MoveClusterRot(p,oldp,parms,pot,istep,parms->npol,parms->debug,t,fproc,my_rank);

                        if(ok>-1) mcount[10]=0;
                        if(ok==1) macc[10]++;
                        mdone[10]++;
                }




/************************************************\
 * 	print di varie ed eventuali		*
\************************************************/

		// print log
		if (iprintlog == parms->nprintlog && ok>-2)
		{
			if(my_rank==0)		/* stampa solo per rank=0 */
				fprintf(parms->flog,"Step = %llu\tE=%lf\tok=%d\tT=%lf\n",istep,p->etot,ok,t);
//			fprintf(stderr,"Step = %d\tE=%lf\t rank %d\ttemp %lf\n",istep,p->etot,my_rank,t);
//			fprintf(fproc,"%d Etot=%lf\tEtot(true)=%lf ok=%d\n",istep,p->etot,
//                              TotalEnergy(p,pot,parms,parms->npol,0,parms->nosidechains,parms->debug,my_rank),ok);
			PrintEnergies_Parallel(stderr,parms->npol,istep,p,my_rank);

			#ifdef DEBUG
			fprintf(stderr,"%llu Etot=%lf\tEtot(true)=%lf ok=%d\n",istep,p->etot,
			TotalEnergy(p,pot,parms,parms->npol,0,parms->nosidechains,parms->debug,my_rank),ok);	// energy update in Move...
			
			if (parms->debug>4 && ok==1)
				PrintStructure(p,parms->npol,stderr,parms->shell);
			#endif
			iprintlog=0;
		}

		// print trajectory
		if (ftrj != NULL && iprinttrj == parms->nprinttrj && ok>-1)
		{
//DEBUG			#ifdef ACTIVE_MPI
//			fprintf(stderr,"printing trajectory MPI: rank %d, nback = %d\n",my_rank,(p)->nback);
//			#endif
			CountContacts(stderr,p,parms,parms->mov);			

			fprintf(stderr,"DEBUG: ok=%d\tTotal Energy %lf\n",ok,TotalEnergy(p,pot,parms,parms->npol,0,parms->nosidechains,parms->debug,my_rank));
			fprintf(stderr,"DEBUG: \t\t\tetot=%lf\n",p->etot);
			sprintf(p->title,"step %llu\tE=%lf",istep,p->etot);
			PrintPDBStream(p,parms->npol,ftrj);
			fflush(ftrj);
			iprinttrj = 0;
			
			

			//fprintf(stderr,"Total Energy:\t%lf\n\n",TotalEnergy(p,pot,parms,parms->npol,0,parms->nosidechains,parms->debug,0));
		}

		// print energies
		if (fe != NULL && iprinte == parms->nprinte && ok>-1)
		{
			PrintEnergies(fe,parms->npol,istep,p);
			fflush(fe);
			iprinte = 0;
		}
//		fprintf(stderr,"ETOT\t%lf\n",p->etot);
		// update shell
		if (parms->shell==1)
		{//DEBUG SHELL
			if (parms->ishell == parms->nshell && ok>-1)
			{
				UpdateShell(p,parms);
				UpdateShell(oldp,parms);					
				parms->ishell=0;
			}
		}		
//ASTEMPERING 	
		#ifdef STEMPERING	
		if (parms->stempering)
                         {
                                 st_o = STempering(p->etot,istep,(parms->p));

				if(!strcmp(parms->op_minim,"none")){
                                if ((parms->p)->st_ttarget_harvest==1 && st_iprint>=(parms->p)->st_printpdb && t <= (parms->p)->st_ttarget + EPSILON && (parms->p)->st_nm==1)
                                 {
                                     sprintf(p->title,"step %llu\tE=%lf",istep,p->etot);
                                        PrintPDBStream(p,parms->npol,(parms->p)->st_pdbf);
                                        nconf++;

                                        if(parms->nconf!=-1)
                                        if(nconf>parms->nconf-1){
                                        fprintf(stderr,"Reached target temperature (Ttarget=%lf).\nENDING SIMULATION\n",(parms->p)->st_ttarget);

					
                                        exit(0);
					
                                        }
                                        

                                        st_iprint = 0;
                                 }
                                 if ((parms->p)->st_ttarget_harvest>0 && ok>-1 && t <= (parms->p)->st_ttarget + EPSILON ) st_iprint ++;

                                 if (st_o == 1) t = (parms->p)->st_temp[(parms->p)->st_itemp];
                         }

		}
          #endif



  #ifdef OPTIMIZEPOT
               if(my_rank==0)
                {

                       if (strcmp(parms->op_minim,"none")  && ok>-1)
			{
			
			#ifdef STEMPERING
                        if (parms->stempering){
                        if ((parms->p)->st_ttarget_harvest==1 && st_iprint>=(parms->p)->st_printpdb && t <= (parms->p)->st_ttarget + EPSILON && (parms->p)->st_nm==1)
                           {

                                        
                                        //#ifdef ACTIVE_MPI

                                        //fprintf(stderr,"> recording polymer %d/%llu\tE=%lf\tT=%lf\n",op->nframes,
                                        //       (parms->nstep-parms->op_wait)/parms->op_deltat,p->etot,parms->T[my_rank]);
                                        //#else
                                        //fprintf(stderr,"> recording polymer %d/%llu\tE=%lf\tT=%lf\n",op->nframes,
                                        //       (parms->nstep-parms->op_wait)/parms->op_deltat,p->etot,parms->T);
                                        //#endif
                                        

                                        op->record = 1;         // fill it,it2,mul
                                        TotalEnergy(p,pot,parms,parms->npol,0,parms->nosidechains,parms->debug,my_rank);      // record energy
                                        #ifdef OP_DEBUG
                                        FILE *fp;
                                        fp = fopen("op_log","a");
                                        for (i=0;i<p->op->ncontacts;i++)
                                                fprintf(fp,"contacts: snap=%5d  %d-%d\tm=%lf\n",op->nframes,
                                        p->op->it1[i],p->op->it2[i],p->op->mul[op->nframes][i]);
                                        fclose(fp);
                                        #endif
                                        op->record = 0;
                                        OP_GetRestrain(op->nframes,p,0,parms->op_input);            // record restrains
                                        #ifdef ACTIVE_MPI       
                                        op->t[op->nframes] = parms->T[my_rank];
                                        #else
                                        op->t[op->nframes] = parms->T;
                                        #endif

                         // record temperature
                            		op->nframes ++;

					 st_iprint = 0;

                                       
                                                        }                               

					if ((parms->p)->st_ttarget_harvest>0 && ok>-1 && t <= (parms->p)->st_ttarget + EPSILON ) st_iprint ++;

                                         if (st_o == 1) t = (parms->p)->st_temp[(parms->p)->st_itemp];

                                        }

                                        
                                else if (op->icount == parms->op_deltat)
                                
				#else
				if (op->icount == parms->op_deltat)
                                

				#endif
                                {        
                                        #ifdef ACTIVE_MPI
                                        
                                        fprintf(stderr,"> recording polymer %d/%llu\tE=%lf\tT=%lf\n",op->nframes,
                                               (parms->nstep-parms->op_wait)/parms->op_deltat,p->etot,parms->T[my_rank]);
                                        #else
                                        fprintf(stderr,"> recording polymer %d/%llu\tE=%lf\tT=%lf\n",op->nframes,
                                               (parms->nstep-parms->op_wait)/parms->op_deltat,p->etot,parms->T);
                                        #endif
                                        

                                        op->record = 1;         // fill it,it2,mul
                                        TotalEnergy(p,pot,parms,parms->npol,0,parms->nosidechains,parms->debug,my_rank);      // record energy
                                        #ifdef OP_DEBUG
                                        FILE *fp;
                                        fp = fopen("op_log","a");
                                        for (i=0;i<p->op->ncontacts;i++)
                                                fprintf(fp,"contacts: snap=%5d  %d-%d\tm=%lf\n",op->nframes,
                                        p->op->it1[i],p->op->it2[i],p->op->mul[op->nframes][i]);
                                        fclose(fp);
                                        #endif
                                        op->record = 0;
                                        OP_GetRestrain(op->nframes,p,0,parms->op_input);            // record restrains
                                        #ifdef ACTIVE_MPI       
                                        op->t[op->nframes] = parms->T[my_rank];
                                        #else
                                        op->t[op->nframes] = parms->T;
                                        #endif

                         // record temperature
                            op->nframes ++;
                                        op->icount=0;
                                }
                                if (istep >= parms->op_wait)
                                        op->icount ++;

                	}
                }

                #endif





                 #ifdef ACTIVE_MPI
                 if(ptempering_count==(parms->nstep_exchange)/2 && parms->ntemp>1)
                 {
			
                 	//Process syncing before exchange

                        MPI_Barrier(MPI_COMM_WORLD);
			
                        if(parms->debug>2) 
				if(my_rank==0)
					fprintf(stderr,"\nstep=%llu\n",istep);
                	ExchangePol(p,replica,oldp,parms,pot,my_rank,parms->ntemp,0,ex_count,ex_acc,Backtype,Sidetype,Rottype,astatus,istep);	
				

		}	
                
		if(ptempering_count==(parms->nstep_exchange) && parms->ntemp>1)
                {

		
                	MPI_Barrier(MPI_COMM_WORLD);
                        if(parms->debug>2) if(my_rank==0)
				fprintf(stderr,"\nstep=%llu\n",istep); 
                        ExchangePol(p,replica,oldp,parms,pot,my_rank,parms->ntemp,1,ex_count,ex_acc,Backtype,Sidetype,Rottype,astatus,istep);
			
			ptempering_count=0;                        
                }
                
		ptempering_count++;
		#endif

		
                




		// advance counters
                for (i=0;i<NMOVES;i++){ 
			mcount[i]++;
			if(mcount[i]>200000000) mcount[i]=0; //reset in order to prevent too big number
			}
                if (ok>-1)
                {
                	istep ++;
                        iprinttrj ++;
                        iprintlog ++;
                        iprinte ++;
                        parms->ishell ++;
                }

	} while (istep<parms->nstep && softexit==0);	
	///////////////////////////////////////// end MC LOOP
		
	// print summary
	#ifdef ACTIVE_MPI		
	fprintf(stderr,"\nRank %d Acceptance:\t%d / %d = %lf\n",my_rank,parms->acc,parms->mov,(double)parms->acc/parms->mov);
	#else
	fprintf(stderr,"\nAcceptance:\t%d / %d = %lf\n",parms->acc,parms->mov,(double)parms->acc/parms->mov);
	#endif
	if(parms->movetype[0]!=-1)
		fprintf(stderr,"\tFlip:\t%d / %d = %lf\n",macc[0],mdone[0],(double)macc[0]/mdone[0]);
	if(parms->movetype[1]!=-1) 
		fprintf(stderr,"\tPivot:\t%d / %d = %lf\n",macc[1],mdone[1],(double)macc[1]/mdone[1]);
	if(parms->movetype[2]!=-1) 
		fprintf(stderr,"\tMPivot:\t%d / %d = %lf\n",macc[2],mdone[2],(double)macc[2]/mdone[2]);
	if(parms->movetype[3]!=-1) 
		fprintf(stderr,"\tSidechain:\t%d / %d = %lf\n",macc[3],mdone[3],(double)macc[3]/mdone[3]);
	if(parms->movetype[4]!=-1) 
		fprintf(stderr,"\tLPivot:\t%d / %d = %lf\n",macc[4],mdone[4],(double)macc[4]/mdone[4]);
	if(parms->movetype[5]!=-1) 
		fprintf(stderr,"\tMFlip:\t%d / %d = %lf\n",macc[5],mdone[5],(double)macc[5]/mdone[5]);
	if(parms->movetype[6]!=-1) 
		fprintf(stderr,"\tCoM shift:\t%d / %d = %lf\n",macc[6],mdone[6],(double)macc[6]/mdone[6]);
	if(parms->movetype[7]!=-1)
		fprintf(stderr,"\tLocal:\t%d / %d = %lf\n",macc[7],mdone[7],(double)macc[7]/mdone[7]); 
	 if(parms->movetype[8]!=-1)
                fprintf(stderr,"\tRotation:\t%d / %d = %lf\n",macc[8],mdone[8],(double)macc[8]/mdone[8]);

         if(parms->movetype[9]!=-1)
                fprintf(stderr,"\tCluster CoM:\t%d / %d = %lf\n",macc[9],mdone[9],(double)macc[9]/mdone[9]);

         if(parms->movetype[10]!=-1)
                fprintf(stderr,"\tCluster Rot:\t%d / %d = %lf\n",macc[10],mdone[10],(double)macc[10]/mdone[10]);


	fprintf(stderr,"\n");
	sprintf(p->title,"final");
	PrintPDBStream(p,parms->npol,ftrj);


}

/*****************************************************************************
 Make a flip in a random backbone atom
	 returns 1 if it could make the move, 0 if nothing was done
 *****************************************************************************/
int MoveBackboneFlip(struct s_polymer *p,struct s_polymer *oldp, struct s_mc_parms *parms, struct s_potential *pot,
		unsigned long long istep, int debug, double t)
{
	int iw,ip,ok=0,out;
	double dw,deltaE=0;

	ip = irand(parms->npol);						// which chain to move
	iw = irand( (p+ip)->nback );						// which backbone to move

	#ifdef DEBUG
	 if (debug>2) fprintf(stderr,"step=%d flip ip = %d iw=%d\n",istep,ip,iw); fflush(stderr);
	#endif

	// old energies
	deltaE = - GetEnergyMonomerRange(p,iw-1,iw+1,ip);			// old energy of iw with the others (iw-1 & iw+1 because I move their sidechains)

	if (!parms->noangpot)
	{
		if (iw>0) deltaE -= (((p+ip)->back)+iw-1)->e_ang ;		// old angular energy
		if (iw<(p+ip)->nback-1) deltaE -= (((p+ip)->back)+iw+1)->e_ang;
	}
	if (!parms->nodihpot)
	{
		if (iw>0) deltaE -= (((p+ip)->back)+iw-1)->e_dih;		// old dihedral energy
		deltaE -= (((p+ip)->back)+iw)->e_dih;
		if (iw<(p+ip)->nback-1) deltaE -= (((p+ip)->back)+iw+1)->e_dih;
		if (iw<(p+ip)->nback-2) deltaE -= (((p+ip)->back)+iw+2)->e_dih;
	}

	// make the move in p
	if (iw==0)
	{	
		ok = MoveHead((p+ip),parms);						// move the head
		ok *= AddSidechain(p,0,1,ip);
	}
	else if (iw==(p+ip)->nback-1)
	{
		ok = MoveTail((p+ip),parms);						// move the tail
		ok *= AddSidechain(p,(p+ip)->nback-2,(p+ip)->nback-1,ip);
	}
	else
	{
		dw = 2. * parms->dw_flip * (0.5 - frand());				// angle to flip
		ok = Flip((p+ip),iw,dw);						// make the flip
		ok *= AddSidechain(p,iw-1,iw+1,ip);
	}

	// check if constrains are violated

	if (parms->a_cloose>0)
	{
		if (iw-2>=0)
			if ( DAbs( Angle( (((p+ip)->back)+iw-2)->pos, (((p+ip)->back)+iw-1)->pos, (((p+ip)->back)+iw)->pos, (p+ip)->tables, &out)
			- (((p+ip)->back)+iw-1)->a_next) > parms->a_cloose || out==0) ok=0;
		if (ok && iw+2<(p+ip)->nback)
			if ( DAbs( Angle( (((p+ip)->back)+iw)->pos, (((p+ip)->back)+iw+1)->pos, (((p+ip)->back)+iw+2)->pos, (p+ip)->tables, &out)
					- (((p+ip)->back)+iw+1)->a_next) > parms->a_cloose || out==0) ok=0;
	}
	if (parms->d_cloose>0)
	{
		if (iw-3>=0 && (((p+ip)->back)+iw-1)->move == 0)
			if ( DAbs( Dihedral( (((p+ip)->back)+iw-3)->pos, (((p+ip)->back)+iw-2)->pos, (((p+ip)->back)+iw-1)->pos, (((p+ip)->back)+iw)->pos, (p+ip)->tables, &out)
						- (((p+ip)->back)+iw-1)->d_next) > parms->d_cloose || out==0) ok=0;
		if (ok && iw+3<(p+ip)->nback && (((p+ip)->back)+iw+2)->move == 0)
			if ( DAbs( Dihedral( (((p+ip)->back)+iw)->pos, (((p+ip)->back)+iw+1)->pos, (((p+ip)->back)+iw+2)->pos, (((p+ip)->back)+iw+3)->pos, (p+ip)->tables, &out)
						- (((p+ip)->back)+iw+2)->d_next) > parms->d_cloose || out==0) ok=0;
	}

	// calculate new energy
	if (ok)
	{
		deltaE += EnergyMonomerRange(p,pot,iw-1,iw+1,ip,parms->npol,parms->shell,1,parms->nosidechains,parms->disentangle,parms->hb);		// new energy of iw with the others
		if (!parms->noangpot)
		{
			if (iw>0) deltaE += EnergyAngles(p,pot,iw-1,ip,1);
			if (iw<(p+ip)->nback-1) deltaE += EnergyAngles(p,pot,iw+1,ip,1);
		}
		if (!parms->nodihpot)
		{
			if (iw>0) deltaE += EnergyDihedrals(p,pot,iw-1,ip,1);
			deltaE += EnergyDihedrals(p,pot,iw,ip,1);
			if (iw<(p+ip)->nback-1) deltaE += EnergyDihedrals(p,pot,iw+1,ip,1);
			if (iw<(p+ip)->nback-2) deltaE += EnergyDihedrals(p,pot,iw+2,ip,1);
		}

	}

	// Metropolis acceptance
	if (ok) ok = Metropolis(deltaE,t,p->tables);

	if (ok == 0) 			// move rejected
	{
		UpdateMonomerRange(oldp,p,iw-1,iw+1,ip,parms->shell);							// return to the old position, contacts, etc.
	}
	else					// move accepted
	{
		UpdateMonomerRange(p,oldp,iw-1,iw+1,ip,parms->shell);							// update oldp
		p->etot += deltaE;
		parms->acc ++;
	}

	#ifdef DEBUG
	if (debug>2) fprintf(stderr,"accept=%d\n",ok); fflush(stderr);
	#endif
	//fprintf(stderr,"** %lf %lf\n",TotalEnergy(p,pot,parms,parms->npol,0,parms->nosidechains,parms->debug),TotalEnergy(oldp,pot,parms,parms->npol,0,parms->nosidechains,parms->debug));
	parms->mov ++;

	if (ok == 0) return 0;
	else return 1;
}

/*****************************************************************************
 Copy an element of a polymer structure
 Watch the copy of the pointer to the rotamers --> do not free.
 *****************************************************************************/
void CopyResiduePositions(struct s_back *from, struct s_back *to)
{
	int k;

	to->irot = from->irot;
	(to->pos).x = (from->pos).x;
	(to->pos).y = (from->pos).y;
	(to->pos).z = (from->pos).z;
	(to->sph).ang = (from->sph).ang;
	(to->sph).dih = (from->sph).dih;
	(to->sph).r = (from->sph).r;

	//sidechain
	for (k=0;k<from->nside;k++)
	{
		(((to->side)+k)->pos).x = (((from->side)+k)->pos).x;
		(((to->side)+k)->pos).y = (((from->side)+k)->pos).y;
		(((to->side)+k)->pos).z = (((from->side)+k)->pos).z;
	}

	// contacts
	for (k=0;k<from->ncontacts;k++)
	{
		(to->contacts)[k] = (from->contacts)[k];
		(to->contacts_p)[k] = (from->contacts_p)[k];
		(to->e)[k] = (from->e)[k];
	}
	to->ncontacts = from->ncontacts;

	//shell
	for (k=0;k<from->nshell;k++)
	{
		(to->shell)[k] = (from->shell)[k];
		(to->shell_p)[k] = (from->shell_p)[k];
	}
	to->nshell = from->nshell;

}



void CopyResiduePositions_NOCONT(struct s_back *from, struct s_back *to)
{
        int k;

        to->irot = from->irot;
        (to->pos).x = (from->pos).x;
        (to->pos).y = (from->pos).y;
        (to->pos).z = (from->pos).z;
        (to->sph).ang = (from->sph).ang;
        (to->sph).dih = (from->sph).dih;
        (to->sph).r = (from->sph).r;
	
	 for (k=0;k<from->nside;k++)
        {
                (((to->side)+k)->pos).x = (((from->side)+k)->pos).x;
                (((to->side)+k)->pos).y = (((from->side)+k)->pos).y;
                (((to->side)+k)->pos).z = (((from->side)+k)->pos).z;
        }


	for (k=0;k<to->nshell;k++)
        {
                (to->shell)[k] = (from->shell)[k];
                (to->shell_p)[k] = (from->shell_p)[k];
        }
        to->nshell = from->nshell;

}

/*****************************************************************************
 Make a pivot move in a random backbone atom (very inefficient!)
	 returns 1 if it could make the move, 0 if nothing was done
 *****************************************************************************/
int MoveBackbonePivot(struct s_polymer *p,struct s_polymer *oldp, struct s_potential *pot, struct s_mc_parms *parms, double t)
{
	int iw,ok=0,half,ip;
	double dw=0,deltaE;
		
	ip = irand(parms->npol);

													// which chain to move
	half = (p+ip)->nback / 2;
	iw = 1 + irand((p+ip)->nback-2);

	if (parms->randdw==1) dw = parms->dw_pivot * (0.5 - frand());						// dihedral to pivot
	else if (parms->randdw==2) dw = parms->dw_pivot * gasdev(&(parms->seed));

	#ifdef DEBUG
		if (parms->debug>2) fprintf(stderr,"pivot ip = %d iw=%d dw=%lf half=%d\n",ip,iw,dw,half);
	#endif

	
	// if iw belongs to the first half, pivot backward (moving dih of iw)
	if (iw < half)
	{

		if ( (((p+ip)->back)+iw+1)->move == 0 ) return 0;

		// calculate old energy
		deltaE = -GetEnergyMonomerRange(p,0,iw,ip);							// two-body energy
		if (!parms->nodihpot) deltaE -= (((p+ip)->back)+iw+1)->e_dih;			// old dihedral energy


		// move
		ok = PivotBackward((p+ip),iw+1,dw,iw,parms);							// moved atoms are in [0,iw-1]; changed dihedral is iw+1
		ok *= AddSidechain(p,0,iw,ip);

		deltaE += EnergyMonomerRange(p,pot,0,iw,ip,parms->npol,0,1,parms->nosidechains,parms->disentangle,parms->hb);			// new energy of iw with the others
		if (!parms->nodihpot) deltaE += EnergyDihedrals(p,pot,iw+1,ip,1);

	}
	// if it belong to the second half, pivot forward (moving dih of iw)
	else
	{
	
		if ( (((p+ip)->back)+iw)->move == 0 ) return 0;

		// calculate old energy
		deltaE = -GetEnergyMonomerRange(p,iw,(p+ip)->nback-1,ip);
		if (!parms->nodihpot) deltaE -= (((p+ip)->back)+iw)->e_dih;				// old dihedral energy

		// move
		ok = PivotForward((p+ip),iw-1,dw,(p+ip)->nback-iw-1,parms);				// moved atoms are in [iw+1,nback-1]; changed dihedral of iw
		ok *= AddSidechain(p,iw,(p+ip)->nback-1,ip);

		deltaE += EnergyMonomerRange(p,pot,iw,(p+ip)->nback-1,ip,parms->npol,0,1,parms->nosidechains,parms->disentangle,parms->hb);
		if (!parms->nodihpot) deltaE += EnergyDihedrals(p,pot,iw,ip,1);
	}

	#ifdef DEBUG
	 int dir = (iw-half)<0 ? -1 : 1;
	 if (parms->debug>2) fprintf(stderr," ip=%d iw=%d dw=%lf direction=%d\n deltaE=%lf\n",ip,iw,dw,dir,deltaE); fflush(stderr);
	#endif

	ok *= Metropolis(deltaE,t,p->tables);

	if (ok == 0) 			// move rejected
	{

		if (iw < half) UpdateMonomerRange(oldp,p,0,iw,ip,0);                               // return to the old position, contacts, etc. (if moved first half)
		else UpdateMonomerRange(oldp,p,iw,(p+ip)->nback-1,ip,0);                           // ... if moved second half
		
	}
	else					// move accepted
	{
	        if (iw < half) UpdateMonomerRange(p,oldp,0,iw,ip,0);                               // update oldp (if moved first half)
	        else UpdateMonomerRange(p,oldp,iw,(p+ip)->nback-1,ip,0);     

		if(parms->shell==1 )		// if shells are active, update them 
		{
			UpdateShell(p,parms);
			CopyShell(p,oldp,parms);
			parms->ishell=0;
		}
		p->etot += deltaE;
		parms->acc ++;
	}

	#ifdef DEBUG
	if (parms->debug>2) fprintf(stderr," accept=%d\n",ok); fflush(stderr);
	#endif

	return ok;
}

/*****************************************************************************
 Make nmul consecutive pivot moves at random
	 returns 1 if it could make the move, 0 if nothing was done
 *****************************************************************************/
int MoveMultiplePivot(struct s_polymer *p, struct s_polymer *oldp, struct s_potential *pot, int nmul, struct s_mc_parms *parms, double t)
{
       int itemp;
//      fprintf(stderr,"\nprocessing protein  with nback=%i\n",(p+0)->nback);
      for(itemp=0;itemp<(p+0)->nback;itemp++)
      {
  //     fprintf(stderr,"\ndih %d =%lf",itemp,(((oldp+0)->back)+itemp)->sph.ang);
      }


	int half,ok=1,m,iw,ip,idir;
	double dw=0,deltaE;

	if(nmul<0) nmul=2+irand(-nmul-1);				//generate random nmul in [2,-mul] if required

	ip = irand(parms->npol);
	half = (p+ip)->nback / 2;
	iw = 1 + irand((p+ip)->nback-2);

	if (iw<half-nmul-1) idir = -1;				// which direction to go
	else if (iw>half+nmul+1) idir=1;
	else idir = 2 * irand(2) - 1;

	#ifdef DEBUG
		if (parms->debug>2) fprintf(stderr,"mpivot ip = %d iw=%d half=%d\n",ip,iw,half);
	#endif

	// if iw belongs to the first half, pivot backward
	if (idir == -1)
	{
		if ( iw+2-nmul < 2 ) nmul = iw;											// if nmul beyond beginning of the chain, decrease it

		deltaE = -GetEnergyMonomerRange(p,0,iw,ip);							// calculate old energy
		if (!parms->nodihpot)
			for (m=0;m<nmul;m++) deltaE -= (((p+ip)->back)+iw-m+1)->e_dih;		// old dihedral energy

		for (m=0;m<nmul;m++)													// apply nmul pivot moves
			if ( (((p+ip)->back)+iw-m+1)->move == 1 )								// check if it you can move it
			{
				if (parms->randdw==1) dw = parms->dw_mpivot * (0.5 - frand());				// dihedral to pivot
				else if (parms->randdw==2) dw = parms->dw_mpivot *gasdev(&(parms->seed));

 				ok *= PivotBackward((p+ip),iw+1-m,dw,iw-m,parms);				// moves in [0,iw-1]
			}
		ok *= AddSidechain(p,0,iw,ip);

//		deltaE += EnergyMonomerRange(p,pot,0,iw,ip,parms->npol,parms->shell,1,parms->nosidechains,parms->disentangle,parms->hb);				// new energy of iw with the others
		deltaE += EnergyMonomerRange(p,pot,0,iw,ip,parms->npol,0,1,parms->nosidechains,parms->disentangle,parms->hb);    


		if (!parms->nodihpot)													// new dihedral energy
			for (m=0;m<nmul;m++) deltaE += EnergyDihedrals(p,pot,iw-m+1,ip,1);
	}
	// if it belong to the second half, pivot forward
	else
	{
		if ( iw-2+nmul > (p+ip)->nback - 3 ) nmul = (p+ip)->nback - 1 - iw;			// if nmul is beyond the end, decrease it

		deltaE = -GetEnergyMonomerRange(p,iw,(p+ip)->nback-1,ip);					// calculate old energy
		if (!parms->nodihpot)
			for (m=0;m<nmul;m++)deltaE -= (((p+ip)->back)+iw+m)->e_dih;				// old dihedral energy

		for (m=0;m<nmul;m++)														// apply nmul pivot moves
			if ( (((p+ip)->back)+iw+m)->move == 1 )								// check if it you can move it
			{
				if (parms->randdw==1) dw = parms->dw_mpivot * (0.5 - frand());					// dihedral to pivot
				else if (parms->randdw==2) dw = parms->dw_mpivot * gasdev(&(parms->seed));

				ok *= PivotForward((p+ip),iw-1+m,dw,(p+ip)->nback-iw-1-m,parms);	// moves in [iw+1,nback-1]
			}
		ok *= AddSidechain(p,iw,(p+ip)->nback-1,ip);

		//shell disattivate: mossa globale
		deltaE += EnergyMonomerRange(p,pot,iw,(p+ip)->nback-1,ip,parms->npol,0,1,parms->nosidechains,parms->disentangle,parms->hb);   
	//	deltaE += EnergyMonomerRange(p,pot,iw,p->nback-1,ip,parms->npol,parms->shell,1,parms->nosidechains,parms->disentangle,parms->hb);			// new energy of iw with the others
		if (!parms->nodihpot)														// new dihedral energy
			for (m=0;m<nmul;m++) deltaE += EnergyDihedrals(p,pot,iw+m,ip,1);
	}

//	fprintf(stderr,"\ndeltaE=%lf\n",deltaE);
      #ifdef DEBUG
	if (parms->debug>2) fprintf(stderr,"deltaE=%lf\n",deltaE); fflush(stderr);
	#endif

	ok *= Metropolis(deltaE,t,p->tables);

	if (ok == 0) 			// move rejected
	{	

	//	if (idir == -1) UpdateMonomerRange(oldp,p,0,iw,ip,parms->shell);				// return to the old position, contacts, etc. (if moved first half)
		if (idir == -1) UpdateMonomerRange(oldp,p,0,iw,ip,0);


		else
			//UpdateMonomerRange(oldp,p,iw,p->nback-1,ip,parms->shell);				// ... if moved second half
			UpdateMonomerRange(oldp,p,iw,(p+ip)->nback-1,ip,0);
			//...if moved second half

			
	}
	else					// move accepted
	{
//		if (idir == -1) UpdateMonomerRange(p,oldp,0,iw,ip,parms->shell);				// update oldp (if moved first half)
		if ( idir == -1) UpdateMonomerRange(p,oldp,0,iw,ip,0);


//		else UpdateMonomerRange(p,oldp,iw,p->nback-1,ip,parms->shell);		
		// ... if moved second half
		 else UpdateMonomerRange(p,oldp,iw,(p+ip)->nback-1,ip,0);

		if(parms->shell==1 )            // if shells are active, update them 
                {
                        UpdateShell(p,parms);
                        CopyShell(p,oldp,parms);
                        parms->ishell=0;
                }


 
		p->etot += deltaE;
		parms->acc ++;
	}

	#ifdef DEBUG
	if (parms->debug>2) fprintf(stderr,"accept=%d\n",ok); fflush(stderr);
	#endif

	parms->mov ++;

	return ok;
}



/*****************************************************************************
 Metropolis acceptance rule
 1 = accept, 0 = reject
 *****************************************************************************/
int Metropolis(double deltaE, double T, struct s_tables *t)
{
	double p;

	if (deltaE<=0) return 1;

	p = FastExp(-deltaE/T,t);
      //fprintf(stderr,"\n p = %lf",p);

    //  fprintf(stderr,"\n %lf < %lf ?",a,p);
	if (frand()<p) return 1;

	return 0;
}


/*****************************************************************************
 Copy the position of the monomer w to another structure, updating also
 contacts, partial energies and shell
 The other structure should be identical to the former except for residue w
 *****************************************************************************/
void UpdateMonomer(struct s_polymer *from, struct s_polymer *to, int w, int n, int shell)
{
	int i,j,l[NCONTMAX2],lp[NCONTMAX2],ntm;//ll[NSHELLMAX],llp[NSHELLMAX];
	// control beginning and end
	if (w<0) w=0;
	if (w>(from+n)->nback-1) w = (from+n)->nback-1;

	ntm=0;
	// Find all residues that are or were in contact with w
	for (i=0;i< (((from+n)->back)+w)->ncontacts;i++)
	{
		l[ntm] = ((((from+n)->back)+w)->contacts)[i];				// which backbone is close to in the from structure
		lp[ntm] = ((((from+n)->back)+w)->contacts_p)[i];			// which chain
		ntm++;
		if (ntm>NCONTMAX2) Error("NCONTMAX2 too small in UpdateMonomerRange");
	}

	for (i=0;i<(((to+n)->back)+w)->ncontacts;i++)
	{
		l[ntm]=((((to+n)->back)+w)->contacts)[i];
		lp[ntm] = ((((to+n)->back)+w)->contacts_p)[i];
		ntm++;
		if (ntm>NCONTMAX2) Error("NCONTMAX2 too small in UpdateMonomerRange");

	}


	// Copy those contacts
	for (i=0;i<ntm;i++)
	{
		(((to+lp[i])->back)+l[i])->ncontacts = (((from+lp[i])->back)+l[i])->ncontacts;
		for (j=0;j<(((from+lp[i])->back)+l[i])->ncontacts;j++)
		{
			((((to+lp[i])->back)+l[i])->contacts)[j] = ((((from+lp[i])->back)+l[i])->contacts)[j];
			((((to+lp[i])->back)+l[i])->contacts_p)[j] = ((((from+lp[i])->back)+l[i])->contacts_p)[j];
			((((to+lp[i])->back)+l[i])->e)[j] = ((((from+lp[i])->back)+l[i])->e)[j];
		}
	}

	(((to+n)->back)+w)->e_ang = (((from+n)->back)+w)->e_ang;
	if (w>0) (((to+n)->back)+w-1)->e_ang = (((from+n)->back)+w-1)->e_ang;
	if (w<(from+n)->nback-1) (((to+n)->back)+w+1)->e_ang = (((from+n)->back)+w+1)->e_ang;

	if (w>1) (((to+n)->back)+w-2)->e_dih = (((from+n)->back)+w-2)->e_dih;
	if (w>0) (((to+n)->back)+w-1)->e_dih = (((from+n)->back)+w-1)->e_dih;
	(((to+n)->back)+w)->e_dih = (((from+n)->back)+w)->e_dih;
	if (w<((from+n)->nback)-1) (((to+n)->back)+w+1)->e_dih = (((from+n)->back)+w+1)->e_dih;
        (((to+n)->back)+w)->irot = (((from+n)->back)+w)->irot;

	if (w<((from+n)->nback)-2) (((to+n)->back)+w+2)->e_dih = (((from+n)->back)+w+2)->e_dih;

	// Finally, copy w itself
	CopyResiduePositions( (((from+n)->back)+w), (((to+n)->back)+w) );

}

void UpdateMonomerRange(struct s_polymer *from, struct s_polymer *to, int wfrom, int wto, int p, int shell)
{
	int i,w,j,l[NCONTMAX2],lp[NSHELLMAX],ntm,k,chk;

	// Control beginning and end
	if (wfrom<0) wfrom=0;
	if (wto>(from+p)->nback-1) wto = (from+p)->nback-1;
	// Find all residues that are or were in contact with moved monomers w

	ntm=0;
	for (w=wfrom;w<=wto;w++)											// w ranges over all changed residues
		for (i=0;i< (((from+p)->back)+w)->ncontacts;i++)			// i ranges over the contacts of w
		{
			chk=0;														// count contacts only once
			for (k=0;k<ntm;k++)
				if (l[k]==((((from+p)->back)+w)->contacts)[i] && lp[k]==((((from+p)->back)+w)->contacts_p)[i] ) chk=1;
			if (chk==0)
			{

				l[ntm] = (((from+p)->back)+w)->contacts[i];				// which backbone is close to in the from structure
				lp[ntm] = (((from+p)->back)+w)->contacts_p[i];			// which chain
				ntm++;
				if (ntm>NCONTMAX2) Error("NCONTMAX2 too small in UpdateMonomerRange");
			}
		}
	
	
	for (w=wfrom;w<=wto;w++)											// w ranges over all changed residues
		for (i=0;i< (((to+p)->back)+w)->ncontacts;i++)
		{
			chk=0;
												// count contacts only once
			for (k=0;k<ntm;k++)
				if (l[k]==((((to+p)->back)+w)->contacts)[i] && lp[k]==((((to+p)->back)+w)->contacts_p)[i] ) chk=1;
			if (chk==0)
			{
				
				l[ntm] = ((((to+p)->back)+w)->contacts)[i];				// same in the to structure
				lp[ntm] = ((((to+p)->back)+w)->contacts_p)[i];
				ntm++;
				if (ntm>NCONTMAX2) Error("NCONTMAX2 too small in UpdateMonomerRange");
			}
		}


	
	// Update angle/dihedral energies from w-1 to w+2
	for (w=wfrom;w<=wto;w++)
	{
		(((to+p)->back)+w)->e_ang = (((from+p)->back)+w)->e_ang;
		if (w>0) (((to+p)->back)+w-1)->e_ang = (((from+p)->back)+w-1)->e_ang;
		if (w<(from+p)->nback-1) (((to+p)->back)+w+1)->e_ang = (((from+p)->back)+w+1)->e_ang;

		if (w>0) (((to+p)->back)+w-1)->e_dih = (((from+p)->back)+w-1)->e_dih;
		(((to+p)->back)+w)->e_dih = (((from+p)->back)+w)->e_dih;
		if (w<((from+p)->nback)-1) (((to+p)->back)+w+1)->e_dih = (((from+p)->back)+w+1)->e_dih;
		if (w<((from+p)->nback)-2) (((to+p)->back)+w+2)->e_dih = (((from+p)->back)+w+2)->e_dih;

		(((to+p)->back)+w)->irot = (((from+p)->back)+w)->irot;
		if (w<((from+p)->nback)-2) (((to+p)->back)+w+2)->e_dih = (((from+p)->back)+w+2)->e_dih;

	}

	// Copy those contacts
	for (i=0;i<ntm;i++)
	{

		        (((to+lp[i])->back)+l[i])->ncontacts = (((from+lp[i])->back)+l[i])->ncontacts;
 	                for (j=0;j<(((from+lp[i])->back)+l[i])->ncontacts;j++)
			{
				((((to+lp[i])->back)+l[i])->contacts)[j] = ((((from+lp[i])->back)+l[i])->contacts)[j];
                                ((((to+lp[i])->back)+l[i])->contacts_p)[j] = ((((from+lp[i])->back)+l[i])->contacts_p)[j];
                                ((((to+lp[i])->back)+l[i])->e)[j] = ((((from+lp[i])->back)+l[i])->e)[j];
			}	


	}
 



//	 Finally, copy w itself
	for (w=wfrom;w<=wto;w++)											// w ranges over all changed residues
		CopyResiduePositions( ((from+p)->back)+w, ((to+p)->back)+w) ;

}






/*****************************************************************************
 Move the sidechain among the possible rotamers
	 returns 1 if it could make the move, 0 if nothing was done,
	 -1 if it couldn't find a sidechain to move
 *****************************************************************************/
int MoveSidechain(struct s_polymer *p,struct s_polymer *oldp, struct s_mc_parms *mc_parms, struct s_potential *pot,
		unsigned long long istep, int debug, double t)
{
	int iw,ok=0,chk=0,ir,ip;
	double deltaE;


	// which chain to move
	ip = irand(mc_parms->npol);

	// select among sidechains which have more than 1 rotamer
	do
	{
		iw = irand((p+ip)->nback);
		chk++;
		if (chk>(p+ip)->nback*3) return -1;
	}
	while ( (((p+ip)->back)+iw)->nside < 1 || (((p+ip)->back)+iw)->nrot ==1 );


	#ifdef DEBUG
	 if (debug>2) fprintf(stderr,"step=%d sidechain iw=%d\n",istep,iw); fflush(stderr);
	#endif

	deltaE = - GetEnergyMonomer(p,ip,iw);								// old energy of iw with the others

	ir = irand( (((p+ip)->back)+iw)->nrot );
	(((p+ip)->back)+iw)->irot = ir;
	ok = AddSidechain(p,iw,iw,ip);

	deltaE += EnergyMonomer(p,pot,iw,ip,mc_parms->npol,1,mc_parms->shell,mc_parms->nosidechains,mc_parms->disentangle,mc_parms->hb);			// new energy of iw with the others
	//deltaE += EnergyMonomer(p,pot,iw,ip,mc_parms->npol,1,0,mc_parms->nosidechains,mc_parms->disentangle,mc_parms->hb); 
	#ifdef DEBUG
	if (debug>2) fprintf(stderr,"deltaE=%lf\n",deltaE); fflush(stderr);
	#endif

	ok *= Metropolis(deltaE,t,p->tables);

	if (ok == 0) 			// move rejected
	{
		UpdateMonomer(oldp,p,iw,ip,mc_parms->shell);
		// return to the old position, contacts, etc.
		//UpdateMonomer(oldp,p,iw,ip,0);
	}
	else					// move accepted
	{
		//fprintf(stderr,"move accepted\n");
		UpdateMonomer(p,oldp,iw,ip,mc_parms->shell);			
		//update oldp
		//UpdateMonomer(p,oldp,iw,ip,0);
		p->etot += deltaE;
		mc_parms->acc ++;

	}

	#ifdef DEBUG
	if (debug>2) fprintf(stderr,"accept=%d\n",ok); fflush(stderr);
	#endif

	mc_parms->mov ++;

	if (ok == 0) return 0;
	else return 1;
}

void CompareStructures(struct s_polymer *a, struct s_polymer *b, int nc, int natoms)
{
	int ic,i,j,irot;


	for (ic=0;ic<nc;ic++)
		for (i=0;i<(a+ic)->nback;i++)
		{
			irot = (((a+ic)->back)+i)->irot;

			if ( (((a+ic)->back)+i)->ia != (((b+ic)->back)+i)->ia ) fprintf(stderr,"*** i=%d ic=%d ia=%d %d\n",i,ic,(((a+ic)->back)+i)->ia,(((b+ic)->back)+i)->ia);
			if ( (((a+ic)->back)+i)->iaa != (((b+ic)->back)+i)->iaa ) fprintf(stderr,"*** i=%d ic=%d iaa=%d %d\n",i,ic,(((a+ic)->back)+i)->iaa,(((b+ic)->back)+i)->iaa);
			if ( (((a+ic)->back)+i)->itype != (((b+ic)->back)+i)->itype ) fprintf(stderr,"*** i=%d ic=%d itype=%d %d\n",i,ic,(((a+ic)->back)+i)->itype,(((b+ic)->back)+i)->itype);
			if ( (((a+ic)->back)+i)->nside != (((b+ic)->back)+i)->nside ) fprintf(stderr,"*** i=%d ic=%d nside=%d %d\n",i,ic,(((a+ic)->back)+i)->nside,(((b+ic)->back)+i)->nside);
			if ( (((a+ic)->back)+i)->move != (((b+ic)->back)+i)->move ) fprintf(stderr,"*** i=%d ic=%d move=%d %d\n",i,ic,(((a+ic)->back)+i)->move,(((b+ic)->back)+i)->move);
			if ( (((a+ic)->back)+i)->irot != (((b+ic)->back)+i)->irot ) fprintf(stderr,"*** i=%d ic=%d irot=%d %d\n",i,ic,(((a+ic)->back)+i)->irot,(((b+ic)->back)+i)->irot);
			if ( (((a+ic)->back)+i)->nrot != (((b+ic)->back)+i)->nrot ) fprintf(stderr,"*** i=%d ic=%d nrot=%d %d\n",i,ic,(((a+ic)->back)+i)->nrot,(((b+ic)->back)+i)->nrot);
			if ( (((a+ic)->back)+i)->ncontacts != (((b+ic)->back)+i)->ncontacts ) fprintf(stderr,"*** i=%d ic=%d ncotacts=%d %d\n",i,ic,(((a+ic)->back)+i)->ncontacts,(((b+ic)->back)+i)->ncontacts);
			if ( (((a+ic)->back)+i)->pos.x != (((b+ic)->back)+i)->pos.x ) fprintf(stderr,"*** i=%d ic=%d pos.x=%lf %lf\n",i,ic,(((a+ic)->back)+i)->pos.x,(((b+ic)->back)+i)->pos.x);
			if ( (((a+ic)->back)+i)->pos.y != (((b+ic)->back)+i)->pos.y ) fprintf(stderr,"*** i=%d ic=%d pos.y=%lf %lf\n",i,ic,(((a+ic)->back)+i)->pos.y,(((b+ic)->back)+i)->pos.y);
			if ( (((a+ic)->back)+i)->pos.z != (((b+ic)->back)+i)->pos.z ) fprintf(stderr,"*** i=%d ic=%d pos.z=%lf %lf\n",i,ic,(((a+ic)->back)+i)->pos.z,(((b+ic)->back)+i)->pos.z);
			if ( (((a+ic)->back)+i)->sph.ang!= (((b+ic)->back)+i)->sph.ang ) fprintf(stderr,"*** i=%d ic=%d ang=%lf %lf\n",i,ic,(((a+ic)->back)+i)->sph.ang,(((b+ic)->back)+i)->sph.ang);
			if ( (((a+ic)->back)+i)->sph.dih!= (((b+ic)->back)+i)->sph.dih ) fprintf(stderr,"*** i=%d ic=%d dih=%lf %lf\n",i,ic,(((a+ic)->back)+i)->sph.dih,(((b+ic)->back)+i)->sph.dih);
			if ( (((a+ic)->back)+i)->sph.r!= (((b+ic)->back)+i)->sph.r ) fprintf(stderr,"*** i=%d ic=%d r=%lf %lf\n",i,ic,(((a+ic)->back)+i)->sph.r,(((b+ic)->back)+i)->sph.r);
			for (j=0;j<(((a+ic)->back)+i)->nside;j++)
			{
				if ( (((((a+ic)->back)+i)->side)+j)->ia != (((((b+ic)->back)+i)->side)+j)->ia ) fprintf(stderr,"*** i=%d ic=%d side %d ia=%d %d\n",i,ic,j,
					(((((a+ic)->back)+i)->side)+j)->ia,(((((b+ic)->back)+i)->side)+j)->ia);
				if ( (((((a+ic)->back)+i)->side)+j)->pos.x != (((((b+ic)->back)+i)->side)+j)->pos.x ) fprintf(stderr,"*** i=%d ic=%d side %d x=%lf %lf\n",i,ic,j,
					(((((a+ic)->back)+i)->side)+j)->pos.x,(((((b+ic)->back)+i)->side)+j)->pos.x);
				if ( (((((a+ic)->back)+i)->side)+j)->pos.y != (((((b+ic)->back)+i)->side)+j)->pos.y ) fprintf(stderr,"*** i=%d ic=%d side %d y=%lf %lf\n",i,ic,j,
					(((((a+ic)->back)+i)->side)+j)->pos.y,(((((b+ic)->back)+i)->side)+j)->pos.y);
				if ( (((((a+ic)->back)+i)->side)+j)->pos.z != (((((b+ic)->back)+i)->side)+j)->pos.z ) fprintf(stderr,"*** i=%d ic=%d side %d z=%lf %lf\n",i,ic,j,
					(((((a+ic)->back)+i)->side)+j)->pos.z,(((((b+ic)->back)+i)->side)+j)->pos.z);
				if ( (((((((a+ic)->back)+i)->side)+j)->rot)+irot)->b1 != (((((((b+ic)->back)+i)->side)+j)->rot)+irot)->b1 ) fprintf(stderr,"*** i=%d ic=%d side %d b1=%d %d\n",
					i,ic,j,(((((((a+ic)->back)+i)->side)+j)->rot)+irot)->b1,(((((((b+ic)->back)+i)->side)+j)->rot)+irot)->b1 );
				if ( (((((((a+ic)->back)+i)->side)+j)->rot)+irot)->b2 != (((((((b+ic)->back)+i)->side)+j)->rot)+irot)->b2 ) fprintf(stderr,"*** i=%d ic=%d side %d b2=%d %d\n",
					i,ic,j,(((((((a+ic)->back)+i)->side)+j)->rot)+irot)->b2,(((((((b+ic)->back)+i)->side)+j)->rot)+irot)->b2 );
				if ( (((((((a+ic)->back)+i)->side)+j)->rot)+irot)->b3 != (((((((b+ic)->back)+i)->side)+j)->rot)+irot)->b3 ) fprintf(stderr,"*** i=%d ic=%d side %d b3=%d %d\n",
					i,ic,j,(((((((a+ic)->back)+i)->side)+j)->rot)+irot)->b3,(((((((b+ic)->back)+i)->side)+j)->rot)+irot)->b3 );
				if ( (((((((a+ic)->back)+i)->side)+j)->rot)+irot)->ang.ang != (((((((b+ic)->back)+i)->side)+j)->rot)+irot)->ang.ang ) fprintf(stderr,"*** i=%d ic=%d side %d ang=%lf %lf\n",
					i,ic,j,(((((((a+ic)->back)+i)->side)+j)->rot)+irot)->ang.ang,(((((((b+ic)->back)+i)->side)+j)->rot)+irot)->ang.ang );
				if ( (((((((a+ic)->back)+i)->side)+j)->rot)+irot)->ang.dih != (((((((b+ic)->back)+i)->side)+j)->rot)+irot)->ang.dih ) fprintf(stderr,"*** i=%d ic=%d side %d dih=%lf %lf\n",
					i,ic,j,(((((((a+ic)->back)+i)->side)+j)->rot)+irot)->ang.dih,(((((((b+ic)->back)+i)->side)+j)->rot)+irot)->ang.dih );
				if ( (((((((a+ic)->back)+i)->side)+j)->rot)+irot)->ang.r != (((((((b+ic)->back)+i)->side)+j)->rot)+irot)->ang.r ) fprintf(stderr,"*** i=%d ic=%d side %d r=%lf %lf\n",
					i,ic,j,(((((((a+ic)->back)+i)->side)+j)->rot)+irot)->ang.r,(((((((b+ic)->back)+i)->side)+j)->rot)+irot)->ang.r );



			}

		}

	for (ic=0;ic<nc;ic++)
		for (i=0;i<natoms;i++)
				{
				//	if ( ((a+ic)->vback[i]) != ((b+ic)->vback[i]) ) fprintf(stderr,"*** ia=%d different address\n",i);
					if ( (*((a+ic)->vback[i])).x != (*((b+ic)->vback[i])).x ) fprintf(stderr,"*** ia=%d x=%lf %lf\n",i,(*((a+ic)->vback[i])).x,(*((b+ic)->vback[i])).x);
					if ( (*((a+ic)->vback[i])).y != (*((b+ic)->vback[i])).y ) fprintf(stderr,"*** ia=%d y=%lf %lf\n",i,(*((a+ic)->vback[i])).y,(*((b+ic)->vback[i])).y);
					if ( (*((a+ic)->vback[i])).z != (*((b+ic)->vback[i])).z ) fprintf(stderr,"*** ia=%d z=%lf %lf\n",i,(*((a+ic)->vback[i])).z,(*((b+ic)->vback[i])).z);
				}
}



/*****************************************************************************
 Make nmul consecutive pivot moves at random and keep fixed  the nmul+1,
 chaniging the bond length
	 returns 1 if it could make the move, 0 if nothing was done
 *****************************************************************************/

int MoveLoosePivot(struct s_polymer *p, struct s_polymer *oldp, struct s_potential *pot, int nmul, struct s_mc_parms *parms, double t)
{
      int half,ok=1,m,iw,ip,idir;
      double dw=0,deltaE,rc2;

      if(nmul<0) nmul=2+irand(-nmul-1);                     //generate random nmul in [2,-mul] if required

      rc2 = parms->r_cloose * parms->r_cloose;

      ip = irand(parms->npol);
      half = (p+ip)->nback / 2;
      iw = 1 + irand((p+ip)->nback-2);

      if (iw<half-nmul-1) idir = -1;                        // which direction to go
      else if (iw>half+nmul+1) idir=1;
      else idir = 2 * irand(2) - 1;

      #ifdef DEBUG
            if (parms->debug>2) fprintf(stderr,"lpivot ip = %d iw=%d half=%d  nmul=%d\n",ip,iw,half,nmul);
      #endif

      if (idir == -1)
      {
            if ( iw-nmul-1 < 2 ) nmul = iw - 1;                                                       // if nmul beyond beginning of the chain, decrease it

            deltaE = -GetEnergyMonomerRange(p,iw-nmul-1,iw,ip);                                 // calculate old energy
            if (!parms->nodihpot)
                  for (m=0;m<nmul+3;m++) deltaE -= (((p+ip)->back)+iw-m+1)->e_dih;  // old dihedral energy in [iw-nmul-1,iw+1]
            if (!parms->noangpot)
            {
                  deltaE -= (((p+ip)->back)+iw-nmul-1)->e_ang ;                                 // old angular energy
                  deltaE -= (((p+ip)->back)+iw-nmul)->e_ang;
            }

            for (m=0;m<nmul;m++)                                                                            // apply nmul pivot moves
                  if ( (((p+ip)->back)+iw-m)->move == 1 )                                             // check if it you can move it
                  {
                        if (parms->randdw==1) dw = parms->dw_lpivot * (0.5 - frand());                            // dihedral to pivot
                        else if (parms->randdw==2) dw = parms->dw_lpivot * gasdev(&(parms->seed));

                        ok *= PivotBackward((p+ip),iw+1-m,dw,nmul-m,parms);                     // moves iw-1-m => overall moves in [iw-nmul,iw-1] (and sidechains of iw-imul-1 and iw)
                  }
            if (iw-nmul-1>=0)
                  if ( Abs(Dist2( (((p+ip)->back)+iw-nmul)->pos, (((p+ip)->back)+iw-nmul-1)->pos ) - (((p+ip)->back)+iw-nmul-1)->d2_next) > rc2 ) ok=0;           // if covalent bond is broken

            if (ok==1)
            {
                  ok = AddSidechain(p,iw-nmul-1,iw,ip);

                  deltaE += EnergyMonomerRange(p,pot,iw-nmul-1,iw,ip,parms->npol,parms->shell,1,
                                                                        parms->nosidechains,parms->disentangle,parms->hb);                      // new energy of iw with the others
                  if (!parms->nodihpot)                                                                                                         // new dihedral energy
                        for (m=0;m<nmul+3;m++) deltaE += EnergyDihedrals(p,pot,iw-m+1,ip,1);
                  if (!parms->noangpot)                                                                                                         // new angle energy
                  {
                        deltaE += EnergyAngles(p,pot,iw-nmul-1,ip,1);
                        deltaE += EnergyAngles(p,pot,iw-nmul,ip,1);
                  }
            }
      }
      else
      {
            if ( iw+2+nmul > (p+ip)->nback - 1 ) nmul = (p+ip)->nback - 3 - iw;                 // if nmul is beyond the end, decrease it

            deltaE = -GetEnergyMonomerRange(p,iw,iw+nmul+1,ip);                                       // calculate old energy
            if (!parms->nodihpot)
                  for (m=0;m<nmul+3;m++) deltaE -= (((p+ip)->back)+iw+m)->e_dih;                // old dihedral energy
            if (!parms->noangpot)
            {
                  deltaE -= (((p+ip)->back)+iw+nmul+1)->e_ang ;                                       // old angular energy
                  deltaE -= (((p+ip)->back)+iw+nmul)->e_ang;
            }

            for (m=0;m<nmul;m++)                                                                                  // apply nmul pivot moves
                  if ( (((p+ip)->back)+iw-1+m)->move == 1 )                                           // check if it you can move it
                  {
                        if (parms->randdw==1) dw = parms->dw_lpivot * (0.5 - frand());                      // dihedral to pivot
                        else if (parms->randdw==2) dw = parms->dw_lpivot * gasdev(&(parms->seed));

                        ok *= PivotForward((p+ip),iw-1+m,dw,nmul-m,parms);                            // moves iw+1+m => overall moves in [iw+1,iw+nmul]
                  }
            if (iw+nmul+1<(p+ip)->nback)
                  if ( Abs(Dist2( (((p+ip)->back)+iw+nmul)->pos, (((p+ip)->back)+iw+nmul+1)->pos ) - (((p+ip)->back)+iw+nmul+1)->d2_next) > rc2 ) ok=0;           // if covalent bond is broken

            if (ok==1)
            {
                  ok = AddSidechain(p,iw,iw+nmul+1,ip);

                  deltaE += EnergyMonomerRange(p,pot,iw,iw+nmul+1,ip,parms->npol,parms->shell,1,
                                                                                    parms->nosidechains,parms->disentangle,parms->hb);                // new energy of iw with the others
                  if (!parms->nodihpot)                                                                                                               // new dihedral energy
                        for (m=0;m<nmul+3;m++) deltaE += EnergyDihedrals(p,pot,iw+m,ip,1);
                  if (!parms->noangpot)
                  {
                        deltaE += EnergyAngles(p,pot,iw+nmul+1,ip,1);                                                                     // new angular energy
                        deltaE += EnergyAngles(p,pot,iw+nmul,ip,1);
                  }
            }
      }

      #ifdef DEBUG
      if (parms->debug>2) fprintf(stderr,"deltaE=%lf\n",deltaE); fflush(stderr);
      #endif

      if (ok==1) ok = Metropolis(deltaE,t,p->tables);

      if (ok == 0)                  // move rejected
      {
            if (idir == -1) UpdateMonomerRange(oldp,p,iw-nmul-1,iw,ip,parms->shell);            // return to the old position, contacts, etc. (if moved first half)
            else UpdateMonomerRange(oldp,p,iw,iw+nmul+1,ip,parms->shell);                       // ... if moved second half
      }
      else                          // move accepted
      {
            if (idir == -1) UpdateMonomerRange(p,oldp,iw-nmul-1,iw,ip,parms->shell);                        // update oldp (if moved first half)
            else UpdateMonomerRange(p,oldp,iw,iw+nmul+1,ip,parms->shell);                       // ... if moved second half
            p->etot += deltaE;
            parms->acc ++;
      }

      #ifdef DEBUG
      if (parms->debug>2) fprintf(stderr,"accept=%d\n",ok); fflush(stderr);
      #endif

      parms->mov ++;

      return ok;
}






void SoftExit()
{
	softexit = 1;
}

/*****************************************************************************
 Make a flip of a whole piece of the chain
	 returns 1 if it could make the move, 0 if nothing was done
 *****************************************************************************/
int MoveMultipleFlip(struct s_polymer *p,struct s_polymer *oldp, struct s_mc_parms *parms, struct s_potential *pot,
		unsigned long long istep, int debug, double t)
{
	int iw1,iw2,ip,ok=0,i,out;
	double dw,deltaE=0;

	ip = irand(parms->npol);						// which chain to move
	iw1 = irand( (p+ip)->nback-6 ) + 3;					// which backbone to move as first
	do { i = irand(2*parms->nmul_mflip) - parms->nmul_mflip; } while ( iw1+i > (p+ip)->nback-4 || iw1+i < 3 || (i>-3 && i<3) );
	iw2 = iw1 + i;								// ... which to move as last
	if (iw2<iw1) { i = iw1; iw1=iw2; iw2=i; }				// iw1 should be lower than iw2 for geometry

	#ifdef DEBUG
	 if (debug>2) fprintf(stderr,"step=%d mflip ip=%d iw1=%d iw2=%d\n",istep,ip,iw1,iw2); fflush(stderr);
	#endif

	// old energies
	deltaE = - GetEnergyMonomerRange(p,iw1-1,iw2+1,ip);			// old energy  (iw1-1 & iw2+1 because I move their sidechains)

	if (!parms->noangpot)
	{
		if (iw1>0) deltaE -= (((p+ip)->back)+iw1-1)->e_ang ;		// old angular energy
		if (iw2<(p+ip)->nback-1) deltaE -= (((p+ip)->back)+iw2+1)->e_ang;
	}
	if (!parms->nodihpot)
	{
		if (iw1>0) deltaE -= (((p+ip)->back)+iw1-1)->e_dih;		// old dihedral energy
		deltaE -= (((p+ip)->back)+iw1)->e_dih;
		if (iw2<(p+ip)->nback-1) deltaE -= (((p+ip)->back)+iw2+1)->e_dih;
		if (iw2<(p+ip)->nback-2) deltaE -= (((p+ip)->back)+iw2+2)->e_dih;
	}

	// make the move in p

	dw = 2. * parms->dw_mflip * (0.5 - frand());				// angle to flip
	ok = FlipFragment((p+ip),iw1,iw2,dw);					// make the flip

	if (ok==0)
	{
		UpdateMonomerRange(oldp,p,iw1-1,iw2+1,ip,parms->shell);		// return to the old position, contacts, etc.
		return 0;
	}

	// check if constrains are violated
	// MODIFICHE: corretto out nell if
	if (parms->a_cloose>0)
	{
		if (iw1-2>=0)
			if ( DAbs( Angle( (((p+ip)->back)+iw1-2)->pos, (((p+ip)->back)+iw1-1)->pos, (((p+ip)->back)+iw1)->pos, (p+ip)->tables, &out)
			- (((p+ip)->back)+iw1-1)->a_next) > parms->a_cloose  || out==0 ) ok=0;
		if (ok && iw2+2<(p+ip)->nback)
			if ( DAbs( Angle( (((p+ip)->back)+iw2)->pos, (((p+ip)->back)+iw2+1)->pos, (((p+ip)->back)+iw2+2)->pos, (p+ip)->tables, &out)
					- (((p+ip)->back)+iw2+1)->a_next) > parms->a_cloose  || out==0 ) ok=0;
	}
	//FINE MODIFICHE
	if (parms->d_cloose>0)
	{
		if (iw1-3>=0 && (((p+ip)->back)+iw1-1)->move == 0)
			if ( DAbs( Dihedral( (((p+ip)->back)+iw1-3)->pos, (((p+ip)->back)+iw1-2)->pos, (((p+ip)->back)+iw1-1)->pos, (((p+ip)->back)+iw1)->pos, (p+ip)->tables, &out)
						- (((p+ip)->back)+iw1-1)->d_next ) > parms->d_cloose || out==0 ) ok=0;
		if (ok && iw2+3<(p+ip)->nback && (((p+ip)->back)+iw2+2)->move == 0)
			if ( DAbs( Dihedral( (((p+ip)->back)+iw2)->pos, (((p+ip)->back)+iw2+1)->pos, (((p+ip)->back)+iw2+2)->pos, (((p+ip)->back)+iw2+3)->pos, (p+ip)->tables, &out)
						- (((p+ip)->back)+iw2+2)->d_next ) > parms->d_cloose || out==0) ok=0;
	}

	if (ok) ok = AddSidechain(p,iw1-1,iw2+1,ip);


	// calculate new energy
	if (ok)
	{
		deltaE += EnergyMonomerRange(p,pot,iw1-1,iw2+1,ip,parms->npol,parms->shell,1,parms->nosidechains,parms->disentangle,parms->hb);		// new energy of iw with the others
		if (!parms->noangpot)
		{
			if (iw1>0) deltaE += EnergyAngles(p,pot,iw1-1,ip,1);
			if (iw2<(p+ip)->nback-1) deltaE += EnergyAngles(p,pot,iw2+1,ip,1);
		}
		if (!parms->nodihpot)
		{
			if (iw1>0) deltaE += EnergyDihedrals(p,pot,iw1-1,ip,1);
			deltaE += EnergyDihedrals(p,pot,iw1,ip,1);
			if (iw2<(p+ip)->nback-1) deltaE += EnergyDihedrals(p,pot,iw2+1,ip,1);
			if (iw2<(p+ip)->nback-2) deltaE += EnergyDihedrals(p,pot,iw2+2,ip,1);
		}

		#ifdef DEBUG
		if (debug>2) fprintf(stderr,"deltaE=%lf\n",deltaE); fflush(stderr);
		#endif
	}
	// Metropolis acceptance
	if (ok) ok = Metropolis(deltaE,t,p->tables);

	if (ok == 0) 			// move rejected
	{
		UpdateMonomerRange(oldp,p,iw1-1,iw2+1,ip,parms->shell);							// return to the old position, contacts, etc.
	}
	else					// move accepted
	{
		UpdateMonomerRange(p,oldp,iw1-1,iw2+1,ip,parms->shell);							// update oldp
		p->etot += deltaE;
		parms->acc ++;

		#ifdef DEBUG
			if (debug>2) fprintf(stderr,"accept=%d\n",ok); fflush(stderr);
		#endif
	}


	//fprintf(stderr,"** %lf %lf\n",TotalEnergy(p,pot,parms,parms->npol,0,parms->nosidechains,parms->debug),TotalEnergy(oldp,pot,parms,parms->npol,0,parms->nosidechains,parms->debug));
	parms->mov ++;

	if (ok == 0) return 0;
	else return 1;
}

/*****************************************************************************
 Moves center of mass of a chain
	 returns 1 if it could make the move, 0 if nothing was done
 *****************************************************************************/
int MoveCoM(struct s_polymer *p,struct s_polymer *oldp, struct s_mc_parms *parms, struct s_potential *pot,
		unsigned long long istep, int debug, double t)
{
	int ip,ok;
	double dx,dy,dz,deltaE;

	ip = irand(parms->npol);
	dx = parms->dx_com * frand() - parms->dx_com / 2;
	dy = parms->dx_com * frand() - parms->dx_com / 2;
	dz = parms->dx_com * frand() - parms->dx_com / 2;

	#ifdef DEBUG
	 if (debug>2) fprintf(stderr,"step=%d com ip=%d\n",istep,ip); fflush(stderr);
	#endif

	// old energies
	deltaE = - GetEnergyMonomerRange(p,0,(p+ip)->nback-1,ip);						// old energy  (iw1-1 & iw2+1 because I move their sidechains)
	
	// move
	DisplaceCoM(p,ip,dx,dy,dz);

	// new energy
	deltaE += EnergyMonomerRange(p,pot,0,(p+ip)->nback-1,ip,parms->npol,0,1,parms->nosidechains,parms->disentangle,parms->hb);

 	#ifdef DEBUG
  	   if (debug>2) fprintf(stderr,"deltaE=%lf\n",deltaE); fflush(stderr);
  	#endif

	// Metropolis acceptance
	ok = Metropolis(deltaE,t,p->tables);

	if (ok == 0) 				// move rejected
	{
		UpdateMonomerRange(oldp,p,0,(p+ip)->nback-1,ip,0);							// return to the old position, contacts, etc.
	}
	else					// move accepted
	{
		UpdateMonomerRange(p,oldp,0,(p+ip)->nback-1,ip,0);							// update oldp
	
                if(parms->shell==1 )            // if shells are active, update them 
                {
                        UpdateShell(p,parms);
                        CopyShell(p,oldp,parms);
                        parms->ishell=0;
                }




		p->etot += deltaE;
		parms->acc ++;

		#ifdef DEBUG
			if (debug>2) fprintf(stderr,"accept=%d\n",ok); fflush(stderr);
		#endif
	}

	parms->mov ++;

	if (ok == 0) return 0;
	else return 1;
}

int MoveRotation(struct s_polymer *p,struct s_polymer *oldp, struct s_mc_parms *parms, struct s_potential *pot,
                unsigned long long istep, int debug, double t)
{
	int ip,ok,idir;
	double dtheta,deltaE;
	ip=irand(parms->npol);

	dtheta= parms->dtheta * (0.5 - frand());


	deltaE = - GetEnergyMonomerRange(p,0,(p+ip)->nback-1,ip);

	idir=irand(3);
	
	if(idir==0)
		RotationX(p,ip,dtheta);
	else if(idir==1)
		RotationY(p,ip,dtheta);
	else
		RotationZ(p,ip,dtheta);
	
	deltaE += EnergyMonomerRange(p,pot,0,(p+ip)->nback-1,ip,parms->npol,0,1,parms->nosidechains,parms->disentangle,parms->hb);

	ok = Metropolis(deltaE,t,p->tables);

        if (ok == 0)                            // move rejected
        {
                UpdateMonomerRange(oldp,p,0,(p+ip)->nback-1,ip,0);                                                     // return to the old position, contacts, etc.
        }
        else                                    // move accepted
        {

                UpdateMonomerRange(p,oldp,0,(p+ip)->nback-1,ip,0); 
       		if(parms->shell==1 )            // if shells are active, update them 
        	{
                	UpdateShell(p,parms);
                	CopyShell(p,oldp,parms);
                        parms->ishell=0;
        	}


       		p->etot += deltaE;
        	parms->acc ++;

        	#ifdef DEBUG
	        	if (debug>2) fprintf(stderr,"accept=%d\n",ok); fflush(stderr);
        	#endif
        }

        parms->mov++;

	return ok;
}


int MoveClusterCoM(struct s_polymer *p,struct s_polymer *oldp, struct s_mc_parms *parms, struct s_potential *pot,
                  int istep, int debug, double t, FILE *fproc, int iproc)
{
      int ip,ok,ic,jc;
      int k,i;//,ind;
      int ctrl, old_ctrl;
      int ncluster,ind_cluster,cont,sum;
      double dx,dy,dz,deltaE;
      int adjacency[NCHAINMAX][NCHAINMAX];
      int cluster[NCHAINMAX][NCHAINMAX];
      int avector[NCHAINMAX];
      int npol_cluster[NCHAINMAX];

      for(ic=0;ic<NCHAINMAX;ic++)
      {
            for(jc=0;jc<NCHAINMAX;jc++)
	    {
            	adjacency[ic][jc]=0;
                cluster[ic][jc]=-1;
            }
      }

      for(ic=0;ic<parms->npol;ic++)
      {
           for (i=0;i<(p+ic)->nback;i++)
	   {
           	for (k=0;k<(((p+ic)->back)+i)->ncontacts;k++)
		{
                        adjacency[ic][((((p+ic)->back)+i)->contacts_p)[k]] = 1;
                }
            }
      }

      cont=0;
      ncluster=0;
      while(1)
      {
            

		ctrl=99;
            for(ic=0;ic<parms->npol;ic++)
                  avector[ic]=adjacency[cont][ic];
            
	    ind_cluster=0;
            old_ctrl=0;


            while( (old_ctrl-ctrl) != 0)
	    {

                  old_ctrl=ctrl;
                  ctrl=0;
                  for(ic=cont;ic<parms->npol;ic++)
		  {

                        sum=0;
                        for(jc=0;jc<parms->npol;jc++)
                              sum += avector[jc]*adjacency[ic][jc];
                        if(sum!=0)
			{
                              cluster[ncluster][ind_cluster] = ic;
                              ind_cluster++;
                              ctrl++;
                              for(jc=0;jc<parms->npol;jc++)
			      {
                                    if(avector[jc]+adjacency[ic][jc]>=1)
				    	avector[jc]=1;

                                    adjacency[ic][jc]=0;
                              }
                        }
                  }
            }

            npol_cluster[ncluster]=ind_cluster;
            ncluster++;
            cont++;

            for(i=cont;i<parms->npol;i++)
	    {
                  if(adjacency[i][i] == 0)
		  	cont++;
                  else 
		  {
                        cont=i;
                        break;
                  }
            }
            if(cont==parms->npol) break;
      }




      ic = irand(ncluster);

      dx = parms->dx_clm * frand() - parms->dx_clm / 2;
      dy = parms->dx_clm * frand() - parms->dx_clm / 2;
      dz = parms->dx_clm * frand() - parms->dx_clm / 2;

      deltaE=0.;
      for(i=0;i<npol_cluster[ic];i++)
      {

            ip = cluster[ic][i];

            #ifdef DEBUG
            if (debug>2) fprintf(fproc,"step=%d com ip=%d\n",istep,ip); fflush(stderr);
            #endif



            deltaE -= GetEnergyMonomerRange(p,0,(p+ip)->nback-1,ip);                                    // old energy  (iw1-1 & iw2+1 because I move their sidechains)
            DisplaceCoM(p,ip,dx,dy,dz);

            //deltaE += EnergyMonomerRange(p,pot,0,(p+ip)->nback,ip,parms->npol,parms->shell,1,parms->nosidechains,parms->disentangle,parms->hb);
	    deltaE += EnergyMonomerRange(p,pot,0,(p+ip)->nback-1,ip,parms->npol,0,1,parms->nosidechains,parms->disentangle,parms->hb);


      }

      #ifdef DEBUG
      if (debug>2) fprintf(fproc,"deltaE=%lf\n",deltaE); fflush(stderr);
      #endif

      ok = Metropolis(deltaE,t,p->tables);

      if (ok == 0)                        // move rejected
      {
            for(i=0;i<npol_cluster[ic];i++)
	    {
                  ip = cluster[ic][i];
                  //UpdateMonomerRange(oldp,p,0,(p+ip)->nback,ip,parms->shell);                                     // return to the old position, contacts, etc.
		  UpdateMonomerRange(oldp,p,0,(p+ip)->nback-1,ip,0);    
            }
      }

      else                          // move accepted
      {
            for(i=0;i<npol_cluster[ic];i++)
	    {
                  ip = cluster[ic][i];
                  //UpdateMonomerRange(p,oldp,0,(p+ip)->nback,ip,parms->shell);                                     // update oldp
                  UpdateMonomerRange(p,oldp,0,(p+ip)->nback-1,ip,0);

                  if(parms->shell==1 )            // if shells are active, update them 
                  {
                        UpdateShell(p,parms);
                        CopyShell(p,oldp,parms);
                        parms->ishell=0;
                  }
  


		  p->etot += deltaE;
                  parms->acc ++;

                  #ifdef DEBUG
                  if (debug>2) fprintf(fproc,"accept=%d\n",ok); fflush(stderr);
                  #endif
            }
      }

      parms->mov ++;

      if (ok == 0) return 0;
      else return 1;

}


int MoveClusterRot(struct s_polymer *p,struct s_polymer *oldp, struct s_mc_parms *parms, struct s_potential *pot,
                  int istep, int npol,int debug, double t, FILE *fproc, int iproc)
{
      int ip,ok,ic,jc;
      int k,i;//,ind;
      int ctrl, old_ctrl;
      int ncluster,ind_cluster,cont,sum;
      //double dx,dy,dz
      double deltaE;
      int adjacency[NCHAINMAX][NCHAINMAX];
      int cluster[NCHAINMAX][NCHAINMAX];
      int avector[NCHAINMAX];
      int npol_cluster[NCHAINMAX];
      double dtheta= parms->dtheta * (0.5 - frand());



      for(ic=0;ic<NCHAINMAX;ic++)
      {
            for(jc=0;jc<NCHAINMAX;jc++)
	    {
                  adjacency[ic][jc]=0;
                  cluster[ic][jc]=-1;
            }
      }

      for(ic=0;ic<parms->npol;ic++)
      {
            for (i=0;i<(p+ic)->nback;i++)
	    {
                  for (k=0;k<(((p+ic)->back)+i)->ncontacts;k++)
		  {
                        adjacency[ic][((((p+ic)->back)+i)->contacts_p)[k]] = 1;
                  }
            }
      }

      cont=0;
      ncluster=0;
      while(1)
      {
            ctrl=99;
            for(ic=0;ic<NCHAINMAX;ic++)
	    {
                  avector[ic]=adjacency[cont][ic];
            }

            ind_cluster=0;
            old_ctrl=0;



            while(old_ctrl-ctrl != 0)
	    {
	
                  old_ctrl=ctrl;
                  ctrl=0;
                  for(ic=cont;ic<parms->npol;ic++)
		  {

                        sum=0;
                        for(jc=0;jc<parms->npol;jc++)
                        	sum += avector[jc]*adjacency[ic][jc];
                        if(sum!=0)
			{


                              cluster[ncluster][ind_cluster] = ic;
                              ind_cluster++;
                              ctrl++;
                              for(jc=0;jc<parms->npol;jc++)
			      {
                                    if(avector[jc]+adjacency[ic][jc]>=1) avector[jc]=1;
                                    adjacency[ic][jc]=0;
                              }
                        }
                  }
            }

            npol_cluster[ncluster]=ind_cluster;
            ncluster++;
            cont++;
            for(i=cont;i<parms->npol;i++)
	    {
                  if(adjacency[i][i] == 0)
		  	cont++;
                  else 
		  {
                        cont=i;
                        break;
                  }
            }

            if(cont==parms->npol) break;
      }



      ic = irand(ncluster);

     deltaE=0.;



      for(i=0;i<npol_cluster[ic];i++)
      {

            ip = cluster[ic][i];

            #ifdef DEBUG
            if (debug>2) fprintf(fproc,"step=%d com ip=%d\n",istep,ip); fflush(stderr);
            #endif

            deltaE -= GetEnergyMonomerRange(p,0,(p+ip)->nback-1,ip);                                    // old energy  (iw1-1 & iw2+1 because I move their sidechains)

//            deltaE += EnergyMonomerRange(p,pot,0,(p+ip)->nback,ip,parms->npol,parms->shell,1,parms->nosidechains,parms->disentangle,parms->hb);
 //	    deltaE += EnergyMonomerRange(p,pot,0,(p+ip)->nback,ip,parms->npol,0,1,parms->nosidechains,parms->disentangle,parms->hb);



      }


	int idir=irand(3);

        if(idir==0)
                RotationClusterX(p,ip,dtheta,ic,cluster,npol_cluster[ic]);
        else if(idir==1)
                RotationClusterY(p,ip,dtheta,ic,cluster,npol_cluster[ic]);
        else
                RotationClusterZ(p,ip,dtheta,ic,cluster,npol_cluster[ic]);


	
	for(i=0;i<npol_cluster[ic];i++)
	{
	ip=cluster[ic][i];
	//deltaE += EnergyMonomerRange(p,pot,0,(p+ip)->nback,ip,parms->npol,parms->shell,1,parms->nosidechains,parms->disentangle,parms->hb);
	deltaE += EnergyMonomerRange(p,pot,0,(p+ip)->nback-1,ip,parms->npol,0,1,parms->nosidechains,parms->disentangle,parms->hb);
	}




      #ifdef DEBUG
      if (debug>2) fprintf(fproc,"deltaE=%lf\n",deltaE); fflush(stderr);
      #endif

      ok = Metropolis(deltaE,t,p->tables);

      if (ok == 0)                        // move rejected
      {
            for(i=0;i<npol_cluster[ic];i++)
	    {
                  ip = cluster[ic][i];
                //  UpdateMonomerRange(oldp,p,0,(p+ip)->nback,ip,parms->shell);                                     // return to the old position, contacts, etc.
                  UpdateMonomerRange(oldp,p,0,(p+ip)->nback-1,ip,0);    
	    }
      }
      else                          // move accepted
      {
            for(i=0;i<npol_cluster[ic];i++)
	    {
                  ip = cluster[ic][i];
 //                 UpdateMonomerRange(p,oldp,0,(p+ip)->nback,ip,parms->shell);                                     // update oldp
		 UpdateMonomerRange(p,oldp,0,(p+ip)->nback-1,ip,0); 
                 if(parms->shell==1 )            // if shells are active, update them 
                 {
                        UpdateShell(p,parms);
                        CopyShell(p,oldp,parms);
                        parms->ishell=0;
                 }
 




                  p->etot += deltaE;
                  parms->acc ++;

                  #ifdef DEBUG
                  if (debug>2) fprintf(fproc,"accept=%d\n",ok); fflush(stderr);
                  #endif
            }
      }

      parms->mov ++;

      if (ok == 0) return 0;
      else return 1;

}

 



/*****************************************************************************
le(old_ctrl-ctrl != 0){
                  old_ctrl=ctrl;
                  ctrl=0;
                  for(ic=cont;ic<parms->npol;ic++){
 Changes temperature to anneal the system 
	status: 0=normal temperature, 1=higher temperature, 2=normal temperature
		but equilibration
	returns 1 to reset the counter
 *****************************************************************************/
/*
void Anneal(struct s_mc_parms *p, double *t, int *counter, int *status, int *ok, int *ishell, int *mcount)
{
	//int i;

	(*counter) ++;
	//for (i=0;i<NMOVES;i++) mcount[i] ++;
	(*ishell) ++;

	// normal temperature
	if ((*status)==0)
	{
		if ((*counter) >= p->anneal_often)
		{
			(*status) = 1;
			(*t) = p->anneal_t;
			(*counter) =0;
			fprintf(stderr,"ANNEALING: increasing temperature to %lf (%d steps)\n",(*t), p->anneal_step);
			return;
		}
	}
	// higher temperature
	else if ((*status)==1)
	{
		(*ok) = -2;
		if ((*counter) >= p->anneal_step)
		{
			(*status) = 2;
			(*t) = p->T;
			(*counter) =0;
			fprintf(stderr,"ANNEALING: annealing temperature to %lf and re-equilibrating (%d steps))\n",(*t),p->anneal_recov);
			return;
		}
	}
	// re-equilibrate
	else if ((*status)==2)
	{
		(*ok) = -2;
		if ((*counter) >= p->anneal_recov)
		{
			(*status) = 0;
			(*counter) =0;
			fprintf(stderr,"ANNEALING: keeping on with actual simulation (T=%lf)\n",(*t));
			return;
		}
	}

}
*/
