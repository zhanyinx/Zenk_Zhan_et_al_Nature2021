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
#include "MPIfunc.h"

#ifdef ACTIVE_MPI

/***********************************************
		MPI Datatypes
***********************************************/

void Create_vector_datatype(MPI_Datatype *Vectortype)
{
	struct vector *x;
	x=calloc(1,sizeof(struct vector));
	MPI_Aint adress[4];
	MPI_Get_address(x, &adress[0]);
        MPI_Get_address(&(*x).x, &adress[1]);
        MPI_Get_address(&(*x).y, &adress[2]);
        MPI_Get_address(&(*x).z, &adress[3]);

	MPI_Datatype type[3]={MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE};

	int blocklen[3]={1,1,1};
        MPI_Aint disp[3];

        int i;
        for(i=0; i<3; i++) {disp[i]=adress[i+1]-adress[0];}

        MPI_Type_create_struct(3,blocklen,disp,type,Vectortype);
        MPI_Type_commit(Vectortype);
        free(x);

        return;


}


void Create_rot_datatype(MPI_Datatype *Rottype)
{
	struct s_rotamers *x;
	x = calloc(1,sizeof(struct s_rotamers));
	if (!x) Error("Cannot allocate struct polymer->back->side->rot");

	MPI_Aint adress[5];
	MPI_Get_address(x, &adress[0]);
	MPI_Get_address(&(*x).b1, &adress[1]);
	MPI_Get_address(&(*x).b2, &adress[2]);
	MPI_Get_address(&(*x).b3, &adress[3]);
	MPI_Get_address(&(*x).ang, &adress[4]);
	
	MPI_Datatype type[4]={MPI_INT, MPI_INT, MPI_INT, MPI_DOUBLE};
	int blocklen[4]={1,1,1,3};
	MPI_Aint disp[4];
	
	int i;
	for(i=0; i<4; i++) {disp[i]=adress[i+1]-adress[0];}
	
	MPI_Type_create_struct(4,blocklen,disp,type,Rottype);
	MPI_Type_commit(Rottype);
	free(x);

	return;
}


void Create_side_datatype(MPI_Datatype *Sidetype)
{
	struct s_side *x;
	x = calloc(1,sizeof(struct s_side));
	if (!x) Error("Cannot allocate struct polymer->back->side");

	MPI_Aint adress[5];
	MPI_Get_address(x, &adress[0]);
	MPI_Get_address(&(*x).ia, &adress[1]);
	MPI_Get_address(&(*x).itype, &adress[2]);
	MPI_Get_address(&(*x).type, &adress[3]);
	MPI_Get_address(&(*x).pos, &adress[4]);
	
	MPI_Datatype type[4]={MPI_INT, MPI_INT, MPI_CHAR, MPI_DOUBLE};
	int blocklen[4]={1,1,5,3};
	MPI_Aint disp[4];
	
	int i;
	for(i=0; i<4; i++) {disp[i]=adress[i+1]-adress[0];}
	
	MPI_Type_create_struct(4,blocklen,disp,type,Sidetype);
	MPI_Type_commit(Sidetype);
	free(x);

	return;
}


void Create_back_datatype(MPI_Datatype *Backtype)
{
	struct s_back *x;
	x = calloc(1,sizeof(struct s_back));
	if (!x) Error("Cannot allocate struct polymer->back");

	MPI_Aint adress[12];
	MPI_Get_address(x, &adress[0]);
	MPI_Get_address(&(*x).ia, &adress[1]);
	MPI_Get_address(&(*x).itype, &adress[2]);
	MPI_Get_address(&(*x).aa, &adress[3]);
	MPI_Get_address(&(*x).type, &adress[4]);
	MPI_Get_address(&(*x).iaa, &adress[5]);
	MPI_Get_address(&(*x).pos, &adress[6]);
	MPI_Get_address(&(*x).nside, &adress[7]);
	MPI_Get_address(&(*x).move, &adress[8]);
	MPI_Get_address(&(*x).irot, &adress[9]);
	MPI_Get_address(&(*x).nrot, &adress[10]);
	MPI_Get_address(&(*x).iapdb,&adress[11]);


	MPI_Datatype type[11]={MPI_INT, MPI_INT, MPI_CHAR, MPI_CHAR, MPI_INT, MPI_DOUBLE, MPI_INT, MPI_INT, MPI_INT, MPI_INT,MPI_INT};
	int blocklen[11]={1,1,5,5,1,3,1,1,1,1,1};
	MPI_Aint disp[11];
	
	int i;
	for(i=0; i<11; i++) {disp[i]=adress[i+1]-adress[0];}
	
	MPI_Type_create_struct(11,blocklen,disp,type,Backtype);
	MPI_Type_commit(Backtype);
	free(x);

	return;
}


void Create_pot_datatype(MPI_Datatype *Pottype)
{
	struct s_potential *x;
	x = calloc(1,sizeof(struct s_potential));
	if (!x) Error("Cannot allocate struct potential");

	MPI_Aint adress[24];
	MPI_Get_address(x, &adress[0]);
	MPI_Get_address(&(*x).dih_periodic, &adress[1]);
	MPI_Get_address(&(*x).dih_tabled, &adress[2]);
	MPI_Get_address(&(*x).g_r0hard, &adress[3]);
	MPI_Get_address(&(*x).g_ehomo, &adress[4]);
	MPI_Get_address(&(*x).g_rhomo, &adress[5]);
	MPI_Get_address(&(*x).g_anglek, &adress[6]);
	MPI_Get_address(&(*x).g_angle0, &adress[7]);
	MPI_Get_address(&(*x).g_dihk1, &adress[8]);
	MPI_Get_address(&(*x).g_dih01, &adress[9]);
	MPI_Get_address(&(*x).g_dihk3, &adress[10]);
	MPI_Get_address(&(*x).g_dih03, &adress[11]);
	MPI_Get_address(&(*x).g_imin, &adress[12]);
	MPI_Get_address(&(*x).g_dihbin, &adress[13]);
	MPI_Get_address(&(*x).g_dihe, &adress[14]);
	MPI_Get_address(&(*x).boxtype, &adress[15]);
	MPI_Get_address(&(*x).boxsize, &adress[16]);
	MPI_Get_address(&(*x).lbox, &adress[17]);
	MPI_Get_address(&(*x).kr2_splice, &adress[18]);
	MPI_Get_address(&(*x).ke_splice, &adress[19]);
	MPI_Get_address(&(*x).splice, &adress[20]);
	MPI_Get_address(&(*x).hc_number, &adress[21]);	
	MPI_Get_address(&(*x).dih_ram, &adress[22]);	
	MPI_Get_address(&(*x).e_dihram, &adress[23]);	

	MPI_Datatype type[23]={MPI_INT, MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_CHAR, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_INT, MPI_INT, MPI_INT, MPI_DOUBLE};
	int blocklen[23]={1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
	MPI_Aint disp[23];
	
	int i;
	for(i=0; i<23; i++) {disp[i]=adress[i+1]-adress[0];}
	
	MPI_Type_create_struct(23,blocklen,disp,type,Pottype);
	MPI_Type_commit(Pottype);
	free(x);

	return;
}


void Create_parms_datatype(MPI_Datatype *Parmstype){

	struct s_mc_parms *x;
	x = calloc(1,sizeof(struct s_mc_parms));
	if (!x) Error("Cannot allocate mc_parms");

	MPI_Aint adress[68];
	MPI_Get_address(x, &adress[0]);
	MPI_Get_address(&(*x).npol, &adress[1]);
	MPI_Get_address(&(*x).nstep, &adress[2]);
	MPI_Get_address(&(*x).seed, &adress[3]);
	MPI_Get_address(&(*x).dw_flip, &adress[4]);
	MPI_Get_address(&(*x).dw_pivot, &adress[5]);
	MPI_Get_address(&(*x).dw_mpivot, &adress[6]);
	MPI_Get_address(&(*x).dw_lpivot, &adress[7]);
	MPI_Get_address(&(*x).dw_mflip, &adress[8]);
	MPI_Get_address(&(*x).fntrj, &adress[9]);
	MPI_Get_address(&(*x).fne, &adress[10]);
	MPI_Get_address(&(*x).flastp, &adress[11]);
	MPI_Get_address(&(*x).fnproc, &adress[12]);
	MPI_Get_address(&(*x).nprinttrj, &adress[13]);
	MPI_Get_address(&(*x).nprintlog, &adress[14]);
	MPI_Get_address(&(*x).nprinte, &adress[15]);
	MPI_Get_address(&(*x).shell, &adress[16]);
	MPI_Get_address(&(*x).nshell, &adress[17]);
	MPI_Get_address(&(*x).r2shell, &adress[18]);
	MPI_Get_address(&(*x).ntemp, &adress[19]);
	MPI_Get_address(&(*x).T, &adress[20]);
	MPI_Get_address(&(*x).randdw, &adress[21]);
	MPI_Get_address(&(*x).debug, &adress[22]);
	MPI_Get_address(&(*x).movetype, &adress[23]);
	MPI_Get_address(&(*x).nmul_mpivot, &adress[24]);
	MPI_Get_address(&(*x).nmul_lpivot, &adress[25]);
	MPI_Get_address(&(*x).nmul_mflip, &adress[26]);
	MPI_Get_address(&(*x).nosidechains, &adress[27]);
	MPI_Get_address(&(*x).noangpot, &adress[28]);
	MPI_Get_address(&(*x).nodihpot, &adress[29]);
	MPI_Get_address(&(*x).nrun, &adress[30]);
	MPI_Get_address(&(*x).always_restart, &adress[31]);
	MPI_Get_address(&(*x).record_native, &adress[32]);
	MPI_Get_address(&(*x).acc, &adress[33]);
	MPI_Get_address(&(*x).mov, &adress[34]);
	MPI_Get_address(&(*x).disentangle, &adress[35]);
	MPI_Get_address(&(*x).stempering, &adress[36]);
	MPI_Get_address(&(*x).dx_com, &adress[37]);
	MPI_Get_address(&(*x).dx_clm, &adress[38]);
	MPI_Get_address(&(*x).r_cloose, &adress[39]);
	MPI_Get_address(&(*x).a_cloose, &adress[40]);
	MPI_Get_address(&(*x).d_cloose, &adress[41]);
	MPI_Get_address(&(*x).hb, &adress[42]);
	MPI_Get_address(&(*x).anneal, &adress[43]);
	MPI_Get_address(&(*x).anneal_often, &adress[44]);
	MPI_Get_address(&(*x).anneal_step, &adress[45]);
	MPI_Get_address(&(*x).anneal_t, &adress[46]);
	MPI_Get_address(&(*x).anneal_recov, &adress[47]);
	#ifdef OPTIMIZEPOT
	MPI_Get_address(&(*x).op_minim, &adress[48]);
	MPI_Get_address(&(*x).op_itermax, &adress[49]);
	MPI_Get_address(&(*x).op_step, &adress[50]);
	MPI_Get_address(&(*x).op_T, &adress[51]);
	MPI_Get_address(&(*x).op_deltat, &adress[52]);
	MPI_Get_address(&(*x).op_stop, &adress[53]);
	MPI_Get_address(&(*x).op_print, &adress[54]);
	MPI_Get_address(&(*x).op_emin, &adress[55]);
	MPI_Get_address(&(*x).op_emax, &adress[56]);
	MPI_Get_address(&(*x).op_wait, &adress[57]);
	MPI_Get_address(&(*x).op_r, &adress[58]);	
	MPI_Get_address(&(*x).op_r0, &adress[59]);
	MPI_Get_address(&(*x).nstep_exchange, &adress[60]);
	MPI_Get_address(&(*x).nmul_local,&adress[61]);
	MPI_Get_address(&(*x).chi2start,&adress[62]);
	MPI_Get_address(&(*x).ishell,&adress[63]);
	MPI_Get_address(&(*x).bgs_a,&adress[64]);
	MPI_Get_address(&(*x).bgs_b,&adress[65]);
	MPI_Get_address(&(*x).dtheta,&adress[66]);

	MPI_Get_address(&(*x).iT_bias,&adress[67]);



	#else
	MPI_Get_address(&(*x).nstep_exchange, &adress[48]);
	#endif
	

	
	#ifdef OPTIMIZEPOT
	MPI_Datatype type[67]={MPI_INT, MPI_INT, MPI_LONG, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_CHAR, MPI_CHAR, MPI_CHAR, MPI_CHAR, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_DOUBLE, MPI_INT, MPI_DOUBLE, MPI_INT, MPI_INT, MPI_INT, MPI_INT,MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT,  MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT,  MPI_DOUBLE,  MPI_DOUBLE,  MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_INT,  MPI_INT,  MPI_INT,  MPI_INT, MPI_DOUBLE, MPI_INT, MPI_CHAR, MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_INT, MPI_DOUBLE, MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_INT, MPI_INT,MPI_INT,MPI_INT,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_INT};
	int blocklen[67]={1,1,1,1,1,1,1,1,50,50,50,50,1,1,1,1,1,1,1,NREPMAX,1,1,NMOVES,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,50,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
	MPI_Aint disp[67];
	
	int i;
	for(i=0; i<67; i++) {disp[i]=adress[i+1]-adress[0];}
	
	MPI_Type_create_struct(67,blocklen,disp,type,Parmstype);
	MPI_Type_commit(Parmstype);
	free(x);
	#else
	MPI_Datatype type[48]={MPI_INT, MPI_INT, MPI_LONG, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_CHAR, MPI_CHAR, MPI_CHAR, MPI_CHAR, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_DOUBLE, MPI_INT, MPI_DOUBLE, MPI_INT, MPI_INT, MPI_INT, MPI_INT,MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT,  MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT,  MPI_DOUBLE,  MPI_DOUBLE,  MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_INT,  MPI_INT,  MPI_INT,  MPI_INT, MPI_DOUBLE, MPI_INT, MPI_INT};
	int blocklen[48]={1,1,1,1,1,1,1,1,50,50,50,50,1,1,1,1,1,1,1,NREPMAX,1,1,NMOVES,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
	MPI_Aint disp[48];
	
	int i;
	for(i=0; i<48; i++) {disp[i]=adress[i+1]-adress[0];}
	
	MPI_Type_create_struct(48,blocklen,disp,type,Parmstype);
	MPI_Type_commit(Parmstype);
	free(x);
	#endif
}



/***********************************************
		Send/Receive ops
***********************************************/

//parameters 
struct s_mc_parms *send_parms(int iproc, int nprocs, MPI_Datatype Parmstype, struct s_mc_parms *parms, MPI_Status astatus)//,double *tmin)
{

	int i;

	if(iproc!=0)
	{
		parms = calloc(1,sizeof(struct s_mc_parms));
		if (!parms) Error("Cannot allocate mc_parms");
	}
	
	if(iproc==0) for(i=1; i<nprocs; i++) MPI_Send(parms, 1, Parmstype, i, 100+i, MPI_COMM_WORLD);
	if(iproc!=0) MPI_Recv(parms, 1, Parmstype, 0, 100+iproc, MPI_COMM_WORLD, &astatus);

	return parms;
}


//number of atoms, types and backbone atoms
void send_struct(int *nback, int iproc, int nprocs, int *nat, int *ntypes, MPI_Status astatus)
{
	int buffer_size, position,i;
	char buffer[buffer_max];
	buffer_size = 3*sizeof(int);
	if(buffer_size>buffer_max)
	{
		fprintf(stderr,"Buffer too small\n");
		MPI_Finalize();
		exit(1);
	}
	
	if(iproc==0) 
	{
		position = 0;	
		MPI_Pack(nat,1,MPI_INT,buffer,buffer_size,&position,MPI_COMM_WORLD);
		MPI_Pack(ntypes,1,MPI_INT,buffer,buffer_size,&position,MPI_COMM_WORLD);
		MPI_Pack(nback,1,MPI_INT,buffer,buffer_size,&position,MPI_COMM_WORLD);
		
		for(i=1;i<nprocs;i++) MPI_Send(buffer,position,MPI_PACKED,i, 200+i, MPI_COMM_WORLD);
	}
	if(iproc!=0)
	{
		MPI_Recv(buffer,buffer_size,MPI_PACKED,0, 200+iproc, MPI_COMM_WORLD, &astatus);
		position=0;
		MPI_Unpack(buffer,buffer_size,&position,nat,1,MPI_INT,MPI_COMM_WORLD);
		MPI_Unpack(buffer,buffer_size,&position,ntypes,1,MPI_INT,MPI_COMM_WORLD);
		MPI_Unpack(buffer,buffer_size,&position,nback,1,MPI_INT,MPI_COMM_WORLD);
	}

	return;
}

//polymer status
struct s_polymer *send_pol(int iproc, int nprocs, int nback, MPI_Datatype Backtype, MPI_Datatype Sidetype,  MPI_Datatype Rottype, struct s_polymer *startp, MPI_Status astatus, int npol, int shell, int nosidechains)
{
	int i, j, k, l,position,buffer_size;

	char buffer[buffer_max];
	(startp+npol)->nback = nback;
	buffer_size=(sizeof(struct s_back)*nback);
//	fprintf(stderr,"send_pol: BUFFER SIZE IS %d\n",buffer_size);
	if(buffer_size>buffer_max)
        {
        	fprintf(stderr,"Buffer too small\n");
                MPI_Finalize();
                exit(1);
        }

	if(iproc==0) 
	{
		position=0;
		for(j=0;j<nback;j++)
			MPI_Pack(((startp+npol)->back)+j,1,Backtype,buffer,buffer_size,&position,MPI_COMM_WORLD);
		for(i=1; i<nprocs; i++)
			MPI_Send(buffer,position,MPI_PACKED,i, 300+100*i, MPI_COMM_WORLD);
	}

	if(iproc!=0)
	{
		MPI_Recv(buffer,buffer_size,MPI_PACKED,0, 300+100*iproc, MPI_COMM_WORLD, &astatus);		
		position=0;
		for(j=0;j<nback;j++)
			MPI_Unpack(buffer,buffer_size,&position,((startp+npol)->back)+j,1,Backtype,MPI_COMM_WORLD);	
	}

	if(!(nosidechains))
        {
		if(iproc==0)
		{
			for(i=1; i<nprocs; i++) for(j=0; j<nback;j++) for(k=0; k<((startp->back)+j)->nside;k++) 
					MPI_Send(((((startp+npol)->back)+j)->side)+k, 1, Sidetype, i, 2000+100*i+10*j+k, MPI_COMM_WORLD);
		}	

	if(iproc!=0) for(j=0; j<nback;j++) for(k=0; k<((startp->back)+j)->nside;k++)
					MPI_Recv(((((startp+npol)->back)+j)->side)+k, 1, Sidetype, 0, 2000+100*iproc+10*j+k, MPI_COMM_WORLD, &astatus);
		if(iproc==0) for(i=1; i<nprocs; i++) for(j=0; j<nback;j++) for(k=0; k<((startp->back)+j)->nside;k++) for(l=0; l<((startp->back)+j)->nrot; l++)
					MPI_Send(((((((startp+npol)->back)+j)->side)+k)->rot)+l, 1, Rottype, i, 10000+1000*i+100*j+10*k+l, MPI_COMM_WORLD);
		if(iproc!=0) for(j=0; j<nback;j++) for(k=0; k<((startp->back)+j)->nside;k++) for(l=0; l<((startp->back)+j)->nrot; l++)
					MPI_Recv(((((((startp+npol)->back)+j)->side)+k)->rot)+l, 1, Rottype, 0, 10000+1000*iproc+100*j+10*k+l , MPI_COMM_WORLD, &astatus);	
	}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	if(iproc!=0)
	{
		for(i=0; i<nback; i++) ((startp+npol)->vback)[(((startp+npol)->back)+i)->ia ] = &(((((startp+npol)->back)+i)->pos));
		for(i=0; i<nback; i++) for(j=0; j<(((startp+npol)->back)+i)->nside; j++) ((startp+npol)->vback)[ (((((startp+npol)->back)+i)->side)+j)->ia ] = &((((((startp+npol)->back)+i)->side)+j)->pos);
	}
		
	return startp;
}
//potential
struct s_potential *send_pot(int nat, int ntypes, int noangpot, int nodihpot, int hb, int iproc, int nprocs, MPI_Datatype Pottype, struct s_potential *u, MPI_Status astatus)
{
	int i;
	
	if(iproc==0) for(i=1; i<nprocs; i++) MPI_Send(u, 1, Pottype, i, 30000+i, MPI_COMM_WORLD);	
	if(iproc!=0) MPI_Recv(u, 1, Pottype, 0, 30000+iproc, MPI_COMM_WORLD, &astatus);	
	
	send_double_matrix(ntypes, ntypes, iproc, u->e, 0);
	send_double_matrix(ntypes, ntypes, iproc, u->r_2, 0);
	send_double_matrix(ntypes, ntypes, iproc, u->r0_2, 0);
	if(!nodihpot){
	send_double_matrix(2, 2, iproc, u->sigma, 0);
	send_int_matrix(2, 2, iproc, u->dih0, 0);
	send_double_matrix(2, NAAMAX, iproc, u->ab_propensity, 0);
	}	

	MPI_Bcast((u->e_ang),nat,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast((u->ang0),nat,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast((u->e_ang),nat,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast((u->ang0),nat,MPI_DOUBLE,0,MPI_COMM_WORLD);
	if(!nodihpot){
	MPI_Bcast((u->e_dih1),nat,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast((u->dih01),nat,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast((u->e_dih3),nat,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast((u->dih03),nat,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast((u->dih_which),nat,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast((u->dih_pa),nat,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast((u->dih_pb),nat,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast((u->dih_f_psi_a),NDIHFMAX,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast((u->dih_f_psi_b),NDIHFMAX,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast((u->dih_f_phi_a),NDIHFMAX,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast((u->dih_f_phi_b),NDIHFMAX,MPI_DOUBLE,0,MPI_COMM_WORLD);
	}
	MPI_Bcast((u->hc_type),ntypes,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast((u->hc_r0),ntypes,MPI_DOUBLE,0,MPI_COMM_WORLD);
	if(hb)
	{
		MPI_Bcast((u->hb),nat,MPI_INT,0,MPI_COMM_WORLD);
		MPI_Bcast((u->hb_iam),nat,MPI_INT,0,MPI_COMM_WORLD);
	}
	return u;

}

void send_double_matrix(int length1, int length2, int iproc, double **m, int source)
{
	int i, j;
	double mm[length1][length2];
	if(iproc==source)
		for(i=0; i<length1; i++) 
			for(j=0; j<length2; j++) mm[i][j]=m[i][j];
	MPI_Bcast(&mm,length1*length2,MPI_DOUBLE,source,MPI_COMM_WORLD);
	if(iproc!=source) 
		for(i=0; i<length1; i++) 
			for(j=0; j<length2; j++) m[i][j]=mm[i][j];
	return;
}

void send_int_matrix(int length1, int length2, int iproc, int **m, int source)
{
	int i, j;
	int mm[length1][length2];
	if(iproc==source)
		for(i=0; i<length1; i++) 
			for(j=0; j<length2; j++) mm[i][j]=m[i][j];
	MPI_Bcast(&mm,length1*length2,MPI_INT,source,MPI_COMM_WORLD);
	if(iproc!=source) 
		for(i=0; i<length1; i++) 
			for(j=0; j<length2; j++) m[i][j]=mm[i][j];
	return;
}


/*****************************************************************************
 Copy a polymer structure from a replica to another 
 
 *****************************************************************************/
int ExchangePol(struct s_polymer *polymer, struct s_polymer *replica, struct s_polymer *oldp, struct s_mc_parms *parms, struct s_potential *u, int iproc, int ntemp, int even, int *ex_count, int *ex_acc, MPI_Datatype Backtype, MPI_Datatype Sidetype, MPI_Datatype Rottype, MPI_Status astatus,unsigned long long istep)
{
	int i,j,k,l, position,a=0;
	int x[ntemp];
	double delta,denomT;
	double E[ntemp];
	int nback = polymer->nback;
	int buffer_back_size=sizeof(struct s_back)*nback;
	int nosidechains=parms->nosidechains;
	int ipolymer;


	E[iproc]=polymer->etot;


	if(iproc==0) CountContacts(stderr,polymer,parms,istep);


	MPI_Barrier(MPI_COMM_WORLD);
	for(i=0;i<ntemp;i++)
		MPI_Bcast(&(E[i]),1,MPI_DOUBLE,i,MPI_COMM_WORLD);  //exchanging energy values between replicas
	MPI_Barrier(MPI_COMM_WORLD);
	//DEBUG
/*	if(iproc==0)
	for(i=0;i<ntemp;i++)
	{
//			fprintf(stderr,"step %llu E[%d]=\t%lf\n",istep,i,E[i]);
	}
*/	for(i=even;i<ntemp-1;i=i+2)
	{
		if(iproc==i) //sender		
		{		
			delta = -(E[i]-E[i+1]);
			denomT = (parms->T[i]*parms->T[i+1])/(parms->T[i+1]-parms->T[i]) ;
			if(Metropolis(delta,denomT,polymer->tables)==1) a=1; 	

			MPI_Send(&a, 1, MPI_INT, i+1, 1000000+i, MPI_COMM_WORLD); 
			
			#ifdef DEBUG
                        if(parms->debug>2) fprintf(stderr,"i=%d\tT(i)=%f\tT(i+1)=%f\tdelta=%f\tacc=%d\n",i,parms->T[i],parms->T[i+1],delta,a);
			#endif			
		}
		if(iproc==i+1) MPI_Recv(&a, 1, MPI_INT, i, 1000000+i, MPI_COMM_WORLD, &astatus); //receiver
		MPI_Barrier(MPI_COMM_WORLD); 
		if(a==1) //confirmed exchange
		{	
			if(iproc==i) fprintf(stderr,"[EX]\t%llu\t%d\t%d\n",istep,i,i+1); 	

		for(ipolymer=0;ipolymer<parms->npol;ipolymer++)
		{
				
			nback=(polymer+ipolymer)->nback;
			buffer_back_size=sizeof(struct s_back)*nback;
			char *buffer_back=malloc(buffer_back_size);

			if(iproc==i)
			{
				position=0;
				for(j=0;j<nback;j++)
        		                MPI_Pack(((polymer+ipolymer)->back)+j,1,Backtype,buffer_back,buffer_back_size,&position,MPI_COMM_WORLD);
	                         MPI_Send(buffer_back,position,MPI_PACKED,i+1, 10000+100*i, MPI_COMM_WORLD);
			}
			if(iproc==i+1)
			{
				MPI_Recv(buffer_back,buffer_back_size,MPI_PACKED,i,10000+100*i, MPI_COMM_WORLD, &astatus);
		                position=0;
                		for(j=0;j<nback;j++)
	               	        	MPI_Unpack(buffer_back,buffer_back_size,&position,((replica+ipolymer)->back)+j,1,Backtype,MPI_COMM_WORLD);

				position=0;
			
				for(j=0;j<nback;j++)
                                        MPI_Pack(((polymer+ipolymer)->back)+j,1,Backtype,buffer_back,buffer_back_size,&position,MPI_COMM_WORLD);
	                        MPI_Send(buffer_back,position,MPI_PACKED,i, 20000+100*i, MPI_COMM_WORLD);

			}
			
			if(iproc==i)
			{
			        MPI_Recv(buffer_back,buffer_back_size,MPI_PACKED,i+1, 20000+100*i, MPI_COMM_WORLD, &astatus);
                                position=0;
                                for(j=0;j<nback;j++)
                                        MPI_Unpack(buffer_back,buffer_back_size,&position,((replica+ipolymer)->back)+j,1,Backtype,MPI_COMM_WORLD);


			}
			
			if(!nosidechains)
			{

				int side_dims[nback];
				int rot_dims[nback];
               			int buffer_side_size,sidechains=0;
				int rotamersl=0;
			        for(j=0;j<nback;j++)
                		{
                        		side_dims[j]=(((polymer+ipolymer)->back)+j)->nside;
					rot_dims[j]=(((polymer+ipolymer)->back)+j)->nrot;
					sidechains+=side_dims[j];
					rotamersl+=rot_dims[j];
                		}

		                buffer_side_size=sizeof(struct s_side)*sidechains;
    				char *buffer_side=malloc(buffer_side_size);
				if(iproc==i)
				{
					position=0;
					for(j=0; j<nback;j++) for(k=0; k<side_dims[j];k++)
						MPI_Pack(((((polymer+ipolymer)->back)+j)->side)+k,1,Sidetype,buffer_side,buffer_side_size,&position,MPI_COMM_WORLD);	
					MPI_Send(buffer_side,position,MPI_PACKED,i+1, 1000000+100*i, MPI_COMM_WORLD);
				}
				if(iproc==i+1) 
				{
					MPI_Recv(buffer_side,buffer_side_size,MPI_PACKED,i,1000000+100*i, MPI_COMM_WORLD, &astatus);
					position=0;
					for(j=0; j<nback;j++) for(k=0; k<side_dims[j];k++) 
						MPI_Unpack(buffer_side,buffer_side_size,&position,((((replica+ipolymer)->back)+j)->side)+k,1,Sidetype,MPI_COMM_WORLD);
				}


				//rotamers	i -> i+1
				if(iproc==i)
				{
//					position=0;
					for(j=0; j<nback;j++) for(k=0; k<side_dims[j];k++) for(l=0; l<rot_dims[j]; l++)
//								MPI_Pack((((((polymer->back)+j)->side)+k)->rot)+l,1,Rottype,buffer_rot,buffer_rot_size,&position,MPI_COMM_WORLD);
//					MPI_Send(buffer_rot,position,MPI_PACKED,i+1, 2000000+1000*i, MPI_COMM_WORLD);
						MPI_Send(((((((polymer+ipolymer)->back)+j)->side)+k)->rot)+l, 1, Rottype, i+1, 2000000+1000*i+100*j+10*k+l, MPI_COMM_WORLD);
				}

				if(iproc==i+1)
				{
//					MPI_Recv(buffer_rot,buffer_rot_size,MPI_PACKED,i,2000000+1000*i, MPI_COMM_WORLD, &astatus);
//					position=0;
					for(j=0; j<nback;j++) for(k=0; k<side_dims[j];k++) for(l=0; l<rot_dims[j]; l++)
//						MPI_Unpack(buffer_rot,buffer_rot_size,&position,(((((replica->back)+j)->side)+k)->rot)+l,1,Rottype,MPI_COMM_WORLD);
						MPI_Recv(((((((replica+ipolymer)->back)+j)->side)+k)->rot)+l, 1, Rottype, i, 2000000+1000*i+100*j+10*k+l , MPI_COMM_WORLD, &astatus);	
				}

			
		//}
				//sidechains	i+i -> i
	//	if(!nosidechains)
	//	{	
				if(iproc==i+1)
				{
					position=0;
                                        for(j=0; j<nback;j++) for(k=0; k<side_dims[j];k++)
                                               MPI_Pack(((((polymer+ipolymer)->back)+j)->side)+k,1,Sidetype,buffer_side,buffer_side_size,&position,MPI_COMM_WORLD);
                                        MPI_Send(buffer_side,position,MPI_PACKED,i, 3000000+100*i, MPI_COMM_WORLD);

//					 for(j=0; j<nback;j++) for(k=0; k<((polymer->back)+j)->nside;k++) 
//						MPI_Send((((polymer->back)+j)->side)+k, 1, Sidetype, i, 3000000+100*i+10*j+k, MPI_COMM_WORLD);
				}
				if(iproc==i)
				{
	  		              	MPI_Recv(buffer_side,buffer_side_size,MPI_PACKED,i+1,3000000+100*i, MPI_COMM_WORLD, &astatus);
                                      	position=0;
                                        for(j=0; j<nback;j++) for(k=0; k<side_dims[j];k++)
                                                MPI_Unpack(buffer_side,buffer_side_size,&position,((((replica+ipolymer)->back)+j)->side)+k,1,Sidetype,MPI_COMM_WORLD);
//					for(j=0; j<nback;j++) for(k=0; k<((replica->back)+j)->nside;k++)
//						MPI_Recv((((replica->back)+j)->side)+k, 1, Sidetype, i+1, 3000000+100*i+10*j+k, MPI_COMM_WORLD, &astatus);
				}

				//rotamers 	i+1 -> i
				if(iproc==i+1)
				{
/*					position=0;
 	                                for(j=0; j<nback;j++) for(k=0; k<((polymer->back)+j)->nside;k++) for(l=0; l<rot_dims[j]; l++)
                                                MPI_Pack((((((polymer->back)+j)->side)+k)->rot)+l,1,Rottype,buffer_rot,buffer_rot_size,&position,MPI_COMM_WORLD);
                                       MPI_Send(buffer_rot,position,MPI_PACKED,i, 4000000+1000*i, MPI_COMM_WORLD);
*/
					 for(j=0; j<nback;j++) for(k=0; k<(((polymer+ipolymer)->back)+j)->nside;k++) for(l=0; l<(((polymer+ipolymer)->back)+j)->nrot; l++)
						MPI_Send(((((((polymer+ipolymer)->back)+j)->side)+k)->rot)+l, 1, Rottype, i, 4000000+1000*i+100*j+10*k+l, MPI_COMM_WORLD);
				}
		
				if(iproc==i)
				{
/*					MPI_Recv(buffer_rot,buffer_rot_size,MPI_PACKED,i+1,4000000+1000*i, MPI_COMM_WORLD, &astatus);
                                        position=0;
                                       for(j=0; j<nback;j++) for(k=0; k<((replica->back)+j)->nside;k++) for(l=0; l<rot_dims[j]; l++)
                                        MPI_Unpack(buffer_rot,buffer_rot_size,&position,(((((replica->back)+j)->side)+k)->rot)+l,1,Rottype,MPI_COMM_WORLD);
*/
					 for(j=0; j<nback;j++) for(k=0; k<(((replica+ipolymer)->back)+j)->nside;k++) for(l=0; l<(((replica+ipolymer)->back)+j)->nrot; l++)
						MPI_Recv(((((((replica+ipolymer)->back)+j)->side)+k)->rot)+l, 1, Rottype, i+1, 4000000+1000*i+100*j+10*k+l , MPI_COMM_WORLD, &astatus);	
				}
		
				free(buffer_side);
			}
			free(buffer_back);		
		}//end of ipolymer loop	

		ex_acc[i]++;	//accepted exchange counter
			x[i]=1; x[i+1]=1;		


		

		}//end of accepted exchange
		// if not accepted 
		else
		{
		x[i]=0; x[i+1]=0;
		}	
		ex_count[i]++; //tried exchange counter
		a=0;
	}
	
	if(x[iproc]==1)
	{
		CopyAllPolymers(replica,polymer,parms->npol,parms->nosidechains,parms->nosidechains);
		CopyAllPolymers(replica,oldp,parms->npol,parms->nosidechains,parms->nosidechains);

		//etot=Totalenergy : polymer and oldp contacts must be recomputed
		polymer->etot = TotalEnergy(polymer,u,parms,parms->npol,1,parms->nosidechains,parms->debug,iproc);
		oldp->etot=TotalEnergy(oldp,u,parms,parms->npol,1,parms->nosidechains,parms->debug,iproc);

		if (parms->shell==1)
		{
			UpdateShell(polymer,parms);
			CopyShell(polymer,oldp,parms);	
			parms->ishell=0;
		}
	}	

//	free(buffer_back);

	return a;
}
#ifdef OPTIMIZEPOT

struct s_optimizepot_input *Allo_op_input(int nrestr)
{
	struct s_optimizepot_input *x;
	x = (struct s_optimizepot_input *) calloc(1,sizeof(struct s_optimizepot_input));
	if (x==NULL) Error("Cannot allocate s_optimizepot_input");
	
	x->ndata = nrestr;
	
	x->expdata = AlloDouble(x->ndata);
	x->sigma = AlloDouble(x->ndata);
	x->datatype = AlloInt(x->ndata);
	x->i1 = AlloInt(x->ndata);
	x->i2 = AlloInt(x->ndata);
	x->i3 = AlloInt(x->ndata);
	x->i4 = AlloInt(x->ndata);

	return x;
}


struct s_optimizepot_input *send_op_input(int nrestr,struct s_optimizepot_input *in)
{
	
	in->ndata = nrestr;
	
	MPI_Bcast((in->expdata),nrestr,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast((in->sigma),nrestr,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast((in->datatype),nrestr,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast((in->i1),nrestr,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast((in->i2),nrestr,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast((in->i3),nrestr,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast((in->i4),nrestr,MPI_INT,0,MPI_COMM_WORLD);


	return in;

}

struct s_optimizepot_input *send_op(int nrestr,struct s_optimizepot_input *in)
{
	
	in->ndata = nrestr;
	
	MPI_Bcast((in->expdata),nrestr,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast((in->sigma),nrestr,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast((in->datatype),nrestr,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast((in->i1),nrestr,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast((in->i2),nrestr,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast((in->i3),nrestr,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast((in->i4),nrestr,MPI_INT,0,MPI_COMM_WORLD);


	return in;

}


#endif

#endif
