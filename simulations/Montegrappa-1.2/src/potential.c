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
 * potential.c
 *
 *  Created on: Sep 23, 2010
 *      Author: guido
 */

#include "montegrappa.h"

double TotalEnergy(struct s_polymer *p, struct s_potential *u, struct s_mc_parms *parms, int npol, int update, int nosidechains, int debug, int iproc)
{
      int i,j,ci,cj,tooclose;
      double etot=0,eang=0,edih=0,ebox=0;
      if (update)
            for (ci=0;ci<npol;ci++)
                  for (i=0;i<(p+ci)->nback;i++)
                  {
                        for (j=0;j<(((p+ci)->back)+i)->ncontacts;j++)  ((((p+ci)->back)+i)->e)[j] = 0;
                        (((p+ci)->back)+i)->ncontacts = 0;
                  }


      for (ci=0;ci<npol;ci++)                                     // loops on polymers
            for (cj=ci;cj<npol;cj++)
                  for (i=0;i<(p+ci)->nback;i++)             // loops on backbone atoms of each polymer
                        for (j=0;j<(p+cj)->nback;j++)
                              if (ci!=cj || i<j)                        // in the same polymer, do not count twice a contact
                              {
                                    tooclose = 0;
                                    if (ci==cj && j-i <= u->g_imin ) tooclose = 1;
                                    etot += EnergyPair(p,u,i,j,ci,cj,update,nosidechains,parms->disentangle,tooclose,parms->hb); 
                             }

      if (debug>0 && iproc==0) fprintf(parms->flog," Epairs = %lf\n",etot);

      if (!parms->noangpot || !parms->nodihpot){
      
      for (ci=0;ci<npol;ci++)           
            for (i=0;i<(p+ci)->nback;i++)
            {
                  if (!parms->noangpot)   eang += EnergyAngles(p,u,i,ci,update);                      // angular potential
                  if (!parms->nodihpot)   edih += EnergyDihedrals(p,u,i,ci,update);             // dihedral potential
                  if (u->boxtype != 'n')                                                  // box potential
                        if (EnergyBox(p,u,i,ci) == 1) ebox = LARGE;
            }

	}

      etot += (eang + edih + ebox);

      if (debug>0 && iproc==0) fprintf(parms->flog," Eang  = %lf\n Edih  = %lf\n Etot  = %lf\n",eang,edih,etot);


      #ifdef OPTIMIZEPOT
       if (strcmp(parms->op_minim,"none"))
       {
            p->op->eold[ p->op->nframes ] = etot;
            p->op->efix[ p->op->nframes ] = eang + edih + ebox;
       }
      #endif


      return etot;
}




/********************************************************************
 Calculate the energy of a pair of residues,
 updating contact lists and partial energies if update==1.
 ********************************************************************/
double EnergyPair(struct s_polymer *p, struct s_potential *u, int i, int j, int ci, int cj, int update, int nosidechains, int disentangle, int tooclose, int hb)
{
	int is,js,a1,a2,ia1,ia2,out;
	double r2,etot=0,e,r02,cos1,cos2,mul;
	
	// within same aa do not interact
	if (ci==cj && (((p+ci)->back)+i)->iaa == (((p+cj)->back)+j)->iaa) return 0;

	// interaction backbone-backbone (i-j)
	if (ci!=cj || (i!=(j+1) && i!=(j-1)))						// if not consecutive (i.e., no covalent interaction)
	{
		a1 = (((p+ci)->back)+i)->itype;
		a2 = (((p+cj)->back)+j)->itype;

		r2 = Dist2( (((p+ci)->back)+i)->pos, (((p+cj)->back)+j)->pos );

		if (tooclose) r02 = u->g_r0hard * u->g_r0hard;
		else r02 = (u->r0_2)[a1][a2];

		if ( r2 < r02 )
		{
			if (!disentangle) return LARGE;
			else etot += LARGE;
		}
		if ( r2 < (u->r_2)[a1][a2] && !tooclose)					// if they are close in space but not in sequence
			if ((u->e)[a1][a2]<-EPSILON || (u->r_2)[a1][a2]>0)			// and their interaction element is defined
			{
				e = (u->e)[a1][a2];
				mul=1.;
				if (u->splice==1)
					if (r2 > (u->r_2)[a1][a2] * u->kr2_splice) mul = u->ke_splice;

				#ifdef OPTIMIZEPOT
				if (p->op->record && mul>0) OP_AddEnergy(p,a1,a2,mul);
				#endif
				
				e *= mul;
				etot += e;
				if (update) AddContact(p,i,j,ci,cj,e);

			}
	}

  	if (nosidechains == 1) return etot;



	for (js=0;js< (((p+cj)->back)+j)->nside;js++)
	{
		// interaction backbone-sidechain (i-js)
		a1 = (((p+ci)->back)+i)->itype;
		a2 = (((((p+cj)->back)+j)->side)+js)->itype;
		r2 = Dist2( (((p+ci)->back)+i)->pos, (((((p+cj)->back)+j)->side)+js)->pos );

		if (tooclose) r02 = u->g_r0hard * u->g_r0hard;
		else r02 = (u->r0_2)[a1][a2];

		if ( r2 < r02 )
		{
			if (!disentangle) return LARGE;
			else etot += LARGE;
		}
		if ( r2 < (u->r_2)[a1][a2] && !tooclose)				// if they are close in space and not in sequence
			if ((u->e)[a1][a2]<-EPSILON || (u->r_2)[a1][a2]>0)		// and their interaction element is defined
			{
				e = (u->e)[a1][a2];
				mul = 1.;
				if (u->splice==1)
					if (r2 > (u->r_2)[a1][a2] * u->kr2_splice) mul = u->ke_splice;
				
				#ifdef OPTIMIZEPOT
				if (p->op->record && mul>0) OP_AddEnergy(p,a1,a2,mul);
				#endif
				
				e *= mul;
				etot += e;
				if (update) AddContact(p,i,j,ci,cj,e);
			}

		for (is=0;is< (((p+ci)->back)+i)->nside;is++)
		{
			// interaction sidechain-sidechain (is-js)
			a1 = (((((p+ci)->back)+i)->side)+is)->itype;
			a2 = (((((p+cj)->back)+j)->side)+js)->itype;

			r2 = Dist2( (((((p+ci)->back)+i)->side)+is)->pos, (((((p+cj)->back)+j)->side)+js)->pos );
			
			if (tooclose) r02 = u->g_r0hard * u->g_r0hard;
			else r02 = (u->r0_2)[a1][a2];

			if ( r2 < r02 )
			{
				//fprintf(stderr,"side %d of %d - side %d of %d : r2<r02\n",is,i,js,j);
				if (!disentangle) return LARGE;
				else etot += LARGE;
			}
			if ( r2 < (u->r_2)[a1][a2] && !tooclose)				// if they are close in space anbd not in sequence
				if ((u->e)[a1][a2]<-EPSILON || (u->r_2)[a1][a2]>0)		// and their interaction element is defined
				{
					e = (u->e)[a1][a2];
					mul = 1;
					if (u->splice==1)
						if (r2 > (u->r_2)[a1][a2] * u->kr2_splice) mul =  u->ke_splice;

					// calculate hydrogen bond (if any)
					if (hb)
					{
						if (u->splice==1)		// if HB, interacts only the first half of splice
							if (r2 > (u->r_2)[a1][a2] * u->kr2_splice) mul = 0;

						ia1 = (((((p+ci)->back)+i)->side)+is)->ia;
						ia2 = (((((p+cj)->back)+j)->side)+js)->ia;
						if ((u->hb)[ia1]*(u->hb)[ia2]==2)
						{
							cos1 = CosAngle( (*(*(((p+ci)->vback)+(u->hb_iam)[ia1]))), (*(*(((p+ci)->vback)+ia1))), 
									(*(*(((p+cj)->vback)+ia2))), p->tables, &out );
							if (out==0) return LARGE;
							cos2 = CosAngle( (*(*(((p+ci)->vback)+ia1))), (*(*(((p+cj)->vback)+ia2))), 
									(*(*(((p+cj)->vback)+(u->hb_iam)[ia2]))), p->tables , &out);
							if (out==0) return LARGE;
							mul = FastSqrt( DAbs(cos1*cos2),p->tables);
						}
					}

					#ifdef OPTIMIZEPOT
					if (p->op->record && mul>0) OP_AddEnergy(p,a1,a2,mul);
					#endif
					if(i==10 && j==43 && r2< 2.) fprintf(stderr,"dist = %lf\n",r2);					
					e *= mul;
					etot += e;
					if (update) AddContact(p,i,j,ci,cj,e);
				}
		}
	}

	// interaction backbone-sidechain (is-j)
	for (is=0;is< (((p+ci)->back)+i)->nside;is++)
	{
			a1 = (((p+cj)->back)+j)->itype;
			a2 = (((((p+ci)->back)+i)->side)+is)->itype;
			r2 = Dist2( (((p+cj)->back)+j)->pos, (((((p+ci)->back)+i)->side)+is)->pos );

			if (tooclose) r02 = u->g_r0hard * u->g_r0hard;
			else r02 = (u->r0_2)[a1][a2];

			if ( r2 < r02 )
			{
				if (!disentangle) return LARGE;
				else etot += LARGE;
			}
			if ( r2 < (u->r_2)[a1][a2] && !tooclose)				// if they are close in space and not in sequence
				if ((u->e)[a1][a2]<-EPSILON || (u->r_2)[a1][a2]>0)		// and their interaction element is defined
				{
					e = (u->e)[a1][a2];
					mul = 1;
					if (u->splice==1)
						if (r2 > (u->r_2)[a1][a2] * u->kr2_splice) mul =  u->ke_splice;
					
					#ifdef OPTIMIZEPOT
					if (p->op->record && mul>0) OP_AddEnergy(p,a1,a2,mul);
					#endif
					
					e *= mul; 
					etot += e;
					if (update) AddContact(p,i,j,ci,cj,e);
				}
	}

	return etot;
}

/********************************************************************
 Calculate the energy of ith residue with all the others,
 updating contact lists and partial energies. (single loop)
 ********************************************************************/
double EnergyMonomer(struct s_polymer *p, struct s_potential *u, int i, int ci, int npol, int update, int shell, int nosidechains,int disentangle,int hb)
{
	int j,cj,k,tooclose;
	double e,etot=0;

/*	// box energy
	if (u->boxtype != 'n')
		if (EnergyBox(p,u,i,ci) == 1) return LARGE;
*/
	// delete contact list and partial energies of i and of its contacting monomers
	if ( update ) ResetContactsMonomer(p,i,ci);


	// box energy
	if (u->boxtype != 'n')
	if (EnergyBox(p,u,i,ci) == 1) return LARGE;


	if (shell == 0)
	{
		for (cj=0;cj<npol;cj++)
			for (j=0;j<(p+cj)->nback;j++)
			{
//				if(i==1 && ci == 1 && j == 7 && cj == 1) fprintf(stderr,"interaction atom %d of %d // %d of %d \n distance %lf",i,ci,j,cj, sqrt(Dist2( (((p+ci)->back)+i)->pos, (((p+cj)->back)+j)->pos )));				
				if ( ci!= cj || i!=j )
				{
					tooclose = 0;
					if (ci == cj && Abs(i-j) <= u->g_imin) tooclose = 1;
					e = EnergyPair(p,u,i,j,ci,cj,update,nosidechains,disentangle,tooclose,hb);
				//	          if(i==136 && ci == 0 && j == 151 && cj == 0)   
				//		  fprintf(stderr,"NOSHELL interaction atom %d of %d // %d of %d E=%lf\tirot1=%d\tirot2=\%d\n",i,ci,j,cj,e,(((p+ci)->back)+i)->irot,(((p+cj)->back)+j)->irot);

					if (e>=LARGE && !disentangle) return LARGE;
					etot += e;
				}
			}
	}

	// use shell
	else
	{
		for (k=0;k<(((p+ci)->back)+i)->nshell;k++)
		{
			j = ((((p+ci)->back)+i)->shell)[k];									// what is the monomer in the shell
			cj = ((((p+ci)->back)+i)->shell_p)[k];	
	//		if(i==26 && ci == 0 && j == 37 && cj == 0) 							// ... what polymers it belongs to
	//		fprintf(stderr,"interaction atom %d of %d // %d of %d \n",i,ci,j,cj);	
			//SSSS
			if( ci!=cj || i!=j)
			{
         			//if(i==26 && ci == 0 && j == 37 && cj == 0)     fprintf(stderr,"interaction atom %d of %d // %d of %d \n",i,ci,j,cj);

				tooclose=0;
				if (ci == cj && Abs(i-j) <= u->g_imin) tooclose = 1;
				e = EnergyPair(p,u,i,j,ci,cj,update,nosidechains,disentangle,tooclose,hb);
				//if(i==136 && ci == 0 && j == 151 && cj == 0)  
				//fprintf(stderr,"SISHELL interaction atom %d of %d // %d of %d E=%lf\tirot1=%d\tirot2=\%d\n",i,ci,j,cj,e,(((p+ci)->back)+i)->irot,(((p+cj)->back)+j)->irot);
 

				if (e>=LARGE && !disentangle) return LARGE;
				etot += e;
			}
		}
	}

	return etot;
}

/********************************************************************
 Calculate the energy of energies of monomers in a range, updating
 contacts and partial energies
 ********************************************************************/
double EnergyMonomerRange(struct s_polymer *p, struct s_potential *u, int from, int to, int ip, int npol, int shell, int update, int nosidechains, int disentangle, int hb)
{
	int i,j,cj,k,tooclose;
	double etot=0,e;

	if (from<0) from=0;
	if (to>(p+ip)->nback - 1) to=(p+ip)->nback - 1;



	// delete contact list and partial energies of i and of its contacting monomers
	if ( update )
		for (i=from;i<=to;i++)
		{
			ResetContactsMonomer(p,i,ip);
		}

	// box energy
	if (u->boxtype != 'n')
        	for (i=from;i<=to;i++)
			if (EnergyBox(p,u,i,ip) == 1) return LARGE;

	// if not shell
	if (shell == 0)
	{
		for (i=from;i<=to;i++)
			for (cj=0;cj<npol;cj++)
				for (j=0;j<(p+cj)->nback;j++)
					if ( ip!=cj || j<from || j>to || i<j )				// if j in [from,to] of the same chain, require i<j to count only once
					{
						//fprintf(stderr,"energypair atomi %d di %d / %d di %d \n",i,ip,j,cj);
						tooclose = 0;
						if (ip == cj && Abs(i-j) <= u->g_imin) tooclose = 1;
						e = EnergyPair(p,u,i,j,ip,cj,update,nosidechains,disentangle,tooclose,hb);
						if (e>=LARGE && !disentangle) return LARGE;
						etot += e;
					}
	}

	// if shell
	else
	{
		for (i=from;i<=to;i++)
			for (k=0;k<(((p+ip)->back)+i)->nshell;k++)
			{
				j = ((((p+ip)->back)+i)->shell)[k];							// what is the monomer in the shell
				cj = ((((p+ip)->back)+i)->shell_p)[k];						// ... what polymers it belongs to
				if ( ip!=cj || j<from || j>to || i<j )					// if j in [from,to] of the same chain, require i<j to count only once
				{
					
					tooclose = 0;
					if (ip == cj && Abs(i-j) <= u->g_imin) tooclose = 1;
					e = EnergyPair(p,u,i,j,ip,cj,update,nosidechains,disentangle,tooclose,hb);
					if (e>=LARGE && !disentangle) return LARGE;
					etot += e;
				}
			}
	}
	return etot;
}

/********************************************************************
 Tell backbones i and j that they are in contact; if they are already in
 contact, update the partial energy
 ********************************************************************/
void AddContact(struct s_polymer *p, int i, int j, int ci, int cj, double e)
{

	int k,w;
	
	// check if j is already among the contacts of i
	w=-1;
	for (k=0;k<(((p+ci)->back)+i)->ncontacts;k++)
		if ( ((((p+ci)->back)+i)->contacts)[k] == j && ((((p+ci)->back)+i)->contacts_p)[k] == cj) w=k;  // w is the element of the contact list
																			// corresponding to j
	if (w==-1) // if not present, add it
	{
		((((p+ci)->back)+i)->contacts)[ (((p+ci)->back)+i)->ncontacts ] = j;			// add contact
		((((p+ci)->back)+i)->contacts_p)[ (((p+ci)->back)+i)->ncontacts ] = cj;			// add contact
		((((p+ci)->back)+i)->e)[ (((p+ci)->back)+i)->ncontacts ] = e;					// add energy
		(((p+ci)->back)+i)->ncontacts ++;
		 if ( (((p+ci)->back)+i)->ncontacts >= NCONTMAX) Error("NCONTMAX too small");
	}
	else	// if already a contact, add the associated energy
	{
		((((p+ci)->back)+i)->e)[w] += e;
	}

	//fprintf(stderr,"added contact atom %d chain %d with atom %d chain %d\n",i,ci,j,cj);
	// check if i is already among the contacts of j
	w=-1;
	for (k=0;k<(((p+cj)->back)+j)->ncontacts;k++)
		if ( ((((p+cj)->back)+j)->contacts)[k] == i && ((((p+cj)->back)+j)->contacts_p)[k] == ci) w=k; // w is the element of the contact list
																		// corresponding to i

	if (w==-1) // if not present, add it
	{
		((((p+cj)->back)+j)->contacts)[ (((p+cj)->back)+j)->ncontacts ] = i;			// add contact
		((((p+cj)->back)+j)->contacts_p)[ (((p+cj)->back)+j)->ncontacts ] = ci;			// add contact
		((((p+cj)->back)+j)->e)[ (((p+cj)->back)+j)->ncontacts ] = e;					// add energy
		(((p+cj)->back)+j)->ncontacts ++;
		if ( (((p+cj)->back)+j)->ncontacts >= NCONTMAX) Error("NCONTMAX too small");
	}
	else	// if already a contact, add the associated energy
	{
		((((p+cj)->back)+j)->e)[w] += e;
	}
	//fprintf(stderr,"added contact atom %d chain %d with atom %d chain %d\n\n\n",i,ci,j,cj);
}

/********************************************************************
 Tell i and j that they are in the same shell
 ********************************************************************/
void AddShell(struct s_polymer *p, int i, int j, int ci, int cj)
{
	int k,ok;

	// check if j is already among the shell of i
	ok=1;
	for (k=0;k<(((p+ci)->back)+i)->nshell;k++)
		if ( ((((p+ci)->back)+i)->shell)[k] == j && ((((p+ci)->back)+i)->shell_p)[k] == cj) ok=0;
	if (ok==1) // if not, add it
	{
		((((p+ci)->back)+i)->shell)[ (((p+ci)->back)+i)->nshell ] = j;
		((((p+ci)->back)+i)->shell_p)[ (((p+ci)->back)+i)->nshell ] = cj;
		(((p+ci)->back)+i)->nshell ++;
		#ifdef DEBUG
			if ( (((p+ci)->back)+i)->ncontacts >= NSHELLMAX) Error("NSHELLMAX too small");
		#endif
	}

	// check if i is already among the shell of j
	ok=1;
	for (k=0;k<(((p+cj)->back)+j)->nshell;k++)
		if ( ((((p+cj)->back)+j)->shell)[k] == i && ((((p+cj)->back)+j)->shell_p)[k] == ci) ok=0;
	if (ok==1) // if not, add it
	{
		((((p+cj)->back)+j)->shell)[ (((p+cj)->back)+j)->nshell ] = i;
		((((p+cj)->back)+j)->shell_p)[ (((p+cj)->back)+j)->nshell ] = ci;
		(((p+cj)->back)+j)->nshell ++;
		#ifdef DEBUG
			if ( (((p+cj)->back)+j)->ncontacts >= NSHELLMAX) Error("NSHELLMAX too small");
		#endif
	}
}

/********************************************************************
 Reset the contacts of i-th monomer, and all references to it of
 its contacting monomers (contact list and energies)
 ********************************************************************/
void ResetContactsMonomer(struct s_polymer *p, int i, int ci)
{
	int k,neigh,l,cn;

	// remove i from its neighbours
	for (k=0;k<(((p+ci)->back)+i)->ncontacts;k++)
	{
		neigh = ((((p+ci)->back)+i)->contacts)[k];					// neigh is the k-th contact of i
		cn = ((((p+ci)->back)+i)->contacts_p)[k];					// ... its chain
		for (l=0;l<(((p+cn)->back)+neigh)->ncontacts;l++)
			if ( ((((p+cn)->back)+neigh)->contacts)[l] == i )		// l is the position of i in the list of neigh
			{								// overwrite l-th with the last (i.e., delete both contact and energy)
				((((p+cn)->back)+neigh)->contacts)[l] = ((((p+cn)->back)+neigh)->contacts)[ (((p+cn)->back)+neigh)->ncontacts - 1 ];
//FRANBOB
				((((p+cn)->back)+neigh)->contacts_p)[l] = ((((p+cn)->back)+neigh)->contacts_p)[ (((p+cn)->back)+neigh)->ncontacts - 1 ];
				((((p+cn)->back)+neigh)->e)[l] = ((((p+cn)->back)+neigh)->e)[ (((p+cn)->back)+neigh)->ncontacts - 1 ];
				//fprintf(stderr,"ncontacts = %d\n",(((p+cn)->back)+neigh)->ncontacts);


				(((p+cn)->back)+neigh)->ncontacts --;
                                //fprintf(stderr,"ncontacts = %d\n",(((p+cn)->back)+neigh)->ncontacts);

				l--;
			}

		//reset its energy
		(((p+ci)->back)+i)->e[k] = 0;
	}

	// reset the neighbour list of i
	(((p+ci)->back)+i)->ncontacts = 0;

}


/********************************************************************
 Read the total interaction energy of a monomer with the others
 from its energy structure
 ********************************************************************/
double GetEnergyMonomer(struct s_polymer *p, int ip, int iw)
{
	double e=0;
	int i;

	for (i=0;i<(((p+ip)->back)+iw)->ncontacts;i++)
	{
		
		e += ((((p+ip)->back)+iw)->e)[i];
	}
	return e;
}

/********************************************************************
 Read the total interaction energy of a set of monomers
 from its energy structure (ip is the id of the chain)
 ********************************************************************/
double GetEnergyMonomerRange(struct s_polymer *p, int from, int to, int ip)
{
	double e=0;
	int i,iw,j,cj;

	for (iw=from;iw<=to;iw++)
		if (iw>=0 && iw<(p+ip)->nback)
			for (i=0;i<(((p+ip)->back)+iw)->ncontacts;i++)
			{
				j = ((((p+ip)->back)+iw)->contacts)[i];
				cj = ((((p+ip)->back)+iw)->contacts_p)[i];
				if (ip!=cj || j<from || j>to || iw<j)				// count only once contacts in [from,to] of the same chain
					e += ((((p+ip)->back)+iw)->e)[i];
			}
	return e;
}

/********************************************************************
 Print contacts and partial energies for debugging purposes
 ********************************************************************/
void PrintContacts(FILE *fp, struct s_polymer *p,int ip, unsigned long long step)
{
	int i,j;

	fprintf(fp,"\n\n------------------------\nSTEP = %llu\n",step);

	for (i=0;i<(p+ip)->nback;i++)
	{
		fprintf(fp,"Back = %d\n\tside =",i);
		for (j=0;j<(((p+ip)->back)+i)->nside;j++)
			fprintf(fp,"%d ",(((((p+ip)->back)+i)->side)+j)->ia);
		fprintf(fp,"\nContacts = ");
		for (j=0;j<(((p+ip)->back)+i)->ncontacts;j++)
			fprintf(fp,"%d ",((((p+ip)->back)+i)->contacts)[j]);
		fprintf(fp,"\nContactsP= ");
	        for (j=0;j<(((p+ip)->back)+i)->ncontacts;j++)
                fprintf(fp,"%d ",((((p+ip)->back)+i)->contacts_p)[j]);
		fprintf(fp,"\n");
		for (j=0;j<(((p+ip)->back)+i)->ncontacts;j++)
			fprintf(fp,"%lf ",((((p+ip)->back)+i)->e)[j]);
		fprintf(fp,"\n");
	}
}


void CountContacts(FILE *fp,struct s_polymer *polymer,struct s_mc_parms *parms,unsigned long long step)
{
                int i,tempncontacts,iback,ipol,icont,matcont[parms->npol];
                for(i=0;i<parms->npol;i++)
                {
				for(ipol=0;ipol<parms->npol;ipol++)
                                        matcont[ipol]=0;
                                tempncontacts=0;
                                fprintf(fp,"> Chain %d\n  %d backbone atoms\n",i,(polymer+i)->nback);
                                for(iback=0;iback<(polymer+i)->nback;iback++)
                                {
                                        tempncontacts+=((polymer+i)->back+iback)->ncontacts;
					for(icont=0;icont<((polymer+i)->back+iback)->ncontacts;icont++)
						matcont[ ((polymer+i)->back+iback)->contacts_p[icont] ]++;
                            
                                
                                }
     
                                fprintf(fp,"  %d contacts ( ",tempncontacts);
				for(ipol=0;ipol<parms->npol;ipol++)
					fprintf(fp," %d ",matcont[ipol]);
				fprintf(fp,")\n");
		}
		

}

/********************************************************************
 Update shells of neighbours
 ********************************************************************/
void UpdateShell(struct s_polymer *p, struct s_mc_parms *parms)
{
	int i,j,is,js,ci,cj;
	double r2;

	/*
	#ifdef DEBUG_SHELL
	for (ci=0;ci<parms->npol;ci++)
                for (i=0;i<(p+ci)->nback;i++)
		if(i==100) fprintf(stderr,"atomo %d, prima -> %d\n",i,(((p+ci)->back)+i)->nshell);
	#endif
	*/
	// reset contacts shells
	for (ci=0;ci<parms->npol;ci++)
		for (i=0;i<(p+ci)->nback;i++)
			(((p+ci)->back)+i)->nshell = 0;


	//  double loop
	for (ci=0;ci<parms->npol;ci++)
		for (cj=ci;cj<parms->npol;cj++)  
			for (i=0;i<(p+ci)->nback;i++)
                                for (j=0;j<(p+cj)->nback;j++)
		{
			
					if (ci!=cj || i<j)
					{
						// interaction backbone-backbone (i-j)
						r2 = Dist2( (((p+ci)->back)+i)->pos, (((p+cj)->back)+j)->pos );
						if ( r2 < parms->r2shell ) AddShell(p,i,j,ci,cj);					// add to shell

						for (js=0;js< (((p+cj)->back)+j)->nside;js++)
						{
							// interaction backbone-sidechain (i-js)
							r2 = Dist2( (((p+ci)->back)+i)->pos, (((((p+cj)->back)+j)->side)+js)->pos );
							if ( r2 < parms->r2shell ) AddShell(p,i,j,ci,cj);			// add to shell

							for (is=0;is< (((p+ci)->back)+i)->nside;is++)
							{
								// interaction sidechain-sidechain (is-js)
								r2 = Dist2( (((((p+ci)->back)+i)->side)+is)->pos, (((((p+cj)->back)+j)->side)+js)->pos );
//								if(i==4&&j==148) fprintf(stderr,"updateshell: i=4,j=148;\tdist2=%lf\n",r2); 
								if ( r2 < parms->r2shell ) AddShell(p,i,j,ci,cj);			// add to shell
							}

						}

						// interaction backbone-sidechain (is-j)
						for (is=0;is< (((p+ci)->back)+i)->nside;is++)
						{
								r2 = Dist2( (((p+cj)->back)+j)->pos, (((((p+ci)->back)+i)->side)+is)->pos );
								if ( r2 < parms->r2shell ) AddShell(p,i,j,ci,cj);			// add to shell
						}
					}
		}
	/*
	#ifdef DEBUG_SHELL
        for (ci=0;ci<parms->npol;ci++)
                for (i=0;i<(p+ci)->nback;i++)
        	if(i==100)        fprintf(stderr,"atomo %d, dopo -> %d\n",i,(((p+ci)->back)+i)->nshell);
        #endif
	
	int pizza;
	fprintf(stderr,"UPDATE SHELL\n");
	for (ci=0;ci<parms->npol;ci++)
                for (i=0;i<(p+ci)->nback;i++)
		{
			fprintf(stderr,"ip= %d \t back = %d\t nshell = %d \n",ci,i,(((p+ci)->back)+i)->nshell);
			for(pizza=0;pizza<(((p+ci)->back)+i)->nshell;pizza++)
				fprintf(stderr,"sh %d chain %d\n",pizza,(((p+ci)->back)+i)->shell_p[pizza]);
		}
	*/
/*
	for(i=0;i<(p+0)->nback;i++)
	{
		fprintf(stderr,"back %d nshell = %d\n",i,((p+0)->back+i)->nshell);
		for(j=0;j<((p+0)->back+i)->nshell;j++)
			fprintf(stderr,"\t %d \t %d \n",((p+0)->back+i)->shell[j],((p+0)->back+i)->shell_p[j]);


	}
*/
	return;
}



void CopyShell(struct s_polymer *from,struct s_polymer *to,struct s_mc_parms *parms)
{

	int i,ci,cs;
	for (ci=0;ci<parms->npol;ci++)
                for (i=0;i<(from+ci)->nback;i++)
		{
                        (((to+ci)->back)+i)->nshell = (((from+ci)->back)+i)->nshell ;
			for(cs=0;cs<(((from+ci)->back)+i)->nshell;cs++)
			{
				(((to+ci)->back)+i)->shell[cs]=(((from+ci)->back)+i)->shell[cs];		
			        (((to+ci)->back)+i)->shell_p[cs]=(((from+ci)->back)+i)->shell_p[cs];	
			}
		}






	return;
}



/********************************************************************
 Angular potential of backbone iw of chain ic
 ********************************************************************/
double EnergyAngles(struct s_polymer *p, struct s_potential *u, int iw, int ic, int update)
{
	int ia,out;
	double ang,e;

	if (iw<1 || iw>(p+ic)->nback-2) return 0;			// angles are not defined for the first and last backbone atom

	ia = (((p+ic)->back)+iw)->ia;

	ang = Angle( (((p+ic)->back)+iw-1)->pos, (((p+ic)->back)+iw)->pos, (((p+ic)->back)+iw+1)->pos, p->tables, &out );
	if (out==0) return LARGE;

	e = (u->e_ang)[ia] * ( ang - (u->ang0)[ia] ) * ( ang - (u->ang0)[ia] );

	if (update)  (((p+ic)->back)+iw)->e_ang = e;

	return e;
}

/********************************************************************
 Dihedral potential of backbone iw of chain ic
 (defined as dihedral among iw-2,iw-1,iw,iw+1)
 ********************************************************************/
double EnergyDihedrals(struct s_polymer *p, struct s_potential *u, int iw, int ic, int update)
{
	int ia,iaa,i, out;
	double e=0,dih01,dih03,dih;

	if (iw<2 || iw>(p+ic)->nback-2) return 0;			// dihedrals are not defined for the first and last two backbone atoms

	ia = (((p+ic)->back)+iw)->ia;
        iaa = (((p+ic)->back)+iw)->iaa - (((p+ic)->back))->iaa;

	dih01 = (u->dih01)[ia];
	dih03 = (u->dih03)[ia];

	dih = Dihedral( (((p+ic)->back)+iw-2)->pos, (((p+ic)->back)+iw-1)->pos, (((p+ic)->back)+iw)->pos, (((p+ic)->back)+iw+1)->pos, p->tables, &out );
	if (out==0) { fprintf(stderr,"WARNING: Error in evaluating dihedral energy"); return LARGE; }

	// periodic kind of potential, i.e. e1*cos(d-d0) + e3*cos(3(d-d0))
	if (u->dih_periodic)
		e = (u->e_dih1)[ia] * ( 1 - FastCos(dih-dih01,p->tables) ) + (u->e_dih3)[ia] * ( 1 - FastCos(3.*(dih-dih03),p->tables) );

	// tabled kind of potential
	if (u->dih_tabled)
	{
		i = (int) ((dih+180.)/u->g_dihbin);
        #ifdef DEBUG
	      if (i<0 || i>NDIHFMAX) { fprintf(stderr,"WARNING: dih=%lf i=%d in EnergyDihedrals\n",dih,i); return LARGE; }
        #endif
		if ((u->dih_which)[ia] == 0)  // is a phi
			e += u->g_dihe * (u->dih_pa)[ia] * (u->dih_f_phi_a)[i] + u->g_dihe * (u->dih_pb)[ia] * (u->dih_f_phi_b)[i];
		else						  // is a psi
			e += u->g_dihe * (u->dih_pa)[ia] * (u->dih_f_psi_a)[i] + u->g_dihe * (u->dih_pb)[ia] * (u->dih_f_psi_b)[i];
	}

	if(u->dih_ram)
        {

	if((ia-1)%3==0)   //di tipo phi
                  for(i=0;i<2;i++)
                        e += - (u->e_dihram) * ( (u->ab_propensity[i][iaa]/u->sigma[i][0]) * exp(-0.5*(((dih - u->dih0[i][0])/u->sigma[i][0])*((dih - u->dih0[i][0])/u->sigma[i][0]))) );
        if((ia+1)%3==0) //di tipo psi
                 for(i=0;i<2;i++)
                        e += - (u->e_dihram) * (u->ab_propensity[i][iaa]/u->sigma[i][1]) * exp(-0.5*(((dih - u->dih0[i][1])/u->sigma[i][1])*((dih - u->dih0[i][1])/u->sigma[i][1])));




	}

	if (update)  (((p+ic)->back)+iw)->e_dih = e;

	return e;
}

/********************************************************************
 Print energy terms
 ********************************************************************/
void PrintEnergies(FILE *fp, int nc, unsigned long long step, struct s_polymer *p)
{
	int i,j;
	double epair=0,eang=0,edih=0;

	for (i=0;i<nc;i++)
		for (j=0;j<(p+i)->nback;j++)
		{
			epair += GetEnergyMonomer(p,i,j) / 2.;
			eang += (((p+i)->back)+j)->e_ang;
			edih += (((p+i)->back)+j)->e_dih;
		}

	fprintf(fp,"%llu\t%9lf\t%9lf\t%9lf\t%9lf\n",step,epair+eang+edih,epair,eang,edih);
}

void PrintEnergies_Parallel(FILE *fp, int nc, unsigned long long step, struct s_polymer *p,int my_rank)
{
        int i,j;
        double epair=0,eang=0,edih=0;

        for (i=0;i<nc;i++)
                for (j=0;j<(p+i)->nback;j++)
                {
                        epair += GetEnergyMonomer(p,i,j)/2;
                        eang += (((p+i)->back)+j)->e_ang;
                        edih += (((p+i)->back)+j)->e_dih;
                }
	
     //   fprintf(fp,"%d\t%llu\t%9lf\t%9lf\t%9lf\t%9lf\n",my_rank,step,epair+eang+edih,epair,eang,edih);
	fprintf(fp,"%d\t%llu\t%9lf\n",my_rank,step,epair+eang+edih);

}


void CopyPotential(struct s_potential *from, struct s_potential *to, int nat, int ntypes)
{
	int i,j;

	for (i=0;i<ntypes;i++)
		for (j=0;j<ntypes;j++)
		{
			to->r0_2[i][j] = from->r0_2[i][j];
			to->r_2[i][j] = from->r_2[i][j];
			to->e[i][j] = from->e[i][j];
		}

	for (i=0;i<nat;i++)
	{
		to->ang0[i] = from->ang0[i];
		to->e_ang[i] = from->e_ang[i];
		to->e_dih1[i] = from->e_dih1[i];
		to->dih01[i] = from->dih01[i];
		to->e_dih3[i] = from->e_dih3[i];
		to->dih03[i] = from->dih03[i];
	}

	to->g_r0hard = from->g_r0hard;
	to->g_ehomo = from->g_ehomo;
	to->g_rhomo = from->g_rhomo;
	to->g_anglek = from->g_anglek;
	to->g_angle0 = from->g_angle0;
	to->g_dihk1 = from->g_dihk1;
	to->g_dih01 = from->g_dih01;
	to->g_dihk3 = from->g_dihk3;
	to->g_dih03 = from->g_dih03;
}

/********************************************************************
 Check overlaps
 	 pr=1	print overlaps
 	 pr=0   only return 1 if something overlaps
 ********************************************************************/

int CheckOverlaps(struct s_polymer *p, struct s_potential *u, struct s_mc_parms *parms, int npol, int nosidechains, int pr, FILE *fproc)
{
      int o=0,tooclose;
      int i,j,ci,cj,a1,a2,is,js,ia1,ia2,iaa1,iaa2;
      char atom1[5],am1[5];
      char atom2[5],am2[5];
      double r2,r02;

      if (pr) fprintf(fproc,"\nChecking overlaps... r0hard = %lf\n",u->g_r0hard);
           // backbone loop
      for (ci=0;ci<npol;ci++)                                     // loops on polymers
            for (cj=ci;cj<npol;cj++)
                  for (i=0;i<(p+ci)->nback;i++)             // loops on backbone atoms of each polymer
                        for (j=0;j<(p+cj)->nback;j++)
                              if (ci!=cj || (i<j && (((p+ci)->back)+i)->iaa != (((p+cj)->back)+j)->iaa))                      // in the same polymer, do not count twice a contact
                              {
                                    tooclose = 0;
                                    if (ci==cj && j-i <= u->g_imin ) tooclose = 1;

                                    // interaction backbone-backbone (i-j)
                                    if (ci!=cj || (i!=j+1 && i!=j-1))                                                   // if not consecutive (i.e., no covalent interaction)
                                    {
                                          a1 = (((p+ci)->back)+i)->itype;
                                          a2 = (((p+cj)->back)+j)->itype;
                                          r2 = Dist2( (((p+ci)->back)+i)->pos, (((p+cj)->back)+j)->pos );
                                          if (tooclose) r02 = u->g_r0hard * u->g_r0hard;
                                          else r02 = (u->r0_2)[a1][a2];
                                          if ( r2 < r02 )
                                          {
                                                strcpy(atom1, (((p+ci)->back)+i)->type);
                                                ia1= (((p+ci)->back)+i)->ia;
                                                strcpy(am1, (((p+ci)->back)+i)->aa);
                                                iaa1= (((p+ci)->back)+i)->iaa;
                                                strcpy(atom2,(((p+cj)->back)+j)->type);
                                                ia2= (((p+cj)->back)+j)->ia;
                                                strcpy(am2,(((p+cj)->back)+j)->aa);
                                                iaa2= (((p+cj)->back)+j)->iaa;

                                                o=1;
                                                if (pr) fprintf(fproc,"WARNING: overlap between backbone %d (chain %d, type %d, atom %s, id %d of amino acid %s %d) and backbone %d (chain %d, type %d, atom %s, id %d of amino acid %s %d): d=%lf d0=%lf\n",
                                                            i,ci,a1,atom1,ia1,am1,iaa1,j,cj,a2,atom2,ia2,am2,iaa2,sqrt(r2),sqrt((u->r0_2)[a1][a2]));
                                          }
                                    }

                                    for (js=0;js< (((p+cj)->back)+j)->nside;js++)
                                    {
                                          // interaction backbone-sidechain (i-js)
                                          a1 = (((p+ci)->back)+i)->itype;
                                          a2 = (((((p+cj)->back)+j)->side)+js)->itype;
                                          r2 = Dist2( (((p+ci)->back)+i)->pos, (((((p+cj)->back)+j)->side)+js)->pos );
                                          if (tooclose) r02 = u->g_r0hard * u->g_r0hard;
                                          else r02 = (u->r0_2)[a1][a2];
                                          if ( r2 < r02 )
                                          {
                                                strcpy(atom1, (((p+ci)->back)+i)->type);
                                                ia1= (((p+ci)->back)+i)->ia;
                                                strcpy(am1, (((p+ci)->back)+i)->aa);
                                                iaa1= (((p+ci)->back)+i)->iaa;
                                                strcpy(atom2, (((((p+cj)->back)+j)->side)+js)->type);
                                                ia2= (((((p+cj)->back)+j)->side)+js)->ia;
                                                strcpy(am2, (((p+cj)->back)+j)->aa);
                                                iaa2= (((p+cj)->back)+j)->iaa;

                                                o=1;
                                                if (pr) fprintf(fproc,"WARNING: overlap between backbone %d (chain %d, type %d, atom %s, id %d of amino acid %s %d) and sidechain %d of backbone %d (chain %d, type %d, atom %s, id %d of amino acid %s %d): d=%lf d0=%lf\n",
                                                                                          i,ci,a1,atom1,ia1,am1,iaa1,js,j,cj,a2,atom2,ia2,am2,iaa2,sqrt(r2),sqrt((u->r0_2)[a1][a2]));
                                          }

                                          for (is=0;is< (((p+ci)->back)+i)->nside;is++)
                                          {
                                                // interaction sidechain-sidechain (is-js)
                                                a1 = (((((p+ci)->back)+i)->side)+is)->itype;
                                                a2 = (((((p+cj)->back)+j)->side)+js)->itype;
                                                r2 = Dist2( (((((p+ci)->back)+i)->side)+is)->pos, (((((p+cj)->back)+j)->side)+js)->pos );
                                                if (tooclose) r02 = u->g_r0hard * u->g_r0hard;
                                                else r02 = (u->r0_2)[a1][a2];
                                                if ( r2 < r02 )
                                                {
                                                      strcpy(atom1, (((((p+ci)->back)+i)->side)+is)->type);
                                                      ia1= (((((p+ci)->back)+i)->side)+is)->ia;
                                                      strcpy(am1, (((p+ci)->back)+i)->aa);
                                                      iaa1= (((p+ci)->back)+i)->iaa;
                                                      strcpy(atom2, (((((p+cj)->back)+j)->side)+js)->type);
                                                      ia2= (((((p+cj)->back)+j)->side)+js)->ia;
                                                      strcpy(am2, (((p+cj)->back)+j)->aa);
                                                      iaa2= (((p+cj)->back)+j)->iaa;

                                                      o=1;
                                                      if (pr)
						      {
						      		fprintf(fproc,"WARNING: overlap between sidechain %d of backbone %d (chain %d, type %d, atom %s, id %d of amino acid %s %d) and sidechain %d of backbone %d (chain %d, type %d, atom %s, id %d of amino acid %s %d): d=%lf d0=%lf\n",is,i,ci,a1,atom1,ia1,am1,iaa1,js,j,cj,a2,atom2,ia2,am2,iaa2,sqrt(r2),sqrt((u->r0_2)[a1][a2]));
						//		fprintf(fproc,"irot 1 = %d \t irot2 = %d \n",(((p+ci)->back)+i)->irot,(((p+cj)->back)+j)->irot);
						      
						//		fprintf(fproc,"OVER: move %d\n",parms->mov);	


							}
                                                }
				}

                                    }

                                    // interaction backbone-sidechain (is-j)
                                    for (is=0;is< (((p+ci)->back)+i)->nside;is++)
                                    {
                                                a1 = (((p+cj)->back)+j)->itype;
                                                a2 = (((((p+ci)->back)+i)->side)+is)->itype;
                                                r2 = Dist2( (((p+cj)->back)+j)->pos, (((((p+ci)->back)+i)->side)+is)->pos );
                                                if (tooclose) r02 = u->g_r0hard * u->g_r0hard;
                                                else r02 = (u->r0_2)[a1][a2];
                                                if ( r2 < r02 )
                                                {
                                                      strcpy(atom1, (((p+cj)->back)+j)->type);
                                                      ia1= (((p+cj)->back)+j)->ia;
                                                      strcpy(am1, (((p+cj)->back)+j)->aa);
                                                      iaa1= (((p+cj)->back)+j)->iaa;
                                                      strcpy(atom2, (((((p+ci)->back)+i)->side)+is)->type);
                                                      ia2= (((((p+ci)->back)+i)->side)+is)->ia;
                                                      strcpy(am2, (((p+ci)->back)+i)->aa);
                                                      iaa2= (((p+ci)->back)+i)->iaa;

                                                      o=1;
                                                      if (pr) fprintf(fproc,"WARNING: overlap between sidechain %d of backbone %d (chain %d, type %d, atom %s, id %d of amino acid %s %d) and backbone %d (chain %d, type %d, atom %s, id %d of amino acid %s %d): d=%lf d0=%lf\n",
                                                                                                                  is,i,ci,a1,atom1,ia1,am1,iaa1,j,cj,a2,atom2,ia2,am2,iaa2,sqrt(r2),sqrt((u->r0_2)[a1][a2]));

						
                                                }
                                    }

                              }
      return o;

}



/********************************************************************
 Energy of chains in a box, returns 1 if any backbone atom is outside
 the box, 0 iw all right.
 'c'=cubic box, 's'=spherical box
 ********************************************************************/
int EnergyBox(struct s_polymer *p, struct s_potential *u, int iw, int ic)
{
	double d2;

	// cubic box
	if (u->boxtype == 'c')
	{
		if ( ((((p+ic)->back)+iw)->pos).x < -u->boxsize || ((((p+ic)->back)+iw)->pos).x > u->boxsize ) return 1;
		if ( ((((p+ic)->back)+iw)->pos).y < -u->boxsize || ((((p+ic)->back)+iw)->pos).y > u->boxsize ) return 1;
		if ( ((((p+ic)->back)+iw)->pos).z < -u->boxsize || ((((p+ic)->back)+iw)->pos).z > u->boxsize ) return 1;
	}

	// spherical box
	else if (u->boxtype == 's')
	{
		d2 = (((((p+ic)->back)+iw)->pos).x * ((((p+ic)->back)+iw)->pos).x) + (((((p+ic)->back)+iw)->pos).y * ((((p+ic)->back)+iw)->pos).y) +
			(((((p+ic)->back)+iw)->pos).z * ((((p+ic)->back)+iw)->pos).z);
		if (d2 > u->boxsize * u->boxsize) return 1;
	}

	return 0;
}

int EnergyBoxPolymer(struct s_polymer *p, struct s_potential *u, int ic)
{
	int i;

	for (i=0;i<(p+ic)->nback;i++)
		if (EnergyBox(p,u,i,ic)==1) return 1;

	return 0;

}
