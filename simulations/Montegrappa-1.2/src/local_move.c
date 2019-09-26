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



struct s_polymer *Allo_Fragment(int npol, int nback,int nangmax, FILE *flog)
{
	
	int ipol;
	struct s_polymer *f;

	f=(struct s_polymer *)calloc(npol,sizeof(struct s_polymer));
	if(!f)      Error("Cannot Allocate fragment");


	for(ipol=0;ipol<npol;ipol++)
	{
		
		(f+ipol)->nback=nback;
		(f+ipol)->back=(struct s_back *)calloc((nback+3),sizeof(struct s_back));
		if(!(f+ipol)->back)     Error("\tAlloFragment(): Cannot Allocate backbone");

                (f+ipol)->A=AlloDoubleMatrix(nangmax,nangmax);
                (f+ipol)->G=AlloDoubleMatrix(nangmax,nangmax);
                (f+ipol)->L=AlloDoubleMatrix(nangmax,nangmax);
                (f+ipol)->Y=AlloDoubleMatrix(nangmax,nangmax);
                (f+ipol)->g_ang=AlloDouble(nangmax);      
                (f+ipol)->d_ang=AlloDouble(nangmax);
	
	}

	return f;

}

double DetTriang( double **Mat, int n )
{

	int i;
	double det=1.0;
	
	for(i=0; i<n; i++)
   		det = det * Mat[i][i];

	return det;

}


void FreeFragment(struct s_polymer *fragment,struct s_mc_parms *parms)
{

	int i;
	for(i=0;i<parms->npol;i++)
        {
	        FreeDoubleMatrix((fragment+i)->A,parms->nmul_local);
                FreeDoubleMatrix((fragment+i)->G,parms->nmul_local);
                FreeDoubleMatrix((fragment+i)->L,parms->nmul_local);
                FreeDoubleMatrix((fragment+i)->Y,parms->nmul_local);
                free((fragment+i)->d_ang );
                free((fragment+i)->g_ang);
		free((fragment+i)->back);	
	}

	free(fragment);

}


void InvertTriang( double **Inv, double **Mat, int n )
{
	int i,k,a;
	
	for(i=0; i<n; i++)
   		for(k=0; k<n; k++)
     			Inv[i][k]=0.0;

	for(i=0; i<n; i++)
   		Inv[i][i] = 1./Mat[i][i];

	for(i=0; i<n; i++)
  		for(a=1; a<=i; a++)
  		{
     			for( k= i-a; k<i; k++)
       				Inv[i][i-a] += (Mat[i][k] * Inv[k][i-a]);
       			Inv[i][i-a] /= -Mat[i][i];
  		}

	return;

}




double Squared_n_Norma( double *vect, int dim )
{
	int i;
	double norm=0.0;

	for(i=0; i<dim; i++)
    		norm += vect[i]*vect[i];
	
	return norm;
}

void TransposedMatOnVect( double **Mat, double *vect, double *risu, int n)
{
       int i,j;

      for(i=0; i<n; i++)
      {
            risu[i] = 0.0;
            for(j=0; j<n; j++)
                  risu[i] += Mat[j][i] * vect[j];
      }

      return;
}

int Gaussian_Angles(double *angles,int n)
{
	int i;
        double r1,r2;

      	for(i=0;i<n;i++)
        {
                r1=(double)rand()/(double)RAND_MAX;
                r2=(double)rand()/(double)RAND_MAX;
                angles[i]=sqrt(-log(r1))*cos(2*PI*r2);
      	}

     	return 0;

}

void MatA( double **A, double **G, int dim,double a,double b)
{

      int i, j;

      for(i=0; i<dim; i++)
            for(j=0; j<dim; j++)
            {
                   
		if(i==j)
                	A[i][i]= a*(1+ b*G[i][j])/2;
                else
                	A[i][j]= a*b*G[i][j]/2;
            
              }
                    return;
}


void Cholesky_2( double **L, double **A, int dim)
{
      int i,j,k;
      double sum;

      for( i=0; i<dim; i++)
            for( j=i; j<dim; j++)
            {
                  sum=A[i][j];
                  for(k=0;k<i;k++)
                        sum += -L[i][k]*L[j][k];
                  if(i == j)
                  {
                        if(sum<=0.0)
                              Error("Cholesky decomposition failed (check the BGS parameters!)");
                        L[i][i]= sqrt(sum);
                  }
                  else
                        L[j][i]= sum/L[i][i];
            }

      return;
}
      

int CopyFragment(struct s_polymer *p,struct s_polymer *f,int iw,int nmul,int natom_fragment,int ip)
{
        int ok=1;
        int i;

	for(i=0;i<natom_fragment;i++)
        {
        	CopyResidueCoordinates( (((p+ip)->back)+iw+i), (((f+ip)->back)+i));
        }

        return ok;
}


void CopyResidueCoordinates(struct s_back *from, struct s_back *to)
{
      (to)->nside=(from)->nside;
      (to->pos).x = (from->pos).x;
      (to->pos).y = (from->pos).y;
      (to->pos).z = (from->pos).z;
      (to)->move=(from)->move;
      strncpy(to->type,from->type,5);
}


int B_Metropolis(double deltaE,double temp,double WN,double WD,struct s_tables *t)
{
	if (deltaE<=0)
        	return 1;
	
	double p;
	p=FastExp(-deltaE/temp,t);
	p=p*(WN/WD);
    
	double random=frand();

	if(random<p)
	{
        	return 1;
      	}

     return 0;
}




      
int Compute_G(struct s_polymer *fragment,struct s_polymer *p,int ip,int iw,int nmul,int nang,struct s_mc_parms *parms)
{
	int ok=1;
	int natom_fragment=nmul+3;
        int i,j,l,i_ang;


      	int k1;
      	int k2;
      	int k3;
      	double deriv1[nang][3],deriv2[nang][3],deriv3[nang][3];
      	double x1,y1,z1,x2,y2,z2,x3,y3,z3;
      	double fact=(1./(LM_DELTAA));


      	CopyFragment(p,fragment,iw,nmul,natom_fragment,ip);


      	k1=nmul-2;
      	k2=nmul-1;
      	k3=nmul;


      	x1=(((fragment+ip)->back)+(k1))->pos.x;
      	y1=(((fragment+ip)->back)+(k1))->pos.y;
	z1=(((fragment+ip)->back)+(k1))->pos.z;
      	x2=(((fragment+ip)->back)+(k2))->pos.x;
      	y2=(((fragment+ip)->back)+(k2))->pos.y;
      	z2=(((fragment+ip)->back)+(k2))->pos.z;
      	x3=(((fragment+ip)->back)+(k3))->pos.x;
      	y3=(((fragment+ip)->back)+(k3))->pos.y;
      	z3=(((fragment+ip)->back)+(k3))->pos.z;
      	i_ang=0;
      

      	for(i=0;i<natom_fragment;i++)
      	{
        	if( (((fragment+ip)->back)+i+1)->move==1 )
            	{
                	PivotForward((fragment+ip),i,LM_DELTAA,natom_fragment-i-2,parms);

                  deriv1[i_ang][0]=(( (((fragment+ip)->back)+(k1))->pos.x )-x1);
                  deriv1[i_ang][1]=( (((fragment+ip)->back)+(k1))->pos.y )-y1;
                  deriv1[i_ang][2]=( (((fragment+ip)->back)+(k1))->pos.z )-z1;
                  deriv2[i_ang][0]=( (((fragment+ip)->back)+(k2))->pos.x )-x2;
                  deriv2[i_ang][1]=( (((fragment+ip)->back)+(k2))->pos.y )-y2;
                  deriv2[i_ang][2]=( (((fragment+ip)->back)+(k2))->pos.z )-z2;
                  deriv3[i_ang][0]=( (((fragment+ip)->back)+(k3))->pos.x )-x3;
                  deriv3[i_ang][1]=( (((fragment+ip)->back)+(k3))->pos.y )-y3;
                  deriv3[i_ang][2]=( (((fragment+ip)->back)+(k3))->pos.z )-z3;

                  deriv1[i_ang][0]=(deriv1[i_ang][0]*fact);
                  deriv1[i_ang][1]=(deriv1[i_ang][1]*fact);
                  deriv1[i_ang][2]=(deriv1[i_ang][2]*fact);
                  deriv2[i_ang][0]=(deriv2[i_ang][0]*fact);
                  deriv2[i_ang][1]=(deriv2[i_ang][1]*fact);
                  deriv2[i_ang][2]=(deriv2[i_ang][2]*fact);
                  deriv3[i_ang][0]=(deriv3[i_ang][0]*fact);
                  deriv3[i_ang][1]=(deriv3[i_ang][1]*fact);
                  deriv3[i_ang][2]=(deriv3[i_ang][2]*fact);

                 CopyFragment(p,fragment,iw,nmul,natom_fragment,ip);

                  i_ang++;
                  if(i_ang==(nang)) break;
            }
      }



     for(i=0;i<nang;i++)
            for(j=0;j<nang;j++)
		(fragment+ip)->G[i][j]=0;
		


      for(i=0;i<nang;i++)
            for(j=0;j<nang;j++)
                  for(l=0;l<3;l++)
      			{                  
			(fragment+ip)->G[i][j]+=(deriv1[i][l]*deriv1[j][l])+(deriv2[i][l]*deriv2[j][l])+(deriv3[i][l]*deriv3[j][l]);

			}

      return ok;
}



int LocalMove(struct s_polymer *p, struct s_polymer *oldp,struct s_polymer *fragment,struct s_potential *pot,int nmul,struct s_mc_parms *parms, double t)
{
	int ok=0,nang=0;
	int iw,ip,iapdbtocheck,i,m;
	double deltaE;
	double detL1,detL2,detA1,detA2,W1,W2,psisquared,e;


	ip=irand(parms->npol);

	iw=(nmul-1)+irand( (p+ip)->nback-nmul+1   );
	//iw=20;
	

	for(i=0;i<nmul;i++)
	{
		nang+=((((p+ip)->back)+iw-nmul+i)+1)->move;
	}

	deltaE = -GetEnergyMonomerRange(p,iw-nmul+1,iw+1,ip);
  

	if(iw==((p+ip)->nback-1)) //pivot forward OK
	{


		for(i=0;i<nang;i++)
			(fragment+ip)->d_ang[i]=parms->dw_mpivot*(0.05-frand());
		if (!parms->nodihpot)
                	for(i=0;i<nmul+3;i++)
			{
				deltaE-=(((p+ip)->back)+iw-nmul+i+1)->e_dih;
			}
		m=0;
		for(i=0;i<nmul;i++)
		{
			if( (((p+ip)->back)+iw+i-nmul+2)->move==1)
                	{
	             		ok*=PivotForward((p+ip),iw-nmul+i+1,(fragment+ip)->d_ang[m],nmul-i-2,parms);
                     		m++;
                     		if(m==nang) break;
                	}

		}
		if(!parms->nosidechains)
		{
			ok*=AddSidechain(p,iw-nmul+1,iw+1,ip);
		}

		deltaE += EnergyMonomerRange(p,pot,iw-nmul+1,iw+1,ip,parms->npol,parms->shell,1,parms->nosidechains,parms->disentangle,parms->hb);
		if (!parms->nodihpot)
                	for(i=0;i<nmul+3;i++)
			{

				deltaE+=EnergyDihedrals(p,pot,iw-nmul+i+1,ip,1);
			}
		
		ok*=Metropolis(deltaE,t,p->tables);

		if(ok==0)
		{
			UpdateMonomerRange(oldp,p,iw-nmul+1,iw+1,ip,parms->shell);
        	}    
		else
		{
			//move accepted

                	UpdateMonomerRange(p,oldp,iw-nmul+1,iw+1,ip,parms->shell);
             		p->etot+=deltaE;
                	parms->acc++;
		}


	} //end of forward
	



		
	else if(iw==nmul-1) //pivot backward OK
	{

		for(i=0;i<nang;i++)
                        (fragment+ip)->d_ang[i]=parms->dw_mpivot*(0.05-frand());
                if (!parms->nodihpot)
                        for(i=0;i<nmul+3;i++)
			{
                                deltaE-=(((p+ip)->back)+i)->e_dih;

			}
	        m=0;

                for(i=0;i<nmul;i++)
                {
                        if( (((p+ip)->back)+iw-i)->move==1)
                        {
		                ok*=PivotBackward((p+ip),iw-i,(fragment+ip)->d_ang[m],nmul-i-2,parms);
                                m++;
                                if(m==nang) break;
                        }

                }
		if(!parms->nosidechains)
			ok*=AddSidechain(p,0,iw+1,ip);
                deltaE += EnergyMonomerRange(p,pot,0,iw+1,ip,parms->npol,parms->shell,1,parms->nosidechains,parms->disentangle,parms->hb);
                if (!parms->nodihpot)
                        for(i=0;i<nmul+3;i++)
			{
		        	deltaE+=EnergyDihedrals(p,pot,i,ip,1);
			}
                ok*=Metropolis(deltaE,t,p->tables);

                if(ok==0)
                {
 
		       UpdateMonomerRange(oldp,p,0,iw+1,ip,parms->shell);
                }
		else
                {
                        UpdateMonomerRange(p,oldp,0,iw+1,ip,parms->shell);

		        p->etot+=deltaE;
                        parms->acc++;
                }


	} //end of backward



	else if(iw!=((p+ip)->nback-2))//pivot local
	{



		if((((p+ip)->back)+iw+1)->iapdb==0)
			iapdbtocheck=1;
		if((((p+ip)->back)+iw+1)->iapdb==1)
			iapdbtocheck=2;
		if((((p+ip)->back)+iw+1)->iapdb==2)
			iapdbtocheck=0;	
		int out;

		Gaussian_Angles((fragment+ip)->g_ang,nang);
		Compute_G(fragment,p,ip,iw-nmul+1,nmul,nang,parms);	


		MatA((fragment+ip)->A,(fragment+ip)->G,nang,parms->bgs_a,parms->bgs_b);
		Cholesky_2((fragment+ip)->L,(fragment+ip)->A,nang);
		psisquared=Squared_n_Norma((fragment+ip)->g_ang,nang);
                e=exp(-psisquared);
		detL1=DetTriang((fragment+ip)->L,nang);
		detA1=detL1*detL1;
		W1=e * sqrt(detA1);
		InvertTriang((fragment+ip)->Y,(fragment+ip)->L,nang);
		TransposedMatOnVect((fragment+ip)->Y,(fragment+ip)->g_ang,(fragment+ip)->d_ang,nang);

		if(!parms->nodihpot)
			for(i=0;i<nmul+3;i++)
			{
				deltaE -= (((p+ip)->back)+iw-nmul+i+1)->e_dih;
			}

			if(!parms->noangpot)
			{		
		
                	deltaE-=(((p+ip)->back)+iw)->e_ang;
  	            	deltaE-=(((p+ip)->back)+iw+1)->e_ang;
			}


		m=0;


		for(i=0;i<nmul;i++)
		{
			if( (((p+ip)->back)+iw+i-nmul+2)->move==1)
			{
				ok*=PivotForward((p+ip),iw-nmul+i+1,(fragment+ip)->d_ang[m],nmul-i-2,parms);
				m++;
				if(m==nang)
				break;
                  	}
	
		}

                                                                   

		
		double rc2=parms->r_cloose*parms->r_cloose;
		double dihedral=Dihedral( (((p+ip)->back)+iw-iapdbtocheck)->pos, (((p+ip)->back)+iw+1-iapdbtocheck)->pos, (((p+ip)->back)+iw+2-iapdbtocheck)->pos, (((p+ip)->back)+iw+3-iapdbtocheck)->pos, p->tables, &out );
//		fprintf(stderr,"Checking Dihedral: backbone atoms %d,%d,%d,%d\n",+iw-iapdbtocheck,+iw+1-iapdbtocheck,+iw+2-iapdbtocheck,+iw+3-iapdbtocheck);
	
		
		if( DAbs (Dist2( (((p+ip)->back)+iw)->pos,(((p+ip)->back)+iw+1)->pos ) -(((p+ip)->back)+iw)->d2_next < rc2 )  &&   DAbs( Angle( (((p+ip)->back)+iw-1)->pos, (((p+ip)->back)+iw)->pos, (((p+ip)->back)+iw+1)->pos, (p+ip)->tables, &out)  - (((p+ip)->back)+iw-1)->a_next) < 0.05 && DAbs(180-DAbs(dihedral)) < DELTAOMEGA)
		{		

			if(!parms->nosidechains)
			{
				ok *= AddSidechain(p,iw-nmul+1,iw+1,ip);
			}
	
			deltaE += EnergyMonomerRange(p,pot,iw-nmul+1,iw+1,ip,parms->npol,parms->shell,1,parms->nosidechains,parms->disentangle,parms->hb);

			if (!parms->nodihpot)
        	        	for(i=0;i<nmul+3;i++)
				{
	                      		deltaE+=EnergyDihedrals(p,pot,iw-nmul+i+1,ip,1);
				}                  

                  	if(!parms->noangpot)
                  	{
				//		fprintf(stderr,"\nCOMPUTING ANGLE ENERGY deltaE=%lf\t",deltaE);
        			deltaE+=EnergyAngles(p,pot,iw,ip,1);
                       		deltaE+=EnergyAngles(p,pot,iw+1,ip,1);
				//		fprintf(stderr,"-> %lf\n",deltaE);
 	              	}

       		

                  	Compute_G(fragment,p,ip,iw-nmul+1,nmul,nang,parms);
 	          	MatA((fragment+ip)->A,(fragment+ip)->G,nang,parms->bgs_a,parms->bgs_b);
			Cholesky_2((fragment+ip)->L,(fragment+ip)->A,nang);
			detL2=DetTriang((fragment+ip)->L,nang);
			detA2=detL2*detL2;
			TransposedMatOnVect((fragment+ip)->L,(fragment+ip)->d_ang,(fragment+ip)->g_ang,nang);
			psisquared=Squared_n_Norma((fragment+ip)->g_ang,nang);
			e=exp(-psisquared);
			W2=e*sqrt(detA2);
			ok*=B_Metropolis(deltaE,t,W2,W1,p->tables);

			if(ok==0) //move rejected
                  	{
      				UpdateMonomerRange(oldp,p,iw-nmul+1,iw+1,ip,parms->shell);
			}
                  	else	//move accepted
                  	{

                	        UpdateMonomerRange(p,oldp,iw-nmul+1,iw+1,ip,parms->shell);
                		
			        p->etot+=deltaE;
       			
		                parms->acc++;
                  	}
			
			
		}//end of loose condition

		else //if not loose pivot
        	{
			ok=0;
	        	UpdateMonomerRange(oldp,p,iw-nmul+1,iw+1,ip,parms->shell);
            	}


	}//end of local pivot


	parms->mov++;

	return ok;

}	

