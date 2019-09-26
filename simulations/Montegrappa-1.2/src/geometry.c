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
 * geometry.c
 *
 *  Created on: Sep 21, 2010
 *      Author: guido
 */

#include "montegrappa.h"

/*****************************************************************************
 Make the flip move of backbone atom iw of angle dw
	 returns 1 if it could make the move
	 returns 0 if nothing could be done
 *****************************************************************************/
int Flip(struct s_polymer *p, int iw, double dw)
{
	struct vector a,b,u,du,newu;
	double invnormb2,coeff,invnormb,norm2b,invnormT,norm2T,T[3][3];

	if (iw<1 || iw>p->nback-2) return 0;		// cannot apply to first and last atom

	// preceding atom
	a.x = (((p->back)+iw-1)->pos).x;
	a.y = (((p->back)+iw-1)->pos).y;
	a.z = (((p->back)+iw-1)->pos).z;

	// following atom
	b.x = (((p->back)+iw+1)->pos).x;
	b.y = (((p->back)+iw+1)->pos).y;
	b.z = (((p->back)+iw+1)->pos).z;

	// atom to be moved around the axis passing by a and b
	u.x = (((p->back)+iw)->pos).x;
	u.y = (((p->back)+iw)->pos).y;
	u.z = (((p->back)+iw)->pos).z;


	// Translate b and u to the system defined by a
	b.x -= a.x;
	b.y -= a.y;
	b.z -= a.z;

	u.x -= a.x;
	u.y -= a.y;
	u.z -= a.z;

	// db = b - u
	du.x = b.x - u.x;
	du.y = b.y - u.y;
	du.z = b.z - u.z;

	// coefficient for the 2nd column of base-change matrix
	norm2b = Norm2(b);
	if (norm2b<EPSILON) return 0;
	invnormb2 = 1. / norm2b;			// inverse of the square norm of b (translated by a)
	if (invnormb2<EPSILON) return 0;		// PARANOID
	coeff = 0.5 * (( Norm2(u) - Norm2(du) ) * invnormb2 + 1.);

	// base-change matrix (first two columns)
	invnormb = FastSqrt(invnormb2,p->tables);

	T[0][0] = b.x * invnormb;		T[0][1] = u.x - coeff * b.x;
	T[1][0] = b.y * invnormb;		T[1][1] = u.y - coeff * b.y;
	T[2][0] = b.z * invnormb;		T[2][1] = u.z - coeff * b.z;

	// normalize 2nd column of the matrix
	norm2T = T[0][1]*T[0][1] + T[1][1]*T[1][1] + T[2][1]*T[2][1];
	if (norm2T<EPSILON) return 0;
	invnormT = 1. / norm2T;

	invnormT = FastSqrt(invnormT,p->tables);
	if (invnormT<EPSILON) return 0;	// PARANOIND
	T[0][1] = T[0][1] * invnormT;
	T[1][1] = T[1][1] * invnormT;
	T[2][1] = T[2][1] * invnormT;

	// the third colum is the vector product of the first two
	T[0][2] = T[1][0] * T[2][1] - T[2][0] * T[1][1];
	T[1][2] = T[2][0] * T[0][1] - T[0][0] * T[2][1];
	T[2][2] = T[0][0] * T[1][1] - T[1][0] * T[0][1];

	// rotate u
	u.x = coeff / invnormb;
	u.y = FastCos(dw,p->tables) / invnormT;
	u.z = FastSin(dw,p->tables) / invnormT;

	MatrixVectorProduct(T,u,&newu);

	// translate back in a
	newu.x += a.x;
	newu.y += a.y;
	newu.z += a.z;

	// actually moves the atom
	(((p->back)+iw)->pos).x = newu.x;
	(((p->back)+iw)->pos).y = newu.y;
	(((p->back)+iw)->pos).z = newu.z;

	return 1;
}

int MoveHead(struct s_polymer *p, struct s_mc_parms *mc_parms)
{
	struct vector c;
	double dt;
	int it;

	// the new origin is the second atom (i.e., +1)
	c.x = (((p->back)+1)->pos).x;
	c.y = (((p->back)+1)->pos).y;
	c.z = (((p->back)+1)->pos).z;

	// shift head atom
	(((p->back)+0)->pos).x -= c.x;
	(((p->back)+0)->pos).y -= c.y;
	(((p->back)+0)->pos).z -= c.z;

	// make a random rotation around a random axis
	dt = 2. * mc_parms->dw_flip * (frand() - 0.5);
	it = irand(3);

	((p->back)+0)->pos = RotateVector(((p->back)+0)->pos,dt,it,p->tables);

	// tranlsate back
	(((p->back)+0)->pos).x += c.x;
	(((p->back)+0)->pos).y += c.y;
	(((p->back)+0)->pos).z += c.z;

	return 1;
}

int MoveTail(struct s_polymer *p, struct s_mc_parms *mc_parms)
{
	struct vector c;
	double dt;
	int it,last;

	// tail atom
	last = p->nback - 1;

	// the new origin is the next-to-tail atom
	c.x = (((p->back)+last-1)->pos).x;
	c.y = (((p->back)+last-1)->pos).y;
	c.z = (((p->back)+last-1)->pos).z;

	// shift tail atom
	(((p->back)+last)->pos).x -= c.x;
	(((p->back)+last)->pos).y -= c.y;
	(((p->back)+last)->pos).z -= c.z;

	// make a random rotation around a random axis
	dt = 2. * mc_parms->dw_flip * (frand() - 0.5);
	it = irand(3);

	((p->back)+last)->pos = RotateVector(((p->back)+last)->pos,dt,it,p->tables);

	// tranlsate back
	(((p->back)+last)->pos).x += c.x;
	(((p->back)+last)->pos).y += c.y;
	(((p->back)+last)->pos).z += c.z;

	return 1;
}

/*****************************************************************************
  Spherical coordinates of D given A,B,C
	  angles in degrees
 *****************************************************************************/
struct angles Cartesian2Spherical( struct vector A, struct vector B,
                                       struct vector C, struct vector D, struct s_tables *tables, int *out)
{
  struct angles X;
  double length=0., dotprod, BC;
  double lengthV, lengthU, lengthprod;
  double temp;
  struct vector V,U;

  *out=0;

  /* Bond length */
  length = (D.x-C.x)*(D.x-C.x) + (D.y-C.y)*(D.y-C.y) + (D.z-C.z)*(D.z-C.z);
  X.r = FastSqrt(length,tables);
  //if (X.r<EPSILON) Error("C and D overlapping in Cartesian2Spherical");
  if (X.r<EPSILON) { fprintf(stderr,"WARNING: C and D overlapping in Cartesian2Spherical\n"); return X; }

  /* Angle beta (BCD) */
  BC = (C.x-B.x)*(C.x-B.x) + (C.y-B.y)*(C.y-B.y) + (C.z-B.z)*(C.z-B.z);
  BC = FastSqrt(BC,tables);
  //if (BC<EPSILON) Error("B and C identical in Cartesian2Spherical");
  if (BC<EPSILON) { fprintf(stderr,"WARNING: B and C identical in Cartesian2Spherical\n");  return X; }

  dotprod = (B.x-C.x)*(D.x-C.x) + (B.y-C.y)*(D.y-C.y) + (B.z-C.z)*(D.z-C.z);

  temp = (dotprod/(BC*X.r));

  //if(temp >=  1)   Error( "temp>1 in Cartesian2Spherical" );
  if(temp>=1) return X; 

  if(temp <= -1)  { X.ang = 180.;   X.dih = 0.0; }
  else
  {
     X.ang = FastAcos(temp,tables);

    /* Dihedral alpha */
     V.x = (B.y - A.y)*(C.z - B.z) - (C.y - B.y)*(B.z - A.z);   // V = AB x BC
     V.y = (C.x - B.x)*(B.z - A.z) - (B.x - A.x)*(C.z - B.z);
     V.z = (B.x - A.x)*(C.y - B.y) - (C.x - B.x)*(B.y - A.y);

     U.x = -(B.y - C.y)*(D.z - C.z) + (D.y - C.y)*(B.z - C.z);  // U = BC x CD
     U.y = -(D.x - C.x)*(B.z - C.z) + (B.x - C.x)*(D.z - C.z);
     U.z = -(B.x - C.x)*(D.y - C.y) + (D.x - C.x)*(B.y - C.y);

     lengthV = V.x*V.x + V.y*V.y + V.z*V.z;
     lengthU = U.x*U.x + U.y*U.y + U.z*U.z;
     lengthprod = FastSqrt(lengthU*lengthV,tables);

      //if(lengthprod < EPSILON)
    	// Error("3 atoms aligned in AbsoluteToInternal\n ");
      if(lengthprod < EPSILON) { fprintf(stderr,"WARNING: 3 atoms aligned in Cartesian2Spherical\n"); return X; }

     dotprod =  V.x * U.x + V.y * U.y + V.z * U.z;
     temp = (dotprod/lengthprod);
     if( temp < 1 && temp > -1 )  X.dih = FastAcos(temp,tables);
     else if (temp >= 1) X.dih = 0.0;
     else X.dih = 180.;

     if ((B.x - A.x)*U.x + (B.y - A.y)*U.y + (B.z - A.z)*U.z  < 0)
       X.dih = -X.dih;
    }

 *out = 1;
 return X;
}

/*****************************************************************************
  Cartesian coordinates of Z given A,B,C
	  angles in degrees
 *****************************************************************************/
struct vector Spherical2Cartesian(struct vector A, struct vector B,
                                  struct vector C, struct angles W, struct s_tables *tables, int *out)
{
  struct vector Z;
  double dumb[3],sb;
  double ZR[3];
  double M[3][3];
  double l_ab,l_bc,l1,l2,inv_bc,inv_ab,inv_l;
  double l_ac;
  int i,j;

  *out = 0;
  for (i=0;i<3;i++) dumb[i]=0.;

  /* From angular internal coords to cartesian internal coords */
  sb = FastSin(W.ang,tables);
  ZR[0] = - W.r * FastCos(W.ang,tables);
  ZR[1] = - W.r * sb * FastCos(W.dih,tables);
  ZR[2] = - W.r * sb * FastSin(W.dih,tables);

  /* Bond lengths */
  l_ab = (B.x-A.x)*(B.x-A.x)+(B.y-A.y)*(B.y-A.y)+(B.z-A.z)*(B.z-A.z);
  l_ab = FastSqrt(l_ab,tables);
  l_bc = (C.x-B.x)*(C.x-B.x)+(C.y-B.y)*(C.y-B.y)+(C.z-B.z)*(C.z-B.z);
  l_bc = FastSqrt(l_bc,tables);
  
  // Check if BA and CB are not aligned
  l_ac = FastSqrt((C.x-A.x)*(C.x-A.x)+(C.y-A.y)*(C.y-A.y)+(C.z-A.z)*(C.z-A.z),tables);
  if (l_ab + l_bc -l_ac < EPSILON) { fprintf(stderr,"WARNING: BA and CB aligned in Spherical2Cartesian\t(move rejected)\n"); return Z; }
   
   // Check if they coincides
   if ( l_ab < EPSILON)  { fprintf(stderr,"WARNING: A and B overlapping in Spherical2Cartesian\t(move rejected)\n");  return Z; }
   if ( l_ac < EPSILON)  { fprintf(stderr,"WARNING: A and C overlapping in Spherical2Cartesian\t(move rejected)\n");  return Z; }
   if ( l_bc < EPSILON)  { fprintf(stderr,"WARNING: B and C overlapping in Spherical2Cartesian\t(move rejected)\n");  return Z; }

  /* Base of internal coords as function of absolute coords */

  /* E_x */
  inv_bc = 1./l_bc;
  M[0][0] = inv_bc * (C.x - B.x);
  M[0][1] = inv_bc * (C.y - B.y);
  M[0][2] = inv_bc * (C.z - B.z);

  /* E_z */
  inv_ab = 1./l_ab;
  M[2][0] = M[0][1] * (B.z - A.z)*inv_ab - M[0][2] * (B.y - A.y)*inv_ab;
  M[2][1] = -M[0][0] * (B.z - A.z)*inv_ab + M[0][2] * (B.x - A.x)*inv_ab;
  M[2][2] = M[0][0] * (B.y - A.y)*inv_ab - M[0][1] * (B.x - A.x)*inv_ab;
  l2 = FastSqrt(M[2][0]*M[2][0]+M[2][1]*M[2][1]+M[2][2]*M[2][2],tables);
  if (l2<EPSILON)  { fprintf(stderr,"WARNING: l2=0 in Spherical2Cartesian\t(move rejected)\n");  return Z; } //PARANOID
  inv_l = 1. / l2;
  if (inv_l<EPSILON)  { fprintf(stderr,"WARNING: inv_l=0 in Spherical2Cartesian\t(move rejected)\n"); return Z; } //PARANOID
  for (i=0;i<3;i++) M[2][i] *= inv_l;

  /* E_y */
  M[1][0] = -M[0][1]*M[2][2] + M[2][1]*M[0][2];
  M[1][1] = M[0][0]*M[2][2] - M[2][0]*M[0][2];
  M[1][2] = -M[0][0]*M[2][1] + M[0][1]*M[2][0];

  l1 = FastSqrt(M[1][0]*M[1][0]+M[1][1]*M[1][1]+M[1][2]*M[1][2],tables);
  if (l1<EPSILON)  { fprintf(stderr,"WARNING: l1=0 in Spherical2Cartesian\t(move rejected)\n"); return Z; } //PARANOID
  inv_l = 1. / l1;
  if (inv_l<EPSILON)  { fprintf(stderr,"WARNING: inv_l=0 Spherical2Cartesian\t(move rejected)\n"); return Z; } //PARANOID
  for (i=0;i<3;i++) M[1][i] *= inv_l;

  /* Absolute coords = Transpose[M] * relative coords */
  for (i=0;i<3;i++)
   for (j=0;j<3;j++)
    dumb[i] += M[j][i] * ZR[j];

  Z.x = C.x + dumb[0];
  Z.y = C.y + dumb[1];
  Z.z = C.z + dumb[2];

  *out = 1;
  return Z;
}

/*****************************************************************************
 Rotate vector v of theta around an axis (w=0,1,2, that is x, y, z)
 *****************************************************************************/
struct vector RotateVector(struct vector v, double theta, int w, struct s_tables *tables)
{
	double M[3][3],st,ct;
	struct vector u;

	st = FastSin(theta,tables);
	ct = FastCos(theta,tables);

	if (w==0)	//rotation on the x axis
	{
		M[0][0] = 1;	M[0][1] = 0;	M[0][2] = 0;
		M[1][0] = 0;	M[1][1] = ct;	M[1][2] = -st;
		M[2][0] = 0;	M[2][1] = st;	M[2][2] = ct;
	}
	else if (w==1)
	{
		M[0][0] = ct;	M[0][1] = 0;	M[0][2] = st;
		M[1][0] = 0;	M[1][1] = 1;	M[1][2] = 0;
		M[2][0] = -st;	M[2][1] = 0;	M[2][2] = ct;
	}
	else if (w==2)
	{
		M[0][0] = ct;	M[0][1] = -st;	M[0][2] = 0;
		M[1][0] = st;	M[1][1] = ct;	M[1][2] = 0;
		M[2][0] = 0;	M[2][1] = 0;	M[2][2] = 1;
	}
	else Error("w should be in [0,2] in RotateVecotor");

	u.x = M[0][0] * v.x + M[0][1] * v.y + M[0][2] * v.z;
	u.y = M[1][0] * v.x + M[1][1] * v.y + M[1][2] * v.z;
	u.z = M[2][0] * v.x + M[2][1] * v.y + M[2][2] * v.z;

	return u;
}

/********************************************************************
 Pivot n atoms after iw of dih dw (moves from iw+2 on)
 affecting dihedral of iw+1
 ********************************************************************/
int PivotForward(struct s_polymer *p, int iw, double dw, int n, struct s_mc_parms *parms)
{
	
	int i;
	struct vector C,V,U,W;
	double norm,inv_norm,M[3][3],cw,sw,ocw;

	#ifdef DEBUG
	if (iw+n+1 >= p->nback) Error("attempting to pivot beyond the end of the chain");
	#endif

	// the pivoting atom
	C.x = (((p->back)+iw)->pos).x;
	C.y = (((p->back)+iw)->pos).y;
	C.z = (((p->back)+iw)->pos).z;

	// the normalized vector around which to pivot
	V.x = (((p->back)+iw+1)->pos).x - C.x;
	V.y = (((p->back)+iw+1)->pos).y - C.y;
	V.z = (((p->back)+iw+1)->pos).z - C.z;

	norm = FastSqrt(Norm2(V),p->tables);
	if (norm<EPSILON) return 0; 		//PARANOID
	inv_norm = 1. / norm;
	if (inv_norm<EPSILON) return 0; 	//PARANOID
	V.x *= inv_norm;
	V.y *= inv_norm;
	V.z *= inv_norm;


	// the rotation matrix
	cw = FastCos(dw,p->tables);
	sw = FastSin(dw,p->tables);
	ocw = 1. - cw;

	M[0][0] = V.x * V.x * ocw + cw;
	M[0][1] = V.x * V.y * ocw + V.z * sw;
	M[0][2] = V.x * V.z * ocw - V.y * sw;

	M[1][0] = V.x * V.y * ocw - V.z * sw;
	M[1][1] = V.y * V.y * ocw + cw;
	M[1][2] = V.y * V.z * ocw + V.x * sw;

	M[2][0] = V.x * V.z * ocw + V.y * sw;
	M[2][1] = V.y * V.z * ocw - V.x * sw;
	M[2][2] = V.z * V.z * ocw + cw;

	// move n atoms from iw+2 to iw+2+n
	for (i=iw+2;i<iw+2+n;i++)
     // for(i=iw;i<iw+n;i++)
	{
//		fprintf(stderr,"\ni = %d \n",i);
            // shift each to be centred in C
		U.x = (((p->back)+i)->pos).x - C.x;
		U.y = (((p->back)+i)->pos).y - C.y;
		U.z = (((p->back)+i)->pos).z - C.z;

		// rotate
		W.x = M[0][0] * U.x + M[1][0] * U.y + M[2][0] * U.z;
		W.y = M[0][1] * U.x + M[1][1] * U.y + M[2][1] * U.z;
		W.z = M[0][2] * U.x + M[1][2] * U.y + M[2][2] * U.z;

		// shift back
		(((p->back)+i)->pos).x = W.x + C.x;
		(((p->back)+i)->pos).y = W.y + C.y;
		(((p->back)+i)->pos).z = W.z + C.z;
	}

	return 1;
}

/*****************************************************************************
 Calculate the x,y,z coordinates of sidechains associated to backbone atoms
 from istart to istop
 *****************************************************************************/
int AddSidechain(struct s_polymer *p, int istart, int istop, int ch)
{
	struct vector b1,b2,b3;
	struct angles ang;
	int i,iside,irot,ib1,ib2,ib3,out;

	for (i=istart;i<=istop;i++)
		for (iside=0;iside<(((p+ch)->back)+i)->nside;iside++)
		{
			irot = (((p+ch)->back)+i)->irot;								// which is the rotamer to insert
			ib1 = (((((((p+ch)->back)+i)->side)+iside)->rot)+irot)->b1;		// which are the atoms which define the dihedrals
			ib2 = (((((((p+ch)->back)+i)->side)+iside)->rot)+irot)->b2;
			ib3 = (((((((p+ch)->back)+i)->side)+iside)->rot)+irot)->b3;
			b1 = *(*(((p+ch)->vback)+ib1));
			b2 = *(*(((p+ch)->vback)+ib2));
			b3 = *(*(((p+ch)->vback)+ib3));
			ang = (((((((p+ch)->back)+i)->side)+iside)->rot)+irot)->ang;
 			(((((p+ch)->back)+i)->side)+iside)->pos = Spherical2Cartesian(b1,b2,b3,ang,p->tables,&out);
			if (out==0) return 0; 
		}

	return 1;
}

/********************************************************************
 Pivot n atoms before iw of dih dw (moves from iw-2 back);
 affecting dihedral of iw
 ********************************************************************/
int PivotBackward(struct s_polymer *p, int iw, double dw, int n, struct s_mc_parms *parms)
{
	int i;
	struct vector C,V,U,W;
	double norm,inv_norm,M[3][3],cw,sw,ocw;

	if (iw-n <0) Error("attempting to pivot beyond the beginning of the chain");

	// the pivoting atom
	C.x = (((p->back)+iw)->pos).x;
	C.y = (((p->back)+iw)->pos).y;
	C.z = (((p->back)+iw)->pos).z;

	// the normalized vector around which to pivot
	V.x = (((p->back)+iw-1)->pos).x - C.x;
	V.y = (((p->back)+iw-1)->pos).y - C.y;
	V.z = (((p->back)+iw-1)->pos).z - C.z;
	norm = FastSqrt(Norm2(V),p->tables);
	if (norm<EPSILON) return 0; 		//PARANOID
	inv_norm = 1. / norm;
	if (inv_norm<EPSILON) return 0; 		//PARANOID
	V.x *= inv_norm;
	V.y *= inv_norm;
	V.z *= inv_norm;

	// the rotation matrix
	cw = FastCos(dw,p->tables);
	sw = FastSin(dw,p->tables);
	ocw = 1. - cw;

	M[0][0] = V.x * V.x * ocw + cw;
	M[0][1] = V.x * V.y * ocw + V.z * sw;
	M[0][2] = V.x * V.z * ocw - V.y * sw;

	M[1][0] = V.x * V.y * ocw - V.z * sw;
	M[1][1] = V.y * V.y * ocw + cw;
	M[1][2] = V.y * V.z * ocw + V.x * sw;

	M[2][0] = V.x * V.z * ocw + V.y * sw;
	M[2][1] = V.y * V.z * ocw - V.x * sw;
	M[2][2] = V.z * V.z * ocw + cw;

	// move n atoms from iw-2 to iw-2-n
	for (i=iw-2;i>iw-2-n;i--)
	{
		// shift each to be centred in C
		U.x = (((p->back)+i)->pos).x - C.x;
		U.y = (((p->back)+i)->pos).y - C.y;
		U.z = (((p->back)+i)->pos).z - C.z;

		// rotate
		W.x = M[0][0] * U.x + M[1][0] * U.y + M[2][0] * U.z;
		W.y = M[0][1] * U.x + M[1][1] * U.y + M[2][1] * U.z;
		W.z = M[0][2] * U.x + M[1][2] * U.y + M[2][2] * U.z;

		// shift back
		(((p->back)+i)->pos).x = W.x + C.x;
		(((p->back)+i)->pos).y = W.y + C.y;
		(((p->back)+i)->pos).z = W.z + C.z;
	}

	return 1;
}

/********************************************************************
 Pivot all chain around atom iw
	 i.e., if iw is beyond the half, around the axis (iw-1) - iw
		   if iw is before the half, around the axis iw - (iw+1)
     putting move=0 at a give iw fixes the rotation around the axis iw - (iw+1)
 ********************************************************************/
int Pivot(struct s_polymer *p, int iw, double dw, int ip, struct s_mc_parms *parms)
{
	int half;
	int ok;

	half = (p+ip)->nback / 2;

	// if iw belongs to the first half, pivot backward
	if (iw < half)
	{
		if ( (((p+ip)->back)+iw)->move == 0 ) return 0;
		ok = PivotBackward((p+ip),iw+1,dw,iw,parms);
		AddSidechain(p,0,iw-1,ip);
	}
	// if it belong to the second half, pivot forward
	else
	{
		if ( (((p+ip)->back)+iw-1)->move == 0 ) return 0;
		ok = PivotForward((p+ip),iw-1,dw,(p+ip)->nback-iw-1,parms);
		AddSidechain(p,iw+1,p->nback-1,ip);
	}

	return ok;
}

/********************************************************************
 Square distance between two vectors
 ********************************************************************/
double Dist2(struct vector a, struct vector b)
{
	return (a.x-b.x)*(a.x-b.x) + (a.y-b.y)*(a.y-b.y) + (a.z-b.z)*(a.z-b.z);
}

/*****************************************************************************
  Angle between B, C and D
 *****************************************************************************/
double Angle( struct vector B, struct vector C, struct vector D, struct s_tables *tables,int *out)
{
  double ang=0;
  double length2=0., dotprod, BC2;
  double temp;

  *out = 0;

  /* Bond length */
  length2 = (D.x-C.x)*(D.x-C.x) + (D.y-C.y)*(D.y-C.y) + (D.z-C.z)*(D.z-C.z);
//  if (length2<EPSILON) { fprintf(stderr,"C and D Angle"); return ang; }

  /* Angle (BCD) */
  BC2 = (C.x-B.x)*(C.x-B.x) + (C.y-B.y)*(C.y-B.y) + (C.z-B.z)*(C.z-B.z);
//  if (BC2<EPSILON) { fprintf(stderr,"BC=0 in Angle"); return ang; }

  dotprod = (B.x-C.x)*(D.x-C.x) + (B.y-C.y)*(D.y-C.y) + (B.z-C.z)*(D.z-C.z);

  temp = dotprod/FastSqrt(BC2*length2,tables);

//  if(temp >=  1) { fprintf(stderr,"temp>1 in Angle"); return ang; }

  if(temp <= -1)  { ang = 180.;  }
  else ang = FastAcos(temp,tables);

  *out = 1;
  return ang;
}

/*****************************************************************************
  Cos of Angle between B, C and D
 *****************************************************************************/
double CosAngle( struct vector B, struct vector C, struct vector D, struct s_tables *tables, int *out)
{
  double length2=0., dotprod, BC2;
  double temp;

  *out = 0;

  /* Bond length */
  length2 = (D.x-C.x)*(D.x-C.x) + (D.y-C.y)*(D.y-C.y) + (D.z-C.z)*(D.z-C.z);
  if (length2<EPSILON) { fprintf(stderr,"C and D overlapping in CosAngle"); return 0; }

  /* Angle (BCD) */
  BC2 = (C.x-B.x)*(C.x-B.x) + (C.y-B.y)*(C.y-B.y) + (C.z-B.z)*(C.z-B.z);

  dotprod = (B.x-C.x)*(D.x-C.x) + (B.y-C.y)*(D.y-C.y) + (B.z-C.z)*(D.z-C.z);

  temp = dotprod/FastSqrt(BC2*length2,tables);

  if(temp >=  1) { fprintf(stderr,"temp>1 in CosAngle"); return 0; }

  *out = 1;

  if(temp <= -1)  return -1;
  else return temp;
}

/*****************************************************************************
  Dihedral between B, C and D
 *****************************************************************************/
double Dihedral( struct vector A, struct vector B,
                                       struct vector C, struct vector D, struct s_tables *tables, int *out)
{
  double dih=0, dotprod;
  double lengthV, lengthU, lengthprod;
  double temp;
  struct vector V,U;

  *out = 0;

  dotprod = (B.x-C.x)*(D.x-C.x) + (B.y-C.y)*(D.y-C.y) + (B.z-C.z)*(D.z-C.z);

  V.x = (B.y - A.y)*(C.z - B.z) - (C.y - B.y)*(B.z - A.z);   // V = AB x BC
  V.y = (C.x - B.x)*(B.z - A.z) - (B.x - A.x)*(C.z - B.z);
  V.z = (B.x - A.x)*(C.y - B.y) - (C.x - B.x)*(B.y - A.y);

  U.x = -(B.y - C.y)*(D.z - C.z) + (D.y - C.y)*(B.z - C.z);  // U = BC x CD
  U.y = -(D.x - C.x)*(B.z - C.z) + (B.x - C.x)*(D.z - C.z);
  U.z = -(B.x - C.x)*(D.y - C.y) + (D.x - C.x)*(B.y - C.y);

  lengthV = V.x*V.x + V.y*V.y + V.z*V.z;
  lengthU = U.x*U.x + U.y*U.y + U.z*U.z;
  lengthprod = FastSqrt(lengthU*lengthV,tables);
   if(lengthprod < EPSILON) { fprintf(stderr,"WARNING: 3 atoms aligned in Dihedral\n "); return dih; }

  dotprod =  V.x * U.x + V.y * U.y + V.z * U.z;
  temp = (dotprod/lengthprod);
  if( temp < 1 && temp > -1 )  dih = FastAcos(temp,tables);
  else if (temp >= 1) dih = 0.0;
  else dih = 180.;

  if ((B.x - A.x)*U.x + (B.y - A.y)*U.y + (B.z - A.z)*U.z  < 0)
	  dih = -dih;

  *out = 1;

  return dih;
}

/*****************************************************************************
  Flip a whole fragment, moving from iftrom to ito around the
  axis defined from (ifrom-1) to (ito+1)
 *****************************************************************************/
int FlipFragment(struct s_polymer *p, int ifrom, int ito, double dw)
{
	int iw;
	struct vector a,b,u,du,newu;
	double normb2,coeff,normb,normT,normT2,T[3][3];

	if (ifrom<2 || ito>p->nback-3 || ifrom >= ito - 1) return 0;		// cannot apply to first and last atom

	// preceding atom
	a.x = (((p->back)+ifrom-1)->pos).x;
	a.y = (((p->back)+ifrom-1)->pos).y;
	a.z = (((p->back)+ifrom-1)->pos).z;

	// following atom
	b.x = (((p->back)+ito+1)->pos).x;
	b.y = (((p->back)+ito+1)->pos).y;
	b.z = (((p->back)+ito+1)->pos).z;

	// Translate b to the system defined by a
	b.x -= a.x;
	b.y -= a.y;
	b.z -= a.z;

	for (iw=ifrom;iw<=ito;iw++)
	{
		// atom to be moved around the axis passing by a and b
		u.x = (((p->back)+iw)->pos).x;
		u.y = (((p->back)+iw)->pos).y;
		u.z = (((p->back)+iw)->pos).z;

		u.x -= a.x;
		u.y -= a.y;
		u.z -= a.z;

		// db = b - u
		du.x = b.x - u.x;
		du.y = b.y - u.y;
		du.z = b.z - u.z;

		// coefficient for the 2nd column of base-change matrix
		normb2 = Norm2(b);
		if (normb2<EPSILON) return 0;				// A and B coincide

		coeff = 0.5 * (( Norm2(u) - Norm2(du) ) / normb2 + 1.);

		// base-change matrix (first two columns)
		normb = FastSqrt(normb2,p->tables);

		T[0][0] = b.x / normb;		T[0][1] = u.x - coeff * b.x;
		T[1][0] = b.y / normb;		T[1][1] = u.y - coeff * b.y;
		T[2][0] = b.z / normb;		T[2][1] = u.z - coeff * b.z;

		// normalize 2nd column of the matrix
		normT2 = T[0][1]*T[0][1] + T[1][1]*T[1][1] + T[2][1]*T[2][1];
		if (normT2<EPSILON) return 0;
		normT = FastSqrt(normT2,p->tables);

		T[0][1] = T[0][1] / normT;
		T[1][1] = T[1][1] / normT;
		T[2][1] = T[2][1] / normT;

		// the third colum is the vector product of the first two
		T[0][2] = T[1][0] * T[2][1] - T[2][0] * T[1][1];
		T[1][2] = T[2][0] * T[0][1] - T[0][0] * T[2][1];
		T[2][2] = T[0][0] * T[1][1] - T[1][0] * T[0][1];

		// rotate u
		u.x = coeff * normb;
		u.y = FastCos(dw,p->tables) * normT;
		u.z = FastSin(dw,p->tables) * normT;

		MatrixVectorProduct(T,u,&newu);

		// translate back in a
		newu.x += a.x;
		newu.y += a.y;
		newu.z += a.z;

		// actually moves the atom
		(((p->back)+iw)->pos).x = newu.x;
		(((p->back)+iw)->pos).y = newu.y;
		(((p->back)+iw)->pos).z = newu.z;
	}

	return 1;
}

/*****************************************************************************
 Copy the whole polymer structure, reconstructing the lookback tables;
 the numeric tables are just copied as pointer
 *****************************************************************************/
void CopyPolymer(struct s_polymer *from, struct s_polymer *to, int cfrom, int cto, int noside, int norot)
{
	int i,j,k;

	(to+cto)->nback = (from+cfrom)->nback;
	(to+cto)->etot = (from+cfrom)->etot;
	strcpy((to+cto)->title,(from+cfrom)->title);
	for (i=0;i<(from+cfrom)->nback;i++)
	{
		(((to+cto)->back)+i)->ia = (((from+cfrom)->back)+i)->ia;
		(((to+cto)->back)+i)->itype = (((from+cfrom)->back)+i)->itype;
		strcpy((((to+cto)->back)+i)->aa,(((from+cfrom)->back)+i)->aa);
		strcpy((((to+cto)->back)+i)->type,(((from+cfrom)->back)+i)->type);
		(((to+cto)->back)+i)->iaa = (((from+cfrom)->back)+i)->iaa;
		((((to+cto)->back)+i)->pos).x = ((((from+cfrom)->back)+i)->pos).x;
		((((to+cto)->back)+i)->pos).y = ((((from+cfrom)->back)+i)->pos).y;
		((((to+cto)->back)+i)->pos).z = ((((from+cfrom)->back)+i)->pos).z;
		((((to+cto)->back)+i)->sph).ang = ((((from+cfrom)->back)+i)->sph).ang;
		((((to+cto)->back)+i)->sph).dih = ((((from+cfrom)->back)+i)->sph).dih;
		((((to+cto)->back)+i)->sph).r = ((((from+cfrom)->back)+i)->sph).r;
		(((to+cto)->back)+i)->nside = (((from+cfrom)->back)+i)->nside;
		(((to+cto)->back)+i)->move = (((from+cfrom)->back)+i)->move;
		(((to+cto)->back)+i)->irot = (((from+cfrom)->back)+i)->irot;
		(((to+cto)->back)+i)->nrot = (((from+cfrom)->back)+i)->nrot;
		(((to+cto)->back)+i)->e_ang = (((from+cfrom)->back)+i)->e_ang;
		(((to+cto)->back)+i)->e_dih = (((from+cfrom)->back)+i)->e_dih;
		(((to+cto)->back)+i)->d2_next = (((from+cfrom)->back)+i)->d2_next;
		(((to+cto)->back)+i)->iapdb = (((from+cfrom)->back)+i)->iapdb;

            (((to+cto)->back)+i)->a_next = (((from+cfrom)->back)+i)->a_next;
		if (!noside)
		for (j=0;j<(((from+cfrom)->back)+i)->nside;j++)
		{
 			strcpy( (((((to+cto)->back)+i)->side)+j)->type, (((((from+cfrom)->back)+i)->side)+j)->type );
 			(((((to+cto)->back)+i)->side)+j)->itype = (((((from+cfrom)->back)+i)->side)+j)->itype;
 			(((((to+cto)->back)+i)->side)+j)->ia = (((((from+cfrom)->back)+i)->side)+j)->ia;
			((((((to+cto)->back)+i)->side)+j)->pos).x = ((((((from+cfrom)->back)+i)->side)+j)->pos).x;
			((((((to+cto)->back)+i)->side)+j)->pos).y = ((((((from+cfrom)->back)+i)->side)+j)->pos).y;
			((((((to+cto)->back)+i)->side)+j)->pos).z = ((((((from+cfrom)->back)+i)->side)+j)->pos).z;
			if (!norot)
			for (k=0;k<(((from+cfrom)->back)+i)->nrot;k++)
			{
				(((((((to+cto)->back)+i)->side)+j)->rot)+k)->b1 = (((((((from+cfrom)->back)+i)->side)+j)->rot)+k)->b1;
				(((((((to+cto)->back)+i)->side)+j)->rot)+k)->b2 = (((((((from+cfrom)->back)+i)->side)+j)->rot)+k)->b2;
				(((((((to+cto)->back)+i)->side)+j)->rot)+k)->b3 = (((((((from+cfrom)->back)+i)->side)+j)->rot)+k)->b3;
				((((((((to+cto)->back)+i)->side)+j)->rot)+k)->ang).ang = ((((((((from+cfrom)->back)+i)->side)+j)->rot)+k)->ang).ang;
				((((((((to+cto)->back)+i)->side)+j)->rot)+k)->ang).dih = ((((((((from+cfrom)->back)+i)->side)+j)->rot)+k)->ang).dih;
				((((((((to+cto)->back)+i)->side)+j)->rot)+k)->ang).r = ((((((((from+cfrom)->back)+i)->side)+j)->rot)+k)->ang).r;
			}
		}
		(((to+cto)->back)+i)->ncontacts = (((from+cfrom)->back)+i)->ncontacts;
		for (j=0;j<(((from+cfrom)->back)+i)->ncontacts;j++)
		{
			*(((((to+cto)->back)+i)->contacts)+j) = *(((((from+cfrom)->back)+i)->contacts)+j);
			*(((((to+cto)->back)+i)->contacts_p)+j) = *(((((from+cfrom)->back)+i)->contacts_p)+j);
			*(((((to+cto)->back)+i)->e)+j) = *(((((from+cfrom)->back)+i)->e)+j);
		}
		(((to+cto)->back)+i)->nshell = (((from+cfrom)->back)+i)->nshell;
		for (j=0;j<(((from+cfrom)->back)+i)->nshell;j++)
		{
				*(((((to+cto)->back)+i)->shell)+j) = *(((((from+cfrom)->back)+i)->shell)+j);
				*(((((to+cto)->back)+i)->shell_p)+j) = *(((((from+cfrom)->back)+i)->shell_p)+j);
		}
	}

	// fills v
	SetLookbackTables(to,cto);

	// the tables are the same, just copies the pointer
	(to+cto)->tables = (from+cfrom)->tables;

	#ifdef OPTIMIZEPOT
	 (to+cto)->op = from->op;
	#endif

}

void CopyAllPolymers(struct s_polymer *from, struct s_polymer *to, int n, int noside, int norot)
{
	int i;

	for (i=0;i<n;i++)
		CopyPolymer(from,to,i,i,noside,norot);
}

void DisplaceCoM(struct s_polymer *p, int ip, double dx, double dy, double dz)
{
	int i,j;

	for (i=0;i<(p+ip)->nback;i++)
	{
		((((p+ip)->back)+i)->pos).x += dx;
		((((p+ip)->back)+i)->pos).y += dy;
		((((p+ip)->back)+i)->pos).z += dz;

		for (j=0;j<(((p+ip)->back)+i)->nside;j++)
		{
			((((((p+ip)->back)+i)->side)+j)->pos).x += dx;
			((((((p+ip)->back)+i)->side)+j)->pos).y += dy;
			((((((p+ip)->back)+i)->side)+j)->pos).z += dz;
		}
	}
}


void	RotationX(struct s_polymer *p,int ip,double dtheta)
{
	/* double rotmat[3][3]= {
				{1,0,0},
				{0,cos(dtheta),-sin(dtheta)},
				{0,sin(dtheta),cos(dtheta)}

	                        };
	*/

        double cmy=0;
        double cmz=0;
        int i,j;
        double fact=1./(p+ip)->nback;
        double y,z,ys,zs;

        for(i=0;i<(p+ip)->nback;i++)
        {
                cmy+=((((p+ip)->back)+i)->pos).y;
                cmz+=((((p+ip)->back)+i)->pos).z;
        }

        cmy=cmy*fact;
        cmz=cmz*fact;

        for (i=0;i<(p+ip)->nback;i++)
        {
                y=((((p+ip)->back)+i)->pos).y-cmy;
                z=((((p+ip)->back)+i)->pos).z-cmz;
                ((((p+ip)->back)+i)->pos).y = cmy+(cos(dtheta)*y)-(sin(dtheta)*z);
                ((((p+ip)->back)+i)->pos).z = cmz+ (sin(dtheta)*y)+(cos(dtheta)*z);
    

                for (j=0;j<(((p+ip)->back)+i)->nside;j++)
                {	
                        ys=((((((p+ip)->back)+i)->side)+j)->pos).y -cmy;
                        zs=((((((p+ip)->back)+i)->side)+j)->pos).z -cmz ;
                        ((((((p+ip)->back)+i)->side)+j)->pos).y = cmy+ (cos(dtheta)*ys)-(sin(dtheta)*zs);
                        ((((((p+ip)->back)+i)->side)+j)->pos).z = cmz+ (sin(dtheta)*ys)+(cos(dtheta)*zs);
                }
        }
}

void    RotationY(struct s_polymer *p,int ip,double dtheta)
{
	/* double rotmat[3][3]= {
                                {cos(dtheta),0,sin(dtheta)},
                                {0,1,0},
                                {-sin(dtheta),0,cos(dtheta)}

                                }; 
	*/
	
	double cmx=0;
	double cmz=0;
	int i,j;
	double fact=1./(p+ip)->nback;
        double x,z,xs,zs;

	for(i=0;i<(p+ip)->nback;i++)
	{
		cmx+=((((p+ip)->back)+i)->pos).x;
		cmz+=((((p+ip)->back)+i)->pos).z;
	}

	cmx=cmx*fact;
	cmz=cmz*fact;

        for (i=0;i<(p+ip)->nback;i++)
        {
		x=((((p+ip)->back)+i)->pos).x-cmx;
		z=((((p+ip)->back)+i)->pos).z-cmz;
                ((((p+ip)->back)+i)->pos).x = cmx+ (cos(dtheta)*x) + (sin(dtheta)*z);
                ((((p+ip)->back)+i)->pos).z = cmz+ (-sin(dtheta)*x)+(cos(dtheta)*z);

                for (j=0;j<(((p+ip)->back)+i)->nside;j++)
                {
			xs=((((((p+ip)->back)+i)->side)+j)->pos).x-cmx ;
	                zs=((((((p+ip)->back)+i)->side)+j)->pos).z -cmz;
                        ((((((p+ip)->back)+i)->side)+j)->pos).x = cmx+(cos(dtheta)*xs) +(sin(dtheta)*zs);
                        ((((((p+ip)->back)+i)->side)+j)->pos).z = cmz+(-sin(dtheta)*xs)+(cos(dtheta)*zs);
                }
        }
}

void    RotationZ(struct s_polymer *p,int ip,double dtheta)
{
	/* double rotmat[3][3]= {
                                {cos(dtheta),-sin(dtheta),0},
                                {sin(dtheta),cos(dtheta),0},
                                {0,0,1}

                             }; 
	*/
	double cmx=0;
        double cmy=0;
        int i,j;
        double fact=1./(p+ip)->nback;
        double x,y,xs,ys;

        for(i=0;i<(p+ip)->nback;i++)
        {
                cmx+=((((p+ip)->back)+i)->pos).x;
                cmy+=((((p+ip)->back)+i)->pos).y;
        }

        cmx=cmx*fact;
        cmy=cmy*fact;
     
        for (i=0;i<(p+ip)->nback;i++)
        {
                x=((((p+ip)->back)+i)->pos).x-cmx;
                y=((((p+ip)->back)+i)->pos).y-cmy;
                ((((p+ip)->back)+i)->pos).x = cmx+(cos(dtheta)*x) - (sin(dtheta)*y);
                ((((p+ip)->back)+i)->pos).y = cmy+(sin(dtheta)*x)+(cos(dtheta)*y);
 
                for (j=0;j<(((p+ip)->back)+i)->nside;j++)
                {
		      	 xs=((((((p+ip)->back)+i)->side)+j)->pos).x -cmx; 
                         ys=((((((p+ip)->back)+i)->side)+j)->pos).y -cmy;
                         ((((((p+ip)->back)+i)->side)+j)->pos).x =cmx+(cos(dtheta)*xs )-( sin(dtheta)*ys); 
                         ((((((p+ip)->back)+i)->side)+j)->pos).y = cmy+(sin(dtheta)*xs)+(cos(dtheta)*ys);
                }

        }
}




void RotationClusterX(struct s_polymer *p,int ip,double dtheta,int icluster,int cluster[NCHAINMAX][NCHAINMAX],int npol_cluster)
{

        double cmy=0;
        double cmz=0;
        int i,j,k;
        double fact=0;
        double y,z,ys,zs;

        for(i=0;i<npol_cluster;i++)
                fact+=(p+cluster[icluster][i])->nback;
        fact=1./fact;



        for(j=0;j<npol_cluster;j++)
        {




	        for(i=0;i<(p+cluster[icluster][j])->nback;i++)
       		 {
                	cmy+=((((p+cluster[icluster][j])->back)+i)->pos).y;
                	cmz+=((((p+cluster[icluster][j])->back)+i)->pos).z;
        	}

	}


        cmy=cmy*fact;
        cmz=cmz*fact;



        for(k=0;k<npol_cluster;k++)
        {




	        for (i=0;i<(p+cluster[icluster][k])->nback;i++)
       		{
                	y=((((p+cluster[icluster][k])->back)+i)->pos).y-cmy;
                	z=((((p+cluster[icluster][k])->back)+i)->pos).z-cmz;
                	((((p+cluster[icluster][k])->back)+i)->pos).y = cmy+(cos(dtheta)*y)-(sin(dtheta)*z);
                	((((p+cluster[icluster][k])->back)+i)->pos).z = cmz+ (sin(dtheta)*y)+(cos(dtheta)*z);
		
	

	                for (j=0;j<(((p+cluster[icluster][k])->back)+i)->nside;j++)
        	        {
                	        ys=((((((p+cluster[icluster][k])->back)+i)->side)+j)->pos).y -cmy;
                       		zs=((((((p+cluster[icluster][k])->back)+i)->side)+j)->pos).z -cmz ;
                        	((((((p+cluster[icluster][k])->back)+i)->side)+j)->pos).y = cmy+ (cos(dtheta)*ys)-(sin(dtheta)*zs);
                        	((((((p+cluster[icluster][k])->back)+i)->side)+j)->pos).z = cmz+ (sin(dtheta)*ys)+(cos(dtheta)*zs);
                	}
        	}

	}



}

void RotationClusterY(struct s_polymer *p,int ip,double dtheta,int icluster,int cluster[NCHAINMAX][NCHAINMAX],int npol_cluster)
{

        double cmx=0;
        double cmz=0;
        int i,j,k;
	double fact=0;
        double x,z,xs,zs;

        for(i=0;i<npol_cluster;i++)
                fact+=(p+cluster[icluster][i])->nback;
        fact=1./fact;

        for(j=0;j<npol_cluster;j++)
        {


	        for(i=0;i<(p+cluster[icluster][j])->nback;i++)
        	{
                	cmx+=((((p+cluster[icluster][j])->back)+i)->pos).x;
                	cmz+=((((p+cluster[icluster][j])->back)+i)->pos).z;
       		}
	}

        cmx=cmx*fact;
        cmz=cmz*fact;


        for(k=0;k<npol_cluster;k++)
        {


	        for (i=0;i<(p+cluster[icluster][k])->nback;i++)
        	{
                	x=((((p+cluster[icluster][k])->back)+i)->pos).x-cmx;
                	z=((((p+cluster[icluster][k])->back)+i)->pos).z-cmz;
                	((((p+cluster[icluster][k])->back)+i)->pos).x = cmx+ (cos(dtheta)*x) + (sin(dtheta)*z);
                	((((p+cluster[icluster][k])->back)+i)->pos).z = cmz+ (-sin(dtheta)*x)+(cos(dtheta)*z);

                	for (j=0;j<(((p+cluster[icluster][k])->back)+i)->nside;j++)
                	{
                        	xs=((((((p+cluster[icluster][k])->back)+i)->side)+j)->pos).x-cmx ;
                        	zs=((((((p+cluster[icluster][k])->back)+i)->side)+j)->pos).z -cmz;
                        	((((((p+cluster[icluster][k])->back)+i)->side)+j)->pos).x = cmx+(cos(dtheta)*xs) +(sin(dtheta)*zs);
                        	((((((p+cluster[icluster][k])->back)+i)->side)+j)->pos).z = cmz+(-sin(dtheta)*xs)+(cos(dtheta)*zs);
                	}
        	}	
	}

}
void RotationClusterZ(struct s_polymer *p,int ip,double dtheta,int icluster,int cluster[NCHAINMAX][NCHAINMAX],int npol_cluster)
{

       /* double rotmat[3][3]= {
	       {cos(dtheta),-sin(dtheta),0},
               {sin(dtheta),cos(dtheta),0},
               {0,0,1}
               }; 
       */

        double cmx=0;
        double cmy=0;
	double fact=0;
        int i,j,k;
	for(i=0;i<npol_cluster;i++)
		fact+=(p+cluster[icluster][i])->nback;
        fact=1./fact;
        double x,y,xs,ys;

	for(j=0;j<npol_cluster;j++)
	{
        	for(i=0;i<(p+cluster[icluster][j])->nback;i++)
        	{
                	cmx+=((((p+cluster[icluster][j])->back)+i)->pos).x;
                	cmy+=((((p+cluster[icluster][j])->back)+i)->pos).y;
        	}
	}

        cmx=cmx*fact;
        cmy=cmy*fact;

	for(k=0;k<npol_cluster;k++)
	{
        	for (i=0;i<(p+cluster[icluster][k])->nback;i++)
        	{	
                	x=((((p+cluster[icluster][k])->back)+i)->pos).x-cmx;
                	y=((((p+cluster[icluster][k])->back)+i)->pos).y-cmy;
                	((((p+cluster[icluster][k])->back)+i)->pos).x = cmx+(cos(dtheta)*x) - (sin(dtheta)*y);
                	((((p+cluster[icluster][k])->back)+i)->pos).y = cmy+(sin(dtheta)*x)+(cos(dtheta)*y);

                	for (j=0;j<(((p+cluster[icluster][k])->back)+i)->nside;j++)
                	{
                         	xs=((((((p+cluster[icluster][k])->back)+i)->side)+j)->pos).x -cmx;
                         	ys=((((((p+cluster[icluster][k])->back)+i)->side)+j)->pos).y -cmy;
                         	((((((p+cluster[icluster][k])->back)+i)->side)+j)->pos).x =cmx+(cos(dtheta)*xs )-( sin(dtheta)*ys);
                         	((((((p+cluster[icluster][k])->back)+i)->side)+j)->pos).y = cmy+(sin(dtheta)*xs)+(cos(dtheta)*ys);
                	}

        	}

	}


}

