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
 * misc.c
 *
 *  Created on: Sep 20, 2010
 *      Author: guido
 */

#include "montegrappa.h"
#include <time.h>

/*****************************************************************************
 Dumb generator of integer random numbers
 *****************************************************************************/
int irand(int r)
{
  int i;

  i=(int)(rand() / (((double)RAND_MAX + 1)/ r)); 
  return i;
}

/*****************************************************************************
 Dumb generator of double random numbers
 *****************************************************************************/
double frand()
{
  double q;

  q=(double) rand()/(RAND_MAX+1.0);
  return q;
}

/*****************************************************************************
 Initialize random number with seed (-1 to use computer time)
 *****************************************************************************/
long Randomize(int n)
{
  int i;

  if (n==-1)
  {
	  n = (time(0) - 760650080) % 3600; 
	  fprintf(stderr,"Random seed = %d\n",n);
  }
 
  for (i=0;i<n;i++) irand(i); 

  return n;
}

/*****************************************************************************
 Fast calculation of sqrt,cos,sin,exp
	 the index should be the closest integer to x, not the integral part,
	 otherwise there is a bias towards smaller values (thus we use HBIN)
 *****************************************************************************/
double FastSqrt(double x, struct s_tables *t)
{
	#ifdef TABLED_F
	int i;

	if (x<0) return sqrt(-1);
	else if (x==0) return 0;
	else if (x<FSQRT_HBIN) return sqrt(x);
	else if (x>=FSQRT_MAX-1) return sqrt(x);

	i = (int) ( (x + FSQRT_HBIN) * FSQRT_BINRC);
	return (t->fast_sqrt)[i];

	#else
	return sqrt(x);
	#endif
}

double FastSin(double x, struct s_tables *t)
{
	#ifdef TABLED_F
	int i;

	// make x between 0 and 360
	while (x<0) { x += 360.; }
	while (x>=360.-FTRIG_HBIN) { x -= 360.; }
	if (x>-FTRIG_HBIN && x<FTRIG_HBIN) return x*PI/180.;

	i = (int) ( (x + FTRIG_HBIN) * FTRIG_BINRC);
	return (t->fast_sin)[i];
	#else
	return sin(x*PI/180.);
	#endif
}

double FastCos(double x, struct s_tables *t)
{
	#ifdef TABLED_F
	int i;

	// make x between 0 and 360
	while (x<0) { x += 360.; }
	while (x>=360.-FTRIG_HBIN) { x -= 360.; }

	if (x>-FTRIG_HBIN && x<FTRIG_HBIN) return 1-0.5*x*x*PI*PI/32400;

	i = (int) ( (x + FTRIG_HBIN) * FTRIG_BINRC);
	return (t->fast_cos)[i];
	#else
	return cos(x*PI/180.);
	#endif
}

double FastAcos(double x, struct s_tables *t)
{
	#ifdef TABLED_F

	int i;

	#ifdef DEBUG
	 if (x<-1 || x>1) { fprintf(stderr,"acos(%lf)\n",x); Error("Argument of FastAcos should be in [-1,1]"); }
    #endif

	if (x>=1.-FACOS_HBIN) return 0;

	i = (int) ( (x + 1.) * FACOS_BINRC);
	return (t->fast_acos)[i];
	#else
	return acos(x)*180./PI;
	#endif
}

double FastExp(double x, struct s_tables *t)
{
	#ifdef TABLED_F
	int i;

	if (x>FEXP_MAX-FEXP_HBIN) return exp(x);
	else if (x<-FEXP_MAX+FEXP_HBIN) return exp(x);
	else if (x>0)
	{
		i = (int) ( (x + FEXP_HBIN) * FEXP_BINRC);
		return (t->fast_expp)[i];
	}
	else
	{
		i = -(int) ( (x + FEXP_HBIN) * FEXP_BINRC);
		return (t->fast_expm)[i];
	}
	#else
	return exp(x);
	#endif
}

/*****************************************************************************
 SQUARE norm of a vector
 *****************************************************************************/
double Norm2(struct vector a)
{
	return a.x*a.x + a.y*a.y + a.z*a.z;
}

/*****************************************************************************
 product of matrix by vector
 *****************************************************************************/
void MatrixVectorProduct(double T[][3], struct vector in, struct vector *out)
{
	out->x = T[0][0] * in.x + T[0][1] * in.y + T[0][2] * in.z;
	out->y = T[1][0] * in.x + T[1][1] * in.y + T[1][2] * in.z;
	out->z = T[2][0] * in.x + T[2][1] * in.y + T[2][2] * in.z;
}

int Abs(int x)
{
	if (x>=0) return x;
	return -x;
}

double DAbs(double x)
{
	if (x>=0) return x;
	return -x;
}

/*****************************************************************************
 Gaussian random number (from numerical recipes)
 *****************************************************************************/
double gasdev(long *idum)
{

        static int iset = 0;
        static double gset;
        double fac, rsq, v1, v2;

        if (*idum < 0)
                iset = 0; /* Reinitialize. */

        if (iset == 0)
                                                        /* We don't have an extra deviate handy, so pick two uniform
                                                         * numbers in the square extending from -1 to +1 in each direction,
                                                         * see if they are in the unit circle, and if they are not, try
                                                         * again. */
        {
          do {
                v1 = 2.0*ran1(idum)- 1.0;
                v2 = 2.0*ran1(idum)- 1.0;
                rsq = v1*v1 + v2*v2;
             } while (rsq >= 1.0 || rsq == 0.0);

          fac = sqrt(-2.0*log(rsq)/rsq);
                                                        /* Now make the Box-Muller transformation to get two normal
                                                         * deviates. Return one and save the other for next time. */
          gset = v1*fac;
          iset = 1; /* Set flag. */

          return v2*fac;
        }

        else     /* We have an extra deviate handy, */
        {
          iset=0; /* so unset the flag, */
          return gset; /* and return it. */
        }
}

float ran1(long *idum)                     // initialize calling with idum<0
{
    int j;
    long k;
    static long iy=0;
    static long iv[NTAB];
    float temp;

      if (*idum <= 0 || !iy)            // initialize
      {
        if (-(*idum) < 1) *idum=1;      // be sure to prevent idum=0
        else *idum = -(*idum);

        for (j=NTAB+7;j>=0;j--)         // Load the shuffle table
        {
            k=(*idum)/IQ;               // after 8 warm-ups
            *idum=IA*(*idum-k*IQ)-IR*k;
            if (*idum < 0) *idum += IM;
            if (j < NTAB) iv[j] = *idum;
        }
        iy=iv[0];
      }

      k=(*idum)/IQ;                         // start here when not initializing
      *idum=IA*(*idum-k*IQ)-IR*k;
      if (*idum < 0) *idum += IM;
      j=iy/NDIV;
      iy=iv[j];
      iv[j] = *idum;
      if ((temp=AM*iy) > RNMX) return RNMX;
      else return temp;
}

double gauss(double av, double sig)
{
  double r1,r2;

  r1 = -log(1-frand());
  r2 =  2*PI*frand();
  r1 =  sqrt(2*r1);
  return sig*r1*cos(r2)+av;
}



/*****************************************************************************
 Copy a double matrix 
 *****************************************************************************/
void CopyDoubleMatrix(double **from, double **to, int n, int m)
{
	int i,j;

	for (i=0;i<n;i++)
		for (j=0;j<m;j++)
			to[i][j] = from[i][j];
}

void CopyVector(struct vector *from, struct vector *to)
{
	to->x = from->x;
	to->y = from->y;
	to->z = from->z;
}
