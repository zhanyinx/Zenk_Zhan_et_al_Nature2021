/*
 * mhistogram.h
 *
 *  Created on: Nov 6, 2009
 *      Author: guido
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "do_mhistogram.h"

#ifndef NTEMPMAX
#define NTEMPMAX 1000
#endif

#define NBINMAX  500
#define NDATAMAX 1000000

#ifndef EPSILON
#define EPSILON 0.000000001
#endif

struct param_s
{
    char title[80];
    int ntemp;
    double ebin;
    double emin;
    int femin;						// 1=force new emin
    int test;						// 1=do a test of the program
    char nfile[NTEMPMAX][200];
    int ncol1[NTEMPMAX];
    int ncol2[NTEMPMAX];
    double temp[NTEMPMAX];
    int debug;						// 0=quiet, 1=general output, 2=detailed output, 3=every piss
    double thresh;					// threshold to consider meaningful a bin of the histogram
    int nohisto;					// 0=read histograms, 1=read row data
    double tminout;					// minimum T to calculate and output thermodynamics
    double tbinout;					// ... its binning
    int ntout;						// ... how many bins
    char nfout[100];				// file name for output thermodynamics
    char nfene[100];				// file name to print dumb energies
    char nfg[100];					// file name to print the density of states
    char nfz[100];					// file name to print the log(z)
    double kb;						// boltzmann constant
    int paranoid;
    int ignoreb;					// ignore temperatures which prevent equation solving
    double deltat;					// ignore temperatures which are closer than deltat
};

struct histo_s
{
    int nbin;
    double emin;
    double emax;
    double ebin;
    double *histo;
};

void Help(FILE *fp);
void Parser(char *fname, struct param_s *p);

// io.c
void Error(char *text);
void ReadHistogram(struct histo_s *h, int itemp, struct param_s *p);
int ReadData(struct param_s *p, double *d, int itemp);
void MakeHistogram(struct histo_s *h, double *d, int itemp, int ndata, struct param_s *p);
void PrintAverageEnergies(FILE *fout, struct histo_s *h, double *t, int ntemp);
void Test(struct param_s *p);

// memory.c
int **AlloIntMat(int l, int m);
double **AlloDoubleMat(int l, int m);
struct histo_s *AlloHisto(int nhisto, int nbin);
double *AlloDouble(int l);
int *AlloInt(int l);

// misc.c
struct histo_s *RebinHistogram(struct histo_s *old, struct param_s *p);
void OrderTemperatures(struct histo_s *h, struct param_s *p, int *n);
void Threshold(struct histo_s *h, struct param_s *p);
int irand(int r);
double frand(void);



