
/* Header file for hdim */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "cx.h"

/* Arrays sizes and limits */
#define INIT 1000    /* Initial space */
#define GROW 1.2     /* Space growth factor */
#define MAXGEN 100
#define MAXIT  400   /* Newton iterations */

/* Constants */
#define FUZZ 1e-11
#define FAC 1e-4

/* Default values */
#define DEBUG 0		/* No debug */
#define PS 0		/* No postscript */
#define EPS 1e-3	/* Markov partition epsilon */
#define MODE SCHOTTKYS 	/* Default mode */
#define ALOW 1		/* Schottky groups */
#define AHIGH 10
#define ASTEP 1
#define NGEN  3
#define VERBOSE 1	/* Reporting on */

/* Input modes */
#define READ      0
#define APOLLO    1
#define SCHOTTKY  2
#define QUAD      3
#define BP        4
#define HECKE     5
#define SERIES    10   /* >= 10 are series */
#define SCHOTTKYS 12
#define QUADS     13
#define BPS       14
#define HECKES    15

/* conceptual constants */
#define ID -1     /* values of cover.g */
#define MAP 0   

/* Structures */
typedef struct 
	{	int nrows,ncols,nentries,*i, *j; 
		double *v;}
	sparse_matrix;

typedef struct
	{	complex c;
		double r;
		int g, t0, t1;}
	tile;

/* External single variables */
extern int acover, apool, debug, maxit, mode, 
	newtons, ncover, ngen, npool, ps, verbose;
extern double alow, ahigh, astep, emin, emax, eps;
extern complex param;
extern FILE *file;

/* External arrays */
extern tile *cover;
extern complex genc[];
extern double *a, genr[];
extern int *pool;
extern sparse_matrix m;

/* Function types */
double  bp_dmap(), klein_dmap(), quad_dmap(),
	(*dmap)();
double cabs2(), critical(), ddelta(), eigenvalue(), 
	pressure(), sum();
complex map_inv(), pre_fp(), quad_fp();
void  bp_map(), klein_map(), quad_map(),(*map)();
int *int_alloc();
double *double_alloc();
