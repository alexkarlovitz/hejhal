
/* Memory management routines */

#include "hdim.h"

int *int_alloc(size)
	int size;
{
	int *p;

	p = (int *) calloc(size,sizeof(int));
	if(p == NULL) space_err();
	return(p);
}

double *double_alloc(size)
	int size;
{
	double *p;

	p = (double *) calloc(size,sizeof(double));
	if(p == NULL) space_err();
	return(p);
}

tile *tile_alloc(size)
	int size;
{
	tile *p;

	p = (tile *) calloc(size,sizeof(tile));
	if(p == NULL) space_err();
	return(p);
}

/* Initialize variables and ponters */
space_init()
{
	cover  = tile_alloc(INIT);
	pool   = int_alloc(INIT);

	acover = INIT;
	apool  = INIT;
	alow   = ALOW;
	ahigh  = AHIGH;
	astep  = ASTEP;
        debug  = DEBUG;
        emin   = EPS;
        emax   = EPS;
        eps    = EPS;
	maxit  = MAXIT;
        mode   = MODE;
	ngen   = NGEN;
        ps     = PS;
        verbose= VERBOSE;
}

/* Contract to space currently used */
pool_reduce()
{
	apool = npool;
	pool = (int *) realloc(pool,apool*sizeof(int));
	if(pool == NULL) space_err();
}

/* Prepare to store into position n */ 
pool_space(n)
	int n;
{
	if(n < apool) return(0);
	apool = GROW * apool;
	pool = (int *) realloc(pool,apool*sizeof(int));
	if(pool == NULL) space_err();
}

/* Prepare to store into position n */ 
cover_space(n)
	int n;
{
	if(n < acover) return(1);
	acover = GROW * acover;
	cover = (tile *) realloc(cover,acover*sizeof(tile));
	if(cover == NULL) space_err();
}

/* Create sparse matrix space */
m_space()
{
	pool_reduce();
	m.nrows = m.ncols = ncover;
	m.nentries = npool;
	m.i = int_alloc(npool);
	m.j = int_alloc(npool);
	m.v = double_alloc(npool);
	a   = double_alloc(ncover);
}

/* Free sparse matrix space */
m_free()
{
	free(m.i);
	free(m.j);
	free(m.v);
	free(a);
}

space_err()
{
	printf("hdim:  Memory capacity exceeded\n");
	exit(1);
}
