
/* Markov partition routines */

#include "hdim.h"

/* Iterated refinement */
finer()
{
	while(refine()) {};
}

/* Refine all tiles larger than eps */
int refine()
{
	int i, nold, s;

	nold=ncover;
	for(s=0; s<nold; s++)
	if(cover[s].r > eps &&
	cover[s].g != ID)
	{	refine_tile(s);
		i=1;
	}
	compose_tiles();
	clean_tiles();
	return(ncover > nold);
}

/* Refine a tile */
refine_tile(s)
	int s;
{
	complex c, cnew;
	double rnew;
	int g, i, t, t0, t1;

	g=cover[s].g;
	c=cover[s].c;
	t0=cover[s].t0;
	t1=cover[s].t1;

	for(i=t0; i<t1; i++)
	{	t=pool[i];
		(*map)(g,cover[t].c,cover[t].r,c,
		&cnew,&rnew);
		pool[i]=ncover;
		add_tile(cnew,rnew,g,&t,1);
	}
	cover[s].g=ID;
}

/* Short-circuit tiles mapping by the identity */
compose_tiles()
{
	int i, j, n, s, t, *u;

	u = int_alloc(ncover);
	if(u==NULL) space_err();
	for(s=0; s<ncover; s++)
	{	n=0;
	for(i=cover[s].t0; i<cover[s].t1; i++)
	{	t=pool[i];
		if(cover[t].g==ID)
		for(j=cover[t].t0; j<cover[t].t1; j++)
			{u[n]=pool[j]; n++;}
		else {u[n] = t; n++;}
	}

	pool_space(n+npool);	/* Request space */
	cover[s].t0 = npool;
	cover[s].t1 = npool+n;
	for(i=0; i<n; i++) pool[npool+i]=u[i];
	npool = npool+n;
	}
	free(u);
}

/* Remove tiles that map by ID */
clean_tiles()
{
	int i, n, *new, s;

	new = int_alloc(ncover);
	n=0;
	for(s=0; s<ncover; s++) 
	if(cover[s].g != ID)
	{	cover[n]=cover[s];
		new[s]=n;
		n++;
	}
	ncover=n;
	for(i=0; i<npool; i++)
		pool[i]=new[pool[i]];
	free(new);
	clean_pool();
}

clean_pool()
{
	int i, n, *new, nold, s, r0, r1;

	new = int_alloc(npool);
	n=0;
	for(s=0; s<ncover; s++)
	{	nold=n;
		for(i=cover[s].t0; i<cover[s].t1; i++)
		{	new[n]=pool[i];
			n++;
		}
		cover[s].t0=nold;
		cover[s].t1=n;
	}
	for(i=0; i<n; i++) pool[i]=new[i];
	npool=n;
	free(new);
}

/* Create new tile, add to cover */
add_tile(c,r,g,im,n)
	complex c;
	double r;
	int g, im[], n;
{
	int i;	

	cover_space(ncover);	/* Request space */
	pool_space(npool+n);

	cover[ncover].c = c;
	cover[ncover].r = r;
	cover[ncover].g = g;
	cover[ncover].t0 = npool;
	cover[ncover].t1 = npool+n;
	ncover++;
	for(i=0; i<n; i++) pool[npool+i] = im[i];
	npool += n;
}

/* Make matrix (m) from tiling (cover,ncover) */

markov()
{
	int g, i, n, s, t;
	complex c;
	double r;

	n = 0;
	for(s=0; s<ncover; s++)
	{	g=cover[s].g;
		for(i=cover[s].t0; i<cover[s].t1; i++)
		{	t = pool[i];
			c = cover[t].c;
			r = cover[t].r;
			m.i[n]=s;
			m.j[n]=t;
			m.v[n]=(*dmap)(g,c);
			n++;
		}
	}
}
