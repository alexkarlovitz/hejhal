
/* Kleinian I-O and Dynamics */

#include "hdim.h"
#define SQRT3 1.732050807568877
#define PI    3.141592653589793

/* I-O Routines */
/* Read reflection generators from stdin */
read_gens()
{
	ngen=0;
	while(EOF != scanf("%lf %lf %lf",
		&genc[ngen].x,
		&genc[ngen].y,
		&genr[ngen])) 
	ngen++;
	if(ngen > MAXGEN || ngen == 0) klein_err();
}

/* n circle Schottky group, each circle making
   a visual angle of param.x degrees from 0.
   Eg:  (3,120) is 3 tangent circles.
*/
schottky_gens()
{
	int i;
	double ang, r, rotate, s;

	if(ngen > MAXGEN || ngen <= 0)
		klein_err();

	ang=M_PI*param.x/360.0;	
	r  =tan(ang);
	s  =1/cos(ang);
	rotate=2*M_PI/ngen;

	for(i=0; i<ngen; i++)
	{	genc[i]=polar(s,i*rotate);
		genr[i] = r;
	}
}

/* Hecke group.  Conjugate to reflection
   in |z|=r and lines Re z = +/- 1.*/
hecke_gens()
{
	complex a;
	double r2, t;

	ngen=3;
	r2 = param.x * param.x;
	a.x=(3-r2)/(3+r2);
	a.y=(2*SQRT3*param.x)/(3+r2);
	genc[0] = polar(2.0,2*PI/3);
	genr[0] = SQRT3;
	genc[1] = polar(2.0,4*PI/3);
	genr[1] = SQRT3;
	genc[2].x = 1.0/a.x;
	genc[2].y = 0.0;
	genr[2] = a.y/a.x;
}

/* Apollonian packing */
apollo_gens()
{
	static complex c[4] = {
	0.0,0.0,
	7.4641016151377545871, 0.0,
	-3.732050807568877294, 6.464101615137754587,
	-3.732050807568877294,-6.464101615137754587};
	static double r[4] = {
	1.0,
	6.4641016151377545871,
	6.4641016151377545871,
	6.4641016151377545871};

	int i;

	ngen=4;
	for(i=0; i<4; i++)
	{	genc[i]=c[i];
		genr[i]=r[i];
	}
}

/* Generate initial cover from generators */
klein_cover()
{
	int i, j, n, t[MAXGEN];

	ncover=npool=0;
	for(i=0; i<ngen; i++)
	{	n=0;
		for(j=0; j<ngen; j++) if(j != i)
		{	t[n]=j;
			n++;
		}
		add_tile(genc[i],genr[i],i,t,n);
	}
}

/* Dynamics Routines */
/* Distance squared */
double cabs2(z,w)
	complex z,w;
{
	return((z.x-w.x)*(z.x-w.x)+
		(z.y-w.y)*(z.y-w.y));
}

/* Invert (c,r) through (d,s) to yield (e,t) */
invert(c,r,d,s,e,t)
	complex c, d, *e;
	double r, s, *t;
{
	double a, b;

	a=cabs2(c,d);
	b=s*s/(a-r*r);
	*t=b*r;
	e->x = d.x + b*(c.x-d.x);
	e->y = d.y + b*(c.y-d.y);
}

/* Apply generator g to circle (c,r)
   yielding (e,t) near point z 
*/
void klein_map(g,c,r,z,e,t)
	int g;
	complex c,*e,z;
	double r,*t;
{
	invert(c,r,genc[g],genr[g],e,t);
}

/* Derivative of dynamics in a ball 
   Bound=0 for center, 1 for upper, -1 for lower */
double klein_dmap(g,c)
	int g;
	complex c;
{
	complex d;
	double deriv, s, t;

	d=genc[g];
	s=genr[g];
	t=cxabs(sub(c,d));
	deriv=(s*s)/(t*t);
	return(deriv);
}

/* Errors */
klein_err()
{	printf("Too many or too few generators\n");
	exit(1);
}
