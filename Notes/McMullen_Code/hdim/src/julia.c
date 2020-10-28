
/* Julia I-O and Dynamics */
/* Map is z^2 + param */

#include "hdim.h"

/* Initialize Julia set Markov partition */
/* with points at angles 1/4 and 3/4     */
julia_cover()
{
	int g, t[2];
	complex z;
	double r;

	t[0] = 0;
	t[1] = 1; 
	g = MAP;    

	ncover=npool=0;
	z = pre_fp();
	z = map_inv(z,z);
	r = 2*cxabs(z);
	add_tile(z,r,g,t,2);   
	z.x = -z.x;
	z.y = -z.y;
	add_tile(z,r,g,t,2);  
}

/* Inverse image of w near z */
complex map_inv(w,z)
	complex w, z;
{
	complex u;
	double dummy;
	
	(*map)(MAP,w,0.0,z,&u,&dummy);
	return(u);
}

/* Find pre-fixed point */
complex pre_fp()
{
	static complex minus_one = {-1.0,0.0};
	complex w, z;

	switch(mode) {
	case BP:
	case BPS:
		return(minus_one);
	case QUAD:
	case QUADS:
		w = quad_fp();
		z = map_inv(w,w);
	if(cxabs(sub(z,w)) < FUZZ)
	{	z.x = -z.x;
		z.y = -z.y;
	}
	return(z);
	}
}

/* Quadratic Map Routines */
/* f[z_]:=z^2 + p */

/* finv[z_]:=Sqrt[z-p] */
void quad_map(g,c,r,z,e,t)
	int g;
	complex c,*e,z;
	double r,*t;
{
	*e = contsqrt(sub(c,param),z);
	*t = r*quad_dmap(MAP,c);
}

/* finv'[z] */
double quad_dmap(g,c)
	int g;
	complex c;
{
	double d;

	d=cxabs(sub(c,param));
	d = 0.5*pow(d,-0.5);
	return(d);
}

/* Find expanding fixed-point */
/* 1/2 +/- Sqrt[1/4-c] */
complex quad_fp()
{
	complex w, z;

	w.x = 0.25 - param.x;
	w.y = param.y;
	w = cxsqrt(w);
	z.x = 0.5 + w.x; 
	z.y = w.y;
	if(cxabs(z) >= 0.5) return(z);
	z.x = 0.5 - w.x; 
	z.y = - w.y;
	return(z);
}

/* Blaschke Dynamics Routines */

/* f[z_]:=((p+2)z^2+(2-p))/((p+2)+(2-p)z^2) */
/* f[1]=1, f'[1]=p */
void bp_map(g,c,r,z,e,t)
	int g;
	complex c,*e,z;
	double r,*t;
{
	static complex two = {2.0,0.0};
	complex frac, num, den;

	num= add(sub(param,two),
		mult(c,add(param,two)));
	den= add(add(param,two),
		mult(c,sub(param,two)));
	frac = cdiv(num,den);
	*e = contsqrt(frac,z);
	*t = r*bp_dmap(MAP,c);
}

/* |finv'[z]| = 4p*|(p+2)z+(p-2)|^(-1/2) *
                   |(p-2)z+(p+2)|^(-3/2)
*/
double bp_dmap(g,c)
	int g;
	complex c;
{
	double d;
	static complex two = {2.0,0.0};

	d = 4*cxabs(param);
	d *= pow(cxabs(add(
		mult(c,add(two,param)),
		sub(param,two)
		)),-0.5);
	d *= pow(cxabs(add(
		mult(c,sub(param,two)),
		add(param,two)
		)),-1.5);
	return(d);
}
