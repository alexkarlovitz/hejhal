#include <math.h>
#include "cx.h"
#define TWO30  1073741824

/* Complex arithmetic and utilities   */

complex add (z, w)
	complex z, w;
{
	complex t;
	t.x = z.x + w.x;
	t.y = z.y + w.y;
	return(t);
}

complex sub (z, w)
	complex z, w;
{
	complex t;
	t.x = z.x - w.x;
	t.y = z.y - w.y;
	return(t);
}

complex mult (z, w)
	complex z, w;
{
	complex t;
	t.x = z.x*w.x - z.y*w.y;
	t.y = z.x*w.y + z.y*w.x;
	return(t);
}

complex cxsqrt(z)
	complex z;
{
	complex w;
	double fabs(), sqrt();
/* Worry about numerical stability */
	if (z.x == 0.0 && z.y == 0.0) return(z); 
	else
	if (z.x > fabs(z.y))
		{
		w.x = sqrt((z.x+sqrt(z.x*z.x+z.y*z.y))/2);
		w.y = z.y/(2*w.x);
		}
	else
		{
		w.y = sqrt((-z.x+sqrt(z.x*z.x+z.y*z.y))/2);
		w.x = z.y/(2*w.y);
		}
	return(w);
}

/* Compute sqrt(z) in the half-plane perpendicular to w. */
complex contsqrt(z,w)
complex z,w;
{
	complex cxsqrt(), t;

	t = cxsqrt(z);
	if (0 > (t.x*w.x + t.y*w.y)) 
		{t.x = -t.x; t.y = -t.y;}
	return(t);
}

double cxabs (z)       /* L 2 norm of z */
	complex z;
{
	return (sqrt (z.x*z.x + z.y*z.y));
}

complex polar (radius, angle) /*Convert to complex. */
	double radius, angle;
{
	complex z;
	z.x = cos (angle) * radius;
	z.y = sin (angle) * radius;
	return(z);
}

double arg(z)
	complex z;
{
	return(atan2(z.y,z.x));
}

complex cdiv (z, w)
        complex z, w;
{
        complex mult(), recip();
        return(mult(z,recip(w)));
}

complex recip (z)
        complex z;
{
        complex w;
        double r;
        r = z.x*z.x + z.y*z.y;
        w.x = z.x / r;
        w.y = -z.y / r;
        return(w);
}


