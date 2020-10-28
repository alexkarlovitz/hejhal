
/* Eigenvalue and pressure routines */

#include "hdim.h" 
#define TARGET 1.0	/* Target eigenvalue */

/* Returns the value of delta s.t. m(i,j)^delta 
   has the target eigenvalue 1.0.

   Frees sparse matrix space before returning 
*/
int eigens;

double critical()
{
	int i, n;
	double d, e;
	static double delta = 1.0;

	n = m.nrows;
	for (i=0; i<n; i++) a[i] = 1.0/n;
	e=pressure(2.0);
	if(e == 0.0) return(0.0); 
	if(e > 1.0)  return(2.0);

	if(debug) report_header();
	for(newtons=0; newtons<maxit; newtons++) 
	{	e = pressure(delta);
		if(debug) report_ev(e,delta);
		if(fabs(e-TARGET) < FUZZ)
			return(delta); 
		d = ddelta(delta);
		delta = delta + (TARGET-e)/d;
		if(delta < 0.0) delta=0.0;
		if(delta > 2.0) delta=2.0;
	}
	printf("Iteration limit of %d reached\n", maxit);
	return(delta);
}

/* Compute eigenvalue of m^delta */ 
double pressure(delta)
	double delta;
{
	int i;
	double e, *q; 

	q = double_alloc(m.nentries);
	for(i=0; i<m.nentries; i++)
	{	q[i]=m.v[i];
		m.v[i]=pow(m.v[i],delta);
	}
	e= eigenvalue();
	for(i=0; i<m.nentries; i++)
		m.v[i]=q[i];
	free(q);
	return(e);
}

/* Compute (d/delta) of m^delta's eigenvalue, 
  where a is a normalized eigenvector. 
*/
double ddelta(delta)
	double delta;
{
	double d, v;
	int i;

	d = 0.0;
	for (i=0; i<m.nentries; i++)
	{	v = m.v[i];
		if (v != 0) 
		d += log(v) * pow(v,delta) * a[m.j[i]];
	}
	return(d);
}

/* Perron-Frobenius eigenvalue and vector of m */

double eigenvalue()
{
	int i, n;
	double e, enew;

	n = m.nrows;
	enew=0;
	eigens=0;
	do {	e=enew;
		eigens++;
		multiply();
		enew = sum(a,n);
		if(enew == 0.0) return(0.0);
		for(i=0; i<n; i++) 
			a[i] *= 1/enew;
	} while (fabs(enew-e) > FUZZ && 
		fabs(e-enew) > FAC*fabs(e-TARGET)) ;
	return(enew);
}

/* Sum of entries in vector a */
/* Split vector in half to reduce error */
double sum(a,n)
	double a[];
	int n;
{
	int i;
	double s;

	if(n < 4)
	{	s = 0.0;
		for (i=0; i<n; i++) s += a[i];
		return(s);
	}
	else
	{	i = n/2;
		return(sum(a,i)+sum(&a[i],n-i));
	}
}

/*  Act on vector a by a square matrix m */
multiply()
{
	int i, n;
	double *b;

	n = m.nrows;
	b = double_alloc(n);
	for (i=0; i<n; i++) b[i] = 0.0;
	for (i=0; i<m.nentries; i++)
		b[m.i[i]] +=
		m.v[i] * a[m.j[i]];
	for (i=0; i<n; i++) a[i] = b[i];
	free(b);
}

report_ev(e,delta)
	double e, delta;
{
	printf("%3d. %.12lf  %.12lf %8d\n",
		newtons,e,delta,eigens);
}

report_header()
{
	printf("\nMarkov partition: %d;",ncover);
	printf("   Transition matrix: %d\n",npool);
	printf("\nStep   Eigenvalue        Delta           M^n\n");
}
