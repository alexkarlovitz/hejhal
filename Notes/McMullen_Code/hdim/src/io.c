
/* Routines to compute and print dimensions */

#include "hdim.h"

/* Main Routine */

dims()
{
	if(mode >= SERIES) 
	{	param.y = 0;
		out_header();
		for(param.x = alow;
		    param.x <= ahigh+astep/2;
		    param.x += astep)
		{gens(); out_dim();}
	} else 
	{	gens();
		out_header();
		for(eps =  emax; 
		    eps >= emin;
		    eps *= 0.5)
		out_dim();
		
	}
}

double dim()
{
	double d;

	init_cover();	/* Initial cover */
	finer();	/* Refine cover */
	ps_report();	/* Report cover */
	m_space();	/* Create space for matrix */
	markov();	/* Transition matrix */
	d=critical();   /* Find critical exponent */
	m_free();	/* Release space */
	return(d);
}

init_cover()
{
	switch(mode) {
	case READ:
	case APOLLO:
	case HECKE:
	case HECKES:
	case SCHOTTKY:
	case SCHOTTKYS:
		klein_cover();
		break;
	case QUAD:
	case QUADS:
	case BP:
	case BPS:
		julia_cover();
		break;
	}
}

/* Set up maps and generators */
gens()
{
	switch(mode) {
	case HECKE:
	case HECKES:
		hecke_gens(); 
		map  = klein_map;
		dmap = klein_dmap;
		break;
	case SCHOTTKY:
	case SCHOTTKYS:
		schottky_gens(); 
		map  = klein_map;
		dmap = klein_dmap;
		break;
	case QUADS:
	case QUAD:
		map  = quad_map;
		dmap = quad_dmap;
		break;
	case APOLLO:
		apollo_gens(); 
		map  = klein_map;
		dmap = klein_dmap;
		break;
	case READ:
		read_gens(); 
		map  = klein_map;
		dmap = klein_dmap;
		break;
	case BPS:
	case BP:
		map = bp_map;
		dmap = bp_dmap;
		break;
	}
}

/* Print head of table */
out_header()
{
	out_title();
	switch(mode) {
	case SCHOTTKYS: 
		printf(" Angle "); break;
	case HECKES: 
		printf(" Radius "); break;
	case QUADS:     
		printf("    C   "); break;
	case BPS:     
		printf(" Lambda"); break;
	default:     
		printf(" Epsilon");
	}

	printf("   Dimension ");
	if(verbose) printf(
	"  Cover  Matrix Steps");
	printf("\n");
}

/* Print title of table */
out_title()
{
	printf("\n");
	switch(mode) {
	case APOLLO:
		printf("   Apollonian gasket");
		break;
	case BP:
		printf("   Blaschke product:");
		printf("  Lambda = %.5lf + %.5lf I",
		param.x,param.y);
		break;
	case BPS:
		printf("   Blaschke products :");
		printf(" Eps %6.2e",eps);
		break;
	case HECKE:
		printf("   Hecke group :");
		printf(" Radius %9.4lf",
			param.x);
		break;
	case HECKES:
		printf("   Hecke Groups :");
		printf(" Eps %6.2e",eps);
		break;
	case QUAD:
		printf("   Julia set for z^2+c :");
		printf(" c = %.5lf + %.5lf I",
		param.x,param.y);
		break;
	case QUADS:
		printf("   Julias sets for z^2+c :");
		printf(" Eps %6.2e",eps);
		break;
	case SCHOTTKY:
		printf("   %d generator Schottky group :",
			ngen);
		printf(" %9.4lf degrees",
			param.x);
		break;
	case SCHOTTKYS:
		printf("   %d Generator Schottky Groups :",
			ngen);
		printf(" Eps %6.2e",eps);
		break;
	case READ:
		printf("   %d Generator Group",
			ngen);
		break;
	}
	printf("\n\n");
}

/* Calculate and output dimension */
out_dim()
{
	double ans;

	ans = dim();
	switch(mode) {
	case SCHOTTKYS:
		printf("%6.1lf ",param.x);
		break;

	case BPS:
	case HECKES:
	case QUADS:
		printf("%8.5lf",param.x);
		break;

	default:
	       	printf("%6.2e",eps);
	}

	printf("  %10.8lf",ans);
	if(verbose)
		printf("%8d%8d%5d",
		ncover,m.nentries,newtons);
	printf("\n");
}

