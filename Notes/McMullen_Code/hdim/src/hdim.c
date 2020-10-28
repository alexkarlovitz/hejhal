/* 
   Calculates Hausdorff dimension of the 
   limit set of an expanding conformal dynamical
   system.

	Global variables
   genc,genr,ngen - generating maps
   cover, ncover - tiling of Markov partition
   m - sparse matrix of derivatives
*/

/* CTM 16 March 1997 */

#include "hdim.h"

/* Declaration of external variables */
int acover, apool, debug, maxit, mode, 
	newtons, ncover, ngen, npool, ps, 
	series, verbose;
double alow, ahigh, astep, emin, emax, eps;
complex param;
FILE *file;

/* Declaration of external arrays */
complex genc [MAXGEN]; 
double  *a, genr [MAXGEN];
sparse_matrix m;
tile *cover;
int *pool; 

/* Declaration of variable functions */
double (*dmap)();
void (*map)();

main(argc,argv)
	int argc;
	char *argv[];
{
	int i;
	double temp;

	space_init();

	for (i=1; i<argc; i++)
	{	if(argv[i][0] != '-') usage();
		switch(argv[i][1]) 
	{
		case 'a':
	mode=APOLLO;
		break;

		case 'b':
	if((i+2) >= argc) usage();
	sscanf(argv[++i],"%lf",&param.x);
	sscanf(argv[++i],"%lf",&param.y);
	mode = BP;
		break;

		case 'B':
	if((i+3) >= argc) usage();
	sscanf(argv[++i],"%lf",&alow);
	sscanf(argv[++i],"%lf",&ahigh);
	sscanf(argv[++i],"%lf",&astep);
	mode=BPS;
		break;

		case 'c':
	if((i+2) >= argc) usage();
	sscanf(argv[++i],"%lf",&param.x);
	sscanf(argv[++i],"%lf",&param.y);
	mode = QUAD;
		break;

		case 'C':
	if((i+3) >= argc) usage();
	sscanf(argv[++i],"%lf",&alow);
	sscanf(argv[++i],"%lf",&ahigh);
	sscanf(argv[++i],"%lf",&astep);
	mode=QUADS;
		break;

		case 'd':
	debug = 1;
		break;

		case 'e':
	if((i+1) >= argc) usage();
	sscanf(argv[++i],"%lf",&eps);
	emax=emin=eps;
		break;

		case 'E':
	if((i+2) >= argc) usage();
	sscanf(argv[++i],"%lf",&emin);
	sscanf(argv[++i],"%lf",&emax);
	if(emin > emax)
	{	temp=emax;
		emax=emin;
		emin=temp;
	}
		break;

		case 'h':
	if((i+1) >= argc) usage();
	sscanf(argv[++i],"%lf",&param.x);
	mode=HECKE;
		break;

		case 'H':
	if((i+3) >= argc) usage();
	sscanf(argv[++i],"%lf",&alow);
	sscanf(argv[++i],"%lf",&ahigh);
	sscanf(argv[++i],"%lf",&astep);
	mode=HECKES;
		break;

		case 'i':
	mode=READ;
		break;

		case 'n':
	if((i+1) >= argc) usage();
	sscanf(argv[++i],"%d",&maxit);
		break;

		case 'p':
	if((i+1) >= argc) usage();
	file = fopen(argv[++i],"w");
	ps=1;
		break;

		case 's':
	if((i+2) >= argc) usage();
	sscanf(argv[++i],"%d",&ngen);
	sscanf(argv[++i],"%lf",&param.x);
	mode=SCHOTTKY;
		break;

		case 'S':
	if((i+4) >= argc) usage();
	sscanf(argv[++i],"%d",&ngen);
	sscanf(argv[++i],"%lf",&alow);
	sscanf(argv[++i],"%lf",&ahigh);
	sscanf(argv[++i],"%lf",&astep);
	mode=SCHOTTKYS;
		break;

		case 'v':
	verbose=!verbose;
		break;

		default:
		usage();
	}
	}

	dims();
}

usage()
{
	fprintf(stderr,"Hausdorff dimension calculator\n");
	fprintf(stderr,"Usage: hdim [options]\n");
	fprintf(stderr,"  -a             Apollonian packing\n");
	fprintf(stderr,"  -b p.x p.y     Julia set of Blaschke product\n");
	fprintf(stderr,"  -B pl ph ps    series of Blaschke products\n");
	fprintf(stderr,"  -c c.x c.y     Julia set of z^2+c\n");
	fprintf(stderr,"  -C cl ch cs    series of real Julia sets\n");
	fprintf(stderr,"  -e eps         circle size epsilon\n");
	fprintf(stderr,"  -E emin emax   range of epsilons\n");
	fprintf(stderr,"  -h r           Hecke group\n");
	fprintf(stderr,"  -H rl rh rs    series of Hecke groups\n");
	fprintf(stderr,"  -i             read circles from stdin\n");
	fprintf(stderr,"  -n maxit       max no. of Newton iterates\n");
	fprintf(stderr,"  -p file        Postcript image of cover\n");
	fprintf(stderr,"  -s n ang       symmetric Schottky group\n");
	fprintf(stderr,"  -S n al ah as  series of Schottky groups\n");
	fprintf(stderr,"  -v             verbose\n");
	fprintf(stderr,"Reads circles (c.x,c.y,r) from stdin\n");
	exit(1);
}
