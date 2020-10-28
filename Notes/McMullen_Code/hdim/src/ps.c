
/* Postscript plotting routines */

#include "hdim.h"

#define PAGEX 540       /* Postscript parameters */
#define PAGEY 720
#define CORNX 36
#define CORNY 36
#define LINEWIDTH 0.2

ps_report()
{
	if(ps) 
	{	ps_open(); 
		ps_plot_cover();
	ps_close();}
}

/* Scale and box */
ps_window(xmin,ymin,xmax,ymax)
	double xmin, ymin, xmax, ymax;
{
	int marginx, marginy;
	double s, sx, sy, tx, ty;

	sx = PAGEX/(xmax-xmin);
	sy = PAGEY/(ymax-ymin);
	s  = (sx < sy) ? sx : sy;
	marginx = (PAGEX-s*(xmax-xmin))/2;
	marginy = (PAGEY-s*(ymax-ymin))/2;

	fprintf(file,"%%%%BoundingBox:  %d %d %d %d\n",
		CORNX+marginx,
		CORNY+marginy,
		CORNX+PAGEX-marginx,
		CORNY+PAGEY-marginy);
	fprintf(file,"%d %d translate\n",
		CORNX+marginx,
		CORNY+marginy);
	fprintf(file,"%.2lf setlinewidth\n",LINEWIDTH);
	fprintf(file,"%%%%\n");
	fprintf(file,"%%%% Change to complex coordinates\n");
	fprintf(file,"%.5lf %.5lf scale\n",s,s);
	fprintf(file,"currentlinewidth %.5lf div setlinewidth\n",s);
	fprintf(file,"%.5lf %.5lf translate\n",-xmin,-ymin); 
	fprintf(file,"newpath\n");
}

ps_circle(t)
	tile t;
{
	fprintf(file,"%.5lf %.5lf %.5lf c\n",
		t.c.x,t.c.y,t.r);
}

ps_open()
{
	int i;

	fprintf(file,"%%!\n");
	fprintf(file,"%%%% Hdim:  Markov covering circles\n");
}

ps_close()
{
	fprintf(file,"showpage\n");
	fflush(file);
}

ps_plot_cover()
{
	int i;
	double bord, xmin, xmax, ymin, ymax;

	if(ncover == 0) return(1);
	xmin = xmax = cover[0].c.x;
	ymin = ymax = cover[0].c.y;
	for(i=0; i<ncover; i++) 
	{	 if(cover[i].c.x-cover[i].r < xmin)
			xmin = cover[i].c.x-cover[i].r;
		 if(cover[i].c.x+cover[i].r > xmax)
			xmax = cover[i].c.x+cover[i].r;
		 if(cover[i].c.y-cover[i].r < ymin)
			ymin = cover[i].c.y-cover[i].r;
		 if(cover[i].c.y+cover[i].r > ymax)
			ymax = cover[i].c.y+cover[i].r;
	}

	fprintf(file,"%%%% Markov partition\n");
	fprintf(file,"%%%%\n");
	bord = (xmax-xmin > ymax-ymin) ? xmax-xmin : ymax-ymin;
	bord = 0.10*bord;
	ps_window(xmin-bord,ymin-bord,xmax+bord,ymax+bord);
	ps_def_circle();
	for(i=0; i<ncover; i++)
		ps_circle(cover[i]);
}

/* Define ps full circle */
ps_def_circle()
{
	fprintf(file,"/c {0 360 arc stroke} def\n");
}
