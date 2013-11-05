/*****************************************************************              
 *  I N S T I T U T   F U E R   E N E R G I E F O R S C H U N G  *              
 *  I E K - 5   P H O T O V O L T A I K                          *              
 *                                                               *              
 *        ########                _   _                          *              
 *     ##########                |_| |_|                         *              
 *    ##########     ##         _ _   _ _     ___ ____ _   _     *              
 *   ##########     ####       | | | | | |   |_ _/ ___| | | |    *              
 *   #########     #####    _  | | | | | |    | | |   | |_| |    *              
 *   #    ###     ######   | |_| | |_| | |___ | | |___|  _  |    *              
 *    #          ######     \___/ \___/|_____|___\____|_| |_|    *              
 *     ##      #######      F o r s c h u n g s z e n t r u m    *              
 *       ##########                                              *              
 *                                                               *              
 *   http://www.fz-juelich.de/iek/iek-5/EN/Home/home_node.html   *              
 *****************************************************************
 *                                                               *
 * Copyright (C) 2013 Dr. Bart E. Pieters                        *
 *                                                               *             
 *****************************************************************/                                                                             
                                                                             
/*  This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see {http://www.gnu.org/licenses/}.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include "main.h"
#include "utils.h"
#include "sourcefield.h"
#include "mod_bessel_01.h"

#define len(a,b) (sqrt((a)*(a)+(b)*(b)))
#define len2(a,b) ((a)*(a)+(b)*(b))
#define veclen(a) (sqrt((a.x)*(a.x)+(a.y)*(a.y)))
#define veclen2(a) ((a.x)*(a.x)+(a.y)*(a.y))
#define MIN(a,b) ((a)<(b) ? (a):(b))
#define MAX(a,b) ((a)<(b) ? (b):(a))

/* maximum dimentions (n*m). The actual dimentions of the field are -n..n and -m..m so the number */
/* of sources is 4*MAX_NXM and the same number for the number of sinks (making a total of 8*MAX_NXM sources+sinks) */
/* MAX_NXM=1e6 seems close enough to infinity for almost all practical purposes */
#define MAX_NXM currentpars.NXM

int VecInRange(vec v, double X1, double X2, double Y1, double Y2)
{
	return (v.x>X1 && v.x<X2 && v.y>Y1 && v.y<Y2);
}
/* compute position of the m,n-th "quasi-origin", i.e. the m,n-th image of the real origin (n,m=0,0) */
static vec quasiOrigin(int n, int m)
{
	vec res={0,0};
	int odd;
	
	if (n<0)
	{
		n=-n;
		odd=n%2;
		res.x=(n+odd)*currentpars.corner1.x - (n-1+(!odd))*currentpars.corner2.x;
	}
	else
	{
		odd=n%2;
		res.x=(n+odd)*currentpars.corner2.x - (n-1+(!odd))*currentpars.corner1.x;
	
	}	
	if (m<0)
	{
		m=-m;
		odd=m%2;
		res.y=(m+odd)*currentpars.corner1.y - (m-1+(!odd))*currentpars.corner2.y;
	}
	else
	{
		odd=m%2;
		res.y=(m+odd)*currentpars.corner2.y - (m-1+(!odd))*currentpars.corner1.y;
	}
	return res;
}

/* compute position of the m,n-th positive contact, coordinate is related to the quasi-origin (see above) */
vec Positive_c(int n, int m)
{
	vec res={0,0};
	int odd;
	
	res=quasiOrigin(n,m);
	
	odd=abs(n%2);
	res.x+=(2*(!odd)-1)*currentpars.contact_p.x;
	
	odd=abs(m%2);
	res.y+=(2*(!odd)-1)*currentpars.contact_p.y;
	
	return res;
}

/* compute position of the m,n-th negative contact, coordinate is related to the quasi-origin (see above)  */
vec Negative_c(int n, int m)
{
	vec res={0,0};
	int odd;
	
	res=quasiOrigin(n,m);
	
	odd=abs(n%2);
	res.x+=(2*(!odd)-1)*currentpars.contact_n.x;
	
	odd=abs(m%2);
	res.y+=(2*(!odd)-1)*currentpars.contact_n.y;
	
	return res;
}

/* compute position of the m,n-th (positive) source (take into account the computed offsets)*/
vec Positive(int n, int m)
{
	vec res={0,0};
	int odd;
	
	res=Positive_c(n, m);
	
	odd=abs(n%2);
	res.x+=(2*(!odd)-1)*currentpars.off_p.x;
	
	odd=abs(m%2);
	res.y+=(2*(!odd)-1)*currentpars.off_p.y;
	
	return res;

}

/* compute position of the m,n-th (negative) sink  (take into account the computed offsets)*/
vec Negative(int n, int m)
{
	vec res={0,0};
	int odd;
	res=Negative_c(n, m);
	
	odd=abs(n%2);
	res.x+=(2*(!odd)-1)*currentpars.off_n.x;
	
	odd=abs(m%2);
	res.y+=(2*(!odd)-1)*currentpars.off_n.y;
	
	return res;
}


static double current_pos(vec Cp, vec Cn, int new)
{
	/* we take a line confined to the lamina at hand and located between the source and sink */
	/* The sum of all currents induced y all sources and sinks is equal to 2pi (assuming q=1)*/
	/* The current is computed to serve as a stop criterium */
	/* When the current is correct there is little current leakage across the edges */
	
	/*          start                                                        
                    o                                                           
                    -\                                                          
                    - \                                                         
                   -   \                                                        
                   -    \                                                       
                  -      \                                                      
                  -       \                                                     
                 - -  alpha\       alpha                                        
                 -   \      \ I = -------                                       
          source O    |    ---->    2 pi                                        
                  --- /       \                                                 
                     ---       \                                                
                        ---     \                                               
                           ---   \                                              
                               ---o                                             
                                  end  
	
	   Take the above illustration.
	   We have a source and a line between the points start and end. As the potentials and currents are computed 
	   with superposition we can simply take that the current from the source flows radially away from the source.
	   This means that if we know alpha we know the portion of the current from the source that crosses the line.
	   
	   The first run we set the points start and end such that the two contacts are on opposite sides of the line. 
	   After that we only have 
	   to take care that:
	   - sources and drains create currents of opposite signs
	   - A source at one side of the line creates a current with opposite sign compared to a source on the other 
	     side of the line.
	*/

	double res=0;	
	vec dse;
	static vec start;
	static vec end;
	static vec dir_t;
	static vec dir_p;
	if (new)
	{
		/* Each new problem we first determine which line to integrate the current through */
		/* line spanned by start and end */
		
		if (fabs(Cp.x-Cn.x)>fabs(Cp.y-Cn.y))
		{
			start.x=(Cp.x+Cn.x)/2;
			end.x=(Cp.x+Cn.x)/2;
			start.y=currentpars.corner1.y;
			end.y=currentpars.corner2.y;
			
			/* current runs from the positive contact to the negative contact at n,m=0,0 */
			/* the computation here gives +1 or -1 such that the direction is correct */
			dir_t.x=0;
			dir_t.y=(2*(double)(start.y-Cp.y>0)-1.0);
			dir_p.x=(2*(double)(start.x-Cp.x>0)-1.0);
			dir_p.y=0;
		}
		else
		{
			start.y=(Cp.y+Cn.y)/2;
			end.y=(Cp.y+Cn.y)/2;
			start.x=currentpars.corner1.x;
			end.x=currentpars.corner2.x;
			dir_t.x=(2*(double)(start.x-Cp.x>0)-1.0);
			dir_t.y=0;
			dir_p.x=0;
			dir_p.y=(2*(double)(start.y-Cp.y>0)-1.0);
		}
		Print(DEBUG,"---------------------------------------------------------------------\n");
		Print(DEBUG,"Integration line (x1,x1) (x2,y2):\n");
		Print(DEBUG,"(%13.6e,%13.6e)\t(%13.6e,%13.6e)\n", start.x,start.y,end.x,end.y);
		Print(DEBUG,"Contact positions (xp,yp) (xn,yn):\n");
		Print(DEBUG,"(%13.6e,%13.6e)\t(%13.6e,%13.6e)\n", Cp.x,Cp.y,Cn.x,Cn.y);
		Print(DEBUG,"---------------------------------------------------------------------\n\n");
	}
	dse.x=(start.x-Cp.x);
	dse.y=(start.y-Cp.y);
	res+=atan((dir_t.x*dse.x+dir_t.y*dse.y)/(dir_p.x*dse.x+dir_p.y*dse.y));
	
	dse.x=(start.x-Cn.x);
	dse.y=(start.y-Cn.y);
	res-=atan((dir_t.x*dse.x+dir_t.y*dse.y)/(dir_p.x*dse.x+dir_p.y*dse.y));
	
	dse.x=(end.x-Cp.x);
	dse.y=(end.y-Cp.y);
	res-=atan((dir_t.x*dse.x+dir_t.y*dse.y)/(dir_p.x*dse.x+dir_p.y*dse.y));
	
	dse.x=(end.x-Cn.x);
	dse.y=(end.y-Cn.y);
	res+=atan((dir_t.x*dse.x+dir_t.y*dse.y)/(dir_p.x*dse.x+dir_p.y*dse.y));
	return res/TWOPI;
}

static double voltage_pos_avg(vec Cp, vec Cn, double *Vn, double *Vp)
#define NANG 10
{
	/* this routine computes the voltage difference induced by the sources and sinks averaged along the contour of the contact*/
	/* this routine is used when computing the potential of the original source and sink. For the images I use the voltage_pos routine */
	/* increments the contact potentials Vn and Vp and returns the delta in voltage due to source and sink Cp and Cn */
	double dVp=1, dVn=1, a;
	int k;
	vec Cc;
	for (k=0;k<NANG;k++)
	{
		a=TWOPI*((double)k)/NANG;
		Cc.x=currentpars.contact_p.x+currentpars.rp*cos(a);
		Cc.y=currentpars.contact_p.y+currentpars.rp*sin(a);	
		dVp*=(len2((Cc.x-Cn.x),Cc.y-Cn.y)/len2((Cc.x-Cp.x),Cc.y-Cp.y));
	}
	dVp=pow(dVp,1.0/(NANG*2.0));
	for (k=0;k<NANG;k++)
	{
		a=TWOPI*((double)k)/NANG;
		Cc.x=currentpars.contact_n.x+currentpars.rn*cos(a);
		Cc.y=currentpars.contact_n.y+currentpars.rn*sin(a);
		dVn*=(len2((Cc.x-Cn.x),Cc.y-Cn.y)/len2((Cc.x-Cp.x),Cc.y-Cp.y));
	}
	dVn=pow(dVn,1.0/(NANG*2.0));
	
	(*Vn)*=dVn;
	(*Vp)*=dVp;
	return log(dVp/dVn);
}
#undef NANG


static double voltage_pos(vec Cp, vec Cn, double *Vn, double *Vp)
{
	/* this routine computes the voltage difference induced by the image sources and sinks */
	/* note that a routine like voltage_pos_avg is more accurate but this routine is only used */
	/* to get a (quite good) estimate of the changes in potentials induced by the images */
	/* and thus is sufficient for the algorithm to function. As we have to monitor the degree to */
	/* which the contacts are equipotential anyway we can determine a good value for the average */
	/* potential in that routine */ 
	/* increments the contact potentials Vn and Vp and returns the delta in voltage due to source and sink Cp and Cn */
	double dVp=1, dVn=1;
	vec Cc;
	Cc.x=currentpars.contact_p.x;
	Cc.y=currentpars.contact_p.y;	
	dVp*=(len2((Cc.x-Cn.x),Cc.y-Cn.y)/len2((Cc.x-Cp.x),Cc.y-Cp.y));
	
	Cc.x=currentpars.contact_n.x;
	Cc.y=currentpars.contact_n.y;
	dVn*=(len2((Cc.x-Cn.x),Cc.y-Cn.y)/len2((Cc.x-Cp.x),Cc.y-Cp.y));
	
	(*Vn)*=dVn;
	(*Vp)*=dVp;
	return log(dVp/dVn);
}

static void Incr_n(int *n, int *m, double *dI, double *dV, double *Vn, double *Vp)
{
	/* This routine increments n and computes the resulting current and voltage change */
	vec Cp, Cn;
	double dv=0, di=0; /* help variables needed for making the code parallel */
	double vvp=1,vvn=1; /* help variables needed for making the code parallel */
	int j;
	(*n)++;
	
	#pragma omp parallel for private(j,Cn,Cp) reduction(+:di, dv) reduction(*:vvn, vvp)
	for(j=-(*m);j<=(*m);j++)
	{
		/* ++ */
		Cp=Positive((*n), j);
		Cn=Negative((*n), j);
		di+=current_pos(Cp, Cn, 0);
		dv+=voltage_pos(Cp, Cn, &vvn, &vvp);
		
		/* -+ */
		Cp=Positive(-(*n), j);
		Cn=Negative(-(*n), j);
		di+=current_pos(Cp, Cn, 0);
		dv+=voltage_pos(Cp, Cn, &vvn, &vvp);
	}
	(*Vp)*=vvp;
	(*Vn)*=vvn;
	(*dI)=di;
	(*dV)=dv;
}

static void Incr_m(int *n, int *m, double *dI, double *dV, double *Vn, double *Vp)
{
	/* This routine increments m and computed the resulting current and voltage change */
	vec Cp, Cn;
	int i;
	double dv=0, di=0; /* help variables needed for making the code parallel */
	double vvp=1,vvn=1; /* help variables needed for making the code parallel */
	(*m)++;
	
	
	#pragma omp parallel for private(i,Cn,Cp) reduction(+:di, dv) reduction(*:vvn, vvp)
	for(i=-(*n);i<=(*n);i++)
	{
		/* ++ */
		Cp=Positive(i, (*m));
		Cn=Negative(i, (*m));
		di+=current_pos(Cp, Cn, 0);
		dv+=voltage_pos(Cp, Cn, &vvn, &vvp);
		
		/* +- */
		Cp=Positive(i, -(*m));
		Cn=Negative(i, -(*m));
		di+=current_pos(Cp, Cn, 0);
		dv+=voltage_pos(Cp, Cn, &vvn, &vvp);
		
	}
	(*Vp)*=vvp;
	(*Vn)*=vvn;
	(*dI)=di;
	(*dV)=dv;
}


void Offsets(void)
/*

             ^
             +
             |
             x
             x
             x
             x
             A
            ===
           //_\\
         _/X.-.X\_
  ______/XX/___\XX\__________________________________09
*/
#define ERR currentpars.Err
#define ERR2 currentpars.Err_Off
#define RP currentpars.rp
#define RN currentpars.rn
#define MAX_ITER currentpars.Off_Iter
{
	
	/* this routine computes offset vectors for the source and sink of the contacts such that the contour of the 
	   circular contact is approximately an equipotential line */ 
	int i,j, iter=0;
	int n,m, fallback=0;
	vec off_p={0,0}, off_n={0,0};
	vec Cp, Cn;
	vec p0, n0;
	double d1, d2, d3, d4, dmin;
	double ra, b;
	double ppx, ppy, pnx, pny;
	
	/* position of the main source and sink, without offsets */
	currentpars.off_p=off_p;
	currentpars.off_n=off_p;
	p0=Positive(0,0);
	n0=Negative(0,0);
	
	
	/* smallest distance in d1, for stop criterium */
	d1=currentpars.contact_p.x-currentpars.rp-currentpars.corner1.x;
	d2=-currentpars.contact_p.x-currentpars.rp+currentpars.corner2.x;	
	d1=MIN(d1,d2);
	
	d2=currentpars.contact_p.y-currentpars.rp-currentpars.corner1.y;
	d1=MIN(d1,d2);
	
	d2=-currentpars.contact_p.y-currentpars.rp+currentpars.corner2.y;
	d1=MIN(d1,d2);
		
	d2=currentpars.contact_n.x-currentpars.rn-currentpars.corner1.x;
	d1=MIN(d1,d2);
		
	d2=-currentpars.contact_n.x-currentpars.rn+currentpars.corner2.x;
	d1=MIN(d1,d2);
		
	d2=currentpars.contact_n.y-currentpars.rn-currentpars.corner1.y;
	d1=MIN(d1,d2);
		
	d2=-currentpars.contact_n.y-currentpars.rn+currentpars.corner2.y;
	d1=MIN(d1,d2);
		
	d2=len(currentpars.contact_n.x-currentpars.contact_p.x,currentpars.contact_n.y-currentpars.contact_p.y)/2;
	d1=MIN(d1,d2);
	
	dmin=d1;

	/*
		n and m are computed according to:
			_  2 r _ 
		N'  =  | ------ |
		         eps Dx  
		
		see Eq. 23 in documentation. (eps is ERR2 in this routine)
	*/
	n=(int)floor(2*MAX(RN,RP)/(ERR2*(currentpars.corner2.x-currentpars.corner1.x))+1.0);
	m=(int)floor(2*MAX(RN,RP)/(ERR2*(currentpars.corner2.y-currentpars.corner1.y))+1.0);


	Print(VERBOSE,"Computing source/sink positional offsets w.r.t the contact positions\n");
	Print(DEBUG,"Field dimensions in Offset calculation (n,m): (%i,%i)\n", n,m);
	Print(VERBOSE,"Offset p (x,y)\t\t\tOffset n (x,y)\t\t\tIter\n");
	Print(VERBOSE,"---------------------------------------------------------------------\n", n,m);
	do
	{
		currentpars.off_p=off_p;
		currentpars.off_n=off_n; 
		Print(DEBUG,"(%13.6e,%13.6e)\t(%13.6e,%13.6e)\t%i\n", currentpars.off_p.x,currentpars.off_p.y,currentpars.off_n.x,currentpars.off_n.y, iter);
		fflush(stdout);
		ppx=1;
		ppy=1;
		pnx=1;
		pny=1;
		#pragma omp parallel for private(i,j,Cn,Cp,d1,d2,d3,d4)  reduction(*:ppx, ppy, pnx, pny) 
		for(i=-n;i<=n;i++)
		for(j=-m;j<=m;j++)
		{
			if (!(i==0 && j==0))
			{
				/* interact both contacts with both image contacts */
				/* get positions */
				Cp=Positive(i,j);
				Cn=Negative(i,j);
				
				/* get lengths^2*/
				d1=len2(p0.x-Cp.x-RP, p0.y-Cp.y);
				d2=len2(p0.x-Cp.x+RP, p0.y-Cp.y);
				d3=len2(p0.x-Cn.x-RP, p0.y-Cn.y);
				d4=len2(p0.x-Cn.x+RP, p0.y-Cn.y);
				
				/* add offsets to new offsetvector */
				ppx*=((d1*d4)/(d2*d3));
				/* get lengths^2 */
				d1=len2(n0.x-Cp.x-RN, n0.y-Cp.y);
				d2=len2(n0.x-Cp.x+RN, n0.y-Cp.y);
				d3=len2(n0.x-Cn.x-RN, n0.y-Cn.y);
				d4=len2(n0.x-Cn.x+RN, n0.y-Cn.y);
				
				/* add offsets to new offsetvector */
				pnx*=((d2*d3)/(d1*d4));
				
				
				/* get lengths^2 */
				d1=len2(p0.x-Cp.x, p0.y-Cp.y-RP);
				d2=len2(p0.x-Cp.x, p0.y-Cp.y+RP);
				d3=len2(p0.x-Cn.x, p0.y-Cn.y-RP);
				d4=len2(p0.x-Cn.x, p0.y-Cn.y+RP);
				
				/* add offsets to new offsetvector */
				ppy*=((d1*d4)/(d2*d3));
				
				/* get lengths^2 */
				d1=len2(n0.x-Cp.x, n0.y-Cp.y-RN);
				d2=len2(n0.x-Cp.x, n0.y-Cp.y+RN);
				d3=len2(n0.x-Cn.x, n0.y-Cn.y-RN);
				d4=len2(n0.x-Cn.x, n0.y-Cn.y+RN);
				
				/* add offsets to new offsetvector */
				pny*=((d2*d3)/(d1*d4));
				
			}
			else
			{
				/* no image contacts, only real ones. i.e. 
				   only interact with the other contact */
				/* get positions */
				Cp=Positive(i,j);
				Cn=Negative(i,j);
				
				/* get lengths^2 */
				d3=len2(p0.x-Cn.x-RP, p0.y-Cn.y);
				d4=len2(p0.x-Cn.x+RP, p0.y-Cn.y);
				
				/* add offsets to new offsetvector */
				ppx*=(d4/d3);
				
				/* get lengths^2 */
				d1=len2(n0.x-Cp.x-RN, n0.y-Cp.y);
				d2=len2(n0.x-Cp.x+RN, n0.y-Cp.y);
				
				/* add offsets to new offsetvector */
				pnx*=(d2/d1);		
				
				/* get lengths^2 */
				d3=len2(p0.x-Cn.x, p0.y-Cn.y-RP);
				d4=len2(p0.x-Cn.x, p0.y-Cn.y+RP);
				
				/* add offsets to new offsetvector */
				ppy*=(d4/d3);
				
				/* get lengths^2 */
				d1=len2(n0.x-Cp.x,n0.y-Cp.y-RN);
				d2=len2(n0.x-Cp.x,n0.y-Cp.y+RN);
				
				/* add offsets to new offsetvector */
				pny*=(d2/d1);		
			}
			
		}
		off_p.x=(1-ppx)/(1+ppx);
		off_p.y=(1-ppy)/(1+ppy);		
		ra=veclen2(off_p);
		if(ra<1)		
			b=RP*(1-sqrt(1-ra))/ra;
		else
		{
			/* If this happens our algorithm fails. */
			/* This basically means there is no position */
			/* for the source which would make v1=v2 and v3=v4 */
			Warning("Cannot compute a good source offset. Try smaller contacts\n");
			fallback=1;
			/* just do something */
			b=0;				
		}
		off_p.x*=b;
		off_p.y*=b;
				
		off_n.x=(1-pnx)/(1+pnx);
		off_n.y=(1-pny)/(1+pny);		
		ra=veclen2(off_n);
		if(ra<1)	
			b=RN*(1-sqrt(1-ra))/ra;
		else
		{
			/* Cannot make approximately equipotential contacts */
			Warning("Cannot compute a good sink offset. Try smaller contacts\n");
			fallback=1;				
			b=0;		
		}
		off_n.x*=b;
		off_n.y*=b;
		d1=veclen(off_p);
		d2=veclen(currentpars.off_p);
		d3=fabs(d1-d2)/dmin;
		d1=veclen(off_n);
		d2=veclen(currentpars.off_n);
		d4=fabs(d1-d2)/dmin;		
		iter++;
	} while (d3+d4>ERR && iter<MAX_ITER && (!fallback));
	currentpars.off_p=off_p;
	currentpars.off_n=off_n;
	Print(VERBOSE,"(%13.6e,%13.6e)\t(%13.6e,%13.6e)\t%i\n", currentpars.off_p.x,currentpars.off_p.y,currentpars.off_n.x,currentpars.off_n.y, iter);
	if(iter==MAX_ITER)
		Warning("Warning, maximum iterations reached in Offsets\n");
	Print(VERBOSE,"---------------------------------------------------------------------\n\n");
}
#undef ERR
#undef ERR2
#undef RN
#undef RP




void ErrorEqPot(int n, int m, double *Vn, double *Vp, double *Emax, double *Erms)
#define NANG currentpars.Nangle
{
	/* this routine computes to what degree the contour of the circular contacts are equipotential lines */  
	/* furthermore it computes the average potential along the contour of the contact and sets Vp and Vn */  
	/* accordingly. (usually the values for Vn and Vp are quite good) */  
	double a;
	vec Cp, Cn, Cc;
	int i,j,k;
	double *E, Vx;
	
	E=malloc((2*NANG+1)*sizeof(double));
	
	(*Emax)=0;
	(*Erms)=0;
	
	for (k=0;k<NANG;k++)
	{
		a=TWOPI*((double)k)/NANG;
		Cc.x=currentpars.contact_p.x+currentpars.rp*cos(a);
		Cc.y=currentpars.contact_p.y+currentpars.rp*sin(a);
		Vx=1;
		#pragma omp parallel for private(i,j,Cn,Cp) reduction(*:Vx) 
		for (i=-n;i<=n;i++)
			for (j=-m;j<=m;j++)
			{
				Cp=Positive(i, j);
				Cn=Negative(i, j);
				Vx*=(len((Cc.x-Cn.x),Cc.y-Cn.y)/len((Cc.x-Cp.x),Cc.y-Cp.y));
			}
		Vx=log(Vx);
		E[k]=Vx/TWOPI;
	}
	for (k=0;k<NANG;k++)
	{
		a=TWOPI*((double)k+0.5)/NANG;
		Cc.x=currentpars.contact_n.x+currentpars.rn*cos(a);
		Cc.y=currentpars.contact_n.y+currentpars.rn*sin(a);
		Vx=1;
		#pragma omp parallel for private(i,j,Cn,Cp) reduction(*:Vx) 
		for (i=-n;i<=n;i++)
			for (j=-m;j<=m;j++)
			{
				Cp=Positive(i, j);
				Cn=Negative(i, j);
				Vx*=(len((Cc.x-Cn.x),Cc.y-Cn.y)/len((Cc.x-Cp.x),Cc.y-Cp.y));
			}
		Vx=log(Vx);
		E[NANG+k]=Vx/TWOPI;
	}
	Vx=0;
	Print(VERBOSE,"Vp         Vp-avg     E            Vn         Vn-avg     E\n");
	Print(VERBOSE,"---------------------------------------------------------------------\n", n,m);
	Print(VERBOSE,"%-11.3e", *Vp);
	for (k=0;k<NANG;k++)
		Vx+=E[k];
	Vx/=NANG;
	Print(VERBOSE,"%-11.3e%-11.3e  ", Vx, ((*Vp)-Vx)/Vx);
	(*Vp)=Vx;
	
	Vx=0;
	Print(VERBOSE,"%-11.3e", *Vn);
	for (k=0;k<NANG;k++)
		Vx+=E[NANG+k];
	Vx/=NANG;
	Print(VERBOSE,"%-11.3e%-11.3e\n", Vx, ((*Vn)-Vx)/Vx);
	Print(VERBOSE,"---------------------------------------------------------------------\n\n");
	(*Vn)=Vx;
	
	for (k=0;k<2*NANG;k++)
	{
		if (k<NANG)
			Vx=(*Vp);
		else
			Vx=(*Vn);
			
		E[k]=E[k]/Vx-1;
		if(E[k]*E[k]>(*Emax))
			(*Emax)=E[k]*E[k];
		(*Erms)+=E[k]*E[k];
	}
	free(E);
	(*Emax)=sqrt((*Emax));
	(*Erms)=sqrt((*Erms)/(2*NANG));
}
#undef NANG

void Dim(double *Vn, double *Vp, double *I, int *n, int *m, double *Emax, double *Erms)
{
	/* This routine increments n and m through the increment functions above */
	/* The routine stops when the required accuracies for the current and voltage are reached */
	/* V, I, n, m are set appropriately */
	double dIm, dIn, dVn, dVm;
	vec Cp, Cn;
	(*n)=0;
	(*m)=0;
	(*I)=0;
	(*Vn)=1;
	(*Vp)=1;
	
	/* compute offsets */
	/* Offsets();*/
		
	Print(VERBOSE,"Computing potentials of the contacts\n");
	Cp=Positive(0, 0);
	Cn=Negative(0, 0);
	(*I)+=current_pos(Cp, Cn, 1);
	voltage_pos_avg(Cp, Cn, Vn, Vp);
	Print(DEBUG,"n\tm\t  exp(Vn) (1A)\t  exp(Vp) (1A)\t  ERR I\n");
	Print(DEBUG,"---------------------------------------------------------------------\n");
	Print(DEBUG,"%i\t%i\t%14.6e\t%14.6e\t%14.6e\n", (*n),(*m),(*Vn), (*Vp),(*I)-1);
	
	Incr_n(n, m, &dIn, &dVn, Vn, Vp);
	(*I)+=dIn;
	/* dIn=fabs(dIn); */
	dVn=fabs(dVn);
	Print(DEBUG,"%i\t%i\t%14.6e\t%14.6e\t%14.6e\n", (*n),(*m),(*Vn), (*Vp),(*I)-1);
	Incr_m(n, m, &dIm, &dVm, Vn, Vp);
	(*I)+=dIm;
	/* dIm=fabs(dIm); */
	dVm=fabs(dVm);
	
	while( (fabs((*I) - 1)>currentpars.Err) && ((*n)*(*m)<MAX_NXM))
	{
		Print(DEBUG,"%i\t%i\t%14.6e\t%14.6e\t%14.6e\n", (*n),(*m),(*Vn), (*Vp),(*I)-1);
		/* the question comes which criterion to use to decide between incrementing n or m */
		/* An obvious method is to decide based on the changes induced by the last increment. */
		/* i.e., increment the one which, with a previous increment, induced the largest change in */
		/* voltage, current or power. All these methods seem to work. Which is best seems situation */
		/* specific, where the current criterion is influenced by the direction of the integration */
		/* line (see current_pos routine), favoring the direction perpendicular to the line. As the */
		/* latter is not really what *should* determine what is to be incremented (or is it?) I */
		/* opted for the voltage. */
		
		if(dVm>dVn)
		{
			Incr_m(n, m, &dIm, &dVm, Vn, Vp);
			(*I)+=dIm;
			/* dIm=fabs(dIm); */
			dVm=fabs(dVm);
		}
		else
		{
			Incr_n(n, m, &dIn, &dVn, Vn, Vp);
			(*I)+=dIn;
			/* dIn=fabs(dIn); */
			dVn=fabs(dVn);
		}
	}
	Print(DEBUG,"---------------------------------------------------------------------\n\n");
	(*Vn)=log((*Vn))/TWOPI;
	(*Vp)=log((*Vp))/TWOPI;
	ErrorEqPot((*n), (*m), Vn, Vp, Emax, Erms);
	
	Print(VERBOSE,"n\tm\t  Vn (1A)\t  Vp (1A)\t  ERR I\n");
	Print(VERBOSE,"---------------------------------------------------------------------\n");
	Print(VERBOSE,"%i\t%i\t%14.6e\t%14.6e\t%14.6e\n",  (*n),(*m),(*Vn), (*Vp),(*I)-1);
	if ((*n)*(*m)>=MAX_NXM)
		Warning("Warning, solution not converged!\n");
	
	Print(VERBOSE,"---------------------------------------------------------------------\n\n");
}

double Potential(vec p, int n, int m, double Vn, double Vp)
{
	/* this routine computes the potential at position p */
	double res=1, dn, dp, dc;
	vec Cp, Cn;
	double *RV=NULL;
	int i,j;
	
	#pragma omp parallel for private(i,j,Cn,Cp,dc) 
	for (i=-n;i<=n;i++)
		for (j=-m;j<=m;j++)
		{
			Cp=Positive_c(i, j);
			dc=len((p.x-Cp.x),p.y-Cp.y);
			if (dc<currentpars.rp)
				RV=(&Vp);
			Cp=Negative_c(i, j);
			dc=len((p.x-Cp.x),p.y-Cp.y);
			if (dc<currentpars.rn)
				RV=(&Vn);		
		
		}
	if (RV)
		return *RV;
	#pragma omp parallel for private(i,j,Cn,Cp,dp,dn) reduction(*:res) 
	for (i=-n;i<=n;i++)
		for (j=-m;j<=m;j++)
		{
			/* ++ */
			Cp=Positive(i, j);
			Cn=Negative(i, j);
			dp=len((p.x-Cp.x),p.y-Cp.y);
			dn=len((p.x-Cn.x),p.y-Cn.y);
			res*=(dn/dp);
		}
	return log(res)/TWOPI;
}

vec Field(vec p, int n, int m)
{
	/* this routine computes the electric field at position p */
	vec res={0,0}, Cp, Cn, v0={0,0};
	vec *RV=NULL;
	double dp, dn ,dc, ex=0, ey=0;
	int i,j;
	
	#pragma omp parallel for private(i,j,Cn,Cp,dc) 
	for (i=-n;i<=n;i++)
		for (j=-m;j<=m;j++)
		{
			Cp=Positive_c(i, j);
			dc=len((p.x-Cp.x),p.y-Cp.y);
			if (dc<currentpars.rp)
				RV=(&v0);
			Cp=Negative_c(i, j);
			dc=len((p.x-Cp.x),p.y-Cp.y);
			if (dc<currentpars.rn)
				RV=(&v0);		
		
		}
	if (RV)
		return *RV;
		
	#pragma omp parallel for private(i,j,Cn,Cp,dp,dn) reduction(+:ex, ey) 
	for (i=-n;i<=n;i++)
		for (j=-m;j<=m;j++)
		{
			/* ++ */
			Cp=Positive(i, j);
			Cn=Negative(i, j);
			dp=len2((p.x-Cp.x),p.y-Cp.y);
			dn=len2((p.x-Cn.x),p.y-Cn.y);
			ex+=((p.x-Cn.x)/dn-(p.x-Cp.x)/dp);
			ey+=((p.y-Cn.y)/dn-(p.y-Cp.y)/dp);
		}
	res.x=ex/TWOPI;
	res.y=ey/TWOPI;
	return res;
}


vec * FieldLine(vec s, int n, int m, double Vn, double Vp, double X1, double X2, double Y1, double Y2, int *len, double DL)
{
	/* follow a fieldline from starting point s (and its potential) */
	/* The starting point (s) should be close to either the source of sink of a source-sink pair. The endpoint then is */
	/* when we get out of the drawing  area or reach the opposite source/sink */
	/* This seems simple enough (and it is) but inherently to the method used here the error accumilates which requires */
	/* small steps, which in turn easily leads to excessive computation times. */
	
	/* to keep it under controll we adapt the step size along the way to the field strength *and the change of */
	/* direction of the field between two succesive points along the field-line*. */ 
	
	/* The parameter DL controlls the step size. */
	/* This algorith sortta worked for me but I make no guarantees for it's accuracy or usefullness. */
	/* You may have to adapt the algorithm to your situation as the routine is only tested for a few cases. */
	   
	vec *res, f1={0,0},f2={0,1};
	
	double psi, dir=1, l, df, Psi_end, d1, d2;
	unsigned int alloc=100;
	
	res=malloc(alloc* sizeof(vec));
	(*len)=0;
	
	psi=Potential(s, n, m, Vn, Vp);
	if(psi-Vn>Vp-psi)
		Psi_end=Vn;
	else
		Psi_end=Vp;
	
	if (psi>Psi_end)
		dir=-1.0;
	f1=Field(s, n, m);
	l=veclen(f1);
	
	/* smallest distance in d1, to scale DL */
	d1=currentpars.contact_p.x-currentpars.rp-currentpars.corner1.x;
	d2=-currentpars.contact_p.x-currentpars.rp+currentpars.corner2.x;
	d1=MIN(d1,d2);
		
	d2=currentpars.contact_p.y-currentpars.rp-currentpars.corner1.y;
	d1=MIN(d1,d2);
		
	d2=-currentpars.contact_p.y-currentpars.rp+currentpars.corner2.y;
	d1=MIN(d1,d2);
		
	d2=currentpars.contact_n.x-currentpars.rn-currentpars.corner1.x;
	d1=MIN(d1,d2);
		
	d2=-currentpars.contact_n.x-currentpars.rn+currentpars.corner2.x;
	d1=MIN(d1,d2);
		
	d2=currentpars.contact_n.y-currentpars.rn-currentpars.corner1.y;
	d1=MIN(d1,d2);
		
	d2=-currentpars.contact_n.y-currentpars.rn+currentpars.corner2.y;
	d1=MIN(d1,d2);
		
	d2=len(currentpars.contact_n.x-currentpars.contact_p.x,currentpars.contact_n.y-currentpars.contact_p.y)/2;
	d1=MIN(d1,d2);
		

	DL*=d1;
	
	Print(DEBUG,"df\tx\ty\tEx\tEy\tdir\n");	
	while(VecInRange(s, X1, X2, Y1, Y2) && dir*(Psi_end-psi)>0)
	{
		if ((*len)==(int)alloc-1)
		{
			alloc+=100;
			res=realloc(res, alloc*sizeof(vec));
		}
		res[(*len)++]=s;
		
		/* f2 is the normalized version of f1 */
		f2=f1;
		f2.x/=l;
		f2.y/=l;
				
		/* compute the field and field-strength at the current position */
		f1=Field(s, n, m);
		l=veclen(f1);
		
		df=l*len(f2.x-f1.x/l,f2.y-f1.y/l);
		
		Print(DEBUG,"%14.6e %14.6e %14.6e %14.6e %14.6e %14.6e %14.6e\n",df,l, s.x, s.y, f1.x, f1.y, dir);
		
		/* empirically the step size DL/((df+1)*l) works pretty well, in my test cases at least */
		s.x+=dir*f1.x*DL/((df+1)*l);
		s.y+=dir*f1.y*DL/((df+1)*l);
			
		psi=Potential(s, n, m, Vn, Vp);
	}
	if ((*len)==(int)alloc-1)
	{
		alloc+=100;
		res=realloc(res, alloc*sizeof(vec));
	}
	res[(*len)++]=s;
	res=realloc(res, ((unsigned int)(*len)+1)*sizeof(vec));
	return res;
}

void InitBessel(int n, int m, double rg, double Rsq, double *rr, double *In, double *Ip, double *Vn, double *Vp, double *Emax, double *Erms)
#define NANG currentpars.Nangle
{
	
	double a, dn, dp;
	vec Cp, Cn, Cc;
	int i,j,k;
	double *E, Vx;
	/* determine spatial scaling factor */
	(*rr)=sqrt(Rsq/rg);
	
	/* compute the current out of the contact with the derivative K0(x).

		d K0 (x)             
		--------  =   - K1(x)
		   d x   
	*/            

	(*Ip)=(*rr)*bessk1((*rr)*currentpars.rp)*currentpars.rp;	
	(*In)=-(*rr)*bessk1((*rr)*currentpars.rn)*currentpars.rn;
		
	E=malloc((2*NANG+1)*sizeof(double));
	
	(*Emax)=0;
	(*Erms)=0;
	
	for (k=0;k<NANG;k++)
	{
		a=TWOPI*((double)k)/NANG;
		Cc.x=currentpars.contact_p.x+currentpars.rp*cos(a);
		Cc.y=currentpars.contact_p.y+currentpars.rp*sin(a);
		Vx=0;
		#pragma omp parallel for private(i,j,Cn,Cp,dp,dn) reduction(+:Vx) 
		for (i=-n;i<=n;i++)
			for (j=-m;j<=m;j++)
			{
				Cp=Positive(i, j);
				Cn=Negative(i, j);
				dp=len((Cc.x-Cp.x),Cc.y-Cp.y);
				dn=len((Cc.x-Cn.x),Cc.y-Cn.y);
				Vx+=(bessk0((*rr)*dp)-bessk0((*rr)*dn));
			}
		E[k]=Vx/TWOPI;
	}
	for (k=0;k<NANG;k++)
	{
		a=TWOPI*((double)k+0.5)/NANG;
		Cc.x=currentpars.contact_n.x+currentpars.rn*cos(a);
		Cc.y=currentpars.contact_n.y+currentpars.rn*sin(a);
		Vx=0;
		#pragma omp parallel for private(i,j,Cn,Cp,dp,dn) reduction(+:Vx) 
		for (i=-n;i<=n;i++)
			for (j=-m;j<=m;j++)
			{
				Cp=Positive(i, j);
				Cn=Negative(i, j);
				dp=len((Cc.x-Cp.x),Cc.y-Cp.y);
				dn=len((Cc.x-Cn.x),Cc.y-Cn.y);
				Vx+=(bessk0((*rr)*dp)-bessk0((*rr)*dn));
			}
		E[NANG+k]=Vx/TWOPI;
	}
	Vx=0;
	Print(VERBOSE,"Vp         Vp-avg     E            Vn         Vn-avg     E\n");
	Print(VERBOSE,"---------------------------------------------------------------------\n", n,m);
	Print(VERBOSE,"%-11.3e", *Vp);
	for (k=0;k<NANG;k++)
		Vx+=E[k];
	Vx/=NANG;
	Print(VERBOSE,"%-11.3e%-11.3e  ", Vx, ((*Vp)-Vx)/Vx);
	(*Vp)=Vx;
	
	Vx=0;
	Print(VERBOSE,"%-11.3e", *Vn);
	for (k=0;k<NANG;k++)
		Vx+=E[NANG+k];
	Vx/=NANG;
	Print(VERBOSE,"%-11.3e%-11.3e\n", Vx, ((*Vn)-Vx)/Vx);
	Print(VERBOSE,"---------------------------------------------------------------------\n\n");
	(*Vn)=Vx;
	
	for (k=0;k<2*NANG;k++)
	{
		if (k<NANG)
			Vx=(*Vp);
		else
			Vx=(*Vn);
			
		E[k]=E[k]/Vx-1;
		if(E[k]*E[k]>(*Emax))
			(*Emax)=E[k]*E[k];
		(*Erms)+=E[k]*E[k];
	}
	free(E);
	(*Emax)=sqrt((*Emax));
	(*Erms)=sqrt((*Erms)/(2*NANG));
	

}
#undef NANG

void Bessel_Offsets(double rr)
#define ERR currentpars.Err
#define ERR2 currentpars.Err_Off
#define RP currentpars.rp
#define RN currentpars.rn
#define MAX_ITER currentpars.Off_Iter
{
	
	/* this routine computes offset vectors for the source and sink of the contacts such that the contour of the 
	   circular contact is approximately an equipotential line */ 
	int i,j, iter=0;
	int n,m, fallback=0;
	vec off_p={0,0}, off_n={0,0};
	vec Cp, Cn;
	vec p0, n0;
	double d1, d2, d3, d4, dmin;
	double ra, b;
	double ppx, ppy, pnx, pny;
	
	/* position of the main source and sink, without offsets */
	currentpars.off_p=off_p;
	currentpars.off_n=off_p;
	p0=Positive(0,0);
	n0=Negative(0,0);
	
	
	/* smallest distance in d1, for stop criterium */
	d1=currentpars.contact_p.x-currentpars.rp-currentpars.corner1.x;
	d2=-currentpars.contact_p.x-currentpars.rp+currentpars.corner2.x;	
	d1=MIN(d1,d2);
	
	d2=currentpars.contact_p.y-currentpars.rp-currentpars.corner1.y;
	d1=MIN(d1,d2);
	
	d2=-currentpars.contact_p.y-currentpars.rp+currentpars.corner2.y;
	d1=MIN(d1,d2);
		
	d2=currentpars.contact_n.x-currentpars.rn-currentpars.corner1.x;
	d1=MIN(d1,d2);
		
	d2=-currentpars.contact_n.x-currentpars.rn+currentpars.corner2.x;
	d1=MIN(d1,d2);
		
	d2=currentpars.contact_n.y-currentpars.rn-currentpars.corner1.y;
	d1=MIN(d1,d2);
		
	d2=-currentpars.contact_n.y-currentpars.rn+currentpars.corner2.y;
	d1=MIN(d1,d2);
		
	d2=len(currentpars.contact_n.x-currentpars.contact_p.x,currentpars.contact_n.y-currentpars.contact_p.y)/2;
	d1=MIN(d1,d2);
	
	dmin=d1;

	/*
		n and m are computed according to:
			_  2 r _ 
		N'  =  | ------ |
		         eps Dx  
		
		see Eq. 23 in documentation. (eps is ERR2 in this routine)
	*/
	n=(int)floor(2*MAX(RN,RP)/(ERR2*(currentpars.corner2.x-currentpars.corner1.x))+1.0);
	m=(int)floor(2*MAX(RN,RP)/(ERR2*(currentpars.corner2.y-currentpars.corner1.y))+1.0);


	Print(VERBOSE,"Computing source/sink positional offsets w.r.t the contact positions\n");
	Print(DEBUG,"Field dimensions in Offset calculation (n,m): (%i,%i)\n", n,m);
	Print(VERBOSE,"Offset p (x,y)\t\t\tOffset n (x,y)\t\t\tIter\n");
	Print(VERBOSE,"---------------------------------------------------------------------\n", n,m);
	do
	{
		currentpars.off_p=off_p;
		currentpars.off_n=off_n; 
		Print(DEBUG,"(%13.6e,%13.6e)\t(%13.6e,%13.6e)\t%i\n", currentpars.off_p.x,currentpars.off_p.y,currentpars.off_n.x,currentpars.off_n.y, iter);
		fflush(stdout);
		ppx=1;
		ppy=1;
		pnx=1;
		pny=1;
		#pragma omp parallel for private(i,j,Cn,Cp,d1,d2,d3,d4)  reduction(*:ppx, ppy, pnx, pny) 
		for(i=-n;i<=n;i++)
		for(j=-m;j<=m;j++)
		{
			if (!(i==0 && j==0))
			{
				/* interact both contacts with both image contacts */
				/* get positions */
				Cp=Positive(i,j);
				Cn=Negative(i,j);
				
				/* get lengths^2*/
				d1=exp(-2*bessk0(rr*len(p0.x-Cp.x-RP, p0.y-Cp.y)));
				d2=exp(-2*bessk0(rr*len(p0.x-Cp.x+RP, p0.y-Cp.y)));
				d3=exp(-2*bessk0(rr*len(p0.x-Cn.x-RP, p0.y-Cn.y)));
				d4=exp(-2*bessk0(rr*len(p0.x-Cn.x+RP, p0.y-Cn.y)));
				
				/* add offsets to new offsetvector */
				ppx*=((d1*d4)/(d2*d3));
				/* get lengths^2 */
				d1=exp(-2*bessk0(rr*len(n0.x-Cp.x-RN, n0.y-Cp.y)));
				d2=exp(-2*bessk0(rr*len(n0.x-Cp.x+RN, n0.y-Cp.y)));
				d3=exp(-2*bessk0(rr*len(n0.x-Cn.x-RN, n0.y-Cn.y)));
				d4=exp(-2*bessk0(rr*len(n0.x-Cn.x+RN, n0.y-Cn.y)));
				
				/* add offsets to new offsetvector */
				pnx*=((d2*d3)/(d1*d4));
				
				
				/* get lengths^2 */
				d1=exp(-2*bessk0(rr*len(p0.x-Cp.x, p0.y-Cp.y-RP)));
				d2=exp(-2*bessk0(rr*len(p0.x-Cp.x, p0.y-Cp.y+RP)));
				d3=exp(-2*bessk0(rr*len(p0.x-Cn.x, p0.y-Cn.y-RP)));
				d4=exp(-2*bessk0(rr*len(p0.x-Cn.x, p0.y-Cn.y+RP)));
				
				/* add offsets to new offsetvector */
				ppy*=((d1*d4)/(d2*d3));
				
				/* get lengths^2 */
				d1=exp(-2*bessk0(rr*len(n0.x-Cp.x, n0.y-Cp.y-RN)));
				d2=exp(-2*bessk0(rr*len(n0.x-Cp.x, n0.y-Cp.y+RN)));
				d3=exp(-2*bessk0(rr*len(n0.x-Cn.x, n0.y-Cn.y-RN)));
				d4=exp(-2*bessk0(rr*len(n0.x-Cn.x, n0.y-Cn.y+RN)));
				
				/* add offsets to new offsetvector */
				pny*=((d2*d3)/(d1*d4));
				
			}
			else
			{
				/* no image contacts, only real ones. i.e. 
				   only interact with the other contact */
				/* get positions */
				Cp=Positive(i,j);
				Cn=Negative(i,j);
				
				/* get lengths^2 */
				d3=exp(-2*bessk0(rr*len(p0.x-Cn.x-RP, p0.y-Cn.y)));
				d4=exp(-2*bessk0(rr*len(p0.x-Cn.x+RP, p0.y-Cn.y)));
				
				/* add offsets to new offsetvector */
				ppx*=(d4/d3);
				
				/* get lengths^2 */
				d1=exp(-2*bessk0(rr*len(n0.x-Cp.x-RN, n0.y-Cp.y)));
				d2=exp(-2*bessk0(rr*len(n0.x-Cp.x+RN, n0.y-Cp.y)));
				
				/* add offsets to new offsetvector */
				pnx*=(d2/d1);		
				
				/* get lengths^2 */
				d3=exp(-2*bessk0(rr*len(p0.x-Cn.x, p0.y-Cn.y-RP)));
				d4=exp(-2*bessk0(rr*len(p0.x-Cn.x, p0.y-Cn.y+RP)));
				
				/* add offsets to new offsetvector */
				ppy*=(d4/d3);
				
				/* get lengths^2 */
				d1=exp(-2*bessk0(rr*len(n0.x-Cp.x,n0.y-Cp.y-RN)));
				d2=exp(-2*bessk0(rr*len(n0.x-Cp.x,n0.y-Cp.y+RN)));
				
				/* add offsets to new offsetvector */
				pny*=(d2/d1);		
			}
			
		}
		off_p.x=(1-ppx)/(1+ppx);
		off_p.y=(1-ppy)/(1+ppy);		
		ra=veclen2(off_p);
		if (ra==0)
			b=RP*0.5;
		else if(ra<1)		
			b=RP*(1-sqrt(1-ra))/ra;
		else
		{
			/* If this happens our algorithm fails. */
			/* This basically means there is no position */
			/* for the source which would make v1=v2 and v3=v4 */
			Warning("Cannot compute a good source offset. Try smaller contacts\n");
			fallback=1;
			/* just do something */
			b=0;				
		}
		off_p.x*=b;
		off_p.y*=b;
				
		off_n.x=(1-pnx)/(1+pnx);
		off_n.y=(1-pny)/(1+pny);		
		ra=veclen2(off_n);
		if (ra==0)
			b=RN*0.5;
		else if(ra<1)	
			b=RN*(1-sqrt(1-ra))/ra;
		else
		{
			/* Cannot make approximately equipotential contacts */
			Warning("Cannot compute a good sink offset. Try smaller contacts\n");
			fallback=1;				
			b=0;		
		}
		off_n.x*=b;
		off_n.y*=b;
		d1=veclen(off_p);
		d2=veclen(currentpars.off_p);
		d3=fabs(d1-d2)/dmin;
		d1=veclen(off_n);
		d2=veclen(currentpars.off_n);
		d4=fabs(d1-d2)/dmin;		
		iter++;
	} while (d3+d4>ERR && iter<MAX_ITER && (!fallback));
	currentpars.off_p=off_p;
	currentpars.off_n=off_n;
	Print(VERBOSE,"(%13.6e,%13.6e)\t(%13.6e,%13.6e)\t%i\n", currentpars.off_p.x,currentpars.off_p.y,currentpars.off_n.x,currentpars.off_n.y, iter);
	if(iter==MAX_ITER)
		Warning("Warning, maximum iterations reached in Offsets\n");
	Print(VERBOSE,"---------------------------------------------------------------------\n\n");
}
#undef ERR
#undef ERR2
#undef RN
#undef RP

double Bessel_Potential(vec p, int n, int m, double rr, double Vn, double Vp)
{
	/* this routine computes the potential at position p */
	/* this time with modified bessel functions and an additional "junction" resistance */
	double res=0, dn, dp, dc;
	vec Cp, Cn;
	double *RV=NULL;
	int i,j;
	
	
	#pragma omp parallel for private(i,j,Cn,Cp,dc) 
	for (i=-n;i<=n;i++)
		for (j=-m;j<=m;j++)
		{
			Cp=Positive_c(i, j);
			dc=len((p.x-Cp.x),p.y-Cp.y);
			if (dc<currentpars.rp)
				RV=(&Vp);
			Cp=Negative_c(i, j);
			dc=len((p.x-Cp.x),p.y-Cp.y);
			if (dc<currentpars.rn)
				RV=(&Vn);		
		
		}
	if (RV)
		return *RV;
	#pragma omp parallel for private(i,j,Cn,Cp,dp,dn) reduction(+:res) 
	for (i=-n;i<=n;i++)
		for (j=-m;j<=m;j++)
		{
			/* ++ */
			Cp=Positive(i, j);
			Cn=Negative(i, j);
			dp=len((p.x-Cp.x),(p.y-Cp.y));
			dn=len((p.x-Cn.x),(p.y-Cn.y));
			res+=(bessk0(rr*dp)-bessk0(rr*dn));
		}
	return res/TWOPI;
}

vec Bessel_Field(vec p, int n, int m, double rr)
{
	/* this routine computes the electric field at position p */
	vec res={0,0}, Cp, Cn, v0={0,0};
	vec *RV=NULL;
	double dp, dn ,dc, ex=0, ey=0, dkn, dkp;
	int i,j;
	
	#pragma omp parallel for private(i,j,Cn,Cp,dc) 
	for (i=-n;i<=n;i++)
		for (j=-m;j<=m;j++)
		{
			Cp=Positive_c(i, j);
			dc=len((p.x-Cp.x),p.y-Cp.y);
			if (dc<currentpars.rp)
				RV=(&v0);
			Cp=Negative_c(i, j);
			dc=len((p.x-Cp.x),p.y-Cp.y);
			if (dc<currentpars.rn)
				RV=(&v0);		
		
		}
	if (RV)
		return *RV;
		
	#pragma omp parallel for private(i,j,Cn,Cp,dp,dn,dkn,dkp) reduction(+:ex, ey) 
	for (i=-n;i<=n;i++)
		for (j=-m;j<=m;j++)
		{
			/* ++ */
			Cp=Positive(i, j);
			Cn=Negative(i, j);
			dp=len((p.x-Cp.x),p.y-Cp.y);
			dn=len((p.x-Cn.x),p.y-Cn.y);
			
			dkn=rr*bessk1(rr*dn);
			dkp=rr*bessk1(rr*dp);
						
			ex+=(dkn*(p.x-Cn.x)/dn-dkp*(p.x-Cp.x)/dp);
			ey+=(dkn*(p.y-Cn.y)/dn-dkp*(p.y-Cp.y)/dp);
		}
	res.x=ex/TWOPI;
	res.y=ey/TWOPI;
	return res;
}


vec * Bessel_FieldLine(vec s, int n, int m, double rr, double Vn, double Vp, double X1, double X2, double Y1, double Y2, int *len, double DL)
{
	/* follow a fieldline from starting point s (and its potential) */
	/* The starting point (s) should be close to either the source of sink of a source-sink pair. The endpoint then is */
	/* when we get out of the drawing  area or reach the opposite source/sink */
	/* This seems simple enough (and it is) but inherently to the method used here the error accumilates which requires */
	/* small steps, which in turn easily leads to excessive computation times. */
	
	/* to keep it under controll we adapt the step size along the way to the field strength *and the change of */
	/* direction of the field between two succesive points along the field-line*. */ 
	
	/* The parameter DL controlls the step size. */
	/* This algorith sortta worked for me but I make no guarantees for it's accuracy or usefullness. */
	/* You may have to adapt the algorithm to your situation as the routine is only tested for a few cases. */
	   
	vec *res, f1={0,0},f2={0,1};
	
	double psi, dir=1, l, df, Psi_end, d1, d2;
	unsigned int alloc=100;
	
	res=malloc(alloc* sizeof(vec));
	(*len)=0;
	
	psi=Bessel_Potential(s, n, m, rr, Vn, Vp);
	if(psi-Vn>Vp-psi)
		Psi_end=Vn;
	else
		Psi_end=Vp;
	
	if (psi>Psi_end)
		dir=-1.0;
	f1=Bessel_Field(s, n, m, rr);
	l=veclen(f1);
	
	/* smallest distance in d1, to scale DL */
	d1=currentpars.contact_p.x-currentpars.rp-currentpars.corner1.x;
	d2=-currentpars.contact_p.x-currentpars.rp+currentpars.corner2.x;
	d1=MIN(d1,d2);
		
	d2=currentpars.contact_p.y-currentpars.rp-currentpars.corner1.y;
	d1=MIN(d1,d2);
		
	d2=-currentpars.contact_p.y-currentpars.rp+currentpars.corner2.y;
	d1=MIN(d1,d2);
		
	d2=currentpars.contact_n.x-currentpars.rn-currentpars.corner1.x;
	d1=MIN(d1,d2);
		
	d2=-currentpars.contact_n.x-currentpars.rn+currentpars.corner2.x;
	d1=MIN(d1,d2);
		
	d2=currentpars.contact_n.y-currentpars.rn-currentpars.corner1.y;
	d1=MIN(d1,d2);
		
	d2=-currentpars.contact_n.y-currentpars.rn+currentpars.corner2.y;
	d1=MIN(d1,d2);
		
	d2=len(currentpars.contact_n.x-currentpars.contact_p.x,currentpars.contact_n.y-currentpars.contact_p.y)/2;
	d1=MIN(d1,d2);
		

	DL*=d1;
	
	Print(DEBUG,"df\tx\ty\tEx\tEy\tdir\n");	
	while(VecInRange(s, X1, X2, Y1, Y2) && dir*(Psi_end-psi)>0)
	{
		if ((*len)==(int)alloc-1)
		{
			alloc+=100;
			res=realloc(res, alloc*sizeof(vec));
		}
		res[(*len)++]=s;
		
		/* f2 is the normalized version of f1 */
		f2=f1;
		f2.x/=l;
		f2.y/=l;
				
		/* compute the field and field-strength at the current position */
		f1=Bessel_Field(s, n, m, rr);
		l=veclen(f1);
		
		df=l*len(f2.x-f1.x/l,f2.y-f1.y/l);
		
		Print(DEBUG,"%14.6e %14.6e %14.6e %14.6e %14.6e %14.6e %14.6e\n",df,l, s.x, s.y, f1.x, f1.y, dir);
		
		/* empirically the step size DL/((df+1)*l) works pretty well, in my test cases at least */
		s.x+=dir*f1.x*DL/((df+1)*l);
		s.y+=dir*f1.y*DL/((df+1)*l);
			
		psi=Bessel_Potential(s, n, m, rr, Vn, Vp);
	}
	if ((*len)==(int)alloc-1)
	{
		alloc+=100;
		res=realloc(res, alloc*sizeof(vec));
	}
	res[(*len)++]=s;
	res=realloc(res, ((unsigned int)(*len)+1)*sizeof(vec));
	return res;
}
