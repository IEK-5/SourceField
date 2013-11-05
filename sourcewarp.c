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
#include "sourcewarp.h"
#include "sourcefield.h"
#include "scd2.h"
#include "mod_bessel_01.h"

#define len(a,b) (sqrt((a)*(a)+(b)*(b)))
#define len2(a,b) ((a)*(a)+(b)*(b))
#define veclen(a) (sqrt((a.x)*(a.x)+(a.y)*(a.y)))
#define veclen2(a) ((a.x)*(a.x)+(a.y)*(a.y))
#define MIN(a,b) ((a)<(b) ? (a):(b))
#define MAX(a,b) ((a)<(b) ? (b):(a))

int InitWarp (double *scale, double *k, double *K)
{
	double a,d;
	a=(currentpars.corner2.y-currentpars.corner1.y);
	d=(currentpars.corner2.x-currentpars.corner1.x);
	(*scale)=2.0*a/d;
	(*k)=ell_mod(*scale);
	(*K)=ell_int(*k);
	(*scale)=2.0*(*K)/d;
	Print(VERBOSE,"----------------------------------------\n");
	Print(VERBOSE,"Scale: %e\nElliptic modulus: %e\nComplete Elliptic Integral: %e\n", *scale, *k, *K);
	Print(VERBOSE,"----------------------------------------\n");
	return 1;
}

vec TransForm (double scale, double K, double k, vec pos)
{
	vec v;
	
	if(sn_i(scale*pos.x, scale*(pos.y-currentpars.corner1.y), K, k, &v.x,&v.y, NULL, NULL, NULL,NULL))
		Error("Transformation failed in TransForm\n");
	Print(DEBUG,"Transform: (%e,%e) --> (%e,%e)\n", pos.x, pos.y, v.x, v.y);
	return v;
}

double Warp_Potential(double scale, double K, double k, vec p, double Vn, double Vp, int eqcontact)
{
	double res, dc;
	vec C,P, Cn, Cp;
	Cn.x=currentpars.contact_n.x;
	Cn.y=currentpars.contact_n.y;
	Cp.x=currentpars.contact_p.x;
	Cp.y=currentpars.contact_p.y;
	
	if (eqcontact)
	{
		dc=len((p.x-Cp.x),p.y-Cp.y);
		if (dc<currentpars.rp)
			return Vp;
		dc=len((p.x-Cn.x),p.y-Cn.y);
		if (dc<currentpars.rn)
			return Vn;
	}
		
	P=TransForm (scale, K, k, p);
	
	C.x=currentpars.contact_n.x+currentpars.off_n.x;
	C.y=currentpars.contact_n.y+currentpars.off_n.y;
	
	C=TransForm (scale, K, k, C);
	res=len2(P.x-C.x,P.y-C.y);
	
	C.y=-C.y;
	res*=len2(P.x-C.x,P.y-C.y);
	
	C.x=currentpars.contact_p.x+currentpars.off_p.x;
	C.y=currentpars.contact_p.y+currentpars.off_p.y;
	C=TransForm (scale, K, k, C);
	res/=len2(P.x-C.x,P.y-C.y);
		
	C.y=-C.y;
	res/=len2(P.x-C.x,P.y-C.y);
	
	return log(sqrt(res))/TWOPI;	
}

vec TransForm_d (double scale, double K, double k, vec pos, double *dudx, double *dudy, double *dvdx, double *dvdy)
{
	vec v;
	
	if(sn_i(scale*pos.x, scale*(pos.y-currentpars.corner1.y), K, k, &v.x,&v.y, dudx, dudy, dvdx, dvdy))
		Error("Transformation failed in TransForm\n");
	return v;
}

vec Warp_Field(double scale, double K, double k, vec p)
{
	double dudx, dudy, dvdx, dvdy;
	vec C,P, E, res;
	
	P=TransForm_d (scale, K, k, p, &dudx, &dudy, &dvdx, &dvdy);
	
	C.x=currentpars.contact_n.x+currentpars.off_n.x;
	C.y=currentpars.contact_n.y+currentpars.off_n.y;
	C=TransForm (scale, K, k, C);
	
	
	E.x=(P.x-C.x)/len2(P.x-C.x,P.y-C.y);
	E.y=(P.y-C.y)/len2(P.x-C.x,P.y-C.y);
	
	C.y=-C.y;
	E.x+=(P.x-C.x)/len2(P.x-C.x,P.y-C.y);
	E.y+=(P.y-C.y)/len2(P.x-C.x,P.y-C.y);
	
	C.x=currentpars.contact_p.x+currentpars.off_p.x;
	C.y=currentpars.contact_p.y+currentpars.off_p.y;
	
	C=TransForm (scale, K, k, C);
	E.x-=(P.x-C.x)/len2(P.x-C.x,P.y-C.y);
	E.y-=(P.y-C.y)/len2(P.x-C.x,P.y-C.y);
		
	C.y=-C.y;
	E.x-=(P.x-C.x)/len2(P.x-C.x,P.y-C.y);
	E.y-=(P.y-C.y)/len2(P.x-C.x,P.y-C.y);
	E.x/=TWOPI;
	E.y/=TWOPI;
	
	res.x=E.x*dudx*scale+E.y*dvdx*scale;
	res.y=E.x*dudy*scale+E.y*dvdy*scale;
	
	return res;	
}
void Warp_Offsets(double scale, double K, double k)
#define ERR currentpars.Err
#define ERR2 currentpars.Err_Off
#define RP currentpars.rp
#define RN currentpars.rn
#define MAX_ITER currentpars.Off_Iter
{
	
	/* this routine computes offset vectors for the source and sink of the contacts such that the contour of the 
	   circular contact is approximately an equipotential line */ 
	int iter=0;
	int fallback=0;
	vec off_p={0,0}, off_n={0,0};
	vec Cp, Cn;
	vec p1, n1, p2, n2, p3, n3, p4, n4;
	double d1, d2, d3, d4, dmin;
	double ra, b;
	double ppx, ppy, pnx, pny;
	
	/* position of the main source and sink, without offsets */
	currentpars.off_p=off_p;
	currentpars.off_n=off_p;
	p1.x=currentpars.contact_p.x+RP;
	p1.y=currentpars.contact_p.y;
	n1.x=currentpars.contact_n.x+RN;
	n1.y=currentpars.contact_n.y;
	
	p2.x=currentpars.contact_p.x-RP;
	p2.y=currentpars.contact_p.y;
	n2.x=currentpars.contact_n.x-RN;
	n2.y=currentpars.contact_n.y;
	
	p3.x=currentpars.contact_p.x;
	p3.y=currentpars.contact_p.y+RP;
	n3.x=currentpars.contact_n.x;
	n3.y=currentpars.contact_n.y+RN;
	
	p4.x=currentpars.contact_p.x;
	p4.y=currentpars.contact_p.y-RP;
	n4.x=currentpars.contact_n.x;
	n4.y=currentpars.contact_n.y-RN;
	
	p1=TransForm (scale, K, k, p1);
	n1=TransForm (scale, K, k, n1);
	p2=TransForm (scale, K, k, p2);
	n2=TransForm (scale, K, k, n2);
	p3=TransForm (scale, K, k, p3);
	n3=TransForm (scale, K, k, n3);
	p4=TransForm (scale, K, k, p4);
	n4=TransForm (scale, K, k, n4);
	
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


	Print(VERBOSE,"Computing source/sink positional offsets w.r.t the contact positions\n");
	Print(VERBOSE,"Offset p (x,y)\t\t\tOffset n (x,y)\t\t\tIter\n");
	Print(VERBOSE,"---------------------------------------------------------------------\n");
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

		/* no image contacts, only real ones. i.e. 
		   only interact with the other contact */
		/* get positions */
		Cp.x=currentpars.contact_p.x+currentpars.off_p.x;
		Cp.y=currentpars.contact_p.y+currentpars.off_p.y;
		Cp=TransForm (scale, K, k, Cp);
		
		Cn.x=currentpars.contact_n.x+currentpars.off_n.x;
		Cn.y=currentpars.contact_n.y+currentpars.off_n.y;
		Cn=TransForm (scale, K, k, Cn);
			
		/* get lengths^2 */
		d3=len2(p2.x-Cn.x, p2.y-Cn.y);
		d4=len2(p1.x-Cn.x, p1.y-Cn.y);
				
		/* add offsets to new offsetvector */
		ppx*=(d4/d3);
		
		/* get lengths^2 */
		d1=len2(n2.x-Cp.x, n2.y-Cp.y);
		d2=len2(n1.x-Cp.x, n1.y-Cp.y);
		
		/* add offsets to new offsetvector */
		pnx*=(d2/d1);		
		
		/* get lengths^2 */
		d3=len2(p4.x-Cn.x, p4.y-Cn.y);
		d4=len2(p3.x-Cn.x, p3.y-Cn.y);
		
		/* add offsets to new offsetvector */
		ppy*=(d4/d3);
		
		/* get lengths^2 */
		d1=len2(n4.x-Cp.x,n4.y-Cp.y);
		d2=len2(n3.x-Cp.x,n3.y-Cp.y);
		
		/* add offsets to new offsetvector */
		pny*=(d2/d1);
		
		/* interact both contacts with both image contacts */
		/* get positions */
		Cp.x=currentpars.contact_p.x+currentpars.off_p.x;
		Cp.y=-currentpars.contact_p.y-currentpars.off_p.y;
		Cn.x=currentpars.contact_n.x+currentpars.off_n.x;
		Cn.y=-currentpars.contact_n.y-currentpars.off_n.y;
		Cp=TransForm (scale, K, k, Cp);
		Cn=TransForm (scale, K, k, Cn);
				
		/* get lengths^2*/
		d1=len2(p2.x-Cp.x, p2.y-Cp.y);
		d2=len2(p1.x-Cp.x, p1.y-Cp.y);
		d3=len2(p2.x-Cn.x, p2.y-Cn.y);
		d4=len2(p1.x-Cn.x, p1.y-Cn.y);
			
		/* add offsets to new offsetvector */
		ppx*=((d1*d4)/(d2*d3));
		/* get lengths^2 */
		d1=len2(n2.x-Cp.x, n2.y-Cp.y);
		d2=len2(n1.x-Cp.x, n1.y-Cp.y);
		d3=len2(n2.x-Cn.x, n2.y-Cn.y);
		d4=len2(n1.x-Cn.x, n1.y-Cn.y);
						
		/* add offsets to new offsetvector */
		pnx*=((d2*d3)/(d1*d4));
		
		
		/* get lengths^2 */
		d1=len2(p4.x-Cp.x, p4.y-Cp.y);
		d2=len2(p3.x-Cp.x, p3.y-Cp.y);
		d3=len2(p4.x-Cn.x, p4.y-Cn.y);
		d4=len2(p3.x-Cn.x, p3.y-Cn.y);
		
		/* add offsets to new offsetvector */
		ppy*=((d1*d4)/(d2*d3));
		
		/* get lengths^2 */
		d1=len2(n4.x-Cp.x, n4.y-Cp.y);
		d2=len2(n3.x-Cp.x, n3.y-Cp.y);
		d3=len2(n4.x-Cn.x, n4.y-Cn.y);
		d4=len2(n3.x-Cn.x, n3.y-Cn.y);
		
		/* add offsets to new offsetvector */
		pny*=((d2*d3)/(d1*d4));
		
				
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
			/* just *do* something */
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
void WarpContactPotential(double scale, double K, double k, double *Vn, double *Vp, double *Emax, double *Erms)
#define NANG currentpars.Nangle
{
	/* this routine computes to what degree the contour of the circular contacts are equipotential lines */  
	/* furthermore it computes the average potential along the contour of the contact and sets Vp and Vn */  
	/* accordingly. (usually the values for Vn and Vp are quite good) */  
	double a;
	vec Cc;
	int kk;
	double *E, Vx=0;
	
	E=malloc((2*NANG+1)*sizeof(double));
	
	(*Emax)=0;
	(*Erms)=0;
	
	for (kk=0;kk<NANG;kk++)
	{
		a=TWOPI*((double)kk)/NANG;
		Cc.x=currentpars.contact_p.x+currentpars.rp*cos(a);
		Cc.y=currentpars.contact_p.y+currentpars.rp*sin(a);
		E[kk]=Warp_Potential(scale, K, k, Cc, 0, 0, 0);
	}
	for (kk=0;kk<NANG;kk++)
	{
		a=TWOPI*((double)kk+0.5)/NANG;
		Cc.x=currentpars.contact_n.x+currentpars.rn*cos(a);
		Cc.y=currentpars.contact_n.y+currentpars.rn*sin(a);
		E[NANG+kk]=Warp_Potential(scale, K, k, Cc, 0, 0, 0);
	}
	Print(VERBOSE,"Vp         	Vn\n");
	Print(VERBOSE,"---------------------------------------------------------------------\n");
	for (kk=0;kk<NANG;kk++)
		Vx+=E[kk];
	Vx/=NANG;
	(*Vp)=Vx;
	Print(VERBOSE,"%-13.3e  ", (*Vp));
	
	Vx=0;
	for (kk=0;kk<NANG;kk++)
		Vx+=E[NANG+kk];
	Vx/=NANG;
	(*Vn)=Vx;
	Print(VERBOSE,"%-13.3e\n", (*Vn));
	Print(VERBOSE,"---------------------------------------------------------------------\n\n");
	
	Vx=(*Vp);
	for (kk=0;kk<2*NANG;kk++)
	{
		if (kk>=NANG)
			Vx=(*Vn);
			
		E[kk]=E[kk]/Vx-1;
		if(E[kk]*E[kk]>(*Emax))
			(*Emax)=E[kk]*E[kk];
		(*Erms)+=E[kk]*E[kk];
	}
	free(E);
	(*Emax)=sqrt((*Emax));
	(*Erms)=sqrt((*Erms)/(2*NANG));
}
#undef NANG
vec * Warp_FieldLine(vec s, double scale, double K, double k, double Vn, double Vp, double X1, double X2, double Y1, double Y2, int *len, double DL)
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
	
	psi=Warp_Potential(scale, K, k, s, Vn, Vp, 1);
	if(psi-Vn>Vp-psi)
		Psi_end=Vn;
	else
		Psi_end=Vp;
	
	if (psi>Psi_end)
		dir=-1.0;
	f1=Warp_Field(scale, K, k, s);
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
		f1=Warp_Field(scale, K, k, s);
		l=veclen(f1);
		
		df=l*len(f2.x-f1.x/l,f2.y-f1.y/l);
		
		Print(DEBUG,"%14.6e %14.6e %14.6e %14.6e %14.6e %14.6e %14.6e\n",df,l, s.x, s.y, f1.x, f1.y, dir);
		
		/* empirically the step size DL/((df+1)*l) works pretty well, in my test cases at least */
		s.x+=dir*f1.x*DL/((df+1)*l);
		s.y+=dir*f1.y*DL/((df+1)*l);
			
		psi=Warp_Potential(scale, K, k, s, Vn, Vp, 1);
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
