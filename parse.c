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
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "main.h"
#include "utils.h"
#include "sourcefield.h"
#include "mod_bessel_01.h"
#include "sourcewarp.h"
#include "parsedef.h"
#include "parse.h"
#define len(a,b) (sqrt((a)*(a)+(b)*(b)))

#define MAXSTRLEN 256

/* parsing utility functions */
static PRSDEF LookupKey (char *name,  const KeyWord * AKeyTable)
{
      unsigned int len;
      len=strlen(name);
      for (; AKeyTable->name; AKeyTable++)
      {
      	    if(strlen(AKeyTable->name)==len)
	    	if (strncmp (name, AKeyTable->name, len) == 0)
			break;
      }
      return AKeyTable->PAR;
}

static char * Begin (char *line)
{
	/* returns a pointer to the beginning of a word, NULL if end is reached */
	while(isspace((*line)) && *line)
		line++; 
	if (*line)
		return line;
	else
		return NULL;
}
static char * End (char *line)
{
	/* returns a pointer to the end of a word */
	while(!isspace((*line)) && *line)
		line++;
	return line;
}


static int CheckPars (void)
#define RP currentpars.rp
#define RN currentpars.rn
#define CXP currentpars.contact_p.x
#define CYP currentpars.contact_p.y
#define CXN currentpars.contact_n.x
#define CYN currentpars.contact_n.y
#define X1 currentpars.corner1.x
#define Y1 currentpars.corner1.y
#define X2 currentpars.corner2.x
#define Y2 currentpars.corner2.y
{
	int res=0;
	if (X1>=X2)
	{
		Warning("x1 should be smaller than x2\n");
		res++;
	}
	if (Y1>=Y2)
	{
		Warning("y1 should be smaller than y2\n");
		res++;
	}
	if ((CXP-RP-X1<=0)||(CXP+RP-X2>=0))
	{
		Warning("Positive Contact: x-position outside field (including radius)\n");
		res++;
	}
	if ((CYP-RP-Y1<=0)||(CYP+RP-Y2>=0))
	{
		Warning("Positive Contact: y-position outside field (including radius)\n");
		res++;
	}
	if ((CXN-RN-X1<=0)||(CXN+RN-X2>=0))
	{
		Warning("Negative Contact: x-position outside field (including radius)\n");
		res++;
	}
	if ((CYN-RN-Y1<=0)||(CYN+RN-Y2>=0))
	{
		Warning("Negative Contact: y-position outside field (including radius)\n");
		res++;
	} 
	if (len(CXP-CXN,CYP-CYN)<=RN+RP)
	{
		Warning("Contact 1 & 2 are on top of each other.\n");
		res++;
	}	
	if (currentpars.Nangle<3)
	{
		Warning("A less than reasonable number of angles. nangle_check is %i\n", currentpars.Nangle);
		res++;
	}
	if (currentpars.NXM<4)
	{
		Warning("A less than reasonable maximum number of images. NxM is %i\n", currentpars.NXM);
		res++;
	}
	if (currentpars.Off_Iter<4)
	{
		Warning("A less than reasonable number of iterations for offset calculations. off_iter is %i\n", currentpars.Off_Iter);
		res++;
	}
	return res;
}
#undef RP
#undef RN
#undef A
#undef X1
#undef Y1
#undef X2
#undef Y2
#undef CXP
#undef CYP
#undef CXN
#undef CYN

void Parse (char *file)
/* some dirty, hard coded parameters for the fieldline computation function */
#define LLEN 50 /* The length of our progress-bar, best to choose it divisable by tics */
#define TICS 5 /* Number of tics*/
#define RR 1.0000001 /* We start just outside the contact */
{
	FILE *f;
	FILE *fout;
	char *line, c;
	/* output filename will be stored here */
	char *outfile=NULL;
	char *mesh=NULL;
	/* pointers to various locations in the read input line, for parsing purposes */
	char *word;
	char *begin, *end;
	/* Current and voltages */
	double Vn, Vp,I;
	double w_scale, w_K, w_k; /* for the sourcewarp calculations */
	double rr, In, Ip; /* for the bessel potantial computations */ 
	/* sheet resistance and injected current */
	double Rsq=1,Ii=1,rg=1;
	/* various data arrays for temporary storage of computed fields and potentials, awayting to be written to a file*/
	double **data;
	double *Xdata;
	double *Ydata;
	/* Various counters, n,m are the dimentions of the source/sink field */
	int n,m, nb, mb, i,j, k, l;
	/* discretization of computation area and specification of number of angles used for fieldlines */
	int Nx, Ny, Nangle=20;
	/* computation area parameters */
	double x1,x2,y1,y2;
	/* Custom computation range parameters */
	double Px1, Px2,Py1, Py2;
	/* Dim indicates whether the problem was solved or not, pdef indicates whether a custom computation area is specified */
	int warp=0, dim=0, bess=0, pdef=0;
	/* various vector variables to fill with temporary things */
	vec E, p;
	/* the fieldline vector array, a series of coordinates making up one fieldline */
	vec *fl;
	/* lengths of the fieldline vector array */
	int fll;
	/* Errors indicating to what degree are the circular contacts really equipotential. */
	double Emax, Erms;
	/* accuracy in field-line computation, see fieldline routine, parameter DL */
	double fl_acc=5e-3;
	PRSDEF key;
	
	if ((f=fopen(file,"r"))==NULL)
		Error("Cannot open file %s\n", file);
	

	line=malloc(MAXSTRLEN*sizeof(char));
	word=malloc(MAXSTRLEN*sizeof(char));
	
    	fgets(line, MAXSTRLEN-1, f);
	
	Px1=Px2=Py1=Py2=0;
	Ny=Nx=40;
	
	while(feof(f)==0)
	{
	
    		k=sscanf(line, " %c", &c);
		if((k!=-1)&&(c!='#'))
		{
			/* read a word, separate at whitespace, check what the word means then read its value or execute */
			begin=Begin(line);
			while(begin)
			{
				end=End(begin); /* the beginning ends here :) */
				word=strncpy(word,begin,end-begin);
				word[end-begin]='\0';
				key=LookupKey (word,  KeyTable);
				switch(key)
				{
					/* parameters */
					case X1:
						/* read next word and process, if no next word trow an error */
						begin=Begin(end); /* the beginning of the end */
						if(!begin)
							goto premature_end;		
						end=End(begin);	 /* we end the beginning */			
						word=strncpy(word,begin,end-begin);
						word[end-begin]='\0';
						currentpars.corner1.x=atof(word);
						dim=0;						
						break;
					case X2:
						/* read next word and process, if no next word trow an error */
						begin=Begin(end); /* the beginning of the end */
						if(!begin)
							goto premature_end;		
						end=End(begin);
						word=strncpy(word,begin,end-begin);
						word[end-begin]='\0';
						currentpars.corner2.x=atof(word);
						dim=0;					
					
						break;
					case Y1:
						/* read next word and process, if no next word trow an error */
						begin=Begin(end); /* the beginning of the end */
						if(!begin)
							goto premature_end;		
						end=End(begin);
						word=strncpy(word,begin,end-begin);
						word[end-begin]='\0';
						currentpars.corner1.y=atof(word);
						dim=0;					
					
						break;
					case Y2:
						/* read next word and process, if no next word trow an error */
						begin=Begin(end); /* the beginning of the end */
						if(!begin)
							goto premature_end;		
						end=End(begin);
						word=strncpy(word,begin,end-begin);
						word[end-begin]='\0';
						currentpars.corner2.y=atof(word);
						dim=0;					
					
						break;
					case CXP:
						/* read next word and process, if no next word trow an error */
						begin=Begin(end); /* the beginning of the end */
						if(!begin)
							goto premature_end;		
						end=End(begin);				
						word=strncpy(word,begin,end-begin);
						word[end-begin]='\0';
						currentpars.contact_p.x=atof(word);
						dim=0;						
						break;
					case CXN:
						/* read next word and process, if no next word trow an error */
						begin=Begin(end); /* the beginning of the end */
						if(!begin)
							goto premature_end;		
						end=End(begin);
						word=strncpy(word,begin,end-begin);
						word[end-begin]='\0';
						currentpars.contact_n.x=atof(word);
						dim=0;					
					
						break;
					case CYP:
						/* read next word and process, if no next word trow an error */
						begin=Begin(end); /* the beginning of the end */
						if(!begin)
							goto premature_end;		
						end=End(begin);
						word=strncpy(word,begin,end-begin);
						word[end-begin]='\0';
						currentpars.contact_p.y=atof(word);
						dim=0;					
					
						break;
					case CYN:
						/* read next word and process, if no next word trow an error */
						begin=Begin(end); /* the beginning of the end */
						if(!begin)
							goto premature_end;		
						end=End(begin);
						word=strncpy(word,begin,end-begin);
						word[end-begin]='\0';
						currentpars.contact_n.y=atof(word);
						dim=0;					
					
						break;
					case OFFSETS:
						Offsets();	
						dim=0;									
						break;
					case R:
						/* read next word and process, if no next word trow an error */
						begin=Begin(end); /* the beginning of the end */
						if(!begin)
							goto premature_end;		
						end=End(begin);
						word=strncpy(word,begin,end-begin);
						word[end-begin]='\0';
						currentpars.rp=atof(word);
						currentpars.rn=currentpars.rp;
						dim=0;					
					
						break;
					case RP:
						/* read next word and process, if no next word trow an error */
						begin=Begin(end); /* the beginning of the end */
						if(!begin)
							goto premature_end;		
						end=End(begin);
						word=strncpy(word,begin,end-begin);
						word[end-begin]='\0';
						currentpars.rp=atof(word);
						dim=0;					
					
						break;
					case RN:
						/* read next word and process, if no next word trow an error */
						begin=Begin(end); /* the beginning of the end */
						if(!begin)
							goto premature_end;		
						end=End(begin);
						word=strncpy(word,begin,end-begin);
						word[end-begin]='\0';
						currentpars.rn=atof(word);
						dim=0;					
					
						break;
					case RSQ:					
						/* read next word and process, if no next word trow an error */
						begin=Begin(end); /* the beginning of the end */
						if(!begin)
							goto premature_end;		
						end=End(begin);
						word=strncpy(word,begin,end-begin);
						word[end-begin]='\0';
						Rsq=atof(word);	
						break;
					case RG:					
						/* read next word and process, if no next word trow an error */
						begin=Begin(end); /* the beginning of the end */
						if(!begin)
							goto premature_end;		
						end=End(begin);
						word=strncpy(word,begin,end-begin);
						word[end-begin]='\0';
						rg=atof(word);	
						break;
					case CURRENT:					
						/* read next word and process, if no next word trow an error */
						begin=Begin(end); /* the beginning of the end */
						if(!begin)
							goto premature_end;		
						end=End(begin);
						word=strncpy(word,begin,end-begin);
						word[end-begin]='\0';
						Ii=atof(word);	
						break;
					case PX1:
						/* read next word and process, if no next word trow an error */
						begin=Begin(end); /* the beginning of the end */
						if(!begin)
							goto premature_end;		
						end=End(begin);
						word=strncpy(word,begin,end-begin);
						word[end-begin]='\0';
						Px1=atof(word);	
						pdef=1;	
					
						break;
					case PX2:
						/* read next word and process, if no next word trow an error */
						begin=Begin(end); /* the beginning of the end */
						if(!begin)
							goto premature_end;		
						end=End(begin);
						word=strncpy(word,begin,end-begin);
						word[end-begin]='\0';
						Px2=atof(word);		
					
						pdef=1;	
						break;
					case PY1:
						/* read next word and process, if no next word trow an error */
						begin=Begin(end); /* the beginning of the end */
						if(!begin)
							goto premature_end;		
						end=End(begin);
						word=strncpy(word,begin,end-begin);
						word[end-begin]='\0';
						Py1=atof(word);	
						pdef=1;		
					
						break;
					case PY2:
						/* read next word and process, if no next word trow an error */
						begin=Begin(end); /* the beginning of the end */
						if(!begin)
							goto premature_end;		
						end=End(begin);
						word=strncpy(word,begin,end-begin);
						word[end-begin]='\0';
						Py2=atof(word);	
						pdef=1;		
					
						break;					
					case NY:
						/* read next word and process, if no next word trow an error */
						begin=Begin(end); /* the beginning of the end */
						if(!begin)
							goto premature_end;		
						end=End(begin);	
						word=strncpy(word,begin,end-begin);
						word[end-begin]='\0';
						Ny=atoi(word);		
					
						break;
					case NX:
						/* read next word and process, if no next word trow an error */
						begin=Begin(end); /* the beginning of the end */
						if(!begin)
							goto premature_end;		
						end=End(begin);
						word=strncpy(word,begin,end-begin);
						word[end-begin]='\0';
						Nx=atoi(word);			
					
						break;
					case NANGLE:
						/* read next word and process, if no next word trow an error */
						begin=Begin(end); /* the beginning of the end */
						if(!begin)
							goto premature_end;		
						end=End(begin);
						word=strncpy(word,begin,end-begin);
						word[end-begin]='\0';
						Nangle=atoi(word);
						break;
					case FL_ACC:
						/* read next word and process, if no next word trow an error */
						begin=Begin(end); /* the beginning of the end */
						if(!begin)
							goto premature_end;		
						end=End(begin);
						word=strncpy(word,begin,end-begin);
						word[end-begin]='\0';
						fl_acc=atof(word);
						break;
					case _QUIET:
						if(!fixverb)
							verbose=QUIET;
						break;
					case _NORMAL:
						if(!fixverb)
							verbose=NORMAL;
						break;
					case _VERBOSE:
						if(!fixverb)
							verbose=VERBOSE;
						break;
					case _DEBUG:
						if(!fixverb)
							verbose=DEBUG;
						break;
					case OUTFILE:
						/* read next word and process, if no next word trow an error */
						begin=Begin(end); /* the beginning of the end */
						if(!begin)
							goto premature_end;		
						end=End(begin);
						if (!outfile)
							outfile=malloc(MAXSTRLEN*sizeof(char));
						outfile=strncpy(outfile,begin,end-begin);
						outfile[end-begin]='\0';					
						break;
					case MESH:
						/* read next word and process, if no next word trow an error */
						begin=Begin(end); /* the beginning of the end */
						if(!begin)
							goto premature_end;		
						end=End(begin);
						if (!mesh)
							mesh=malloc(MAXSTRLEN*sizeof(char));
						mesh=strncpy(mesh,begin,end-begin);
						mesh[end-begin]='\0';
						pdef=2;						
						break;
					case ERR:
						/* read next word and process, if no next word trow an error */
						begin=Begin(end); /* the beginning of the end */
						if(!begin)
							goto premature_end;		
						end=End(begin);
						word=strncpy(word,begin,end-begin);
						word[end-begin]='\0';
						currentpars.Err=atof(word);	
						dim=0;				
					
						break;
					case ERR_OFF:
						/* read next word and process, if no next word trow an error */
						begin=Begin(end); /* the beginning of the end */
						if(!begin)
							goto premature_end;		
						end=End(begin);
						word=strncpy(word,begin,end-begin);
						word[end-begin]='\0';
						currentpars.Err_Off=atof(word);	
						dim=0;				
					
						break;
					case NANGLE_CHECK:
						/* read next word and process, if no next word trow an error */
						begin=Begin(end); /* the beginning of the end */
						if(!begin)
							goto premature_end;		
						end=End(begin);
						word=strncpy(word,begin,end-begin);
						word[end-begin]='\0';
						currentpars.Nangle=atoi(word);	
						dim=0;				
					
						break;
					case OFF_ITER:
						/* read next word and process, if no next word trow an error */
						begin=Begin(end); /* the beginning of the end */
						if(!begin)
							goto premature_end;		
						end=End(begin);
						word=strncpy(word,begin,end-begin);
						word[end-begin]='\0';
						currentpars.Off_Iter=atoi(word);	
						dim=0;						
						break;
					case NXM:
						/* read next word and process, if no next word trow an error */
						begin=Begin(end); /* the beginning of the end */
						if(!begin)
							goto premature_end;		
						end=End(begin);
						word=strncpy(word,begin,end-begin);
						word[end-begin]='\0';
						currentpars.NXM=atoi(word);	
						dim=0;						
						break;
					/* execute calculation */
					case RESIST:
						if(!dim) /* we need to solve it all */
						{
							if (CheckPars())
								break;
							Dim(&Vn, &Vp, &I, &n, &m, &Emax, &Erms);
							dim=1;	
							bess=0;		
							warp=0;			
						}	
						Print(QUIET, "%.12g",Rsq*(Vp-Vn)/I);
						Print(NORMAL, " Ohm\nMax EQP Error %6.3f %%\tRMS EQP Error %6.3f %%", 100*Emax, 100*Erms);
						Print(QUIET, "\n");
						break;
					case TWOPOINT:
						if(!dim) /* we need to solve it all */
						{
							if (CheckPars())
								break;
							Dim(&Vn, &Vp, &I, &n, &m, &Emax, &Erms);
							dim=1;	
							bess=0;		
							warp=0;				
						}	
						Print(QUIET, "%.12g",Ii*Rsq*(Vp-Vn)/I);
						Print(NORMAL, " V\nMax EQP Error %6.3f %%\tRMS EQP Error %6.3f %%", 100*Emax, 100*Erms);
						Print(QUIET, "\n");
						break;
					case FOURPOINT:
					{
						double v1, v2;
						if(!dim) /* we need to solve it all */
						{
							if (CheckPars())
								break;
							Dim(&Vn, &Vp, &I, &n, &m, &Emax, &Erms);
							dim=1;	
							bess=0;	
							warp=0;			
						}
						p.x=Px1;
						p.y=Py1;
						v1=Potential(p,n,m, Vn, Vp);
						p.x=Px2;
						p.y=Py2;
						v2=Potential(p,n,m, Vn, Vp);
						Print(QUIET, "%.12g",Ii*Rsq*(v2-v1)/I);
						Print(NORMAL, " V");
						Print(QUIET, "\n");
						break;
					}
					case FIELD:
					{
						int Nm;
						if(!outfile)
							Error("No output file specified\n");
						if ((fout=fopen(outfile,"w"))==NULL)
							Error("Cannot open %s for saving data\n", outfile);
						
						if(!dim) /* we need to solve it all */
						{
							if (CheckPars())
								break;
							Dim(&Vn, &Vp, &I, &n, &m, &Emax, &Erms);
							Print(NORMAL, "Max EQP Error %6.3f %%\tRMS EQP Error %6.3f %%\n\n", 100*Emax, 100*Erms);
							dim=1;
							bess=0;	
							warp=0;	
						}
						
						/* generate the computation mesh */						
						if(!pdef)
							Nm=MakeMesh(&Xdata, &Ydata, currentpars.corner1.x, currentpars.corner2.x, currentpars.corner1.y, currentpars.corner2.y, Nx, Ny);
						else if (!mesh)
							Nm=MakeMesh(&Xdata, &Ydata, Px1, Px2, Py1, Py2, Nx, Ny);
						else
						{
							/* read the mesh */
							Nm=ReadMesh(mesh, &Xdata, &Ydata);
						}								
						data=malloc(3*sizeof(double *));
						data[0]=malloc((Nm+2)*sizeof(double));
						data[1]=malloc((Nm+2)*sizeof(double));
						
						for (i=0;i<=Nm;i++)
						{
							p.x=Xdata[i];
							p.y=Ydata[i];
							E=Field(p,n,m);
							data[0][i]=Ii*Rsq*E.x;	
							data[1][i]=Ii*Rsq*E.y;
						}
						SaveMesh(fout, data, Xdata, Ydata, Nm, 2);
						fclose(fout);
						
						free(data[0]);
						free(data[1]);						
						free(data);
						free(Xdata);
						free(Ydata);
						
						break;
					}
					case POTENTIAL:
					{
						int Nm;
						if(!outfile)
							Error("No output file specified\n");
						
						if ((fout=fopen(outfile,"w"))==NULL)						
							Error("Cannot open %s for saving data\n", outfile);
						
						if(!dim) /* we need to solve it all */
						{
							if (CheckPars())
								break;
							Dim(&Vn, &Vp, &I, &n, &m, &Emax, &Erms);
							Print(NORMAL, "Max EQP Error %6.3f %%\tRMS EQP Error %6.3f %%\n\n", 100*Emax, 100*Erms);
							dim=1;	
							bess=0;
							warp=0;
						}		
						
						/* generate the computation mesh */						
						if(!pdef)
							Nm=MakeMesh(&Xdata, &Ydata, currentpars.corner1.x, currentpars.corner2.x, currentpars.corner1.y, currentpars.corner2.y, Nx, Ny);
						else if (!mesh)
							Nm=MakeMesh(&Xdata, &Ydata, Px1, Px2, Py1, Py2, Nx, Ny);
						else
						{
							/* read the mesh */
							Nm=ReadMesh(mesh, &Xdata, &Ydata);
						}		
							
						data=malloc(2*sizeof(double *));
						data[0]=malloc((Nm+2)*sizeof(double));
						
						for (i=0;i<=Nm;i++)
						{
							p.x=Xdata[i];
							p.y=Ydata[i];
							data[0][i]=Ii*Rsq*Potential(p,n,m, Vn, Vp);
						}
						SaveMesh(fout, data, Xdata, Ydata, Nm, 1);
						fclose(fout);
						
						free(data[0]);					
						free(data);
						free(Xdata);
						free(Ydata);
						
						break;
					}
					case FIELDLINE:
					{	
						int pc=-100/LLEN, TL;					
						if(!outfile)
							Error("No output file specified\n");
							
						if ((fout=fopen(outfile,"w"))==NULL)
							Error("Cannot open %s for saving data\n", outfile);
						
						if(!dim) /* we need to solve it all */
						{
							if (CheckPars())
								break;	
							Dim(&Vn, &Vp, &I, &n, &m, &Emax, &Erms);
							Print(NORMAL, "Max EQP Error %6.3f %%\tRMS EQP Error %6.3f %%\n\n", 100*Emax, 100*Erms);
							dim=1;	
							bess=0;
							warp=0;
						}		
							
						if(!pdef || pdef==2)
						{
							x1=currentpars.corner1.x;
							x2=currentpars.corner2.x;
							y1=currentpars.corner1.y;
							y2=currentpars.corner2.y;	
						}
						else
						{
							x1=Px1;
							x2=Px2;
							y1=Py1;
							y2=Py2;	
						
						}
						/* compute the total amount of work */
						Print(NORMAL,"Computing fieldlines from:\n");
						TL=0;
						for (i=-n;i<=n;i++)
							for (j=-m;j<=m;j++)
							{
								p=Positive_c(i,j);
								if(VecInRange(p,x1,x2,y1,y2))
								{
									TL+=Nangle;
									Print(NORMAL," Source to Sink (%2i %2i)\n", i, j);
								}
								else if (VecInRange(p=Negative_c(i,j),x1,x2,y1,y2))
								{
									TL+=Nangle;
									Print(NORMAL," Sink to Source (%2i %2i)\n", i, j);
								}
							}
						
								
						Print(NORMAL,"Computing fieldlines. This may take some time!\n");
						Print(NORMAL,"For your entertainment I have an utterly useless\n");
						Print(NORMAL,"progress-bar, creeping slowly too 100%%. Enjoy!\n");
						l=0;
						while (l<LLEN)
						{							
							l++;
							if (l==1)
							{
								Print(NORMAL,"0 %%");
								l+=2;
							}
							else if (l==LLEN)
							{
								Print(NORMAL,"100 %%");
								l+=4;
							}
							else
								Print(NORMAL,"-");	
						}
						Print(NORMAL,"\n");		
						l=0;
						for (i=-n;i<=n;i++)
							for (j=-m;j<=m;j++)
							{
								p=Positive_c(i,j);
								if(VecInRange(p,x1,x2,y1,y2))
								{						
									fflush(stdout);
									for (k=0;k<Nangle;k++)
									{
										E.x=p.x+RR*currentpars.rp*cos(TWOPI*((double)k+0.5)/(double)Nangle);
										E.y=p.y+RR*currentpars.rp*sin(TWOPI*((double)k+0.5)/(double)Nangle);
										fl=FieldLine(E, n, m, Vn, Vp, x1, x2, y1, y2, &fll, fl_acc);
										SaveVecs(fout, fl, fll);
										free(fl);
										l++;
										pc=ProgressBar((100*l)/TL, pc, LLEN, TICS);
											
									}
								}
								else if (VecInRange(p=Negative_c(i,j),x1,x2,y1,y2))
								{
									fflush(stdout);
									for (k=0;k<Nangle;k++)
									{
										fflush(stdout);
										E.x=p.x+RR*currentpars.rn*cos(TWOPI*((double)k+0.5)/(double)Nangle);
										E.y=p.y+RR*currentpars.rn*sin(TWOPI*((double)k+0.5)/(double)Nangle);
										fl=FieldLine(E, n, m, Vn, Vp, x1, x2, y1, y2, &fll, fl_acc);
										SaveVecs(fout, fl, fll);
										free(fl);
										l++;
										pc=ProgressBar((100*l)/TL, pc, LLEN, TICS);
									}
								
								}
							
							}
									
						break;
					}
					case WARP_OFFSETS:
						if(!warp) /* we need to solve it all */
						{
							if (CheckPars())
								break;	
							warp=InitWarp(&w_scale, &w_k, &w_K);						
							WarpContactPotential(w_scale, w_K, w_k, &Vn, &Vp, &Emax, &Erms);
							dim=0;	
							bess=0;
							warp=1;				
						}
						Warp_Offsets(w_scale, w_K, w_k);				
						break;
					case WARP_RESIST:
						if(!warp) /* we need to solve it all */
						{
							if (CheckPars())
								break;	
							warp=InitWarp(&w_scale, &w_k, &w_K);							
							WarpContactPotential(w_scale, w_K, w_k, &Vn, &Vp, &Emax, &Erms);
							dim=0;	
							bess=0;
							warp=1;				
						}
						Print(QUIET, "%.12g",Rsq*(Vp-Vn));
						Print(NORMAL, " Ohm\nMax EQP Error %6.3f %%\tRMS EQP Error %6.3f %%", 100*Emax, 100*Erms);
						Print(QUIET, "\n");
						break;
					case WARP_TWOPOINT:
						if(!warp) /* we need to solve it all */
						{
							if (CheckPars())
								break;	
							warp=InitWarp(&w_scale, &w_k, &w_K);							
							WarpContactPotential(w_scale, w_K, w_k, &Vn, &Vp, &Emax, &Erms);
							dim=0;	
							bess=0;
							warp=1;				
						}
						Print(QUIET, "%.12g",Ii*Rsq*(Vp-Vn));
						Print(NORMAL, " V\nMax EQP Error %6.3f %%\tRMS EQP Error %6.3f %%", 100*Emax, 100*Erms);
						Print(QUIET, "\n");
						break;
					case WARP_FOURPOINT:
					{
						double v1, v2;
						if(!warp) /* we need to solve it all */
						{
							if (CheckPars())
								break;	
							warp=InitWarp(&w_scale, &w_k, &w_K);							
							WarpContactPotential(w_scale, w_K, w_k, &Vn, &Vp, &Emax, &Erms);
							Print(NORMAL, "Max EQP Error %6.3f %%\tRMS EQP Error %6.3f %%\n\n", 100*Emax, 100*Erms);
							dim=0;	
							bess=0;
							warp=1;				
						}
						p.x=Px1;
						p.y=Py1;
						v1=Warp_Potential(w_scale, w_K, w_k, p, Vn, Vp,1);
						p.x=Px2;
						p.y=Py2;
						v2=Warp_Potential(w_scale, w_K, w_k, p, Vn, Vp,1);
						Print(QUIET, "%.12g",Ii*Rsq*(v2-v1));
						Print(NORMAL, " V");
						Print(QUIET, "\n");
						break;
					}
					case WARP_FIELD:
					{
						int Nm;
						if(!outfile)
							Error("No output file specified\n");
						if ((fout=fopen(outfile,"w"))==NULL)
							Error("Cannot open %s for saving data\n", outfile);
						
						if(!warp) /* we need to solve it all */
						{
							if (CheckPars())
								break;	
							warp=InitWarp(&w_scale, &w_k, &w_K);							
							WarpContactPotential(w_scale, w_K, w_k, &Vn, &Vp, &Emax, &Erms);
							Print(NORMAL, "Max EQP Error %6.3f %%\tRMS EQP Error %6.3f %%\n\n", 100*Emax, 100*Erms);
							dim=0;	
							bess=0;
							warp=1;				
						}
						
						/* generate the computation mesh */						
						if(!pdef)
							Nm=MakeMesh(&Xdata, &Ydata, currentpars.corner1.x, currentpars.corner2.x, currentpars.corner1.y, currentpars.corner2.y, Nx, Ny);
						else if (!mesh)
							Nm=MakeMesh(&Xdata, &Ydata, Px1, Px2, Py1, Py2, Nx, Ny);
						else
						{
							/* read the mesh */
							Nm=ReadMesh(mesh, &Xdata, &Ydata);
						}								
						data=malloc(3*sizeof(double *));
						data[0]=malloc((Nm+2)*sizeof(double));
						data[1]=malloc((Nm+2)*sizeof(double));
						
						for (i=0;i<=Nm;i++)
						{
							p.x=Xdata[i];
							p.y=Ydata[i];
							E=Warp_Field(w_scale, w_K, w_k, p);
							data[0][i]=Ii*Rsq*E.x;	
							data[1][i]=Ii*Rsq*E.y;
						}
						SaveMesh(fout, data, Xdata, Ydata, Nm, 2);
						fclose(fout);
						
						free(data[0]);
						free(data[1]);						
						free(data);
						free(Xdata);
						free(Ydata);
						
						break;
					}
					case WARP_POTENTIAL:
					{
						int Nm;
						if(!outfile)
							Error("No output file specified\n");
						
						if ((fout=fopen(outfile,"w"))==NULL)						
							Error("Cannot open %s for saving data\n", outfile);
						
						if(!warp) /* we need to solve it all */
						{
							if (CheckPars())
								break;	
							warp=InitWarp(&w_scale, &w_k, &w_K);							
							WarpContactPotential(w_scale, w_K, w_k, &Vn, &Vp, &Emax, &Erms);
							Print(NORMAL, "Max EQP Error %6.3f %%\tRMS EQP Error %6.3f %%\n\n", 100*Emax, 100*Erms);
							dim=0;	
							bess=0;
							warp=1;				
						}	
						
						/* generate the computation mesh */						
						if(!pdef)
							Nm=MakeMesh(&Xdata, &Ydata, currentpars.corner1.x, currentpars.corner2.x, currentpars.corner1.y, currentpars.corner2.y, Nx, Ny);
						else if (!mesh)
							Nm=MakeMesh(&Xdata, &Ydata, Px1, Px2, Py1, Py2, Nx, Ny);
						else
						{
							/* read the mesh */
							Nm=ReadMesh(mesh, &Xdata, &Ydata);
						}		
							
						data=malloc(2*sizeof(double *));
						data[0]=malloc((Nm+2)*sizeof(double));
						
						for (i=0;i<=Nm;i++)
						{
							p.x=Xdata[i];
							p.y=Ydata[i];
							data[0][i]=Ii*Rsq*Warp_Potential(w_scale, w_K, w_k, p, Vn, Vp, 1);
						}
						SaveMesh(fout, data, Xdata, Ydata, Nm, 1);
						fclose(fout);
						
						free(data[0]);					
						free(data);
						free(Xdata);
						free(Ydata);
						
						break;
					}
					case WARP_FIELDLINE:
					{	
						int pc=-100/LLEN, TL;					
						if(!outfile)
							Error("No output file specified\n");
							
						if ((fout=fopen(outfile,"w"))==NULL)
							Error("Cannot open %s for saving data\n", outfile);
						
						if(!warp) /* we need to solve it all */
						{
							if (CheckPars())
								break;	
							warp=InitWarp(&w_scale, &w_k, &w_K);							
							WarpContactPotential(w_scale, w_K, w_k, &Vn, &Vp, &Emax, &Erms);
							Print(NORMAL, "Max EQP Error %6.3f %%\tRMS EQP Error %6.3f %%\n\n", 100*Emax, 100*Erms);
							dim=0;	
							bess=0;
							warp=1;				
						}
							
						if(!pdef || pdef==2)
						{
							x1=currentpars.corner1.x;
							x2=currentpars.corner2.x;
							y1=currentpars.corner1.y;
							y2=currentpars.corner2.y;	
						}
						else
						{
							x1=Px1;
							x2=Px2;
							y1=Py1;
							y2=Py2;	
						
						}
						/* with the warp version we conly compute field lines within the rectangle */
						Print(NORMAL,"Computing fieldlines from:\n");
						TL=0;
						if (VecInRange(currentpars.contact_p,x1,x2,y1,y2))
						{
							TL+=Nangle;
							Print(NORMAL," Source to Sink\n");
						}
						else if (VecInRange(currentpars.contact_n,x1,x2,y1,y2))
						{
							TL+=Nangle;
							Print(NORMAL," Sink to Source\n");
						}
								
						Print(NORMAL,"Computing fieldlines.\n");
						l=0;
						while (l<LLEN)
						{							
							l++;
							if (l==1)
							{
								Print(NORMAL,"0 %%");
								l+=2;
							}
							else if (l==LLEN)
							{
								Print(NORMAL,"100 %%");
								l+=4;
							}
							else
								Print(NORMAL,"-");	
						}
						Print(NORMAL,"\n");		
						l=0;
						if (VecInRange(currentpars.contact_p,x1,x2,y1,y2))
						{						
							fflush(stdout);
							for (k=0;k<Nangle;k++)
							{
								p=currentpars.contact_p;
								E.x=p.x+RR*currentpars.rp*cos(TWOPI*((double)k+0.5)/(double)Nangle);
								E.y=p.y+RR*currentpars.rp*sin(TWOPI*((double)k+0.5)/(double)Nangle);
								fl=Warp_FieldLine(E,w_scale, w_K, w_k,  Vn, Vp, x1, x2, y1, y2, &fll, fl_acc);
								SaveVecs(fout, fl, fll);
								free(fl);
								l++;
								pc=ProgressBar((100*l)/TL, pc, LLEN, TICS);
									
							}
						}
						else if (VecInRange(currentpars.contact_n,x1,x2,y1,y2))
						{
							fflush(stdout);
							for (k=0;k<Nangle;k++)
							{
								fflush(stdout);
								p=currentpars.contact_n;
								E.x=p.x+RR*currentpars.rn*cos(TWOPI*((double)k+0.5)/(double)Nangle);
								E.y=p.y+RR*currentpars.rn*sin(TWOPI*((double)k+0.5)/(double)Nangle);
								fl=Warp_FieldLine(E,w_scale, w_K, w_k,  Vn, Vp, x1, x2, y1, y2, &fll, fl_acc);
								SaveVecs(fout, fl, fll);
								free(fl);
								l++;
								pc=ProgressBar((100*l)/TL, pc, LLEN, TICS);
							}
						}
									
						break;
					}
					case BESS_FAST:
					{
						double dx, dy, sum_b, sum, bE, bd;
						if(!dim) /* we need to solve it all */
						{
							if (CheckPars())
								break;
							Dim(&Vn, &Vp, &I, &n, &m, &Emax, &Erms);
							dim=1;	
							bess=0;
							warp=0;
						}
						
						dx=(currentpars.corner2.x-currentpars.corner1.x);
						dy=(currentpars.corner2.y-currentpars.corner1.y);
						rr=sqrt(Rsq/rg);
						bE=1/currentpars.Err;
						
						sum=0;
						for (k=1;k<n;k++)
						{
							sum+=bessk0(rr*(double)(k+1)*dx)-bessk0(rr*(double)k*dx);
						}
						sum_b=0;
						bE=1/currentpars.Err;
						nb=1;
						while ((sum_b/sum < bE)&&(nb<n))
						{
							bd=bessk0(rr*(double)(nb+1)*dx)-bessk0(rr*(double)nb*dx);
							sum_b+=bd;
							sum-=bd;
							nb++;
						}
						
						sum=0;
						for (k=1;k<m;k++)
						{
							sum+=bessk0(rr*(double)(k+1)*dy)-bessk0(rr*(double)k*dy);
						}
						sum_b=0;
						mb=1;
						while ((sum_b/sum < bE)&&(mb<m))
						{
							bd=bessk0(rr*(double)(mb+1)*dy)-bessk0(rr*(double)mb*dy);
							sum_b+=bd;
							sum-=bd;
							mb++;
						}
						
						if(!bess) /* we need to solve it all */
						{
							if (CheckPars())
								break;
								
							InitBessel(nb, mb, rg, Rsq, &rr, &In, &Ip,  &Vn, &Vp, &Emax, &Erms);
							Print(NORMAL, "Vn %.12g V\tVp %.12g V\n", Vn, Vp);
							Print(NORMAL, "In %.12g A\tIp %.12g A\n", In*Ii, Ip*Ii);
							Print(NORMAL, "Max EQP Error %6.3f %%\tRMS EQP Error %6.3f %%\n\n", 100*Emax, 100*Erms);
							bess=1;
						}
						Print(NORMAL, "Fast Bessel enabled\n(n,m) was (%i,%i), now (%i,%i)\n", n, m, nb, mb);	
						break;
					}
					case BESS_OFFSETS:
						if(!dim) /* we need to solve it all */
						{
							if (CheckPars())
								break;
							Dim(&Vn, &Vp, &I, &n, &m, &Emax, &Erms);
							dim=1;	
							bess=0;	
							warp=0;
						}
						Bessel_Offsets(sqrt(Rsq/rg));
						bess=0;						
						break;
					case BESS_RESIST:
						if(!dim) /* we need to solve it all */
						{
							if (CheckPars())
								break;
							Dim(&Vn, &Vp, &I, &n, &m, &Emax, &Erms);
							dim=1;	
							bess=0;
							warp=0;
						}
						if(!bess) /* we need to solve it all */
						{
							if (CheckPars())
								break;
							InitBessel(n, m, rg, Rsq, &rr, &In, &Ip,  &Vn, &Vp, &Emax, &Erms);
							Print(NORMAL, "Vn %.12g V\tVp %.12g V\n", Vn, Vp);
							Print(NORMAL, "In %.12g A\tIp %.12g A\n", In*Ii, Ip*Ii);
							Print(NORMAL, "Max EQP Error %6.3f %%\tRMS EQP Error %6.3f %%\n\n", 100*Emax, 100*Erms);
							bess=1;
							nb=n;
							mb=m;
						}
						Print(NORMAL, "Rn ");
						Print(QUIET, "%.12g",Rsq*fabs(Vn/(In*Ii)));
						Print(NORMAL, " Ohm");
						Print(QUIET, "\n");
						Print(NORMAL, "Rp ");
						Print(QUIET, "%.12g",Rsq*fabs(Vp/(Ip*Ii)));
						Print(NORMAL, " Ohm");
						Print(QUIET, "\n");
						break;
					case BESS_FIELD:
					{
						int Nm;
						if(!outfile)
							Error("No output file specified\n");
						if ((fout=fopen(outfile,"w"))==NULL)
							Error("Cannot open %s for saving data\n", outfile);
						
						if(!dim) /* we need to solve it all */
						{
							if (CheckPars())
								break;
							Dim(&Vn, &Vp, &I, &n, &m, &Emax, &Erms);
							dim=1;	
							bess=0;
							warp=0;
						}
						if(!bess) /* we need to solve it all */
						{
							if (CheckPars())
								break;
							InitBessel(n, m, rg, Rsq, &rr, &In, &Ip,  &Vn, &Vp, &Emax, &Erms);
							Print(NORMAL, "Vn %.12g V\tVp %.12g V\n", Vn, Vp);
							Print(NORMAL, "In %.12g A\tIp %.12g A\n", In*Ii, Ip*Ii);
							Print(NORMAL, "Max EQP Error %6.3f %%\tRMS EQP Error %6.3f %%\n\n", 100*Emax, 100*Erms);
							bess=1;
							nb=n;
							mb=m;
						}
						
						/* generate the computation mesh */						
						if(!pdef)
							Nm=MakeMesh(&Xdata, &Ydata, currentpars.corner1.x, currentpars.corner2.x, currentpars.corner1.y, currentpars.corner2.y, Nx, Ny);
						else if (!mesh)
							Nm=MakeMesh(&Xdata, &Ydata, Px1, Px2, Py1, Py2, Nx, Ny);
						else
						{
							/* read the mesh */
							Nm=ReadMesh(mesh, &Xdata, &Ydata);
						}		
							
						data=malloc(3*sizeof(double *));
						data[0]=malloc((Nm+2)*sizeof(double));
						data[1]=malloc((Nm+2)*sizeof(double));
						
						for (i=0;i<=Nm;i++)
						{
							p.x=Xdata[i];
							p.y=Ydata[i];
							E=Bessel_Field(p,nb,mb,rr);
							data[0][i]=Ii*Rsq*E.x;	
							data[1][i]=Ii*Rsq*E.y;
						}
						SaveMesh(fout, data, Xdata, Ydata, Nm, 2);
						fclose(fout);
						
						free(data[0]);
						free(data[1]);						
						free(data);
						free(Xdata);
						free(Ydata);
						
						break;
					}
					case BESS_POTENTIAL:
					{
						int Nm;
						if(!outfile)
							Error("No output file specified\n");
						
						if ((fout=fopen(outfile,"w"))==NULL)						
							Error("Cannot open %s for saving data\n", outfile);
						
						if(!dim) /* we need to solve it all */
						{
							if (CheckPars())
								break;
							Dim(&Vn, &Vp, &I, &n, &m, &Emax, &Erms);
							dim=1;	
							bess=0;
							warp=0;
						}
						if(!bess) /* we need to solve it all */
						{
							if (CheckPars())
								break;
							InitBessel(n, m, rg, Rsq, &rr, &In, &Ip,  &Vn, &Vp, &Emax, &Erms);
							Print(NORMAL, "Vn %.12g V\tVp %.12g V\n", Vn, Vp);
							Print(NORMAL, "In %.12g A\tIp %.12g A\n", In*Ii, Ip*Ii);
							Print(NORMAL, "Max EQP Error %6.3f %%\tRMS EQP Error %6.3f %%\n\n", 100*Emax, 100*Erms);
							bess=1;
							nb=n;
							mb=m;
						}
							
							
						/* generate the computation mesh */						
						if(!pdef)
							Nm=MakeMesh(&Xdata, &Ydata, currentpars.corner1.x, currentpars.corner2.x, currentpars.corner1.y, currentpars.corner2.y, Nx, Ny);
						else if (!mesh)
							Nm=MakeMesh(&Xdata, &Ydata, Px1, Px2, Py1, Py2, Nx, Ny);
						else
						{
							/* read the mesh */
							Nm=ReadMesh(mesh, &Xdata, &Ydata);
						}		
							
						data=malloc(2*sizeof(double *));
						data[0]=malloc((Nm+2)*sizeof(double));
						
						for (i=0;i<=Nm;i++)
						{
							p.x=Xdata[i];
							p.y=Ydata[i];
							data[0][i]=Ii*Rsq*Bessel_Potential(p,nb,mb, rr, Vn, Vp);
						}
						SaveMesh(fout, data, Xdata, Ydata, Nm, 1);
						fclose(fout);
						
						free(data[0]);					
						free(data);
						free(Xdata);
						free(Ydata);
						
						break;
					}
					case BESS_FIELDLINE:
					{	
						int pc=-100/LLEN, TL;				
						if(!outfile)
							Error("No output file specified\n");
							
						if ((fout=fopen(outfile,"w"))==NULL)
							Error("Cannot open %s for saving data\n", outfile);
						
						
						if(!dim) /* we need to solve it all */
						{
							if (CheckPars())
								break;
							Dim(&Vn, &Vp, &I, &n, &m, &Emax, &Erms);
							dim=1;	
							bess=0;
							warp=0;
						}
						if(!bess) /* we need to solve it all */
						{
							if (CheckPars())
								break;
							InitBessel(n, m, rg, Rsq, &rr, &In, &Ip,  &Vn, &Vp, &Emax, &Erms);
							Print(NORMAL, "Vn %.12g V\tVp %.12g V\n", Vn, Vp);
							Print(NORMAL, "In %.12g A\tIp %.12g A\n", In*Ii, Ip*Ii);
							Print(NORMAL, "Max EQP Error %6.3f %%\tRMS EQP Error %6.3f %%\n\n", 100*Emax, 100*Erms);
							bess=1;
							nb=n;
							mb=m;
						}	
									
						if(!pdef || pdef==2)
						{
							x1=currentpars.corner1.x;
							x2=currentpars.corner2.x;
							y1=currentpars.corner1.y;
							y2=currentpars.corner2.y;	
						}
						else
						{
							x1=Px1;
							x2=Px2;
							y1=Py1;
							y2=Py2;	
						
						}
						/* compute the total amount of work */
						Print(NORMAL,"Computing fieldlines from:\n");
						TL=0;
						for (i=-n;i<=n;i++)
							for (j=-m;j<=m;j++)
							{
								p=Positive_c(i,j);
								if(VecInRange(p,x1,x2,y1,y2))
								{
									TL+=Nangle;
									Print(NORMAL," Source to Sink (%2i %2i)\n", i, j);
								}
								else if (VecInRange(p=Negative_c(i,j),x1,x2,y1,y2))
								{
									TL+=Nangle;
									Print(NORMAL," Sink to Source (%2i %2i)\n", i, j);
								}
							}
						
								
						Print(NORMAL,"Computing fieldlines. This may take some time!\n");
						Print(NORMAL,"For your entertainment I have an utterly useless\n");
						Print(NORMAL,"progress-bar, creeping slowly too 100%%. Enjoy!\n");
						l=0;
						while (l<LLEN)
						{							
							l++;
							if (l==1)
							{
								Print(NORMAL,"0 %%");
								l+=2;
							}
							else if (l==LLEN)
							{
								Print(NORMAL,"100 %%");
								l+=4;
							}
							else
								Print(NORMAL,"-");	
						}
						Print(NORMAL,"\n");		
						l=0;
						for (i=-n;i<=n;i++)
							for (j=-m;j<=m;j++)
							{
								p=Positive_c(i,j);
								if(VecInRange(p,x1,x2,y1,y2))
								{						
									fflush(stdout);
									for (k=0;k<Nangle;k++)
									{
										E.x=p.x+RR*currentpars.rp*cos(TWOPI*((double)k+0.5)/(double)Nangle);
										E.y=p.y+RR*currentpars.rp*sin(TWOPI*((double)k+0.5)/(double)Nangle);
										fl=Bessel_FieldLine(E, nb, mb, rr, Vn, Vp, x1, x2, y1, y2, &fll, fl_acc);
										SaveVecs(fout, fl, fll);
										free(fl);
										l++;
										pc=ProgressBar((100*l)/TL, pc, LLEN, TICS);
											
									}
								}
								else if (VecInRange(p=Negative_c(i,j),x1,x2,y1,y2))
								{
									fflush(stdout);
									for (k=0;k<Nangle;k++)
									{
										fflush(stdout);
										E.x=p.x+RR*currentpars.rn*cos(TWOPI*((double)k+0.5)/(double)Nangle);
										E.y=p.y+RR*currentpars.rn*sin(TWOPI*((double)k+0.5)/(double)Nangle);
										fl=Bessel_FieldLine(E, nb, mb, rr, Vn, Vp, x1, x2, y1, y2, &fll, fl_acc);
										SaveVecs(fout, fl, fll);
										free(fl);
										l++;
										pc=ProgressBar((100*l)/TL, pc, LLEN, TICS);
									}
								
								}
							
							}
									
						break;
					}
					case CONTACT_COORD:
						if(!outfile)
							Error("No output file specified\n");
							
						if ((fout=fopen(outfile,"w"))==NULL)
							Error("Cannot open %s for saving data\n", outfile);
						
						if(!dim) /* we need to solve it all */
						{
							if (CheckPars())
								break;	
							Dim(&Vn, &Vp, &I, &n, &m, &Emax, &Erms);
							Print(NORMAL, "Max EQP Error %6.3f %%\tRMS EQP Error %6.3f %%\n\n", 100*Emax, 100*Erms);
							dim=1;	
							bess=0;
							warp=0;
						}
							
						if(!pdef || pdef==2)
						{
							x1=currentpars.corner1.x;
							x2=currentpars.corner2.x;
							y1=currentpars.corner1.y;
							y2=currentpars.corner2.y;	
						}
						else
						{
							x1=Px1;
							x2=Px2;
							y1=Py1;
							y2=Py2;	
						
						}
						fl=malloc(2*sizeof(vec));
						
						for (i=-n;i<=n;i++)
							for (j=-m;j<=m;j++)
							{	
								fll=0;
								p=Positive_c(i,j);
								if(VecInRange(p,x1,x2,y1,y2))
									fl[fll++]=p;
									
								p=Negative_c(i,j);
								if(VecInRange(p,x1,x2,y1,y2))
									fl[fll++]=p;
								if (fll)	
									SaveVecs(fout, fl, fll);							
							}
						free(fl);
									
						break;
					case SOURCE_COORD:
						if(!outfile)
							Error("No output file specified\n");
							
						if ((fout=fopen(outfile,"w"))==NULL)
							Error("Cannot open %s for saving data\n", outfile);
						
						if(!dim) /* we need to solve it all */
						{
							if (CheckPars())
								break;	
							Dim(&Vn, &Vp, &I, &n, &m, &Emax, &Erms);
							Print(NORMAL, "Max EQP Error %6.3f %%\tRMS EQP Error %6.3f %%\n\n", 100*Emax, 100*Erms);
							dim=1;	
							bess=0;
							warp=0;
						}
							
						if(!pdef || pdef==2)
						{
							x1=currentpars.corner1.x;
							x2=currentpars.corner2.x;
							y1=currentpars.corner1.y;
							y2=currentpars.corner2.y;	
						}
						else
						{
							x1=Px1;
							x2=Px2;
							y1=Py1;
							y2=Py2;	
						
						}
						
						fl=malloc(2*sizeof(vec));
						
						for (i=-n;i<=n;i++)
							for (j=-m;j<=m;j++)
							{	
								fll=0;
								p=Positive(i,j);
								if(VecInRange(p,x1,x2,y1,y2))
									fl[fll++]=p;
									
								p=Negative(i,j);
								if(VecInRange(p,x1,x2,y1,y2))
									fl[fll++]=p;
									
								if (fll)	
									SaveVecs(fout, fl, fll);						
							}
						free(fl);
									
						break;
					case NONE:
						Warning("Warning: Word \"%s\" is not recognized\n", word);			
					default:
						break;
premature_end:
						Error("Premature end of input\n");
						exit(1);
				
				}
				begin=Begin(end); /* a new beginning starts at the end */
			}
			
		}
    		fgets(line, MAXSTRLEN-1, f);
	}
	free(line);
	free(word);
	if(outfile)
		free(outfile);
	if (mesh)
		free(mesh);
	return;
}
#undef LLEN
#undef RR
