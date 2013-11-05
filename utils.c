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
#include <stdarg.h>
#include "main.h"
#include "utils.h"
#define MAXSTRLEN 256

void Print(VERB_LEVEL v_level, const char *format_str, ...)
{
      	va_list ap;
      	va_start (ap, format_str);
	if (v_level<=verbose)
		vprintf(format_str, ap);
}

void Error( const char *format_str, ...)
{
      	va_list ap;
      	va_start (ap, format_str);
	vfprintf(stderr,format_str, ap); 
	exit(1);
}


void Warning( const char *format_str, ...)
{
      	va_list ap;
      	va_start (ap, format_str);
	vfprintf(stderr,format_str, ap); 
}

int ProgressBar(int pcn, int pco, int len, int tics)
/* pcn: new percentage complete */
/* pco: old percentage complete */
/* len: length of the progress bar */
/* tics: number of tics */
{
	int pc_n=(pcn*len)/100;
	int pc=pco, tic;
	int i;
	pco=(pco*len)/100;
	if (pco==len)
		return pc;
	
	tic=len/tics;
	
	for (i=pco+1;i<=pc_n;i++)
	{
		if(i%tic == 0)
		{
			if (i==0)
				printf("|");
			else
			{
				printf("\b\b\b\b");
				printf("|");
			}
		}
		else
		{
			printf("\b\b\b\b");
			printf("=");
		}
		pc=pcn;
		printf("%3i%%",pc);
	}
	if (pcn>pc)
	{
		printf("\b\b\b\b");
		printf("%3i%%",pcn);
		pc=pcn;	
	}
		
	if (pc_n==len)
		printf("\n");	
	fflush(stdout);
	
	return pc;
}

int MakeMesh(double **x, double **y, double x1, double x2, double y1, double y2, int Nx, int Ny)
{
	int i,j;
	double xs, ys;
	
	if (Nx==0)
		xs=0;
	else
		xs=(x2-x1)/(double)Nx;
	
	if (Ny==0)
		ys=0;
	else	
		ys=(y2-y1)/(double)Ny;
		
	(*x)=malloc((Nx*Ny+Nx+Ny+1)*sizeof(double));	
	(*y)=malloc((Nx*Ny+Nx+Ny+1)*sizeof(double));	
	
	for (i=0;i<=Nx;i++)
		for (j=0;j<=Ny;j++)
		{
			(*x)[i*(Ny+1)+j]=x1+i*xs;
			(*y)[i*(Ny+1)+j]=y1+j*ys;
		}
	return Nx*Ny+Nx+Ny;
}

int ReadMesh(char *fn, double **x, double **y)
{
	FILE *f;
	char c, *line;
	int k, Nm=0, Nalloc=100;
	
	(*x)=malloc((Nalloc+1)*sizeof(double));	
	(*y)=malloc((Nalloc+1)*sizeof(double));	
	
	if ((f=fopen(fn,"r"))==NULL)
		Error("Cannot open %s for reading\n", fn);
	
	line=malloc(MAXSTRLEN*sizeof(char));
    	fgets(line, MAXSTRLEN-1, f);
	
	while(feof(f)==0)
	{
		
    		k=sscanf(line, " %c", &c);
		if((k!=-1)&&(c!='#'))
		{
			k=sscanf(line, " %le %le", (*x)+Nm,(*y)+Nm);
			if(k!=-1)
			{
				Nm++;
				if (Nalloc==Nm)
				{
					Nalloc+=100;
					(*x)=realloc((*x), (Nalloc+1)*sizeof(double));	
					(*y)=realloc((*y), (Nalloc+1)*sizeof(double));					
				}
			}
			
		}
    		fgets(line, MAXSTRLEN-1, f);
	}
	free(line);
	(*x)=realloc((*x), (Nm+1)*sizeof(double));	
	(*y)=realloc((*y), (Nm+1)*sizeof(double));
	fclose(f);	
	return Nm-1;
}

void SaveMesh(FILE *f, double **Data, double *x, double *y, int Nm, int Nd)
{
	int i, j;
	for (i=0;i<=Nm;i++)
	{
		fprintf(f,"%.12e\t%.12e", x[i], y[i]);
		for (j=0;j<Nd;j++)
			fprintf(f,"\t%.12e", Data[j][i]);
		fprintf(f,"\n");
	}	
}

void SaveArray(FILE *f, double **Data, double *x, double *y, int Nx, int Ny)
{
	int i, j;
	for (i=0;i<=Nx;i++)
	{
		for (j=0;j<=Ny;j++)
			fprintf(f,"%.12e\t%.12e\t%.12e\n", x[i], y[j], Data[i][j]);
		fprintf(f,"\n");
	}	
}
void SaveVecs(FILE *f, vec *v, int Nv)
{
	int i;
	for (i=0;i<Nv;i++)
		fprintf(f,"%.12e\t%.12e\n", v[i].x, v[i].y);
	fprintf(f,"\n");
	fflush(f); /* keep file current for loooooong calculations and impatient users */
}

void SaveXY(FILE *f, double *x, double *y, int N)
{
	int i;
	for (i=0;i<=N;i++)
		fprintf(f,"%.12e\t%.12e\n", x[i], y[i]);
}


void PrintHeader()
{
	if (NORMAL>verbose)
		return;
	printf(" ____                           _____ _      _     _ \n"); 
	printf("/ ___|  ___  _   _ _ __ ___ ___|  ___(_) ___| | __| |\n"); 
	printf("\\___ \\ / _ \\| | | | '__/ __/ _ \\ |_  | |/ _ \\ |/ _` |\n"); 
	printf(" ___) | (_) | |_| | | | (_|  __/  _| | |  __/ | (_| |\n"); 
	printf("|____/ \\___/ \\__,_|_|  \\___\\___|_|   |_|\\___|_|\\__,_|\n"); 
	printf("SourceField Version %s   %s   B.E. Pieters \n", VERSION, __DATE__); 
	printf("IEK-5 Photovoltaik Forschungszentrum Juelich, Germany\n");  
	printf("\n");
}

void Disclamer()
{                                            
	printf("\n");
	printf("DISCLAMER:\n");
	printf("No warranties, either express or implied, are hereby  \n");
	printf("given. All software is supplied as is, without        \n"); 
	printf("guarantee. The user assumes all responsibility for    \n"); 
	printf("damages resulting from the use of this software,      \n");
	printf("including, but not limited to, attacks from angry     \n");
	printf("little animals (like kittens or bunnies).             \n");
	printf("\n");
}
