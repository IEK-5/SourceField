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
#include "parse.h"
#include "main.h"
#include "utils.h"
#include "defaults.h"
#include "keywords.h"


int main(int argc, char **argv)
{
	int f=1;
	if(argc==2||argc==3)
	{
		if (argc==3)
		{
			fixverb=1;
			if (strncmp(argv[1],"-q",2)==0)
				verbose=QUIET;
			else if (strncmp(argv[1],"-v",2)==0)
				verbose=VERBOSE;
			else if (strncmp(argv[1],"-db",3)==0)
				verbose=DEBUG;
			else if (strncmp(argv[1],"-n",2)==0)
				verbose=NORMAL;
			else
			{
				PrintHeader();
				fprintf(stderr,"USAGE:\n%s [verbose-level] <inputfile>\n", argv[0]);
				fprintf(stderr,"optional verbose-level argument can be:\n");
				fprintf(stderr,"       -q    -   quiet\n");
				fprintf(stderr,"       -n    -   normal\n");
				fprintf(stderr,"       -v    -   verbose\n");
				fprintf(stderr,"       -db   -   debug\n");
				return 1;
			}
			f++;				
		}
		if (strncmp(argv[1],"-h",2)==0)
		{
			int i;
			for (i=0;i<NumKWT;i++)
				printf("%s\n", KWTARR[i]);
			return 0;
		}
		PrintHeader();
		Parse (argv[f]);
	}
	else
	{
		PrintHeader();
		fprintf(stderr,"USAGE:\n%s [verbose-level] <inputfile>\n", argv[0]);
		fprintf(stderr,"optional verbose-level argument can be:\n");
		fprintf(stderr,"       -q        -   quiet\n");
		fprintf(stderr,"       -n        -   normal\n");
		fprintf(stderr,"       -v        -   verbose\n");
		fprintf(stderr,"       -db       -   debug\n");
		Disclamer();
		return 1;
	}
	return 0;	
}
