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
 *                    http://www.fz-juelich.de/ief/ief5/index    *              
 *****************************************************************
 *                                                               *
 * Dr. Bart E. Pieters 2011                                      *
 *                                                               *             
 *****************************************************************/                                                                             

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
