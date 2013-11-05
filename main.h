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
                                                                           
/* 2pi @ 16 significant digits */
#define TWOPI 6.283185307179587

/* a 2D vector struct */
typedef struct {
	double 	x,y;
} vec;

/* the geometry of our structure */ 
/* corner1 and corner2:		specify the rectanglular lamina. Both coordinates of    */
/*				corner1 are nagative and both of corner2 are positive   */
/*				(i.e., the origine lies within the rectangle).   	*/
/* contact_p and contact_n:	Specify the centers of the contacts (p=source, n=sink)	*/
/* off_p and off_n:		The offsets of the actual source/sink positions with 	*/
/*				respect to the center of the contacts (specified above).*/
/*				These parameters are computed such that the circumfence */
/*				of the conatcts are approximately equipotential-lines.	*/
/* rn, rp:			Radii of the circular contacts				*/
/* Rsq:				Ohm square value of the lamina				*/
/* Err:				Maximal alowable error in current confinement. The 	*/
/*				program uses the currenmt confinement as a stop 	*/
/*				criterium. This means that source/sink pairs are added  */
/*				untill less than the error fraction of the current 	*/
/*				leaks across the lamina boundaries.			*/
/* Err_Off:			Maximal alowable error offset calculation	 	*/
/* NXM:				Maximal number of source sink pairs/4			*/
typedef struct {
	vec corner1, corner2;
	vec contact_p, contact_n;
	vec off_p, off_n;
	double 	rp, rn,
		Err,
		Err_Off;
	int NXM;
	int Off_Iter;
	int Nangle;
} Pars;

/* verbose level */ 
typedef enum {QUIET=0, NORMAL, VERBOSE, DEBUG} VERB_LEVEL;

/* initialization in utils.c */ 
extern VERB_LEVEL verbose;
extern int fixverb;
/* initialization in sourcefield.c */ 
extern Pars currentpars;
