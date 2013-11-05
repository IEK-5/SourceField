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
Pars currentpars={ 
	{-1, -1},		/* corner1 - lower left corner 				*/
	{1, 1}, 		/* corner2 - upper right corner 			*/
	{-0.5, 0.0}, 		/* contact_p - the p contact 				*/
	{0.5, 0.0}, 		/* contact_n - the n contact 				*/ 
	{0.0, 0.0},  		/* off_p - offset for the p contact 			*/
	{0.0, 0.0},   		/* off_n - offset for the n contact 			*/
	1e-4, 1e-4, 1e-4,   	/* rp, rn - radii, Rsq - sheet resistance and err - error	*/ 
	1e-2,			/* Err_Off */
	1000000,		/* NxM */
	100,			/* Off_Iter*/
	40			/* Nangle*/
};
VERB_LEVEL verbose=NORMAL;
int fixverb=0;
