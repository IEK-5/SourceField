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
#include <math.h>
#include <stdlib.h>


/*  polynomial coefficients given by Abramowitz and Stegun
    Abramowitz, M., and Stegun, I.A. 1964, Handbook of Mathematical Functions, Applied Mathe-
    matics Series, Volume 55 (Washington: National Bureau of Standards; reprinted 1968 by
    Dover Publications, New York).
    
     C.W.CLENSHAW, CHEBYSHEV SERIES FOR MATHEMATICAL FUNCTIONS,
     MATHEMATICAL TABLES, VOL.5, 1962.
*/

double bessi0(double x)
/* Returns the modified Bessel function I0 (x) for any real x. */
{
	double ax,ans;
	static double P[7]={1.0, 3.5156229, 3.0899424, 1.2067492, 
	                    0.2659732, 0.360768e-1, 0.45813e-2};
	static double Q[9]={0.39894228, 0.1328592e-1, 0.225319e-2, 
	                   -0.157565e-2, 0.916281e-2, -0.2057706e-1, 
	                   0.2635537e-1, -0.1647633e-1, 0.392377e-2};
	double y;
	if ((ax=fabs(x)) < 3.75) 
	{
		y=x/3.75;
		y*=y;
		ans=P[0]+y*(P[1]+y*(P[2]+y*(P[3]+y*(P[4]+y*(P[5]+y*P[6])))));
			
	} 
	else 
	{
		y=3.75/ax;
		ans=(exp(ax)/sqrt(ax))*(Q[0]+y*(Q[1]+y*(Q[2]+y*(Q[3]+y*(Q[4]+y*(Q[5]+y*(Q[6]+y*(Q[7]+y*Q[8]))))))));
	}
	return ans;
}


double bessi1(double x)
/* Returns the modified Bessel function I1 (x) for any real x. */
{
	double ax,ans;
	static double P[7]={0.5, 0.87890594, 0.51498869, 0.15084934,
      			    0.2658733e-1, 0.301532e-2, 0.32411e-3};
	static double Q[9]={0.39894228, -0.3988024e-1, -0.362018e-2,
			    0.163801e-2, -0.1031555e-1, 0.2282967e-1,
			   -0.2895312e-1, 0.1787654e-1, -0.420059e-2};
	double y;
	if ((ax=fabs(x)) < 3.75) 
	{
		y=x/3.75;
		y*=y;
		ans=x*(P[0]+y*(P[1]+y*(P[2]+y*(P[3]+y*(P[4]+y*(P[5]+y*P[6]))))));
			
	} 
	else 
	{
		y=3.75/ax;
		ans=(exp(ax)/sqrt(ax))*(Q[0]+y*(Q[1]+y*(Q[2]+y*(Q[3]+y*(Q[4]+y*(Q[5]+y*(Q[6]+y*(Q[7]+y*Q[8]))))))));
	}
	return ans;
}

double bessk0(double x)
/* Returns the modified Bessel function K0 (x) for any real x. */
{
	double y,ans;
	static double P[7]={-0.57721566, 0.42278420, 0.23069756, 
	                    0.3488590e-1, 0.262698e-2, 0.10750e-3, 0.74e-5};
	static double Q[7]={1.25331414, -0.7832358e-1, 0.2189568e-1, 
	                   -0.1062446e-1, 0.587872e-2, -0.251540e-2, 0.53208e-3};

	if (x <= 2.0) 
	{        
		y=x*x/4.0;
		ans=(-log(x/2.0)*bessi0(x))+(P[0]+y*(P[1]+y*(P[2]+y*(P[3]+y*(P[4]+y*(P[5]+y*P[6]))))));
	} 
	else
	{
		y=2.0/x;
		ans=(exp(-x)/sqrt(x))*(Q[0]+y*(Q[1]+y*(Q[2]+y*(Q[3]+y*(Q[4]+y*(Q[5]+y*Q[6]))))));
	}
	return ans;
}

double bessk1(double x)  
/* Returns the modified Bessel function K1 (x) for any real x. */
{
	double y,ans;
	static double P[7]={1.0, 0.15443144, -0.67278579, -0.18156897,
	                    -0.1919402e-1, -0.110404e-2, -0.4686e-4};
	static double Q[7]={1.25331414, 0.23498619, -0.3655620e-1, 
	                   0.1504268e-1, -0.780353e-2, 0.325614e-2, -0.68245e-3};

	if (x <= 2.0) 
	{
		y=x*x/4.0;
		ans=log(x/2.0)*bessi1(x)+(1.0/x)*(P[0]+y*(P[1]+y*(P[2]+y*(P[3]+y*(P[4]+y*(P[5]+y*P[6]))))));
	} 
	else
	{
		y=2.0/x;
		ans=(exp(-x)/sqrt(x))*(Q[0]+y*(Q[1]+y*(Q[2]+y*(Q[3]+y*(Q[4]+y*(Q[5]+y*Q[6]))))));
	}
	return ans;
}

