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

vec Positive(int n, int m);
vec Negative(int n, int m);
vec Positive_c(int n, int m);
vec Negative_c(int n, int m);

void Offsets(void);
void Dim(double *Vn, double *Vp, double *I, int *n, int *m, double *Emax, double *Erms);
double Potential(vec p, int n, int m, double Vn, double Vp);
vec Field(vec p, int n, int m);
int VecInRange(vec v, double X1, double X2, double Y1, double Y2);
vec * FieldLine(vec s, int n, int m, double Vn , double Vp, double X1, double X2, double Y1, double Y2, int *len, double DL);
void InitBessel(int n, int m, double rg, double Rsq, double *rr, double *In, double *Ip, double *Vn, double *Vp, double *Emax, double *Erms);
void Bessel_Offsets(double rr);
double Bessel_Potential(vec p, int n, int m, double rr, double Vn, double Vp);
vec Bessel_Field(vec p, int n, int m, double rr);
vec * Bessel_FieldLine(vec s, int n, int m, double rr, double Vn, double Vp, double X1, double X2, double Y1, double Y2, int *len, double DL);
