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

int InitWarp (double *scale, double *k, double *K);
vec TransForm (double scale, double K, double k, vec pos);
double Warp_Potential(double scale, double K, double k, vec p, double Vn, double Vp, int eqcontact);
vec Warp_Field(double scale, double K, double k, vec p);
void Warp_Offsets(double scale, double K, double k);
void WarpContactPotential(double scale, double K, double k, double *Vn, double *Vp, double *Emax, double *Erms);
vec * Warp_FieldLine(vec s, double scale, double K, double k, double Vn, double Vp, double X1, double X2, double Y1, double Y2, int *len, double DL);
