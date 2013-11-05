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


void Print(VERB_LEVEL  v_level,  const char *format_str, ...);
void Error( const char *format_str, ...);
void Warning( const char *format_str, ...);

int ProgressBar(int pcn, int pco, int len, int tics);

int MakeMesh(double **x, double **y, double x1, double x2, double y1, double y2, int Nx, int Ny);
int ReadMesh(char *fn, double **x, double **y);

void SaveMesh(FILE *f, double **Data, double *x, double *y, int Nm, int Nd);
void SaveArray(FILE *f, double **Data, double *x, double *y, int Nx, int Ny);
void SaveVecs(FILE *f, vec *v, int Nv);
void SaveXY(FILE *f, double *x, double *y, int N);

void PrintHeader(void);
void Disclamer(void);
