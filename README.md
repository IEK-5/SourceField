SourceField
===========
 *****************************************************************              
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
 * Dr. Bart E. Pieters 2013                                      *
 *                                                               *             
 *****************************************************************                                                                             


This document describes how to use the compile the SourceField program. The program can be used to evaluate 4-point
measurements of the sheet resistance or, in conjunction with the modulesim scripts, it can be used to compute the 
impact of local defects (shunts) on the performance of a thin-film solar cell in a module. The model is described in:
B. E. Pieters, and U. Rau, "A new 2D model for the electrical potential in a cell stripe in thin-film solar modules 
including local defects", Prog. Photovolt: Res. Appl. (2013), DOI: 10.1002/pip.2436

Binaries are available under:
https://github.com/IEK-5/SourceField/releases/latest

COMPILING:
The program uses a simple make-script to build. You may have to adapt this script to your system. If all is 
succesfull you should end up with a SourceField executable.

USAGE:
./SourceField [-h] [verbose-level] <inputfile>
optionale -h option returns a list of commands with a description
optional verbose-level argument can be:
       -q        -   quiet
       -n        -   normal
       -v        -   verbose
       -db       -   debug   
In the input file one can specify a sequence of commands to the program.
