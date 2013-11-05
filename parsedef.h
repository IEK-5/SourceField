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

/* parse keys and keywords */

/* enummeration of keys */
typedef enum {
	X1, X2,
	Y1, Y2,
	CXP, CXN,
	CYP, CYN,
	R, 
	RP, 
	RN, 
	RSQ,
	RG,
	CURRENT, 
	PX1, PX2, 
	PY1, PY2, 
	NY, NX,
	OFFSETS,
	OUTFILE, 
	MESH, 
	ERR,
	ERR_OFF,
	NANGLE_CHECK,
	NXM,
	OFF_ITER,
	CONTACT_COORD,
	SOURCE_COORD,
	RESIST,
	TWOPOINT,
	FOURPOINT,
	FIELD,
	POTENTIAL,
	FIELDLINE,
	WARP_RESIST,
	WARP_TWOPOINT,
	WARP_FOURPOINT,
	WARP_OFFSETS,
	WARP_FIELD,
	WARP_FIELDLINE,
	WARP_POTENTIAL,
	BESS_FAST,
	BESS_OFFSETS,
	BESS_RESIST,
	BESS_FIELD,
	BESS_POTENTIAL,
	BESS_FIELDLINE,
	NANGLE,
	FL_ACC,
	_QUIET,
	_NORMAL,
	_VERBOSE,
	_DEBUG,
	NONE
} PRSDEF;

/* Keyword to key mapping struct */
typedef struct {
	const char *name;
	PRSDEF PAR;
} KeyWord;


/* The keyword table mapping keywords to keys */
const KeyWord KeyTable[] =
{
      	{"x1", X1},
	{"x2", X2},
	{"y1", Y1},
	{"y2", Y2},
	{"cxp",CXP},
	{"cxn",CXN},
	{"cyp",CYP},
	{"cyn",CYN},
	{"r",  R},
	{"rp",  RP},
	{"rn",  RN},
	{"Rsq",RSQ},
	{"rg",RG},
	{"I",CURRENT},
	{"px1",PX1},
	{"px2",PX2},
	{"py1",PY1},
	{"py2",PY2},
	{"ny", NY},
	{"nx", NX},
	{"offsets",OFFSETS},
	{"file",OUTFILE},
	{"mesh",MESH},
	{"err", ERR},
	{"err_off", ERR_OFF},
	{"nangle_check", NANGLE_CHECK},
	{"NxM", NXM},
	{"off_iter", OFF_ITER},
	{"contact_coord",CONTACT_COORD},
	{"source_coord",SOURCE_COORD},
	{"resistance",RESIST},
	{"2point",TWOPOINT},
	{"4point",FOURPOINT},
	{"field",FIELD},
	{"potential",POTENTIAL},
	{"fieldline",FIELDLINE},
	{"warp_offsets",WARP_OFFSETS},
	{"warp_resistance",WARP_RESIST},
	{"warp_2point",WARP_TWOPOINT},
	{"warp_4point",WARP_FOURPOINT},
	{"warp_field",WARP_FIELD},
	{"warp_fieldline",WARP_FIELDLINE},
	{"warp_potential",WARP_POTENTIAL},	
	{"bess_fast",BESS_FAST},
	{"bess_offsets",BESS_OFFSETS},
	{"bess_resistance",BESS_RESIST},
	{"bess_field",BESS_FIELD},
	{"bess_potential",BESS_POTENTIAL},
	{"bess_fieldline",BESS_FIELDLINE},
	{"nangle", NANGLE},
	{"fl_acc",FL_ACC},
	{"out_quiet",_QUIET},
	{"out_normal",_NORMAL},
	{"out_verbose",_VERBOSE},
	{"out_debug",_DEBUG},
      	{NULL, NONE}
};
