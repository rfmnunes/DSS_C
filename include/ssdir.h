#ifndef SSDIR_INCLUDED
#define SSDIR_INCLUDED

#include "boost_include.h"
//#include <time.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
//
#include <fstream>
#include <vector>
#include <iostream>


//#define _REPRODUCTIBLE 	// if active does the exact path as the sequential, but is slower
							// if commented is faster, and geostatistically accurate, but  simultions 
							// with more than one processor are not reproductible, even with same seed.
#include "log.h"
#include "logging.h"

#include "zones.h"

#include "krige.h"

#include "simulation.h"




/*void     Readdata	(KrigePars& krige, GSLibGridPars& grid_def, 
					GSLibGrid<float>** secondary_grid,GSLibGrid<float> **LVM, GSLibGrid<float>** cc_grid,	
					float nosvalue, BlocksPars &blocks, Zones_Pars &Zones, int headerflag,  int filetype);
					*/
int      Readparm	(registry *reg, Cl_mtrng** outer_RNG);


#endif