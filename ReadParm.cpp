//
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//
//    OBJECTIVES:
//      The input parameters and data are read in from their files. Some quick
//      error checking is performed and the statistics of all the variables
//      being considered are written to standard output.
//
//    INPUT VARIABLES:
//      filepar - name of the parameter datafile
//
//    OUTPUT VARIABLES:
//
//
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//

//#include "SSdir.h"
//#include "krige.h"
#include "registry.h"
#include "math_random.h"
#include "utils.h"
#include "log.h"

int Readparm (registry *reg, Cl_mtrng** outer_RNG)
{
	reg_key *k;
	double p;

	long seed;
	k = get_key(reg,  ("GENERAL"),  ("SEED"));
	if (k){
		seed = get_long(k);
		//seed_random(seed);
		
		*outer_RNG = new Cl_mtrng(seed);
		
		p = (*outer_RNG)->random_double();

		(void)p;
	}else{
	
		//generate a time seed
		long s;
		time_t seconds;
		s = time ( &seconds ); /* get CPU seconds since 01/01/1970 */
				
		seed = abs(((s*181)*((83)*359))%104729);

		*outer_RNG = new Cl_mtrng(seed);
		
		p = (*outer_RNG)->random_double();

		(void)p;

	}
	logger << " Seed: " << utilities::make_string(seed) <<"\n";

	
	// end of general


	return 1; // exit ok
}

