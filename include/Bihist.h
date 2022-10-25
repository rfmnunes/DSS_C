#ifndef BIHIST_INCLUDED
#define BIHIST_INCLUDED

#include "boost_include.h"
#include "hardata.h" 
#include "math_random.h"
#include "log.h"

class BihistogramPars{

public:

	std::string		bihistfl, // bihistogram
					auxbihistfl; // bihistogram auxiliary

	GSLibGrid<float>* auxiliary_grid;

	int		bihflag;

	long	*n_points_in_bihst_class;
	//double	*curr_mean_in_bihst_class;  //bihist control means correction
	//double	*curr_var_in_bihst_class;  //bihist control var correction
	double *Mnext_in_bihst_class;
	double *M_in_bihst_class;
	double *S_in_bihst_class;

	int		imclas;
	int		ntbi; // number of bihist points
	
	float	*orig_distr1, // bihist main var
			*orig_distr2;  // bihist acessory var
	
private:
	double	*target_mean_bihst_class;
	double	*target_var_bihst_class;
	
	float dclaslocal; // espacamento entre classes 

public:

	// constructor
	//BihistogramPars(registry_type *, int, misc&, GSLibGridPars&, int );
	//
	~BihistogramPars();

	int get_bihist_pars(registry_type *, int);

	int load_bihist_aux(GSLibGridPars& , float , int, int);

	int Loadbihist(misc& );

	/*Distribution*/void Makelocaldist(misc& , float , double *, double *, int , int ,Cl_mtrng&, Distribution &bih_local_distr);

	float bihistlimits(double& , double, double, double, double);

	int	get_current_bihist_class(const float);
}; // end class bihist 


#endif