#ifndef SIMULATION_INCLUDED
#define SIMULATION_INCLUDED




#include "zones.h"
#include "krige.h"
#include "log.h"
class Cl_Simulation
{

public:

	unsigned int numprocess;

	long *nnodesim; //+
	//double *zmean;
	//double *zvariance;


	double* M_curr;//+
	double* Mnext_curr;//+
	double* S_curr;//+


	GSLibGridPars *grid_def;

	GSLibGrid<unsigned int> *order;
	GSLibGrid<float> *sim;

	Cl_mtrng *outer_RNG; 
	unsigned long *outer_seeds; //+

	boost::mutex *mut1;

	bool *hEvent_aux;//+
	boost::condition_variable *hEvent;//+

	// for the synch of the threads, know in which node the threads are.
	long thread_in[8];

	// multiple grid for genpath - normally unused

	int nmult, mults;


	//Constructor
	Cl_Simulation (GSLibGridPars&, Cl_mtrng&, int);
	
	Cl_Simulation(const Cl_Simulation & in);

	// destructor 
	~Cl_Simulation();

	void Genpath (GSLibGridPars&, int, int);
	
	GSLibGrid<float>* Simulate( int ,misc&, Zones_Pars&, KrigePars&,
								BlocksPars&, pseudoharddata &,SearchPars &);

	void SimulNode (long long, long long,KrigePars&,
				BlocksPars&, Zones_Pars&,  misc&, pseudoharddata&,SearchPars& );

};

#endif