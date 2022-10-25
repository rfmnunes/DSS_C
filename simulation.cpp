// 
// Simulates a Grid
//

#include "simulation.h"

#define _CRTDBG_MAP_ALLOC
#include<iostream>
#include <crtdbg.h>


#ifdef _DEBUG
	#define DEBUG_NEW new(_NORMAL_BLOCK, __FILE__, __LINE__)
#else
	#define DEBUG_NEW new 
#endif


//#define _REPRODUCTIBLE 	// if active does the exact path as the sequential, but is slower


Cl_Simulation::Cl_Simulation(GSLibGridPars& grid_pars, Cl_mtrng& RNG, int n_procs )
{
	this->grid_def =  &grid_pars;
	this->outer_RNG = &RNG;
	this->numprocess = n_procs;
	this->mut1 = new boost::mutex; 

	//this->order = new GSLibGrid<long>(*grid_def);
};

Cl_Simulation::Cl_Simulation(const Cl_Simulation & in)
{
	this->grid_def = in.grid_def;
	this->outer_RNG = in.outer_RNG;
	this->numprocess = in.numprocess;
	this->mut1 = new boost::mutex;

	//this->order = new GSLibGrid<long>(*grid_def);
};

Cl_Simulation::~Cl_Simulation() 
{  
	//delete mut1;
	//long *nnodesim;
	//double *zmean;
	//double *zvariance;


	//double* M_curr;
	//double* Mnext_curr;
	//double* S_curr;

	//GSLibGridPars *grid_def;

	//GSLibGrid<unsigned int> *order;
	//GSLibGrid<float> *sim;

	//Cl_mtrng *outer_RNG;
	//unsigned long *outer_seeds; // seeds

	//boost::mutex *mut1;

	//bool *hEvent_aux;
	//boost::condition_variable *hEvent;

	//delete mut1;
//	delete[] sim;

//	//	delete grid_def;
////	delete outer_RNG;
////
////	delete secondary;
////	delete cc;
////	delete LVM;
////
};



void Cl_Simulation ::  Genpath (GSLibGridPars &grid_def, int nmult, int mults)
{

	unsigned int  ix,iy,iz;
	long long nn;
	long long temp2;
	long long swape;

		
	unsigned int nx= grid_def.get_nx();
	unsigned int ny= grid_def.get_ny();
	unsigned int nz= grid_def.get_nz();
	long long nxyz=grid_def.get_nxyz();

	// The Do loop assigns values to the arrays order and sim. order is now (1,2,...,nxyz-1,nxyz)
	// and sim has random numbers in the [0;1] subset */
	
	for (long long i=1; i<=nxyz; ++i)
	{
		order->grid[i] = i;
	}


	for (long long i=1; i<=nxyz; ++i)
	{
		temp2=  (long long) (outer_RNG->random_double() *(nxyz+1-i)) + i;

		swape = order->grid[i];
		order->grid[i]= order->grid[temp2];
		order->grid[temp2]=swape;
	}


	// TO DO SEE IF THIS CAN BE DONE
	// The multiple grid search works with multiples of 4 (yes, that is somewhat arbitrary):
	if (mults == 1) {
		nn=0;
		for (iz=1; iz<=nz; iz=iz+nmult) {
			for (iy=1; iy<=ny; iy=iy+nmult) {
				for (ix=1; ix<=nx; ix=ix+nmult) {
//					index = ix + (iy-1)*nx + (iz-1)*nxy;
					nn=nn+1;
				}
			}
		}
	}

	return;
};


GSLibGrid<float>* Cl_Simulation :: Simulate(int isim, misc& utils,
											Zones_Pars& Zones, KrigePars& krige, 
											BlocksPars& blocks, pseudoharddata& pseudo,
											SearchPars &search_pars)
{

	boost::thread_group tgroup;

	boost::posix_time::ptime SimStartTime = boost::posix_time::microsec_clock::local_time();

	logger << "Preparing simulation\n";

	// create a sim grid
	sim= new GSLibGrid<float>(*grid_def);


	long long nxyz=grid_def->get_nxyz();

	// put nullvalues on grid 
	for (long long i=1; i<=nxyz; ++i)
		sim->grid[i]=utils.nosvalue;


	int currzone;

	order= new GSLibGrid<unsigned int>(*grid_def);

	// Random path for this realization
		// first get the pseudo hard
	if(pseudo.usepseudo ==1) 
		pseudo.get_pseudodata_random_path(*grid_def, *outer_RNG,utils, *order);

	else// genpath will prepare order
		Genpath (*grid_def, nmult, mults);


	// initialize some stuff for zones
	
	nnodesim = new long[Zones.n_zones];
	M_curr = new double[Zones.n_zones];
	Mnext_curr = new double[Zones.n_zones];
	S_curr = new double[Zones.n_zones];

	for (unsigned int i = 0; i < Zones.n_zones; i++) {
		nnodesim[i] = 0;
		M_curr[i] = 0;
		Mnext_curr[i] = 0;
		S_curr[i] = 0;
	}

	// put the harddata on the grid
	if (search_pars.sstrat == 1)
	{ // TODO change the dtonode, remove a parameter
		for (unsigned int i = 0; i < Zones.n_zones; i++)
			Zones.zone[i]->harddata->Dtonode ( *grid_def, utils, sim->grid, *(Zones).zone_grid, i, Zones.n_zones);

		for (long long l=0; l<grid_def->get_nxyz(); ++l)
		{
			int id = (int)sim->grid[l+1];
			
			if (id > 0)
			{
				currzone = Zones.get_zone_from_index(l+1);
			
				if (currzone != -1) 
				{
					sim->grid[l + 1] = Zones.zone[currzone]->harddata->point[id].value;
					//point[id].value;

					++nnodesim[currzone];
					Mnext_curr[currzone] += (sim->grid[l + 1] - M_curr[currzone]) / nnodesim[currzone];
					S_curr[currzone] += (sim->grid[l + 1] - M_curr[currzone])*(sim->grid[l + 1] - Mnext_curr[currzone]);
					M_curr[currzone] = Mnext_curr[currzone];
					//M_curr[currzone] += sim->grid[l + 1];
				}
				//Mnext_curr[i] = 0;
				//S_curr[i] = 0;

				//Mnext_curr[currzone] += (newsim - M_curr[currzone]) / nnodesim[currzone];
				//S_curr[currzone] += (newsim - M_curr[currzone])*(newsim - Mnext_curr[currzone]);
				//M_curr[currzone] = Mnext_curr[currzone];
			}
		}

		// pass all zones
		//for (int currzone = 0; currzone < Zones.n_zones; currzone++)
		//	M_curr[currzone] = M_curr[currzone] / nnodesim[currzone];
	}


	// initialize some stuff for zones

	//nnodesim= new long[Zones.n_zones];

	//M_curr = new double[Zones.n_zones];
	//Mnext_curr = new double[Zones.n_zones];
	//S_curr = new double[Zones.n_zones];
			
	for (unsigned int i = 0; i < Zones.n_zones; i++){
		nnodesim[i]=0;
		//M_curr[i] = 0;
		//Mnext_curr[i] = 0;
		//S_curr[i] = 0;
	}

	// now we enter the parallel bits, we do it in parts to save memory 
	// in the coordination of the threads. 

	long long nodemax=0;

	logger << "Simulating\n";

	do // while (nodemax <= grid_def.nxyz);
	{
		// prepare stuff for paralelization
		hEvent= new boost::condition_variable[MAXNODES_THREAD];
		hEvent_aux= new bool[MAXNODES_THREAD];
		outer_seeds = new unsigned long[MAXNODES_THREAD];
	
		for(int i=0; i<MAXNODES_THREAD; ++i)	hEvent_aux[i] = false;

		for(int i=0; i<MAXNODES_THREAD; ++i)	outer_seeds[i] = outer_RNG->random_long();


		logger<<" .progress. "<<utilities::make_string(100.*((float)nodemax/(float)grid_def->get_nxyz()))<<"% \n";
		
		//Eigen::initParallel();

		for(unsigned int i=0; i<numprocess; ++i)
		{
			//auto binded_fun = boost::bind(&Cl_Simulation::SimulNode, this, i + 1, nodemax, krige, blocks,Zones, utils, pseudo, search_pars);
			//boost::thread *in_t = new boost::thread(boost::bind(&Cl_Simulation::SimulNode, this, i + 1, nodemax, krige, blocks, Zones, utils, pseudo, search_pars));

			tgroup.add_thread(new boost::thread(boost::bind(&Cl_Simulation::SimulNode, this, i + 1, nodemax, boost::ref(krige), blocks, boost::ref(Zones), utils, pseudo, search_pars)));
		}

		// wait for threads to finish, update nodemax
		tgroup.join_all();
		nodemax = nodemax + MAXNODES_THREAD;

		// delete what we used in this cycle
		
		delete[] hEvent;
		delete[] hEvent_aux;
		delete[] outer_seeds;

	} while (nodemax <=nxyz);

	
	delete [] nnodesim;
	delete[] M_curr;
	delete[] Mnext_curr;
	delete[] S_curr;

	//delete [] zmean;
	//delete [] zvariance;

	//order->headerflag = utils.headerflag;
	//order->null_data = utils.nosvalue;
	//order->n_vars = 1;
	//order->var_name = new std::string[order->n_vars];
	//order->var_name[0] = Zones.zone[0].harddata.colsname[Zones.zone[0].harddata.ivrl];
	//order->title = new std::string;
	//*order->title = order->var_name[0] + "_Simulation_path" + utilities::make_string(isim);


	//switch (utils.filetype)
	//{
	//case GEOEAS: // case 0 is an ascii file
	//	order->filename = new std::string;
	//	*order->filename = utils.outfl + "_" + utilities::make_string(isim) + "_path.out";
	//	order->write_grid_ascii();
	//	break;
	//case SGEMS: // case 1 is sgems binary
	//	order->filename = new std::string;
	//	*order->filename = utils.outfl + "_" + utilities::make_string(isim) + "_path.sgems";
	//	order->write_grid_binary_sgems();
	//	break;
	//}


	delete order;

	boost::posix_time::ptime SimEndTime = boost::posix_time::microsec_clock::local_time();

	boost::posix_time::time_duration msdiff = SimEndTime - SimStartTime;
	
	logger <<  ("-------------------------------------------------------------------------------\n");
	logger <<   " Simulation Elapsed time: "  << msdiff.total_seconds() << " seconds              \n" ;
	logger <<  ("-------------------------------------------------------------------------------\n");


	return sim;

}; 




double keep_distribution_limits(pseudoharddata &pseudo, int currentpseudo, HarddataPars &current_hardata, int in, double* M_curr, int currzone, double newsim, double simval)
{

	if ((pseudo.usepseudo == 1&&  /*currentpseudo !=-1 */in <= pseudo.n_pseudo_hard))
	{
		if (pseudo.pseudo_corr == 1)
		{
			newsim = newsim - (float)(M_curr[currzone] - pseudo.pseudo_distr.at(currentpseudo).get_distribution_average());
		}

		if (newsim < (*(pseudo.pseudo_distr[currentpseudo].orig_distr))[1])
		{
			newsim = /*(float)*/(*(pseudo.pseudo_distr[currentpseudo].orig_distr))[1];
		}
		if (newsim > /*(float)*/(*pseudo.pseudo_distr[currentpseudo].orig_distr)[pseudo.pseudo_distr[currentpseudo].n_data]) 
		{
			newsim = /*(float)*/  (*pseudo.pseudo_distr[currentpseudo].orig_distr)[pseudo.pseudo_distr[currentpseudo].n_data];

		}
	}
	else
	{
		if (newsim <= current_hardata.tmin)
		{
			//if (simval >= current_hardata.tmin)
			//	newsim = simval;
			//else
				newsim = (float)current_hardata.tmin;
		}
		if (newsim >= current_hardata.tmax)
		{
			//if (simval <= current_hardata.tmax)
			//	newsim = simval;
			//else
				newsim = (float)current_hardata.tmax;
		}
	}
	return newsim;
}

void Cl_Simulation ::SimulNode (long long start, long long nodemax, KrigePars &krige,
				BlocksPars &blocks, Zones_Pars &zones, misc &utils, pseudoharddata& pseudo,
				SearchPars &search_pars)
{	

	unsigned int     k;
	unsigned int currentpseudo=-1;

	//float   gauss_krig_mean;
	double simval = 0;
	double   newsim;

	float auxiliary=-999.25; // to avoid something strange with pointers
//	float cpdev;

	int currzone; // stores current zone 

	boost::mutex critical_sim;
	//boost::mutex critical_sim2;

	float dclaslocal;
	int current_bihst_class;

	// Loop over the nodes
	//	boost::posix_time::time_duration msdiff;

	//if (fmod(1000.*in,grid_def->get_nxyz()) == 0.0){
	//logger<<" .progress. "<<utilities::make_string(100.*((float)nodemax/(float)grid_def->get_nxyz()))<<"\n";
	//}

	//pseudo.updated_distribution = new Distribution(zones.zone[0]->harddata->harddata_distr);
	



	for (long long in = start + nodemax; in <= nodemax + MAXNODES_THREAD; in = in + numprocess) {
		if (in > grid_def->get_nxyz()) return;



#ifdef _REPRODUCTIBLE
		critical_sim.lock();
		thread_in[start - 1] = in;
		critical_sim.unlock();
#endif
		Nodes_search nodes;
		// get the new seed for this node
		Cl_mtrng inner_RNG(outer_seeds[in - nodemax - 1]);


		//get the current zone
		currzone = zones.get_zone_from_index(order->grid[in]);


		//VariogramPars* current_variogram = &zones.zone[currzone].variogram;

		//if (fmod(1000.*in,grid_def->get_nxyz()) == 0.0){
		//	logger<<" .progress. "<<utilities::make_string(100.*(float)in/(float)grid_def->get_nxyz())<<"\n";
		//}

		if (currzone == -1) //mask
		{
			hEvent_aux[in - nodemax - 1] = true;
			hEvent[in - nodemax - 1].notify_one();
			continue;
		}
		
		if (zones.zone[currzone]->bihist.bihflag == 1)
		{
			auxiliary = zones.zone[currzone]->bihist.auxiliary_grid->grid[order->grid[in]];

		}
		if (
			((float)fabs(sim->grid[order->grid[in]] - utils.nosvalue) > EPSLON)|| // already has value
			((krige.secondary_grid != NULL)&& (krige.secondary_grid->grid[order->grid[in]] == utils.nosvalue))|| // null on secondary
			((zones.zone[currzone]->bihist.bihflag == 1)&& (auxiliary == utils.nosvalue)) //no aux point in bihist
			)
		{
			hEvent_aux[in - nodemax - 1] = true;
			hEvent[in - nodemax - 1].notify_one();
			continue;
		}

		//	critical_sim.lock();
		//std::shared_ptr<VariogramPars> current_variogram;
		//current_variogram = std::make_shared<VariogramPars>(zones.zone[currzone]->variogram);

		// search the conditioning data

		//critical_sim.lock();


		//search_results *srch_res = search_pars.search_data(*mut1, hEvent, hEvent_aux, nodemax, *current_variogram,
		//	*grid_def, zones.zone[currzone].harddata, nodes,
		//	utils, sim->grid, order->grid, in, numprocess, thread_in, *current_variogram->search_radius);

		std::shared_ptr<search_results> srch_res=(search_pars.search_data(*mut1, hEvent, hEvent_aux, nodemax, zones.zone[currzone]->variogram,
			*grid_def, zones.zone[currzone]->harddata, nodes,
			utils, sim->grid, order->grid, in, numprocess, thread_in, *zones.zone[currzone]->variogram->search_radius));

		//critical_sim.unlock();


		// get one without blocks, maybe
		// get block krige stuff
		//we now have the nearby nodes. lets get the near blocks


		// search the conditioning data

		long *incloseblocks;
		

		if (blocks.blocksflag == 1) 
		{		
			//long *incloseblocks2;
			//float *close_blk_cov;
			
			incloseblocks = new long[blocks.n_blocks + 1]();
			//incloseblocks2 = new long[blocks.n_blocks + 1]();
			//close_blk_cov = new float[blocks.n_blocks + 1]();


			boost::scoped_array<long> incloseblocks2(new long[blocks.n_blocks + 1]);
			boost::scoped_array<float> close_blk_cov(new float[blocks.n_blocks + 1]);


			float c[1];
			//unsigned int n;

			float tol = (float)pow(10, -9);

			blocks.ncloseblocks = 0;

			// copia o array de covariancias bloco ponto e guarda que bloco e e a sua cov.
			for (unsigned int n = 1; n < blocks.n_blocks + 1; n++) {

				//usa o incloseblocks
				incloseblocks2[n] = n;
				close_blk_cov[n] = blocks.block2point_covt[n - 1][order->grid[in]]; // ver se +1 e mesmo assim

			}

			// agora ordena o array de covariancias arrastando o incloseblocks

			utils.sortit(1, blocks.n_blocks, &close_blk_cov[0], 1, (float *)&incloseblocks2[0], c, c, c, c, c, c);

			// agora guarda os maxblocks de maior covarincia 

			//blocks.ncloseblocks = 0
			for (unsigned int n = 0; n<blocks.n_blocks; n++) {

				if (fabs(close_blk_cov[blocks.n_blocks - n]) > tol) {

					incloseblocks[n + 1] = incloseblocks2[blocks.n_blocks - n];
					++blocks.ncloseblocks;
				}
			}

			if (blocks.ncloseblocks > blocks.maxblocks)
				blocks.ncloseblocks = blocks.maxblocks;

			
			//delete[] incloseblocks2;
			//delete[] close_blk_cov;

		}
		else
		{
			incloseblocks = 0;
		}
		//critical_sim.unlock();


		double cmean = 0;
		double cstdev = 0;

		double current_local_average;
		//float global_avg = zones.zone[currzone]->harddata->harddata_distr.get_distribution_average();

		//HarddataPars* current_hardata = &zones.zone[currzone]->harddata;

		//std::shared_ptr<HarddataPars> current_hardata(zones.zone[currzone]->harddata);

		if (pseudo.usepseudo == 1 && in <= pseudo.n_pseudo_hard)
		{
			if (pseudo.has_pseudo->grid[order->grid[in]] != 0)
			{
				currentpseudo = pseudo.has_pseudo->grid[order->grid[in]] - 1;
				current_local_average = zones.zone[currzone]->harddata->vmedexp;
				//current_local_average = pseudo.pseudo_distr[currentpseudo].get_distribution_average();
		//		//global_avg = zones.zone[currzone]->harddata->harddata_distr.get_distribution_average(); // pseudo.updated_distribution->get_distribution_average();
			}
		}
		else
		{ 
			//break;
			current_local_average = zones.zone[currzone]->harddata->vmedexp;
			//global_avg = current_hardata->vmedexp;
		}

		//current_local_average = zones.zone[currzone]->harddata->vmedexp;

		int calc_aux;
		if (blocks.blocksflag == 1)
		{
			calc_aux= (int)nodes.close_node_index->size() + blocks.ncloseblocks;
		}
		else
			calc_aux = srch_res->n_close_data + (int)nodes.close_node_index->size();

		//float * rhs_avg;

		boost::scoped_array<float> rhs_avg(new float[calc_aux]);
		


		// now get the local average of the found points
		if (pseudo.usepseudo == 1 /*&& in <= pseudo.n_pseudo_hard*/)
		{ // we have locals, lets try to find them
			for (int i = 0; i < nodes.close_node_index->size(); ++i)
			{
				long curr_point = nodes.close_node_index->at(i);
				// check if this point is a psewudosample
				if (pseudo.has_pseudo->grid[curr_point] != 0)
				{// this is a pseudowell
					//rhs_avg[i] = pseudo.pseudo_distr[pseudo.has_pseudo->grid[curr_point]-1].get_distribution_average();
					//rhs_avg[i] = sim->grid[curr_point];
					rhs_avg[i] = zones.zone[currzone]->harddata->vmedexp;
				}
				else
				{
					//rhs_avg[i] = sim->grid[curr_point];
					rhs_avg[i] = zones.zone[currzone]->harddata->vmedexp;
				}
			}
		}
		else{
			// all points have the global average
			for (int i = 0; i < srch_res->n_close_data+ nodes.close_node_index->size(); ++i)
			{
				rhs_avg[i] = zones.zone[currzone]->harddata->vmedexp;
			}
			if (blocks.blocksflag == 1)
			{
				for (int i = 0; i < nodes.close_node_index->size()+blocks.ncloseblocks; ++i)
				{
					rhs_avg[i] = zones.zone[currzone]->harddata->vmedexp;
				}
			}		
		}


		krig_results krgres(zones.zone[currzone]->variogram->cbb);

		switch (krige.ktype)
		{
		case 0:
			//devia ser algo do género:
			// krige.simple_krg(found_samples_cov, found_blocks_cov ) e retornar aprenas um krg_result

			krige.Simple_krg(blocks, *zones.zone[currzone]->harddata, *srch_res, *zones.zone[currzone]->variogram, *grid_def,  nodes,  order->grid[in], krgres, incloseblocks, current_local_average, rhs_avg);
			break;
		case 1:
			krige.Ordinary_krg(blocks, *zones.zone[currzone]->harddata, *srch_res, *zones.zone[currzone]->variogram, *grid_def, nodes, order->grid[in], krgres, incloseblocks);
			break;
		case 2:
			krige.LVM_krg(blocks, *zones.zone[currzone]->harddata, *srch_res, *zones.zone[currzone]->variogram, *grid_def, nodes, sim->grid, order->grid[in], krgres, incloseblocks);
			break;
		case 3:
			krige.Ext_drift_krg(blocks, *zones.zone[currzone]->harddata, *srch_res, *zones.zone[currzone]->variogram, *grid_def, nodes, sim->grid, order->grid[in], krige.secondary_grid->grid, krgres, incloseblocks);
			break;
		case 4:
/*O_Krig*/  krige.O_ColocCosim_krg(blocks, *zones.zone[currzone]->harddata, *srch_res, *zones.zone[currzone]->variogram, *grid_def, nodes, sim->grid, order->grid[in], krgres, incloseblocks,zones.zone[currzone]->med_exp_sec);
			//krige.ColocCosim_krg(blocks, *zones.zone[currzone]->harddata, *srch_res, *zones.zone[currzone]->variogram, *grid_def, nodes, sim->grid, order->grid[in], krgres, incloseblocks);
			break;
		case 5:
/*O_Krig*/	krige.O_LocalCC_ColocCosim_krg(blocks, *zones.zone[currzone]->harddata, *srch_res, *zones.zone[currzone]->variogram, *grid_def, nodes, sim->grid, order->grid[in], krige.cc_grid->grid[order->grid[in]], krgres, incloseblocks, zones.zone[currzone]->med_exp_sec);
			//krige.LocalCC_ColocCosim_krg(blocks, *current_hardata, *srch_res, *current_variogram, *grid_def, nodes, sim->grid, order->grid[in], krige.cc_grid->grid[order->grid[in]], krgres, incloseblocks);
			break;
		case 6:
			krige.LVM_ColocCosim_krg(blocks, *zones.zone[currzone]->harddata, *srch_res, *zones.zone[currzone]->variogram, *grid_def, nodes, sim->grid, order->grid[in], krige.colocorr, krgres, incloseblocks);
			break;
		case 7:
			krige.LVM_ColocCosim_krg(blocks, *zones.zone[currzone]->harddata, *srch_res, *zones.zone[currzone]->variogram, *grid_def, nodes, sim->grid, order->grid[in], krige.cc_grid->grid[order->grid[in]], krgres, incloseblocks);
			break;
		}

 		cmean = krgres.cmean;
		cstdev = krgres.cstdev;


		//delete[] rhs_avg;
		//if (cmean < current_hardata.tmin)
		//{
		//	cmean = current_hardata.tmin;
		//}
		//if (cmean > current_hardata.tmax)
		//{
		//	cmean = current_hardata.tmax;
		//}


		if (blocks.blocksflag == 1) 
		{
			delete[] incloseblocks;
		}
		//else
		//	delete incloseblocks;


		//average correction on pseudo
		if (pseudo.pseudo_corr == 1 && pseudo.usepseudo == 1 && nnodesim[currzone] <= pseudo.n_pseudo_hard && nnodesim[currzone] >1) //changes in avg and vr on datatonode can break this.
		{
			//	newsim= newsim+(zmean[currzone]-pseudo.pseudo_distr.at(currentpseudo).orig_distr[pseudo.pseudo_distr.at(currentpseudo).n_data/2]);
			//newsim = newsim + (M_curr[currzone] - pseudo.pseudo_distr.at(currentpseudo).get_distribution_average());
			//cmean = utils.average_correction(cmean, M_curr[currzone], pseudo.pseudo_distr.at(currentpseudo).get_distribution_average());
			double avg = pseudo.pseudo_distr[currentpseudo].get_distribution_average();
			double median = pseudo.pseudo_distr[currentpseudo].get_distribution_median();
			//cmean = utils.average_correction(cmean, M_curr[currzone], median);
			cmean = utils.average_correction(cmean, M_curr[currzone], avg);

		}

		if ((pseudo.usepseudo == 0 && nnodesim[currzone] > 1 && zones.zone[currzone]->bihist.bihflag == 0 )
			|| (pseudo.usepseudo == 1 && (nnodesim[currzone] > pseudo.n_pseudo_hard)||(in > pseudo.n_pseudo_hard)))
		{
			if (utils.icmean == 1)
			{
				cmean = utils.average_correction(cmean, M_curr[currzone], zones.zone[currzone]->harddata->vmedexp);
			}
			if (utils.icvar == 1)
			{
				cstdev = utils.variance_correction(cstdev, S_curr[currzone]/ std::max((int)(nnodesim[currzone] - 1), 1), zones.zone[currzone]->harddata->vvarexp);
			}
		}

	
		// what distribution are we using?
		Distribution* current_distribution;
		//std::unique_ptr<Distribution> current_distribution;
		
		if (zones.zone[currzone]->bihist.bihflag == 1)
		{
			current_distribution = new Distribution();
			zones.zone[currzone]->bihist.Makelocaldist(utils, auxiliary, &cmean, &cstdev, utils.icmean, utils.icvar, inner_RNG, *current_distribution );
		}
		else if (pseudo.usepseudo == 1 && in <= pseudo.n_pseudo_hard)
		{
			//if (pseudo.has_pseudo->grid[order->grid[in]]!=0)
			current_distribution = &pseudo.pseudo_distr.at(currentpseudo);
			//else
			//{
			//	current_distribution = &(zones.zone[currzone]->harddata->harddata_distr);
			//}
			//current_distribution = std::make_unique<Distribution>(&pseudo.pseudo_distr.at(currentpseudo))
		}
		else
		{
			current_distribution = &(zones.zone[currzone]->harddata->harddata_distr);
		}

		//simval = current_distribution.simulate_value_in_distribution(cmean, cstdev, currharddata, inner_RNG, utils.ntry);

		float gauss_krig_mean = (float)current_distribution->get_limits(cmean);


		double p;
		float bias_corr_avg = 0.0;
		float bias_corr_var = 0.0;

	
		k = 0;
		while (k < utils.ntry)
		{
			// Simulate a random gaussian value
			p = inner_RNG.random_double();

			double xp=Distribution::gauinv(p);
			
			//xp = Distribution::NormalCDFInverse(p);

			xp = xp * cstdev/*sqrt(cstdev)*/ + gauss_krig_mean;


			simval = current_distribution->backtr(xp
												, zones.zone[currzone]->harddata->ltail
												, zones.zone[currzone]->harddata->ltpar
												, zones.zone[currzone]->harddata->utail
												, zones.zone[currzone]->harddata->utpar);


			bias_corr_avg += simval;
			++k;

		}
		bias_corr_avg = bias_corr_avg / utils.ntry;
		//welford
		//bias_corr_var = S/(k);


		// Reter o valor simulado (SDSIM)
		newsim = simval;

		if (utils.ntry > 1) 
		{
			newsim += (float)(cmean - bias_corr_avg);
		}

		// LITO BIHIST
				
		if (zones.zone[currzone]->bihist.bihflag == 1)
		{
			newsim = zones.zone[currzone]->bihist.bihistlimits(newsim,
				zones.zone[currzone]->harddata->harddata_distr.zmin,
				zones.zone[currzone]->harddata->harddata_distr.zmax,
				current_distribution->zmin, current_distribution->zmax);
		}


		//also limit on a pseudo correction flag

		
		// now lets keep the limits of the distribution
		newsim = keep_distribution_limits(pseudo, currentpseudo, *zones.zone[currzone]->harddata, in, M_curr, currzone, newsim, simval);
		

		critical_sim.lock();	
		{
			//boost::lock_guard<boost::mutex> lk(m1);
			// try only when no point dist 
			++nnodesim[currzone];

			{
				// Condicionamento as medias locais
				
				Mnext_curr[currzone] += (newsim - M_curr[currzone]) / nnodesim[currzone];
				S_curr[currzone] += (newsim - M_curr[currzone])*(newsim - Mnext_curr[currzone]);
				M_curr[currzone] = Mnext_curr[currzone];

			}
			// bihist local mean correction

		
			if ( zones.zone[currzone]->bihist.bihflag==1)
			{
			
				//dclaslocal = (zones.zone[currzone]->bihist.orig_distr1[zones.zone[currzone]->bihist.ntbi]-zones.zone[currzone]->bihist.orig_distr1[1])/zones.zone[currzone]->bihist.imclas;
				//current_bihst_class = (int) ((auxiliary-zones.zone[currzone]->bihist.orig_distr1[1])/dclaslocal)+1;
		
				current_bihst_class = zones.zone[currzone]->bihist.get_current_bihist_class(auxiliary);
				//if (current_bihst_class<1 ) current_bihst_class =1;
		
				//if (current_bihst_class>(*zones).bihist[currzone].imclas ) current_bihst_class =(*zones).bihist[currzone].imclas;
				zones.zone[currzone]->bihist.n_points_in_bihst_class[current_bihst_class]++;

				zones.zone[currzone]->bihist.Mnext_in_bihst_class[current_bihst_class] += (newsim - zones.zone[currzone]->bihist.M_in_bihst_class[current_bihst_class]) / zones.zone[currzone]->bihist.n_points_in_bihst_class[current_bihst_class];
				zones.zone[currzone]->bihist.S_in_bihst_class[current_bihst_class] += (newsim - zones.zone[currzone]->bihist.M_in_bihst_class[current_bihst_class])*(newsim - zones.zone[currzone]->bihist.Mnext_in_bihst_class[current_bihst_class]);
				zones.zone[currzone]->bihist.M_in_bihst_class[current_bihst_class] = zones.zone[currzone]->bihist.Mnext_in_bihst_class[current_bihst_class];


				//zones.zone[currzone].bihist.Mnext_in_bihst_class[current_bihst_class] += newsim;
				//zones.zone[currzone].bihist.curr_var_in_bihst_class[current_bihst_class] += (float)pow((float) newsim, (int) 2);
		
			}
		
		}
		critical_sim.unlock();

		//sim->grid[order->grid[in]] = current_distribution->zmax-current_distribution->zmin;
		sim->grid[order->grid[in]] = newsim;

		//sim->grid[order->grid[in]] = krgres.cmean;
		{
			boost::lock_guard<boost::mutex> lk(*mut1);
			hEvent_aux[in-nodemax-1]=true;
			hEvent[in-nodemax-1].notify_all();
		}	

		//add the simulated point to harddata if we are using pseudowells, and recalculate the stats. 
		//if ((pseudo.usepseudo == 1 && in <= pseudo.n_pseudo_hard))
		//{
		//	float c[1];

		//	//pseudo.updated_distribution

		//	//Distribution *hd_dist = &zones.zone[currzone].harddata.harddata_distr;

		//	pseudo.updated_distribution->orig_distr.insert(pseudo.updated_distribution->orig_distr.end(),newsim);

		//	utils.sortit(1, (int)pseudo.updated_distribution->orig_distr.size()-1, &pseudo.updated_distribution->orig_distr[0], 0, c, c, c, c, c, c, c);

		//	pseudo.updated_distribution->update_gaussian_distribution();

		//	//zones.zone[currzone].harddata.vmedexp = pseudo.updated_distribution->get_distribution_average();
		//	//zones.zone[currzone].harddata.vvarexp = pseudo.updated_distribution->get_distribution_variance();

		//	//update limits
		//	
		//	if (pseudo.updated_distribution->orig_distr[1] <= zones.zone[currzone].harddata.tmin)
		//		zones.zone[currzone].harddata.tmin = pseudo.updated_distribution->orig_distr[1];
		//	if (pseudo.updated_distribution->orig_distr[pseudo.updated_distribution->orig_distr.size() - 1] >= zones.zone[currzone].harddata.tmax)
		//		zones.zone[currzone].harddata.tmax = pseudo.updated_distribution->orig_distr[pseudo.updated_distribution->orig_distr.size() - 1];
		//	//int ntr; // n after trimmed
		//	//int nt; // n data

		//}

		//sim->grid[order->grid[in]] = krgres.cmean;

		if (zones.zone[currzone]->bihist.bihflag == 1)
			delete current_distribution;

		//delete srch_res;
		//delete current_hardata;
	}

	//delete pseudo.updated_distribution;
	//logger <<   " krige Elapsed time: "  << msdiff.total_seconds() << " seconds              \n" ;


	return;
}
