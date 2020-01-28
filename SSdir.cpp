//
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//            Direct Sequential Simulation and Cosimulation (DSSIM and CoDSSIM)
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//

#include "ssdir.h"

// this block of code here is to catch memory leaks
//#define _CRTDBG_MAP_ALLOC  
//#include <stdlib.h>  
//#include <crtdbg.h>  
//
////#ifdef _DEBUG
////#define DBG_NEW new ( _NORMAL_BLOCK , __FILE__ , __LINE__ )
////// Replace _NORMAL_BLOCK with _CLIENT_BLOCK if you want the
////// allocations to be of _CLIENT_BLOCK type
////#else
////#define DBG_NEW new
////#endif
//
//#ifdef _DEBUG
//#define DEBUG_NEW new(_NORMAL_BLOCK, __FILE__, __LINE__)
//#define new DEBUG_NEW
//#endif

#ifdef LICENSED_VERSION
	#include "licensing.h"
#endif

clss_log g_log;
//src::logger lg;
Logger logger;
int main (int argc, char *argv[])
{
	// also for debug memory leaks
	_CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF);

	unsigned long   i;
	int    isim;

	char default_registry_file[] = "ssdir.par";
	char *reg_file;
	registry  *r;
	int numprocess; //number of processors to use
	int max_proc;

	GSLibGridPars *grid_def;
	
	BlocksPars blocks;

	misc utils;

	Cl_mtrng *outer_RNG;

	int opcao_escreve=1;

	
	max_proc=64/*MAX_THREADS*/;
	//if (argc >= 2) strcpy (filepar,argv[1]);
	if (argc >= 3) max_proc=atoi(argv[2]);
	if (argc >= 4) opcao_escreve=atoi(argv[3]);

	// Determine how many processors are in the system
	numprocess=boost::thread::hardware_concurrency();//SystemInfo.dwNumberOfProcessors;

	if (numprocess==0){
		return -1; // cannot find number of processors, indicate manualy
	}

	if (numprocess > max_proc && max_proc>0)
		numprocess=max_proc;
	

	//numprocess=1;


	//let's check license
	//std::cout << "Passone";
#ifdef LICENSED_VERSION
	char* licfeature = "Interwell_DSS";
	char* licversion = "20180801.0";

	for (int i = 1; i <= numprocess; ++i)
	{
		//std::cout << "Passtwo";
		license_begin(licfeature, licversion);
		//std::cout << "Passone";
	}
#endif
/*
	Zones_Pars Zones;*/
	//= new Zones_Pars;
	

	// Begin execution
	std::cout <<  "-------------------------------------------------------------------------------\n";
	std::cout <<  "                      Direct Sequential Simulation, SSDir                      \n";
	std::cout <<  "                      Parallel mode: " << numprocess << " cores\n";
	std::cout <<  "-------------------------------------------------------------------------------\n";


	// Begin elapsed time counter
	boost::posix_time::ptime StartTime = boost::posix_time::microsec_clock::local_time();


	// Read the parameters

	if (argc >= 2)  /* evaluate if a registry file (or more) was passed as parameter */
		reg_file = argv[1];
	else     /* switch to default file */
	{
		std::cout << "No registry file passed as argument! Trying default file...\n";
		reg_file = default_registry_file;
	}

	r = new_registry(reg_file);

	if (r)  /* attempts to load registry file */
		std::cout << "Registry loaded...\n";		/* evaluate if the first registry file was loaded successfully */
	else
		return 1;
	// first prepare the log
	g_log.open_log(r);
		
	logger.set_file(g_log.dbgfl.c_str());
	logger << " Debugging file (1/2/3): " << g_log.dbgfl << " "<< g_log.idbg << "\n";


	// And output...
    logger << "Log prepared!\n";
	
	//if ( !krige.Read_krige_pars(r))
	//	//logger << "Kriging parameters loaded...\n"; 
	//	logger << "Kriging parameters loaded...\n";
	//	//g_log.log_string( "Kriging parameters loaded...\n");
	//else 
	//	return 1;

	pseudoharddata* pseudo = new pseudoharddata();
	
	if ( !(pseudo->read_pseudo_parameters(r)))
		logger << "Pseudo Harddata parameters loaded...\n";
		//g_log.log_string("Pseudo Harddata parameters loaded...\n");
	else 
		return 1;

	if (!blocks.read_blocks_parameters(r))
		logger << "Blocks parameters loaded...\n";
		//g_log.log_string("Blocks parameters loaded...\n");
	else 
		return 1;

	if (!utils.read_miscelaneous_parameters(r))
		logger << "Miscelaneous parameters loaded...\n";
	else 
		return 1;

	//SuperBlocksPars superblocks;
	SearchPars search_pars;
	search_pars.superblocks= new SuperBlocksPars;
	search_pars.get_search_pars(r);

	if (Readparm (r, &outer_RNG))
		logger << "Parameters loaded...\n";
	else
		return 1;


	//now we know how many zones do we have. allocate hardata 

//	harddata = new HarddataPars[n_zones];
	
	
	// get parameters
	grid_def = new GSLibGridPars(r);


	
	// define the dimension agnostic grid
	// first prepare inputs
	//std::vector<int> vector_n_bl;
	//std::vector<float> vector_o_bl;
	//std::vector<float> vector_s_bl;

	//std::vector<int> vector_c_bl;
	//
	//vector_n_bl.push_back(grid_def->get_nx());
	//vector_o_bl.push_back(grid_def->get_ox());
	//vector_s_bl.push_back(grid_def->get_xsiz());

	//vector_n_bl.push_back(grid_def->get_ny());
	//vector_o_bl.push_back(grid_def->get_oy());
	//vector_s_bl.push_back(grid_def->get_ysiz());

	//vector_n_bl.push_back(grid_def->get_nz());
	//vector_o_bl.push_back(grid_def->get_oz());
	//vector_s_bl.push_back(grid_def->get_zsiz());

	//vector_c_bl.push_back(grid_def->get_nx());
	//vector_c_bl.push_back(grid_def->get_nx()*grid_def->get_ny());
	//vector_c_bl.push_back(grid_def->get_nx()*grid_def->get_ny()*grid_def->get_nz());

	//CoordinatesClass<int> in_n_blocks(&vector_n_bl);
	//CoordinatesClass<float> in_origin_blocks(&vector_o_bl);
	//CoordinatesClass<float> in_size_blocks(&vector_s_bl);
	//
	//CoordinatesClass<int> in_cumulative_blocks(&vector_c_bl);


	//GSLibGridPars_generic *grid_def_MD = new GSLibGridPars_generic(in_n_blocks, in_origin_blocks, in_size_blocks, in_cumulative_blocks);


	//CoordinatesClass<float> test = grid_def_MD->get_xyz_from_index_f(123);


		
	// now allocate and read zone-dependent parameters - also harddata, bihist, variogram...

	Zones_Pars Zones(r, "");

	//if (!Zones.get_zone_pars(r,""))
	//	logger << "Zoning Parameters loaded...\n";
	//else{
	//	logger <<"Zoning Parameters failed to load!\n";
	//	return 1;
	//}

	Zones.read_zone_data(utils,*grid_def);

	KrigePars krige(r, *grid_def, utils, Zones);
	// Read the data
	// pode haver problema com deriva externa...


	// we no longer need the parameters
	delete_registry(r);
	
	//Readdata (krige, *grid_def, &secondary_grid, &LVM, &cc_grid,
	//		utils.nosvalue, blocks, Zones,utils.headerflag,utils.filetype );
	//// for each zone read the hardata and bihist 
	

	// allocate nd dependent stuff
	//search_pars.inclose = new long[Zones.nd_g];



	// Setup the rotation/anisotropy matrices that are needed for the variogram and search
	
	//for (i = 1; i <= Zones.n_zones; i++){
	//	Zones.variogram[i-1].prepare_rotation(); // TODO PUT ON ZONES 
	//}

	// Set up the super block search
	if (search_pars.sstrat == 0) {
		logger << " INITIALIZE SUPER-BLOCKS \n\n";
		search_pars.superblocks->setsupr (utils, *grid_def,Zones.zone[0]->harddata);// TODO  remove hardcoded 0
		search_pars.superblocks->picksupr ( Zones.zone[0]->variogram); // TODO SEE THIS VARIOGRAM
	}

	// Set up the covariance table for spiral search
	
	for (i = 0; i < Zones.n_zones; i++){
			
		Zones.zone[i]->variogram->Ctable (*grid_def,utils); // TODO PUT ON ZONES
		
	}

	if ( blocks.blocksflag == 1) {
		// prepare block kriging if in use
		// we put the blocks in a harddata object  (for now)


		// now set a search for block centroids

		blocks.get_blocks_covtable(Zones.zone[0]->variogram);// TODO GET THIS WORKING WITH BLOCKS
		blocks.get_block_to_point_covtable( Zones.zone[0]->variogram, *grid_def );

	}// blocks


	// pseudo-hard data stuff
	if(pseudo->usepseudo ==1)
	{	
		//switch (utils.filetype)
		//{
		//case GEOEAS: // case 0 is an ascii file
			pseudo->read_pseudo_harddata(utils, *grid_def);
		//	break;
		//case SGEMS: // case 1 is sgems binary
		//	pseudo->read_pseudo_harddata_binary(utils, *grid_def);
		//	break;
		//}
	}
	
	//   // Now, simulate images
	
	Cl_Simulation *Simulator_engine= new Cl_Simulation ( *grid_def, *outer_RNG, numprocess);

	logger << " BEGIN SIMULATION\n";
	for (isim=1; isim<=utils.nsim; ++isim) {
		logger << " Working on realization number " << isim << " \n";

		// call a simulation function, return the simulated grid
		Simulator_engine->Simulate(isim, utils, Zones, krige, blocks, *pseudo, search_pars);

		 
		//Simulator_engine->sim->filename  = utils.outfl + "_" + utilities::make_string(isim) + ".out" ;
		Simulator_engine->sim->headerflag=utils.headerflag;
		Simulator_engine->sim->null_data=utils.nosvalue;
		Simulator_engine->sim->n_vars=1;		
		Simulator_engine->sim->var_name=new std::string[Simulator_engine->sim->n_vars];
		Simulator_engine->sim->var_name[0]= Zones.zone[0]->harddata->file_header->colsname[Zones.zone[0]->harddata->ivrl];
		Simulator_engine->sim->title = new std::string;
		*Simulator_engine->sim->title=Simulator_engine->sim->var_name[0]+"_Simulation_"+ utilities::make_string(isim);


		switch(utils.filetype)
		{
		case GEOEAS: // case 0 is an ascii file
			Simulator_engine->sim->filename = new std::string;
			*Simulator_engine->sim->filename  = utils.outfl + "_" + utilities::make_string(isim) + ".out" ;
			Simulator_engine->sim->write_grid_ascii();
			break;
		case SGEMS: // case 1 is sgems binary
			Simulator_engine->sim->filename = new std::string;
			*Simulator_engine->sim->filename  = utils.outfl + "_" + utilities::make_string(isim) + ".sgems" ;
			Simulator_engine->sim->write_grid_binary_sgems();
			break;
		}


        delete Simulator_engine->sim;


	}
	delete Simulator_engine;
	boost::posix_time::ptime EndTime = boost::posix_time::microsec_clock::local_time();

	boost::posix_time::time_duration msdiff = EndTime - StartTime;
	logger <<  ("-------------------------------------------------------------------------------\n");
	logger << " Elapsed time: "  << msdiff.total_seconds() << " seconds               Press any key to continue...   \n" ;
	logger <<  ("-------------------------------------------------------------------------------\n");
	// if (opcao_escreve == 1)
	  //getch();


	// deallocations
	//int j;
	//delete[] Zones->variogram->ixnode;
 //   delete[] Zones->variogram->iynode;
 //   delete[] Zones->variogram->iznode;
 //   for (i = 0; i <= (2*Zones->variogram->nctx+2); ++i) {
 //       for (j = 0; j <= (2*Zones->variogram->ncty+2); ++j) {
 //           delete[] Zones->variogram->covtab[i][j];
 //       }
 //       delete[] Zones->variogram->covtab[i];
 //   }
 //   delete[] Zones->variogram->covtab;
 //   for (i = 0; i <= Zones->variogram->nst; ++i) {
 //       for (j = 0; j < 3; ++j) {
 //           delete[] Zones->variogram->rotmat[i][j];
 //       }
 //       delete[] Zones->variogram->rotmat[i];
 //   }
 //   delete[] Zones->variogram->rotmat;
 //   delete[] Zones->variogram->it;
 //   delete[] Zones->variogram->cc;
 //   delete[] Zones->variogram->aa;
 //   delete[] Zones->variogram->ang1;
 //   delete[] Zones->variogram->ang2;
 //   delete[] Zones->variogram->ang3;
 //   delete[] Zones->variogram->anis1;
 //   delete[] Zones->variogram->anis2;
 //   delete[] Zones->variogram;
	//if (Zones->bihist->bihflag ==1){
 ////   delete[] Zones->bihist->auxiliary_grid->grid;
 //   //delete Zones->bihist->auxiliary_grid;
 //   //delete[] Zones->bihist->vrtr1;
 //   //delete[] Zones->bihist->vrtr2;
 //   //delete[] Zones->bihist->n_points_in_bihst_class;
 //   //delete[] Zones->bihist->curr_mean_in_bihst_class;
 //   //delete[] Zones->bihist->curr_var_in_bihst_class;
 //   delete[] Zones->bihist;
	//}
 //   delete Zones->zone_grid;
 //   delete[] Zones->harddata;
    //delete Zones;
	//if (secondary_grid!=NULL){

 //  // delete secondary_grid;
	//}.

    //delete outer_RNG;
    //delete grid_def;
	
    //delete[] search_pars.inclose;


	return 0;
}
