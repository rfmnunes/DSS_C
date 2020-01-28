#ifndef ZONES_INCLUDED
#define ZONES_INCLUDED

#include "variogram.h"
#ifndef BIHIST_INCLUDED
#include "bihist.h"
#endif


class Zone 
{
public:
	//VariogramPars variogram;
	std::shared_ptr<VariogramPars> variogram;
	
	BihistogramPars bihist;
	
	//HarddataPars harddata;	
	std::shared_ptr<HarddataPars> harddata;

	double med_exp_sec;
	double var_exp_sec;

	Zone()
	{

	};

	Zone(registry *reg, std::string sufix, int i)
	{ // includes the variogram, bhist and hardata


	  // now allocate and read zone-dependent parameters - harddata, bihist, variogram...

		//harddata.reset(new HarddataPars(reg, "HARDDATA" + sufix, i));
		harddata = std::make_shared<HarddataPars>(reg, "HARDDATA" + sufix, i);

		//if (!harddata.get_harddata_pars(reg, "HARDDATA" + sufix, i))
		//	logger << "Hard Data Parameters loaded...\n";
		//else
		//{
		//	logger << "Hard Data Parameters failed to load!\n";
		//	exit(-99);
		//}
		if (!bihist.get_bihist_pars(reg, i))
			logger << "Bihistogram Parameters loaded...\n";
		else
		{
			logger << "Bihistogram Parameters failed to load!\n";
			exit(-99);
		}

		// variogram loading is a good place to start putting ewxceptions and error handling

		//variogram.reset(reg, "VARIOGRAM" + sufix + "Z", i);

		//variogram.reset(new VariogramPars(reg, ("VARIOGRAM" + sufix + "Z"), i));
		variogram=std::make_shared<VariogramPars>(reg, ("VARIOGRAM" + sufix + "Z"), i);
		// end zone parameters


	};


	//now made by constructor
	//int get_zone_pars(registry *reg, std::string sufix, int i)
	//{ // includes the variogram, bhist and hardata
	//	// now allocate and read zone-dependent parameters - harddata, bihist, variogram...
	//	//harddata.reset(new HarddataPars(reg, "HARDDATA" + sufix, i));
	//	harddata = std::make_shared<HarddataPars>(reg, "HARDDATA" + sufix, i);
	//	//if (!harddata.get_harddata_pars(reg, "HARDDATA" + sufix, i))
	//	//	logger << "Hard Data Parameters loaded...\n";
	//	//else
	//	//{
	//	//	logger << "Hard Data Parameters failed to load!\n";
	//	//	return 1;
	//	//}
	//	if (!bihist.get_bihist_pars(reg, i))
	//		logger << "Bihistogram Parameters loaded...\n";
	//	else
	//	{
	//		logger << "Bihistogram Parameters failed to load!\n";
	//		return 1;
	//	}
	//	// variogram loading is a good place to start putting ewxceptions and error handling
	//	//variogram.reset(new VariogramPars(reg, ("VARIOGRAM" + sufix + "Z"), i));
	//	variogram = std::make_shared<VariogramPars>(reg, ("VARIOGRAM" + sufix + "Z"), i);
	//	//variogram(reg, "VARIOGRAM" + sufix + "Z", i);
	//	//if (!variogram->get_variogram_pars(reg, "VARIOGRAM" + sufix + "Z", i))
	//	//	logger << "Variogram Parameters loaded...\n";
	//	//else
	//	//{
	//	//	logger << "Variogram Parameters failed to load!\n";
	//	//	return 1;
	//	//}
	//	// end zone parameters
	//	return 0;
	//};
	//


	int read_zone_data(misc& utils, GSLibGridPars& grid_def)
	{


		harddata->read_harddata(utils);

		if (bihist.bihflag == 1)
		{
			if (!bihist.Loadbihist(utils))
				return 1;

			//if (i == 1)
			//{
			if (!bihist.load_bihist_aux(grid_def, utils.nosvalue, utils.headerflag, utils.filetype))
				return 1;
			//}
			//else
			//{
			//	(bihist[i - 1]).auxiliary_grid = (bihist[0]).auxiliary_grid;
			//}
		} // bihist flag


		//get_global_samples_avg_var();

		return 0;
	};


};




// change zones to a vector of zones, instead of a class with internal vectors 
class Zones_Pars{
public:

	unsigned int n_zones;
	GSLibGrid<int> *zone_grid;
	//VariogramPars *variogram;
	//BihistogramPars *bihist;
	//HarddataPars *harddata;

	Zone* *zone;

	// total average and variance
	double vmedexp_g;
	double vvarexp_g;
	
	//double *med_exp_sec;
	//double *var_exp_sec;

	int nd_g;

	std::string zonefl;

	Zones_Pars()
	{
	};
	
	// another constructor
	Zones_Pars(registry *reg, std::string sufix)
	{
	 // includes the variogram, bhist and hardata

		reg_key *k;


		k = get_key(reg, "ZONES", "NZONES");
		if (k)
			n_zones = get_int(k);
		else 
			exit(-99);

		if (n_zones<1)
		{
			logger << "We must have at least one zone! Setting NZONES to 1...\n";
			n_zones = 1;
		}

		if (n_zones>1)
		{

			if ((k = get_key(reg, (char*) "ZONES", "ZONESFILE")) != NULL)
				zonefl = get_string(k);

			boost::algorithm::trim(zonefl);

			logger << " Zones file: " << zonefl << "\n";
		}

		// now allocate and read zone-dependent parameters - harddata, bihist, variogram...

		zone = new Zone*[n_zones];

		//zone = new *Zone[n_zones];
		
		//zone= new Zone()[n_zones]

		for (unsigned int i = 0; i < n_zones; ++i)
		{
			zone[i] = new Zone(reg, sufix, i + 1);
			//if (!zone[i].get_zone_pars(reg, sufix, i + 1))
			//	logger << "Zone loaded\n";
			//else
			//{
			//	logger << "Zone load failed\n";
			//	exit(-1);
			//}

		}


	};
	
	//copy constructor
	Zones_Pars(const Zones_Pars& zones_to_copy)
	{
		n_zones = zones_to_copy.n_zones;
		zone_grid = zones_to_copy.zone_grid;
		//variogram = zones_to_copy.variogram;
		//bihist = zones_to_copy.bihist;
		//harddata = zones_to_copy.harddata;
		zone = zones_to_copy.zone;


		//med_exp_sec= zones_to_copy.med_exp_sec;
		//var_exp_sec= zones_to_copy.var_exp_sec;

		vmedexp_g = zones_to_copy.vmedexp_g;
		vvarexp_g = zones_to_copy.vvarexp_g;
		nd_g = zones_to_copy.nd_g;

		zonefl = zones_to_copy.zonefl;

	};
	
	//~Zones_Pars()
	//{

	//	delete[] variogram;

	//	delete[] harddata;

	//	delete[]	med_exp_sec;
	//	delete[]	var_exp_sec;
	//	
	//	if (bihist->bihflag ==1){delete[] bihist;}
	//	//delete zone_grid;
	//
	//};



	inline int get_zone_from_index(long index)
	{
		return (int)zone_grid->grid[index];
	};

	//prepare zones

	/*now done on constructor
	int get_zone_pars(registry *reg, std::string sufix)
	{ // includes the variogram, bhist and hardata pars
		
		reg_key *k;

				
		k = get_key(reg,  "ZONES",  "NZONES");
		if (k)
			n_zones = get_int(k);
		else return -1;
	
		if (n_zones<1)
		{
			logger << "We must have at least one zone! Setting NZONES to 1...\n";
			n_zones = 1;
		}

		if (n_zones>1)
		{

			if ((k = get_key(reg, (char*) "ZONES",  "ZONESFILE")) != NULL)
			zonefl= get_string(k);

			boost::algorithm::trim(zonefl);

			logger << " Zones file: "<< zonefl <<"\n";
		}
		
// now allocate and read zone-dependent parameters - harddata, bihist, variogram...

		//harddata=	new HarddataPars[n_zones];
		//bihist =	new BihistogramPars[n_zones];
		//variogram = new VariogramPars[n_zones];

		//med_exp_sec = new double[n_zones];
		//var_exp_sec = new double[n_zones];
		zone = new Zone*[n_zones];

		for (unsigned int i = 0; i < n_zones; ++i)
		{
			zone[i] = new Zone(reg, sufix, i + 1);
			//if (!zone[i].get_zone_pars(reg, sufix, i+1))
			//	logger << "Zone loaded\n";
			//else
			//{
			//	logger << "Zone load failed\n";
			//	exit(-1);
			//}
			//if (!harddata[i-1].get_harddata_pars(reg,"HARDDATA" + sufix,i))
			//	logger << "Hard Data Parameters loaded...\n";
			//else
			//{
			//	logger <<"Hard Data Parameters failed to load!\n";
			//	return 1;
			//}
			//if (!bihist[i-1].get_bihist_pars(reg,i))
			//	logger <<"Bihistogram Parameters loaded...\n";
			//else
			//{
			//	logger <<"Bihistogram Parameters failed to load!\n";
			//	return 1;
			//}

			//if (!variogram[i-1]. get_variogram_pars(reg, "VARIOGRAM" + sufix + "Z",i))
			//	logger <<"Variogram Parameters loaded...\n";
			//else
			//{
			//	logger <<"Variogram Parameters failed to load!\n";
			//	return 1;
			//}
		}

		// end zone parameters
		return 0;
	
	};
	*/

	int read_zone_data(misc& utils, GSLibGridPars& grid_def)
	{

		if (n_zones>1)	
			read_zone_grid(grid_def, utils,  zonefl);
		else 
		{ // we are not using zones, we are always on 1 zone 

			zone_grid= new GSLibGrid<int>(grid_def);

			for(unsigned int l =0 ; l<= grid_def.get_nxyz(); ++l)
				zone_grid->grid[l]=0; //zones start on 0 

		} 

		
		for (unsigned int i = 0; i < n_zones; ++i)
		{
			if (!zone[i]->read_zone_data(utils, grid_def))
				logger << "zone data loaded\n";
			else
			{
				logger << "zone data failed to load\n";
				exit(0);
			}
		
			//harddata[i-1].read_harddata(utils);
		
			//if ( bihist[i-1].bihflag==1)
			//{ 
			//	if (!bihist[i-1].Loadbihist(utils)) 
			//		return 1;
			//	
			//	if (i==1)
			//	{
			//		if (!bihist[i-1].load_bihist_aux(grid_def, utils.nosvalue, utils.headerflag,utils.filetype))
			//			return 1;
			//	}
			//	else
			//	{
			//		(bihist[i-1]).auxiliary_grid = (bihist[0]).auxiliary_grid;
			//	}
			//} // bihist flag
		
		}


		get_global_samples_avg_var();
	
		return 0;
	};



	int get_global_samples_avg_var()
	{
	
		double avg=0;
		double var=0;

		nd_g=0;
		
		for (unsigned int i = 0; i < n_zones; i++)
		{

			avg += zone[i]->harddata->vmedexp * zone[i]->harddata->harddata_distr.n_data;
			var += (zone[i]->harddata->vvarexp + zone[i]->harddata->vmedexp *zone[i]->harddata->vmedexp)* zone[i]->harddata->harddata_distr.n_data ;

			nd_g += zone[i]->harddata->harddata_distr.n_data;
		}

		vmedexp_g= avg/nd_g;
		vvarexp_g = (var/nd_g) - vmedexp_g*vmedexp_g;

		return 0;

	};


	//void get_zonal_avg_var(Zones_Pars & Zones);
	template <class templt> void get_zonal_avg_var(GSLibGrid<templt>& cube_in)
	{

		double av = 0;
		int nxyz = cube_in.get_nxyz();
		int currzone;

		long *n_points_zone = new long[n_zones];
		double *Med_next = new double[n_zones];

		for (unsigned int i = 0; i < n_zones; i++)
		{
			zone[i]->med_exp_sec = 0;
			zone[i]->var_exp_sec = 0;
			n_points_zone[i] = 0;
		}


		for (long i = 1; i <= nxyz; ++i)
		{
			currzone = get_zone_from_index(i);

			//accumulate on correct place
			if ((cube_in.grid[i] != cube_in.null_data) && currzone != -1)
			{
				++n_points_zone[currzone];
				Med_next[currzone] = zone[currzone]->med_exp_sec + (cube_in.grid[i] - zone[currzone]->med_exp_sec) / n_points_zone[currzone];
				zone[currzone]->var_exp_sec = zone[currzone]->var_exp_sec + (cube_in.grid[i] - zone[currzone]->med_exp_sec)*(cube_in.grid[i] - Med_next[currzone]);
				zone[currzone]->med_exp_sec = Med_next[currzone];

			}

		}


		for (unsigned int i = 0; i < n_zones; i++)
		{
			zone[i]->var_exp_sec /= n_points_zone[i];
		}

		return;

	};



	private:
	int read_zone_grid(GSLibGridPars& grid_def, misc& utils, std::string zonefl)
	{
	
		// prepare zoning data grid object
		zone_grid= new GSLibGrid<int>(grid_def);

		zone_grid->filename = new std::string;
		*zone_grid->filename=zonefl; // get filename

		zone_grid->headerflag=utils.headerflag; // check if we use header

		zone_grid->null_data=utils.nosvalue; // set nulls

		logger << " Reading zones datafile \n";

		//zone_grid->read_grid();//reading grid, puts on object.grid

		switch (utils.filetype)
		{
		case GEOEAS: // case 0 is an ascii file
			zone_grid->read_grid();//reading grid, puts on object.grid
			break;
		case SGEMS: // case 1 is sgems binary
			zone_grid->read_grid_binary_sgems();
			break;
		}

		zone_grid ->get_average();
		
		return 0;	
	
	};

};


#endif