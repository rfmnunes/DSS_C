#ifndef BIHIST_INCLUDED
	#include "bihist.h"
#endif

	// copy constructor
	//BihistogramPars::BihistogramPars(registry_type *reg, int sufix,misc& utils, GSLibGridPars& grid_def, int headerflag ){
	//// just call the other functions
	//	get_bihist_pars(reg,sufix);
	//	load_bihist_aux(grid_def,utils.nosvalue,headerflag);
	//	Loadbihist(utils);
	//};

//BihistogramPars::BihistogramPars(registry_type *, int, misc &, GSLibGridPars &, int)
//{
//}

	BihistogramPars::~BihistogramPars()
	{
		delete[] orig_distr1;
		delete[] orig_distr2;
		delete[] n_points_in_bihst_class;
		//delete[] curr_mean_in_bihst_class;
		//delete[] curr_var_in_bihst_class;
		delete[] M_in_bihst_class;
		delete[] S_in_bihst_class;
		delete[] Mnext_in_bihst_class;
		delete[] auxiliary_grid;
	};

	int  BihistogramPars::get_bihist_pars(registry_type *reg, int sufix)
		// sufix tells which are we loading
	{ 

		reg_key *k;
		int i;

		//parse the bihistogram section
				
		std::string group_id;
		group_id ="BIHIST" + utilities::make_string(sufix);

		//buffer = NULL;
		k = get_key(reg, (char*)group_id.c_str(), "USEBIHIST");
		if (k)
			bihflag = get_int(k);
		else return -1;
		logger << " Using Bihistogram? (0)No /(1)Yes: " << utilities::make_string(bihflag)<<"\n";

		if (bihflag==1){
			
			if ((k = get_key(reg, (char*)group_id.c_str(), "BIHISTFILE")) != NULL)
				bihistfl= get_string(k);
			boost::algorithm::trim(bihistfl);
						
			logger<<" Bihistogram file: " << bihistfl << "\n";


			k = get_key(reg, (char*)group_id.c_str(), "NCLASSES");
			if (k)
				imclas = get_int(k);
			else return -1;
			logger << " Number of classes for bihistogram: "<<utilities::make_string(imclas) <<"\n";

			if ((k = get_key(reg, (char*)group_id.c_str(), "AUXILIARYFILE")) != NULL)
				auxbihistfl= get_string(k);
			boost::algorithm::trim(auxbihistfl);
			
			
			logger<<" Auxiliary bihistogram file: " << auxbihistfl << "\n";

			//allocate bihist stuff

			target_mean_bihst_class= new double[imclas+1+1];//+1 put on 03-06-13
			target_var_bihst_class= new double[imclas+1+1];
			n_points_in_bihst_class=new long[imclas+1+1];
			//curr_mean_in_bihst_class=new double[imclas+1+1];
			//curr_var_in_bihst_class=new double[imclas+1+1];
			M_in_bihst_class = new double[imclas + 1 + 1];
			S_in_bihst_class = new double[imclas + 1 + 1];
			Mnext_in_bihst_class = new double[imclas + 1 + 1];

			// and initialize stuff:

			for(i = 0; i<=imclas+1; i++){

				target_mean_bihst_class[i]=0;
				target_var_bihst_class[i]=0;
				n_points_in_bihst_class[i]=0;
				//curr_mean_in_bihst_class[i]=0;
				//curr_var_in_bihst_class[i]=0;
				M_in_bihst_class[i] = 0;
				S_in_bihst_class[i] = 0;
				Mnext_in_bihst_class[i] = 0;
			}

		}
		// end bihistogram


		return 0;

	};

	int BihistogramPars::load_bihist_aux(GSLibGridPars& grid_def, float nosvalue, int headerflag,int filetype)
	{

		auxiliary_grid = new GSLibGrid<float>(grid_def);


		// get filename
		auxiliary_grid->filename = new std::string;
		*auxiliary_grid->filename=auxbihistfl;

		auxiliary_grid->headerflag=headerflag; // check if we use header

		auxiliary_grid->null_data=nosvalue; // set nulls

		// see if it exists
		if (!boost::filesystem::exists(bihistfl.c_str()))
		{
			logger << "Bihistogram auxiliary file " << *auxiliary_grid->filename << " does not exist!";
			exit(-9);
		}

		logger << " Reading bihist auxiliary datafile \n";

		switch (filetype)
		{
		case GEOEAS: // case 0 is an ascii file
			auxiliary_grid->read_grid();//reading grid, puts on object.grid
			break;
		case SGEMS: // case 1 is sgems binary
			auxiliary_grid->read_grid_binary_sgems();
			break;
		}

		auxiliary_grid->get_average();
		auxiliary_grid->get_variance();

		
		logger << "   acceptable data  = " << (*auxiliary_grid).valid_points << "\n";
		logger << "  weighted average  = " << (*auxiliary_grid).average << "\n";
		logger << " weighted variance  = " << (*auxiliary_grid).variance << "\n\n";

		return 1;

	};

	int BihistogramPars::Loadbihist(misc& utils)
	{
		if (!boost::filesystem::exists(bihistfl.c_str())) 
		{
			logger << "Bihistogram file " << bihistfl.c_str() << " does not exist!";
			exit(-9);
		}
		//Read the bidistribution: two columns of data. The first one is ordered

		std::ifstream bhdata(bihistfl.c_str());

		std::vector<float> *buffer;
 		buffer = new std::vector<float>[2];


		// open bihist file, check if exists
		//if (bhdata==NULL){
		//	logger << " Bihistogram file doesnt exist! \n";
		//	return 0;
		//}

		logger << " Reading bihistogram file \n";


		// read line while it exists
		//while ( (bhdata >> buffer[0].push_back()) && 
		//		(bhdata >> buffer[1][ntbi+1]) && 
		//		!(bhdata.eof())) ntbi++ ; 

		float _one, _two;

		while (bhdata >> _one && bhdata >> _two) {
			buffer[0].push_back(_one);
			buffer[1].push_back(_two);
		};
		bhdata.close();

		// now allocate bihist stuff and put data there

		orig_distr1 = new float[buffer[0].size() + 1];
		orig_distr2 = new float[buffer[0].size() + 1];

		ntbi = (unsigned int) buffer[0].size();

		orig_distr1[0] = (float)-999.99;
		orig_distr2[0] = (float)-999.99;

		for (int j = 1; j <= ntbi; ++j) {
			orig_distr1[j] = buffer[0][j - 1];
			orig_distr2[j] = buffer[1][j - 1];
		}


		float c[1];
		utils.sortit(1, ntbi, orig_distr1, 1, orig_distr2, c, c, c, c, c, c);


		//dealloc buffer

		//TO DO CHECK THIS DEALLOC

		//for (int j=0; j<2; ++j) delete [] buffer[j];
		//delete[] buffer;


		//// CALCULATE target averages and varianes per class


		//float dclaslocal = (orig_distr1[ntbi] - orig_distr1[1]) / imclas;
	
		//int nclasses = 0;
		//int n_p_classes = 0;
		//for (int i = 1; i < ntbi;++i)
		//{
		//	if (orig_distr1[i] <= orig_distr1[1] + dclaslocal*nclasses) 
		//	{
		//		target_mean_bihst_class[nclasses] += orig_distr1[i];
		//		target_var_bihst_class[nclasses] += pow(orig_distr1[i], (int)2);
		//		++n_p_classes;

		//	}

		//	else 
		//	{
		//		target_mean_bihst_class[nclasses] /= n_p_classes;
		//		target_var_bihst_class[nclasses] = target_var_bihst_class[nclasses] / n_p_classes - pow(target_mean_bihst_class[nclasses], 2);
		//		if (target_var_bihst_class[nclasses] < 0) target_var_bihst_class[nclasses] = 0;
		//		++nclasses;
		//		n_p_classes = 0;

		//		target_mean_bihst_class[nclasses] += orig_distr1[i];
		//		target_var_bihst_class[nclasses] += pow(orig_distr1[i], (int)2);
		//		++n_p_classes;
		//	
		//	}			
		//
		//}

		//target_mean_bihst_class[nclasses] /= n_p_classes;
		//target_var_bihst_class[nclasses] -= pow(target_mean_bihst_class[nclasses], 2) / n_p_classes;


		//

		logger << " Bihistogram file loaded \n";

		dclaslocal = (orig_distr1[ntbi] - orig_distr1[1]) / imclas;
		
		return 1;
	};

	int	BihistogramPars::get_current_bihist_class(const float aux)
	{
		//float dclaslocal = (orig_distr1[ntbi] - orig_distr1[1]) / imclas; // espacamento entre classes
		int current_bihst_class = (int)((aux - orig_distr1[1]) / dclaslocal) + 1;
		if (current_bihst_class < 1)
			current_bihst_class = 1;
		
		else if (current_bihst_class > imclas + 1)
			current_bihst_class = imclas + 1;
		
		return current_bihst_class;
	};

	/*Distribution*/void BihistogramPars::Makelocaldist(misc& utils, float auxiliary, double *cmean, double *cstdev, int icmean, int icvar,Cl_mtrng& inner_RNG, Distribution &bih_local_distr){

		//create local distribution from bihist law

		int kkinit,
			kkend;

		boost::mutex critical_sim;
	
		float b[1];

		double  vrg;
		double w;


		//we have the cmean and cstdev, lets see in which class we are

		//float cpdev;
		//Distribution bih_local_distr;


		//float dclaslocal = (orig_distr1[ntbi]-orig_distr1[1])/imclas; // espacamento entre classes
		//int current_bihst_class = (int) ((auxiliary-orig_distr1[1])/dclaslocal)+1;
		//if (current_bihst_class<1 )
		//{ 
		//	current_bihst_class =1;
		//}
		//else if (current_bihst_class>imclas +1 ) 
		//{
		//	current_bihst_class =imclas +1;
		//} 
	
		int current_bihst_class = get_current_bihist_class(auxiliary);
		//critical_sim.lock(); 
		

		//if (n_points_in_bihst_class[current_bihst_class]>0) 
		//{
		//	if (icmean == 1)
		//	{

		//		//*cmean-= (float)((curr_mean_in_bihst_class[current_bihst_class]/n_points_in_bihst_class[current_bihst_class])-target_mean_bihst_class[current_bihst_class]);
		//		*cmean = utils.average_correction(*cmean, M_in_bihst_class[current_bihst_class], target_mean_bihst_class[current_bihst_class]);

		//	}
		//	if (icvar == 1) 
		//	{
		//		//if (target_var_bihst_class[current_bihst_class]>0)
		//		*cstdev = utils.variance_correction(*cstdev, S_in_bihst_class[current_bihst_class]/ n_points_in_bihst_class[current_bihst_class], target_var_bihst_class[current_bihst_class]);
		//	}

		//}

			
	

		//critical_sim.unlock();
	//	}


		// weird acorni call, TO DO
		int imclas2= ntbi/imclas;

		double ran_1 = inner_RNG.random_double();

		int randnum =(int) (ran_1*imclas2 /*+ imclas2*/);

		bih_local_distr.n_data=0;

		for (int k=1; k<=ntbi; ++k)
		{

			
			if (auxiliary<=orig_distr1[k]) // find the nearest in the bihist
			{

				kkinit = k-randnum;
				kkend= k + randnum;
				if (kkinit<1) 
					kkinit=1;
				if (kkend>ntbi) 
					kkend=ntbi;

				// check limits 
				for (int j = kkinit; j <= k; j++)
				{
					if (((auxiliary - dclaslocal) > orig_distr1[j]))
						kkinit = j;
					else
						break;
				}

				for (int j = k; j <= kkend; j++)
				{
					if (((auxiliary + dclaslocal) < orig_distr1[j]))
					{
						kkend = j;
						break;
					}
				}
				//if ((auxiliary - dclaslocal) > orig_distr1[kkinit]) // has less than one class of difference

				bih_local_distr.zmin= 9999999;
				bih_local_distr.zmax=-9999999;

				//long	sum_l_11 = 0;
				//double zmed_l_11 = 0;
				//double zvar_l_11 = 0;

				double zMnext = 0;
				double zS = 0;
				double zM = 0;

				bih_local_distr.orig_distr->push_back(-999.25);
				bih_local_distr.gauss_distr->push_back(-999.25);

				for (int j=kkinit; j<=kkend;j++)
				{

					//if (((auxiliary - dclaslocal) > orig_distr1[j])
					//	||
					//	((auxiliary + dclaslocal) < orig_distr1[j]))
					//	continue;


					++bih_local_distr.n_data;		



					//curr_val_bih=orig_distr2[j]
					bih_local_distr.orig_distr->push_back(orig_distr2[j]);
					
					if (orig_distr2[j] > bih_local_distr.zmax)
						bih_local_distr.zmax = orig_distr2[j];
					if (orig_distr2[j] < bih_local_distr.zmin)
						bih_local_distr.zmin = orig_distr2[j];

					// lito bihist - update max and min
					//if ((*bih_local_distr.orig_distr)[bih_local_distr.n_data]>bih_local_distr.zmax) 
					//	bih_local_distr.zmax=(*bih_local_distr.orig_distr)[bih_local_distr.n_data];
					//if ((*bih_local_distr.orig_distr)[bih_local_distr.n_data]<bih_local_distr.zmin) 
					//	bih_local_distr.zmin=(*bih_local_distr.orig_distr)[bih_local_distr.n_data];

					zMnext += (orig_distr2[j] - zM) / bih_local_distr.n_data;
					zS += (orig_distr2[j] - zM)*(orig_distr2[j] - zMnext);
					zM = zMnext;


					//zmed_l_11 += bih_local_distr.orig_distr[bih_local_distr.n_data];
					//zvar_l_11 += pow((float) bih_local_distr.orig_distr[bih_local_distr.n_data], (int) 2);
					//++sum_l_11;
										
				}
				
				if(bih_local_distr.n_data >0)
				{
				//target_mean_bihst_class[current_bihst_class]=zmed_l_11/sum_l_11;
				//target_var_bihst_class[current_bihst_class]=(zvar_l_11-pow(zmed_l_11,2)/sum_l_11)/sum_l_11;
					target_mean_bihst_class[current_bihst_class] = zM;
					target_var_bihst_class[current_bihst_class]=zS/ bih_local_distr.n_data;
				}
				//double var = bih_local_distr.get_distribution_variance();
				break;

			}//if (auxiliary<=orig_distr1[k])

		}//for (int k=1; k<=ntbi; ++k)

		if (n_points_in_bihst_class[current_bihst_class]>0)
		{
			if (icmean == 1)
			{
				//*cmean-= (float)((curr_mean_in_bihst_class[current_bihst_class]/n_points_in_bihst_class[current_bihst_class])-target_mean_bihst_class[current_bihst_class]);
				*cmean = utils.average_correction(*cmean, M_in_bihst_class[current_bihst_class], target_mean_bihst_class[current_bihst_class]);

			}
			if (icvar == 1)
			{
				if (target_var_bihst_class[current_bihst_class]>0 && S_in_bihst_class[current_bihst_class]>0)
					*cstdev = utils.variance_correction(*cstdev, S_in_bihst_class[current_bihst_class] / std::max((int)(n_points_in_bihst_class[current_bihst_class] - 1), 1), target_var_bihst_class[current_bihst_class]);
			}

		}


		int istart=1;
		int iend= bih_local_distr.n_data;

		if (iend==0) 
		{
			std::cout << "LOCAL HISTOGRAM FROM BIHIST HAS ZERO POINTS. AUXILIARY VALUES MUST BE CONTAINED IN BIHISTOGRAM!";
		}
		
		utils.sortit(istart,iend,&(*bih_local_distr.orig_distr)[0],0,b,b,b,b,b,b,b);

		// compute cumulatives and write transform table

		float cp = 0.0;
		
		for(int j=istart;j<iend+1;j++)
		{

			cp = cp+1;
			w = (double) ((cp-0.5)/ bih_local_distr.n_data);

			vrg = Distribution::gauinv(w);

			// reset the weight to the normal scores value
			bih_local_distr.gauss_distr->push_back((float) vrg);


		}

		return /*bih_local_distr*/;

	};

	float BihistogramPars::bihistlimits(double& newsim, double zmin, double zmax , double zmin_local,double zmax_local ){

		float quartz=(zmax-zmin)/imclas;// class size

		//if (newsim<zmin) newsim=zmin;

		if (newsim<= zmin_local)
		{
			newsim = zmin_local;

			float dzmin=zmin_local-zmin; //diff from local to global min
			
			if (dzmin>0)
			{
				float rzmin=(zmin_local- newsim)/dzmin;
				newsim=zmin_local-rzmin*quartz;
				if(newsim<zmin) newsim=zmin;
			}
		}


		//if (newsim>zmax)newsim=zmax;

		if (newsim>= zmax_local)
		{
			newsim = zmax_local;
			float dzmax=zmax-zmax_local;

			if (dzmax>0)
			{
				float rzmax=(newsim-zmax_local)/dzmax;
				newsim=zmax_local-rzmax*quartz;
				if(newsim>zmax) newsim=zmax;
			}

		}

		return newsim;

	};

