#ifndef POINTDISTRIBUTION_INCLUDED
#define POINTDISTRIBUTION_INCLUDED

#include "utils.h"

class pseudoharddata{
public:
	int usepseudo;
	std::string pseudofl;
	int pseudo_corr;

	unsigned int n_pseudo_hard;
//	unsigned int n_distr_pseudo;


	std::vector<double> pseudo_X;
	std::vector<double> pseudo_Y;
	std::vector<double> pseudo_Z;
	

	long  *ind_pseudo;
		
	std::vector<Distribution> pseudo_distr;


	//Distribution* updated_distribution;

	bool did_pseudo;

	GSLibGrid<int> *has_pseudo;

	//empty constructor
	pseudoharddata()
	{
		usepseudo = 0;
		n_pseudo_hard = 0;
	};


	//copy constructor
	pseudoharddata (const pseudoharddata &in_pseudo)
	{
		usepseudo = in_pseudo.usepseudo;
		if (usepseudo == 1) 
		{
			pseudofl = in_pseudo.pseudofl;
			pseudo_corr = in_pseudo.pseudo_corr;
			n_pseudo_hard = in_pseudo.n_pseudo_hard;

			pseudo_X = in_pseudo.pseudo_X;
			pseudo_Y = in_pseudo.pseudo_Y;
			pseudo_Z = in_pseudo.pseudo_Z;

			ind_pseudo = new long[n_pseudo_hard];
			for (int i = 0; i < n_pseudo_hard; ++i)
				ind_pseudo[i] = in_pseudo.ind_pseudo[i];

			pseudo_distr = in_pseudo.pseudo_distr;

			did_pseudo = in_pseudo.did_pseudo;

			has_pseudo = new GSLibGrid<int>(*in_pseudo.has_pseudo);
		}
	};

	// assignment
	pseudoharddata& operator=(const pseudoharddata &in_pseudo)
	{
		usepseudo = in_pseudo.usepseudo;
		pseudofl = in_pseudo.pseudofl;
		pseudo_corr = in_pseudo.pseudo_corr;
		n_pseudo_hard = in_pseudo.n_pseudo_hard;

		pseudo_X = in_pseudo.pseudo_X;
		pseudo_Y = in_pseudo.pseudo_Y;
		pseudo_Z = in_pseudo.pseudo_Z;

		ind_pseudo = new long[n_pseudo_hard];
		for (unsigned int i = 0; i < n_pseudo_hard; ++i)
			ind_pseudo[i] = in_pseudo.ind_pseudo[i];

		pseudo_distr = in_pseudo.pseudo_distr;

		did_pseudo = in_pseudo.did_pseudo;

		has_pseudo = new GSLibGrid<int>(*in_pseudo.has_pseudo);
	
	};


	//destructor
	~pseudoharddata()
	{
		if (usepseudo == 1)
		{
			delete has_pseudo;
			delete[] ind_pseudo;
		}
	};

		
	int read_pseudo_parameters(registry *reg)
	{	
	
		reg_key *k;
		
		// pseudohard data
		k = get_key(reg,  "PSEUDOHARD",  "USEPSEUDO");
		if (k)
			usepseudo = get_int(k);
		else return -1;

		if (usepseudo==1)
		{	
			if ((k = get_key(reg,  "PSEUDOHARD",  "PSEUDOFILE")) != NULL)
			pseudofl = get_string(k);
			boost::algorithm::trim(pseudofl);
			std::cout<<" Pseudo hard data file: "<< pseudofl <<"\n";
		}

		k = get_key(reg,  "PSEUDOHARD",  "PSEUDOCORR");
		if (k)
			pseudo_corr = get_int(k);
		else return -1;


		return 0;
	}

	void read_pseudo_harddata(misc& utils, GSLibGridPars &grid_def){

		static const double epslon = 1.0E-20;
		//unsigned int i,j; // cycle control

		//int ret; // return value

		double twt=0,vrg;
		double cp,oldcp,w;
		std::ifstream fpseudo; // file with pseudo hard data
		float  c[1]; // auxiliar to sort

		double buff;

		// open file and all that jazz
		std::cout << "Loading point distributions";
		//fpseudo.open("test.txt", std::ifstream::in);
		fpseudo.open((const char*)pseudofl.c_str(), std::ifstream::in);
		if (!fpseudo.is_open())
		{
			std::cout << "file open failed!";
			exit(99);
		}
			
		// read the definitions
		//fpseudo >> n_pseudo_hard;
		//fpseudo.is_open();
		std::string line;



		// read the data


		int i = 0 ;
		//std::getline(fpseudo, line);
		while(std::getline(fpseudo, line))
		{

			// insert on 0 to mantain fortran style arrays

			Distribution distr_buff;


			distr_buff.orig_distr->insert(distr_buff.orig_distr->end(),-999.25);
			distr_buff.gauss_distr->insert(distr_buff.gauss_distr->end(),-999.25);
			distr_buff.n_data=0;

			twt=0;

			std::istringstream iss(line);
			
			float bufferx, buffery, bufferz;

			iss >> bufferx >> buffery >> bufferz;

			pseudo_X.push_back(bufferx);
			pseudo_Y.push_back(buffery);
			pseudo_Z.push_back(bufferz);


			while (iss>>buff)
			{

				distr_buff.orig_distr->push_back(buff);
				distr_buff.gauss_distr->push_back(1);
				++distr_buff.n_data;

				twt += 1; //distr_buff.gauss_distr.at(1/*pseudo_distr[i].gauss_distr.end()*/);
				//++j;
			}
			// now sort data by value
			utils.sortit(1,distr_buff.n_data,&(*distr_buff.orig_distr)[0],1,&(*distr_buff.gauss_distr)[0],c,c,c,c,c,c);

			twt   = (double)std::max(twt,epslon);
			oldcp = 0.0;
			cp    = 0.0;

			for (unsigned int j=0; j< distr_buff.n_data; ++j) 
			{

				cp =  cp + (double)((*distr_buff.gauss_distr)[j+1]/twt);
				w  = (cp + oldcp)*0.5;

				vrg = Distribution::gauinv(w);

				oldcp = cp;

				// Now, reset the weight to the normal scores value:
				(*distr_buff.gauss_distr)[j+1] = vrg; 
			}

			// put max and mins here
			distr_buff.zmin=(*distr_buff.orig_distr)[1];
			distr_buff.zmax=(*distr_buff.orig_distr)[distr_buff.n_data];


			pseudo_distr.push_back(distr_buff);
			++i;
			//delete *distr_buff;
		}

		n_pseudo_hard=i;

		ind_pseudo = new long[n_pseudo_hard];


		has_pseudo = new GSLibGrid<int>(grid_def);
		has_pseudo->zero_the_grid();

		for (unsigned int i = 0; i < n_pseudo_hard; ++i) {

			//psudodata index
			cl_coord<double> pseudo_coord(pseudo_X.at(i), pseudo_Y.at(i), pseudo_Z.at(i));
			ind_pseudo[i] = grid_def.get_index_from_xyz_f(pseudo_coord);
			has_pseudo->grid[(int)ind_pseudo[i]] = i+1;
		}

		//set did pseudo to FALSe

		did_pseudo= false;

	};

	void read_pseudo_harddata_binary(misc& utils, GSLibGridPars &grid_def)
	{

		static const double epslon = 1.0E-20;
		//unsigned int i,j; // cycle control

		//int ret; // return value

		double twt = 0, vrg;
		double cp, oldcp, w;
		std::ifstream fpseudo; // file with pseudo hard data
		float  c[1]; // auxiliar to sort

//		double buff;

		// open file and all that jazz

		fpseudo.open(pseudofl.c_str(),std::ios::binary);

		if (!fpseudo.is_open()) {
			logger << pseudofl << " failed to open!";
			return;
		}

		// read the definitions
		//fpseudo >> n_pseudo_hard;

		//std::string line;

		//float *buffer = new float[sizeof(float)*3];//3 collumns of floats coordinates

		double buffer[3];

		// read the data


		int i = 0;



		//fdata.read((char*)&nx, sizeof(unsigned int));

		//double doub_x;
		///*doub_x << */fpseudo >> doub_x;

		while (fpseudo.read((char*)&buffer, sizeof(double)*3))
		{
			// we have a new point, let's findthe rest

			// insert on 0 to mantain fortran style arrays

			Distribution distr_buff;


			distr_buff.orig_distr->push_back(-999.25);
			distr_buff.gauss_distr->push_back(-999.25);
			distr_buff.n_data = 0;

			twt = 0;


			//std::istringstream iss(line);

			//float bufferx, buffery, bufferz;

			//iss >> bufferx >> buffery >> bufferz;


			pseudo_X.push_back(buffer[0]);
			pseudo_Y.push_back(buffer[1]);
			pseudo_Z.push_back(buffer[2]);
			
			int n_pnts;
			fpseudo.read((char*)&distr_buff.n_data, sizeof(int));

			for (int i = 0; i < distr_buff.n_data; ++i)
			{
				float val;
				fpseudo.read((char*)&val, sizeof(float));
				distr_buff.orig_distr->push_back(val);
				distr_buff.gauss_distr->push_back(1);
				//++distr_buff.n_data;
			}

			//float val;
			//fpseudo.read((char*)&val, sizeof(float));

			//while (val!= 0x0A)
			//{



			//	twt = twt + distr_buff.gauss_distr->at(1/*pseudo_distr[i].gauss_distr.end()*/);
			//	//++j;
			//	fpseudo.read((char*)&val, sizeof(float));
			//}

			
			//fpseudo.read((char *)&val, sizeof(char));
			//while (val != 0x0A)
			//{ 
			//return a byte and reread

			//}


			//int n_points_in_distr;

			//fpseudo.read((char*)&n_points_in_distr, sizeof(int));

			//float distrib_val;

			//for (int k = 0; k < distr_buff.n_data; ++k) {

			//	fpseudo.read((char*)&distrib_val, sizeof(float));
			//
			//	distr_buff.orig_distr.insert(distr_buff.orig_distr.end(), distrib_val);
			//	distr_buff.gauss_distr.insert(distr_buff.gauss_distr.end(), 1);
			//	++distr_buff.n_data;

			//	twt = twt + distr_buff.gauss_distr.at(1);
			//	
			//}
			// now sort data by value
			utils.sortit(1, distr_buff.n_data, &(*distr_buff.orig_distr)[0], 1, &(*distr_buff.gauss_distr)[0], c, c, c, c, c, c);

			twt = (double)std::max(twt, epslon);
			oldcp = 0.0;
			cp = 0.0;

			for (unsigned int j = 1; j <= distr_buff.n_data; ++j)
			{

				cp = cp + (double)((*distr_buff.gauss_distr)[j] / twt);
				w = (cp + oldcp)*0.5;

				vrg = Distribution::gauinv(w);

				oldcp = cp;

				// Now, reset the weight to the normal scores value:
				(*distr_buff.gauss_distr)[j] = vrg;
			}

			// put max and mins here
			distr_buff.zmin = (*distr_buff.orig_distr).at(1);
			distr_buff.zmax = (*distr_buff.orig_distr).at(distr_buff.n_data);


			pseudo_distr.insert(pseudo_distr.end(), distr_buff);
			++i;
			//delete *distr_buff;
		}

		n_pseudo_hard = i;



		ind_pseudo = new long[n_pseudo_hard];


		has_pseudo = new GSLibGrid<int>(grid_def);
		has_pseudo->zero_the_grid();

		for (i = 0; i < n_pseudo_hard; ++i) {

			//psudodata index
			cl_coord<double> pseudo_coord(pseudo_X.at(i), pseudo_Y.at(i), pseudo_Z.at(i));
			ind_pseudo[i] = grid_def.get_index_from_xyz_f(pseudo_coord);
			has_pseudo->grid[(long)ind_pseudo[i]] = i + 1;
		}


		//set did pseudo to FALSe

		did_pseudo = false;

	};

	void get_pseudodata_random_path(GSLibGridPars grid_def,Cl_mtrng &outer_RNG, misc &utils,GSLibGrid<unsigned int> &order){

		unsigned int i; // cycle control

		long swape;


		unsigned long nxyz=grid_def.get_nxyz();
		//for(i = 0; i < n_pseudo_hard; ++i) {

		//	//psudodata index
		//	cl_coord<double> pseudo_coord(pseudo_X.at(i), pseudo_Y.at(i),pseudo_Z.at(i));
		//	ind_pseudo[i]=grid_def.get_index_from_xyz_f(pseudo_coord);
		//	
		//}

		for (i=1; i<=nxyz; ++i) {
			order.grid[i] = i;
		}

		//boost::scoped_array<float> aux_rands(new float[n_pseudo_hard]);
		//boost::scoped_array<int> aux_position(new int[n_pseudo_hard]);
		//// make an auxiliary vector with rands
		//// and other with ints

		//for (int i = 0; i < n_pseudo_hard; ++i) 
		//{
		//	aux_rands[i] = outer_RNG.random_float();
		//	aux_position[i] = i;
		//}
		//float c[1];
		//utils.sortit(0, n_pseudo_hard-1, &aux_rands[0], 1, &aux_position[0], c, c, c, c, c, c);
		////now we have the internal order of pseudo points
		////lets add them to the order vector
		//
		////indpseudo has the index of the pseudo dists 

		//for (int i = 0; i < n_pseudo_hard; ++i)
		//{
		//	order.grid[ind_pseudo[aux_position[i]]] = i + 1;
		//	order.grid[i + 1] = ind_pseudo[aux_position[i]];
		//}

		//unsigned long nxyz;

		//nxyz=grid_def.get_nxyz();

		////aux order
		//for (i=1; i<=nxyz; ++i) {
		//	order.grid[i] = i;
		//}
		
	// if the pseudo go in ordered, it always works


		// copy the array
		int *auxpseudo= new int[n_pseudo_hard];
		for (i = 0; i < n_pseudo_hard; i++)
		{
			auxpseudo[i] = ind_pseudo[i];
		}


		// ordeer the array 
		float c[1];
		utils.sortit(0, n_pseudo_hard-1, auxpseudo, 0, c, c, c, c, c, c, c);

		// now swap the ordered 
		for (i = 0; i < n_pseudo_hard; i++)
		{

			swape = order.grid[auxpseudo[i]];
			order.grid[auxpseudo[i]] = order.grid[i + 1];
			order.grid[i + 1] = swape;

		}



		 //now we do auxiliary random vectors  and sort

		//float *temp = new float[nxyz+1];
		boost::scoped_array<float> temp(new float[nxyz + 1]);

		for (i = 1; i <= nxyz; i++)
		{
			temp[i]=outer_RNG.random_float();
		}

		
		utils.sortit(1,n_pseudo_hard,&temp[0], 1,&order.grid[0], c,c,c,c,c,c );
		
		utils.sortit(n_pseudo_hard+1,nxyz,&temp[0], 1,&order.grid[0], c,c,c,c,c,c );
		

		//delete [] temp;

		return;

	} 



};


#endif