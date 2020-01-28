#ifndef HARDDATA_INCLUDED
#define HARDDATA_INCLUDED

// KNOWN BUG - EXTRA LINES IN HARDDATA FILE DOUBLE THE LÇAST POINT: 

#include <vector>
#include <iostream>
#include <boost/algorithm/string/trim.hpp>


#include "utils.h"
#include "gslib_grid.h"
#include "math_random.h"
#include "log.h"

#include "distribution.h"
#include "point_distribution.h"

#include "coordinates.h"



template <class T> class cl_data_point
{
public:

	cl_coord<double> xyz; 

	T value;

	cl_data_point<T>()
	{
		xyz;
		value;
	};

	cl_data_point(cl_coord<double> xyz_in, T value_in)
	{
		xyz= xyz_in;
		value = value_in;
	};

	cl_data_point operator=(const cl_data_point& b)
	{         
		//cl_data_point copied;
		this->xyz= b.xyz;
		this->value=b.value;
        return *this;
	};

};

class GSLIBHeader 
{	
public:
	char *name;
	int  nvars;
	char **colsname;

	GSLIBHeader(std::ifstream& fin)
	{ //constructs from file
		name = new char[512];
		fin >> name;


		fin >> nvars;

		colsname = new char*[nvars + 1];

		for (unsigned int l = 1; l <= nvars; ++l) {

			(colsname)[l] = new char[256];
			fin >> (colsname)[l];
		}

	};



	~GSLIBHeader()
	{
	
		for (unsigned int l = 1; l <= nvars; ++l) {
			delete [] colsname[l];
		}
		delete[] colsname;

		delete[] name;
	}

	GSLIBHeader& operator= (const GSLIBHeader& in_header)
	{
		name = in_header.name;
		nvars = in_header.nvars;

		colsname = new char*[nvars + 1];

		for (unsigned int l = 1; l <= nvars; ++l) {

			(colsname)[l] = new char[256];
			(colsname)[l] = in_header.colsname[l];
		}

	
	}

};



class HarddataPars{
public:

	// input data file name
	std::string datafl;

	//int nd;

	std::vector< cl_data_point<double>> point;

	std::vector<float> wt;
	std::vector<float> vsec;

	Distribution harddata_distr;

	double vmedexp;
	double vvarexp;

	//int ntr; // n after trimmed
	int nt; // n data

	//double zmin, zmax;
	int ltail,utail;
	double utpar,ltpar;

	// harddata pars 
	HarddataPars() {};

	unsigned int nvari; //number of collumns in harddata file

	unsigned int ixl, iyl, izl; // coordinate columns
	int ivrl, iwt; // columns for variable and weight

	double tmin, tmax; // trimming limits 

	int itrans;
	// transformation file name
	std::string transfl;

	GSLIBHeader *file_header;
	//char **colsname; // header columns names

	//std::vector<std::vector<char>> colsname;

	// constructor

	HarddataPars(registry_type *reg, std::string keyword, int sufix) { 
		// sufix tells wich are we loading

		reg_key *k;

		std::string group_id;

		group_id = keyword + utilities::make_string(sufix);

		if ((k = get_key(reg, (char*)group_id.c_str(), "DATAFILE")) != NULL)
			datafl = get_string(k);

		boost::algorithm::trim(datafl);

		std::cout << " Data file: " << datafl << "\n";

		k = get_key(reg, (char*)group_id.c_str(), "COLUMNS");
		if (k)
			nvari = get_int(k);
		else 
			exit(-99);

		k = get_key(reg, (char*)group_id.c_str(), "XCOLUMN");
		if (k)
			ixl = get_int(k);
		else exit(-99);
		if (ixl < 0 || ixl > nvari) {
			std::cout << " Invalid column for XCOLUMN ";
			exit(-99);
		}

		k = get_key(reg, (char*)group_id.c_str(), ("YCOLUMN"));
		if (k)
			iyl = get_int(k);
		else exit(-99);
		if (iyl < 0 || iyl > nvari) {
			std::cout << " Invalid column for YCOLUMN ";
			exit(-99);
		}

		k = get_key(reg, (char*)group_id.c_str(), ("ZCOLUMN"));
		if (k)
			izl = get_int(k);
		else exit(-99);
		if (izl < 0 || izl > nvari) {
			std::cout << " Invalid column for ZCOLUMN ";
			exit(-99);
		}

		k = get_key(reg, (char*)group_id.c_str(), ("VARCOLUMN"));
		if (k)
			ivrl = get_int(k);
		else exit(-99);
		if (ivrl < 0 || ivrl > nvari) {
			std::cout << " Invalid column for VARCOLUMN ";
			exit(-99);
		}

		k = get_key(reg, (char*)group_id.c_str(), ("WTCOLUMN"));
		if (k)
			iwt = get_int(k);
		else exit(-99);
		if (iwt < 0 || iwt > nvari) {
			std::cout << " Invalid column for WTCOLUMN ";
			exit(-99);
		}

		k = get_key(reg, (char*)group_id.c_str(), ("MINVAL"));
		if (k)
			tmin = get_double(k);
		else exit(-99);

		k = get_key(reg, (char*)group_id.c_str(), ("MAXVAL"));
		if (k)
			tmax = get_double(k);
		else exit(-99);

		std::cout << " Trimming  limits: " << tmin << " " << tmax << "\n";

		k = get_key(reg, (char*)group_id.c_str(), ("USETRANS"));
		if (k)
			itrans = get_int(k);
		else exit(-99);

		if ((k = get_key(reg, (char*)group_id.c_str(), ("TRANSFILE"))) != NULL)
		{
			transfl = get_string(k);
			boost::algorithm::trim(transfl);
		}

		std::cout << " Transformation file (1/0):" << transfl << " " << itrans << "\n";


		// now load global information - maybe we must Zonate this as well. TO DO TODO 

		k = get_key(reg, (char*)keyword.c_str(), ("ZMIN"));
		if (k)
			harddata_distr.zmin = get_double(k);
		else exit(-99);

		k = get_key(reg, (char*)keyword.c_str(), ("ZMAX"));
		if (k)
			harddata_distr.zmax = get_double(k);
		else exit(-99);

		std::cout << " Data limits (distribution tails): " << harddata_distr.zmin << " " << harddata_distr.zmax << "\n";


		k = get_key(reg, (char*)keyword.c_str(), ("LTAIL"));
		if (k)
			ltail = get_int(k);
		else exit(-99);

		k = get_key(reg, (char*)keyword.c_str(), ("LTPAR"));
		if (k)
			ltpar = get_double(k);
		else exit(-99);

		k = get_key(reg, (char*)keyword.c_str(), ("UTAIL"));
		if (k)
			utail = get_int(k);
		else exit(-99);

		k = get_key(reg, (char*)keyword.c_str(), ("UTPAR"));
		if (k)
			utpar = get_double(k);
		else exit(-99);

		std::cout << " Lower tail (interpolation type, min.): " << ltail << " " << ltpar << "\n";
		std::cout << " Upper tail (interpolation type, min.): " << utail << " " << utpar << "\n";

		// end harddata


	};

	~HarddataPars()
	{
		vsec.clear();
		wt.clear();
		point.clear();

		//delete transfl;
		//delete datafl;
		delete file_header;
	}

	//removed as pars are now loaded by constructor
/*
	// pars read
	int get_harddata_pars(registry_type *reg, std::string keyword,int sufix){ // sufix tells wich are we loading

		reg_key *k; 

		std::string group_id;

		group_id= keyword+ utilities::make_string(sufix);

		if ((k = get_key(reg, (char*) group_id.c_str(),  "DATAFILE")) != NULL)
			datafl = get_string(k);
		
		boost::algorithm::trim(datafl);
				
		std::cout << " Data file: "<< datafl << "\n";

		k = get_key(reg, (char*)group_id.c_str(),  "COLUMNS");
		if (k)
			nvari = get_int(k);
		else return -1;

		k = get_key(reg, (char*)group_id.c_str(),  "XCOLUMN");
		if (k)
			ixl = get_int(k);
		else return -1; 
		if (ixl < 0 || ixl > nvari) {
			std::cout << " Invalid column for XCOLUMN ";
			return -1;
		}

		k = get_key(reg,(char*)group_id.c_str(),  ("YCOLUMN"));
		if (k)
			iyl = get_int(k);
		else return -1;
		if (iyl < 0 || iyl > nvari) {
			std::cout << " Invalid column for YCOLUMN ";
			return -1;
		}

		k = get_key(reg, (char*)group_id.c_str(),  ("ZCOLUMN"));
		if (k)
			izl = get_int(k);
		else return -1;
		if (izl < 0 || izl > nvari) {
			std::cout << " Invalid column for ZCOLUMN ";
			return -1;
		}

		k = get_key(reg, (char*)group_id.c_str(),  ("VARCOLUMN"));
		if (k)
			ivrl = get_int(k);
		else return -1;
		if (ivrl < 0 || ivrl > nvari) {
			std::cout << " Invalid column for VARCOLUMN ";
			return -1;
		}

		k = get_key(reg, (char*)group_id.c_str(),  ("WTCOLUMN"));
		if (k)
			iwt = get_int(k);
		else return -1;
		if (iwt < 0 || iwt > nvari) {
			std::cout << " Invalid column for WTCOLUMN ";
			return -1;
		}

		k = get_key(reg, (char*)group_id.c_str(),  ("MINVAL"));
		if (k)
			tmin = get_double(k);
		else return -1;

		k = get_key(reg, (char*)group_id.c_str(),  ("MAXVAL"));
		if (k)
			tmax = get_double(k);
		else return -1;

		std::cout << " Trimming  limits: "<< tmin << " " << tmax << "\n";

		k = get_key(reg, (char*)group_id.c_str(),  ("USETRANS"));
		if (k)
			itrans = get_int(k);
		else return -1;

		if ((k = get_key(reg, (char*)group_id.c_str(),  ("TRANSFILE"))) != NULL)
		{
			transfl = get_string(k);
			boost::algorithm::trim(transfl);
		}

		std::cout << " Transformation file (1/0):" << transfl <<" "<< itrans << "\n";


		// now load global information - maybe we must Zonate this as well. TO DO TODO 

		k = get_key(reg,  (char*)keyword.c_str(),  ("ZMIN"));
		if (k)
			harddata_distr.zmin = get_double(k);
		else return -1;

		k = get_key(reg,  (char*)keyword.c_str(),  ("ZMAX"));
		if (k)
			harddata_distr.zmax = get_double(k);
		else return -1;

		std::cout << " Data limits (distribution tails): "<< harddata_distr.zmin << " "<< harddata_distr.zmax << "\n";


		k = get_key(reg,  (char*)keyword.c_str(),  ("LTAIL"));
		if (k)
			ltail = get_int(k);
		else return -1;

		k = get_key(reg,  (char*)keyword.c_str(),  ("LTPAR"));
		if (k)
			ltpar = get_double(k);
		else return -1;

		k = get_key(reg,  (char*)keyword.c_str(),  ("UTAIL"));
		if (k)
			utail = get_int(k);
		else return -1;

		k = get_key(reg,  (char*)keyword.c_str(),  ("UTPAR"));
		if (k)
			utpar = get_double(k);
		else return -1;

		std::cout << " Lower tail (interpolation type, min.): "<< ltail << " "<< ltpar << "\n";
		std::cout << " Upper tail (interpolation type, min.): "<< utail << " "<< utpar << "\n";

		// end harddata


		return 0;

	};
*/
	
int read_harddata(misc& utils)
	{
		unsigned int ntr;

		static const double epslon = 1.0E-20;
		// local variables
				
		std::ofstream ofdata;
		std::ifstream fdata;

		//int r; // return value 

		unsigned int j; //cycle control

		double cp,oldcp,w;
		float  vrg,c[1];
		double twt;
		double av, ss; // average and 


		long   icolvr,icolwt;

		std::string tmpfl;
		//char *title= new char[512];//;// 512 characters max for header title
		
		//double *var = new double[nvari];

		// end local variables
		std::cout << "READ HARD DATAFILE \n";


		fdata.open(datafl.c_str());

		if (fdata.fail())
		{ 
			std::cout << "Failed to open datafile \n";
			exit(-1);
		}
		

		fdata.close();

		harddata_distr.n_data = 0;
		av = (double)0.0;
		ss = (double)0.0;

		if(itrans == 1)
		{
			std::cout << " Setting up transformation table \n";
			
			tmpfl=datafl;
			icolvr = ivrl;
			icolwt = iwt;

			// Open up the file with reference distribution:
			fdata.open(tmpfl.c_str());//fdata=fopen (tmpfl,"r");
			if (fdata.fail()) {
				std::cout << "ERROR: " << tmpfl << " does not exist! \n";
				std::cout << "     this file is needed! \n";
				//				getch();
				exit (0);
			}

			// Now, read in the actual data:
			nt  = 0;
			ntr = 0;
			twt = 0.0;

			if (utils.headerflag == 1)
				file_header = new GSLIBHeader(fdata);
				
				//if(!ReadGSLibHeader(fdata, title ,&nvari/*, &colsname*/))
				//{
				//	std::cout << " ERROR reading gslib header! \n";
				//	exit(0);
				//};

			// insert on 0 to mantain fortran style arrays
			harddata_distr.orig_distr->insert(harddata_distr.orig_distr->end(),-999.25); 
			harddata_distr.gauss_distr->insert(harddata_distr.gauss_distr->end(),-999.25);

			double *var = new double[file_header->nvars];
			while (!(fdata.eof())) {

				if (fdata.eof()) 
					break;
				
				for (j=0; j<file_header->nvars; ++j)
					fdata >> var[j];

				// Trim this data?
				if (var[icolvr-1] < tmin || var[icolvr-1] > tmax) {
					++nt;
					continue;
				}
				++ntr;

				if (icolvr > file_header->nvars || icolwt > file_header->nvars) {
					std::cout << " ERROR: too few columns in ref data \n";
					exit (0);
				}

				// Keep this data: Assign the data value and coordinate location:
				harddata_distr.orig_distr->insert(harddata_distr.orig_distr->end(),(var[icolvr-1]));
				if(icolwt <= 0) {
					harddata_distr.gauss_distr->insert(harddata_distr.gauss_distr->end(),1.0);

				}
				else {
					harddata_distr.gauss_distr->insert(harddata_distr.gauss_distr->end(),var[icolwt-1]);
				}
				if ((*harddata_distr.gauss_distr)[ntr] <= 0.0) {
					ntr = ntr - 1;
					nt  = nt  + 1;
					continue;
				}
				twt = twt + (*harddata_distr.gauss_distr)[ntr];

				// Go back for another datum:
			}

			delete[] var;

			fdata.close();//fclose (fdata);
			if(ntr <= 1) {
				std::cout << "ERROR: too few data for transformation \n";
				exit (0);
			}

			// Write transformation table:
			ofdata.open(transfl.c_str());

			//ofdata=fopen (transfl,"w");

			// Sort data by value:
			utils.sortit(1,ntr,&(*harddata_distr.orig_distr)[0],1,&(*harddata_distr.gauss_distr)[0],c,c,c,c,c,c);

			// now update the min and max, to deal with the limits problem 

			//if(tmin < harddata_distr.orig_distr[1]+harddata_distr.orig_distr[1]*0.01) tmin = harddata_distr.orig_distr[1]+harddata_distr.orig_distr[1]*0.01;
			//if(tmax > harddata_distr.orig_distr[ntr]-harddata_distr.orig_distr[ntr]*0.01)tmax = harddata_distr.orig_distr[ntr]-harddata_distr.orig_distr[ntr]*0.01;
			//if(harddata_distr.zmin < harddata_distr.orig_distr[1]+harddata_distr.orig_distr[1]*0.01) harddata_distr.zmin = harddata_distr.orig_distr[1]+harddata_distr.orig_distr[1]*0.01;
			//if(harddata_distr.zmax > harddata_distr.orig_distr[ntr]-harddata_distr.orig_distr[ntr]*0.01)harddata_distr.zmax = harddata_distr.orig_distr[ntr]-harddata_distr.orig_distr[ntr]*0.01;

			//harddata_distr.orig_distr[1]=harddata_distr.orig_distr[1]-harddata_distr.orig_distr[1]*0.03;
			//harddata_distr.orig_distr[ntr]=harddata_distr.orig_distr[ntr]+harddata_distr.orig_distr[ntr]*0.03;
			// Compute the cumulative probabilities and write transformation table
			twt   = (double)std::max(twt,epslon);
			oldcp = 0.0;
			cp    = 0.0;
			for (j=1; j<=ntr; ++j) 
			{
				cp =  cp + (double)((*harddata_distr.gauss_distr)[j]/twt);
				w  = (cp + oldcp)*0.5;

				vrg = (float)Distribution::gauinv(w);
				//	vrg = utils.nosvalue;
				//}

				ofdata << (*harddata_distr.orig_distr)[j] << vrg;

				oldcp = cp;

				// Now, reset the weight to the normal scores value:
				(*harddata_distr.gauss_distr)[j] = vrg;
			}

			ofdata.close();
			//fclose(ofdata);
		}

		// Now, read the data if the file exists:

		fdata.open(datafl.c_str()); 
		if (fdata.fail()) 
			return 1;

		std::cout << " Reading input datafile \n";

		if (ixl > file_header->nvars || iyl > file_header->nvars || izl > file_header->nvars || ivrl > file_header->nvars || iwt> file_header->nvars) {
			std::cout << " ERROR: you have asked for a column number \n";
			std::cout << "        greater than available in file \n";
			exit (0);
		}

		// Read all the data until the end of the file:
		harddata_distr.n_data  = 0;
		nt  = 0;
		twt = 0.0;

		// insert on 0 position to keep fortran indexes

		wt.insert(wt.end(),-999.25);
		//vsec.insert(vsec.end(),-999.25);

		cl_coord<double> dummy_xyz(-999.25,-999.25,-999.25);
		double dummy_float = -999.25;
		point.insert(point.end(),cl_data_point<double>(dummy_xyz,dummy_float) ) ;

		
		if (utils.headerflag==1) 
			file_header = new GSLIBHeader(fdata);
			//if(!ReadGSLibHeader(fdata, title ,&nvari/*, &colsname*/))
			//{
			//	std::cout << " ERROR reading gslib header! \n";
			//	exit(0);
			//};

		//double M_old = 0;
		//double M_curr = 0;
		//double S_curr = 0;

		double wSum = 0;
		double wSum2 = 0;
		double mean = 0;
		double meanOld = 0;
		double S = 0;

		double *var = new double[file_header->nvars];
		//while (!(fdata.eof())) 
		//{
		while (fdata >> var[0] >> var[1] >> var[2] >> var[3])
		{
			//if (fdata.eof()) 
			//	break;


			//for (j=0; j<file_header->nvars; ++j)
			//	fdata >> var[j];


			if (var[ivrl-1] < tmin || var[ivrl-1] > tmax) 
			{
				nt = nt + 1;
				continue;
			}

			++harddata_distr.n_data;


			// Acceptable data, assign the value, X, Y, Z coordinates, and weight:
			cl_coord<double> xyz_point(var[ixl-1],var[iyl-1],var[izl-1]);
			point.insert(point.end(),cl_data_point<double>(xyz_point,var[ivrl-1]) ) ;

			if(iwt <= 0) 
				wt.insert(wt.end(),1.0);
			else 	
				wt.insert(wt.end(), var[iwt-1]);
			

			//if(isecvr <= 0) {
			//	vsec.insert(vsec.end(),1.0);
			//}
			//else {
			//	vsec.insert(vsec.end(),var[isecvr-1]);
			//}

			//twt += wt[harddata_distr.n_data];
			//av  += var[ivrl-1]*wt[harddata_distr.n_data];
			//ss  += var[ivrl-1]*var[ivrl-1]*wt[harddata_distr.n_data];
		
			////no weights welford
			//M_old = M_curr;
			//M_curr += (var[ivrl - 1] - M_curr) / harddata_distr.n_data;
			//S_curr += (var[ivrl - 1] - M_curr)*(var[ivrl - 1] - M_old);


			// weighted welford, as per wikipedia
			wSum += wt[harddata_distr.n_data];
			wSum2 += std::pow(wt[harddata_distr.n_data] , 2);
			meanOld = mean;
			mean = meanOld + (wt[harddata_distr.n_data] / wSum) * (var[ivrl - 1] - meanOld);
			S += wt[harddata_distr.n_data] * (var[ivrl - 1] - meanOld) * (var[ivrl - 1] - mean);
		
		}

		delete[] var;
		fdata.close();

		// Compute the averages and variances as an error check for the user:
		//av /= (double)std::max(twt,epslon);
		//ss = (ss / (double)std::max(twt,epslon)) - av * av;
		//vmedexp = (double)av;
		//vvarexp = (double)ss;

		//double popvar = S / (wSum-1);
		vmedexp = mean;
		vvarexp = S / (wSum - 1);


		//delete[] title;
		//delete[] var;

		
		std::cout << "   acceptable data  = " << harddata_distr.n_data << "\n";
		std::cout << "    number trimmed  = " << nt << "\n";
		std::cout << "  weighted average  = " << mean << "\n";
		std::cout << " weighted variance  = " << vvarexp << "\n\n";


		return 0;
	};

	void Dtonode (GSLibGridPars& grid_def,misc& utils, boost::shared_array<float> &sim, GSLibGrid<int>& zones, int current_zone, int n_zones)
	{
		
		long  ind;
		
		int warning=0;

		// Check all samples
		long nexdata=0;
		for (unsigned long id=1; id<point.size(); ++id) 
		{
			// Check if samples are within study area
			//point[id].xyz.x < grid_def.get_ox() ||  
				//point[id].xyz.y < grid_def.get_oy() || 
				//point[id].xyz.z < grid_def.get_oz() ||
				//point[id].xyz.x > grid_def.get_ox()+grid_def.get_nx()*grid_def.get_xsiz() ||  
				//point[id].xyz.y > grid_def.get_oy()+grid_def.get_ny()*grid_def.get_ysiz() || 
				//point[id].xyz.z > grid_def.get_oz()+grid_def.get_nz()*grid_def.get_zsiz() 
				//)
			if (!grid_def.is_xyz_inside(point[id].xyz))
			{
				++nexdata;
				continue;
			}

			cl_coord<int> sample_indexes=grid_def.get_ijk_form_xyz_f ( point[id].xyz);

			ind = grid_def.get_index_from_xyz_f ( point[id].xyz);

			cl_coord<double> xxyyzz= grid_def.get_xyz_from_ijk_f(sample_indexes);

			double dist = (double)fabs(xxyyzz.x- point[id].xyz.x) + (double)fabs(xxyyzz.y- point[id].xyz.y) + (double)fabs(xxyyzz.z- point[id].xyz.z);
			// Assign this data to the node (unless there is a closer data)
			if (sim[ind] != utils.nosvalue) 
			{
				int id2 = (int)sim[ind];
				

				double dist2 = (double)fabs(xxyyzz.x- point[id2].xyz.x) + (double)fabs(xxyyzz.y- point[id2].xyz.y) + (double)fabs(xxyyzz.z- point[id2].xyz.z);
				
				if (dist <= dist2) sim[ind] = (float)id;

				if ((current_zone == zones.grid[ind]) || (n_zones <= 1))
					sim[ind] = (float)id;

				if (warning == 0) 
				{
					warning=1;
					logger << " WARNING: data values at \n" <<  "  "<< utilities::make_string(point[id].xyz.x)<< ","<< utilities::make_string(point[id].xyz.y)<< ","<< utilities::make_string(point[id].xyz.z) << " - "<< utilities::make_string(point[id2].xyz.x)<< ","<< utilities::make_string(point[id2].xyz.y)<< ","<< utilities::make_string(point[id2].xyz.z)<< "\n";
				}else
					logger <<"  " << utilities::make_string(point[id].xyz.x)<< ","<< utilities::make_string(point[id].xyz.y)<< ","<< utilities::make_string(point[id].xyz.z) << " - "<< utilities::make_string(point[id2].xyz.x)<< ","<< utilities::make_string(point[id2].xyz.y)<< ","<< utilities::make_string(point[id2].xyz.z)<<"\n";
			
			}else
				if ((current_zone == zones.grid[ind])||(n_zones<=1))
					sim[ind] = (float)id; 

			//sim[ind]=vr[id];
			
				//sim[ind]=point[id].value;
		}
		
		


		if (warning == 1)
		{
			logger << " are both assigned to same node - taking closest: \n";
			//logger << ;
		}

		if (nexdata > 0)
			logger <<" Samples outside area: "<< utilities::make_string(nexdata) <<"\n\n";

		// Now, enter data values into the simulated grid
		//for (long i=1; i<=grid_def.get_nxyz(); ++i) {
		//	int id = (int)sim[i];
		//	
		//	if (id > 0)	sim[i] = point[id].value;

		//}

		return;
	};


	float check_limits(float in_value) 
	{
	};

private:

	//int ReadGSLibHeader(std::ifstream& fin, char *name, unsigned int *nvars)
	//{
	//	
	//	fin >> name;
	//	fin >> *nvars;


	//	colsname = new char*[*nvars+1];

	//	for (unsigned int l = 1; l<= *nvars;++l){

	//		//colsname.push_back()
	//		(colsname)[l] = new char[256];
	//		fin >> (colsname)[l];

	//	}

	//	return 1; // 1 ok
	//};

};


#endif