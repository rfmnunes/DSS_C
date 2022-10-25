#ifndef VARIOGRAM_INCLUDED
#define VARIOGRAM_INCLUDED

#define PI 3.141592654
//#define DEG2RAD    PI/180.0
#define TINY                    0.00001
#define EPSLON                 1.0E-20

#include "gslib_grid.h"
#include "utils.h"
#include "log.h"

#include "coordinates.h"



class variogram_structure 
{
public:

	int	    model;					//model type
	float   covariance_proportion; //covariance
	float	max_amplitude;			// amplitudes
	double	ang1, ang2, ang3;	//angles
	double	anis1, anis2;			//anisotropies

	double   **rotmat; // rotation matrix for this structure

	// populate structure #i

	variogram_structure()
	{
		this->model = 0;
		this->covariance_proportion = 0;
		this->max_amplitude = 0;
		this->ang1 = 0;
		this->ang2 = 0;
		this->ang3 = 0;
		this->anis1 = 0;
		this->anis2 = 0;

		this->rotmat = new double*[3];
		for (int j = 0; j < 3; ++j)
			this->rotmat[j] = new double[3];

	
	};

	variogram_structure(int model,float covariance_proportion,
		float max_amplitude,double ang1,double ang2,double ang3,
		double anis1,double anis2)
	{
		this->model=model;					
		this->covariance_proportion= covariance_proportion; 
		this->max_amplitude=max_amplitude;
		this->ang1=ang1;
		this->ang2=ang2;
		this->ang3=ang3;	
		this->anis1=anis1;
		this->anis2=anis2;			

		rotmat = new double*[3];
		for (int j = 0; j < 3; ++j) 
		{
			rotmat[j] = new double[3];
		}

		Setrot();
			
	};


	//copy
	variogram_structure(const variogram_structure& in_structure)
	{
		this->model = in_structure.model;
		this->covariance_proportion = in_structure.covariance_proportion;
		this->max_amplitude = in_structure.max_amplitude;
		this->ang1 = in_structure.ang1;
		this->ang2 = in_structure.ang2;
		this->ang3 = in_structure.ang3;
		this->anis1 = in_structure.anis1;
		this->anis2 = in_structure.anis2;
	
		this->rotmat = new double*[3];
		for (int j = 0; j < 3; ++j)
			this->rotmat[j] = new double[3];
			
		for (int j = 0; j < 3; ++j)
			for (int i = 0; i < 3;++i )
				this->rotmat[j][i] = in_structure.rotmat[j][i];
		

	};

	//assignment
	variogram_structure& operator= (const variogram_structure& in_structure)
	{
		this->model = in_structure.model;
		this->covariance_proportion = in_structure.covariance_proportion;
		this->max_amplitude = in_structure.max_amplitude;
		this->ang1 = in_structure.ang1;
		this->ang2 = in_structure.ang2;
		this->ang3 = in_structure.ang3;
		this->anis1 = in_structure.anis1;
		this->anis2 = in_structure.anis2;

		this->rotmat = new double*[3];
		for (int j = 0; j < 3; ++j)
			this->rotmat[j] = new double[3];

		for (int j = 0; j < 3; ++j)
			for (int i = 0; i < 3; ++i)
				this->rotmat[j][i] = in_structure.rotmat[j][i];

		return *this;

	};


	~variogram_structure()
	{
		
		for (int j = 0; j < 3; ++j) 
			delete[] rotmat[j];
		
		delete[] rotmat;

	};


	int get_structure(registry_type *reg, std::string group_id, int i) 
	{
		// populates the object with data
		
		float aa1, aa2; // auxiliary vars
		reg_key *k;


		std::string variogram = group_id + "S" + utilities::make_string(i);
	
		logger << "Loading " << variogram.c_str() << "\n";

		if ((k = get_key(reg, (char*)variogram.c_str(), (char*)("TYPE"))) == NULL)
			return 1;
		model = get_int(k);

		if (model < 1 || model > 5) 
		{
			logger << " Selected model is NOT allowed. Choose a different model and restart ";
			exit(-1);
			//return 1;
		}

		logger << "Variogram model: " << model << "\n";

		if ((k = get_key(reg, (char*)variogram.c_str(), (char*)("COV"))) == NULL)
			return 1;	
		covariance_proportion = get_float(k);

		logger << "Covariance proportion: " << covariance_proportion << "\n";


		if ((k = get_key(reg, (char*)variogram.c_str(), (char*)("ANG1"))) == NULL)
			return 1;
		ang1 = get_float(k);
		logger << "Angle one: " << ang1 << "\n";

		if ((k = get_key(reg, (char*)variogram.c_str(), (char*)("ANG2"))) == NULL)
			return 1;
		ang2 = get_float(k);
		logger << "Angle two: " << ang2 << "\n";

		if ((k = get_key(reg, (char*)variogram.c_str(), (char*)("ANG3"))) == NULL)
			return 1;
		ang3 = get_float(k);
		logger << "Angle three: " << ang3 << "\n";

		if ((k = get_key(reg, (char*)variogram.c_str(), (char*)("AA"))) == NULL)
			return 1;
		max_amplitude = get_float(k);
		logger << "Amplitude one: " << max_amplitude << "\n";

		if ((k = get_key(reg, (char*)variogram.c_str(), (char*)("AA1"))) == NULL)
			return 1;
		aa1 = get_float(k);
		logger << "Amplitude two: " << aa1 << "\n";


		if ((k = get_key(reg, (char*)variogram.c_str(), (char*)("AA2"))) == NULL)
			return 1;
		aa2 = get_float(k);
		logger << "Amplitude three: " << aa2 << "\n";


		anis1 = aa1 / (float)std::max((float)(max_amplitude), (float)1.0E-20);
		anis2 = aa2 / (float)std::max((float)(max_amplitude), (float)1.0E-20);

		
		//printf(" it,cc,ang[1,2,3]: %d %6.2f %6.2f %6.2f %6.2f \n", it[i], cc[i], ang1[i], ang2[i], ang3[i]);
		//printf("         a1 a2 a3: %4.0f %4.0f %4.0f \n", aa[i], aa1, aa2);


		// now prepare the rotation matrix


		rotmat = new double*[3];
		for (int j = 0; j < 3; ++j) 
		{
			rotmat[j] = new double[3];
		}

		Setrot();
		
		return 0;
	//	logger << " INITIALIZE ROTATION AND ANISOTROPY \n\n";

		//Setrot(0, sang1, sang2, sang3, sanis1, sanis2);
		//for (unsigned int is = 1; is <= nst; ++is)

	};

	
	void Setrot()
	{
		static const double pi = 3.14159265358979323846;
		static const double deg2rad = pi / 180.0;
		static const double epslon = 1.0E-20;
		
		double alpha;

		//if (ang1 >= 0 && ang1 < 270)
		//{
		//	alpha = (90.0 - ang1) * deg2rad;
		//}
		//else
		//{
		//	alpha = (450.0 - ang1) * deg2rad;
		//}

		alpha =	ang1 * deg2rad;
		//double beta =	-1*ang2 * deg2rad; 
		double beta = ang2 * deg2rad;
		double theta =	ang3 * deg2rad;

		// Get the required sines and cosines:
		double sina = (double)sin(alpha);
		double sinb = (double)sin(beta);
		double sint = (double)sin(theta);
		double cosa = (double)cos(alpha);
		double cosb = (double)cos(beta);
		double cost = (double)cos(theta);

		// Construct the rotation matrix in the required memory:
		double afac1 = 1.0 / std::max(anis1, epslon);
		double afac2 = 1.0 / std::max(anis2, epslon);

		rotmat[0][0] = (double)(cosb * cosa);
		rotmat[0][1] = (double)(cosb * sina);
		rotmat[0][2] = (double)(-sinb);
		rotmat[1][0] = (double)(afac1*(-cost*sina + sint*sinb*cosa));
		rotmat[1][1] = (double)(afac1*(cost*cosa + sint*sinb*sina));
		rotmat[1][2] = (double)(afac1*(sint * cosb));
		rotmat[2][0] = (double)(afac2*(sint*sina + cost*sinb*cosa));
		rotmat[2][1] = (double)(afac2*(-sint*cosa + cost*sinb*sina));
		rotmat[2][2] = (double)(afac2*(cost * cosb));

		return;
	};


	double	sqdist(cl_coord<double> &point1, cl_coord<double> &point2)
	{

		// Compute component distance vectors and the squared distance:
		double sqdista = 0.0;

		cl_coord<double> temp = point1 - point2;

		for (int i = 0; i < 3; ++i)
			sqdista += pow(rotmat[i][0] * (double)temp.x + rotmat[i][1] * (double)temp.y +	rotmat[i][2] 
							* (double)temp.z, 2.);

		return sqdista;
	};


	double get_covariance(cl_coord<double> &point1, cl_coord<double> &point2 )
	{

		double cova = 0.0;

		// Compute the appropriate distance:
		double h = sqrt(sqdist(point1, point2));
		double hr = h / max_amplitude;

		switch (model)
		{

		case 1:            // Spherical Variogram Model?
			if (hr < (float)1.0)
				cova = (covariance_proportion*(1. - hr*(1.5 - .5*hr*hr)));
			break;
		case 2:            // Exponential Variogram Model?
			cova = (covariance_proportion*exp(-3.0*hr));
			break;
		case 3:            // Gaussian Variogram Model ?
			cova = (covariance_proportion*exp(-3.0*pow(hr, 2)));
			break;

		case 4: // Hole effect model
			cova = covariance_proportion*pow(cos(hr*PI), 2);
			break;
		//case 5:
		//	// Dampened hole effect
		//	cova = covariance_proportion*(1 - exp(-3 * h / 320))*pow(cos(hr*PI), 2);

		}

		return cova;
	};


	//double get_covariance_and_sqdist(cl_coord<double> &point1, cl_coord<double> &point2, double &sqdistan)
	//{

	//	double cova = 0.0;

	//	// Compute the appropriate distance:
	//	sqdistan = sqdist(point1, point2);
	//	
	//	double h = sqrt(sqdistan);
	//	double hr = h / max_amplitude;

	//	switch (model)
	//	{

	//	case 1:            // Spherical Variogram Model?
	//		if (hr < (float)1.0)
	//			cova = (covariance_proportion*(1. - hr*(1.5 - .5*hr*hr)));
	//		break;
	//	case 2:            // Exponential Variogram Model?
	//		cova = (covariance_proportion*exp(-3.0*hr));
	//		break;
	//	case 3:            // Gaussian Variogram Model ?
	//		cova = (covariance_proportion*exp(-3.0*pow(hr, 2)));
	//		break;

	//	case 4: // Hole effect model
	//		cova = covariance_proportion*pow(cos(hr*PI), 2);
	//		break;
	//	case 5:
	//		// Dampened hole effect
	//		cova = covariance_proportion*(1 - exp(-3 * h / 320))*pow(cos(hr*PI), 2);

	//	}

	//	return cova;
	//};


};



//class covariance_table
//{
//	float    ***covtab;
//
//	//constrctors
//	covariance_table();
//
//	//copy
//	covariance_table(const covariance_table& in_covt) {};
//
//	//assignment
//	covariance_table& operator=(const covariance_table& in_covt) {};
//
//
//
//	void	Ctable_notext(GSLibGridPars& grid, misc& utils)
//	{
//		long ncts; // max number of points in the covariance table
//		double hsqd;
//	
//		int ic, jc, kc, *temp2;
//		float *temp1;
//		float c[1];
//		long loc;
//	
//		//logger << " BUILDING COVARIANCE TABLE \n";
//	
//	
//		// Size of the look-up table:
//		nctxyz.x = std::min((unsigned int)((MAXCT.x - 1) / 2), (grid.get_nx() - 1));
//		nctxyz.y = std::min((unsigned int)((MAXCT.y - 1) / 2), (grid.get_ny() - 1));
//		nctxyz.z = std::min((unsigned int)((MAXCT.z - 1) / 2), (grid.get_nz() - 1));
//	
//		ncts = (2 * nctxyz.x + 1)*(2 * nctxyz.y + 1)*(2 * nctxyz.z + 1);
//	
//	
//		// Allocalte size for grid variables
//		temp1 = new float[ncts];
//		temp2 = new int[ncts];
//	
//	
//		// NOTE: If dynamically allocating memory, and if there is no shortage
//		//    it would a good idea to go at least as far as the radius and
//		//    twice that far if you wanted to be sure that all covariances
//		//    in the left hand covariance matrix are within the table look-up.
//	
//		// Initialize the covariance subroutine and cbb at the same time:
//		cl_coord<double> zerocoord(0, 0, 0);
//		cbb = (float)covariance(zerocoord, zerocoord);
//	
//		// Now, set up the table and keep track of the node offsets that are within the search radius:
//		nlooku = 0;
//	
//	
//		//float xx, yy, zz;
//		cl_coord<double> xxyyzz;
//	
//		float x_siz = grid.get_xsiz();
//		float y_siz = grid.get_ysiz();
//		float z_siz = grid.get_zsiz();
//	
//		int ctab_size_x = (nctxyz.x * 2 + 1);
//		int ctab_size_xy = (nctxyz.x * 2 + 1)*(nctxyz.y * 2 + 1);
//	
//		for (int i = -nctxyz.x; i <= nctxyz.x; ++i)
//		{
//			xxyyzz.x = i * x_siz;
//			ic = nctxyz.x + i + 1;
//			for (int j = -nctxyz.y; j <= nctxyz.y; ++j) {
//				xxyyzz.y = j * y_siz;
//				jc = nctxyz.y + j + 1;
//	
//				int precalc_patial_position = ic + (jc - 1)*ctab_size_x;
//				for (int k = -nctxyz.z; k <= nctxyz.z; ++k) {
//					xxyyzz.z = k * z_siz;
//					kc = nctxyz.z + k + 1;
//	
//					//cl_coord<float> xxyyzz(xx,yy,zz);
//					//covtab[ic][jc][kc] = (float)covariance(zerocoord, xxyyzz);
//	
//	
//					covtab[ic][jc][kc] = covariance_and_sqdist(zerocoord, xxyyzz, hsqd);
//	
//	
//	
//					//hsqd = search_radius->sqdist(zerocoord, xxyyzz);
//	
//					if ((float)hsqd <= radsqd)
//					{
//						// We want to search by closest variogram distance (and use the
//						// anisotropic Euclidean distance to break ties:
//						temp1[nlooku] = (float)-(covtab[ic][jc][kc] - TINY * hsqd);
//	
//						//temp2[nlooku] = (int) (ic + (jc-1)*MAXCT.x + (kc-1)*MAXCT.x*MAXCT.y);
//						temp2[nlooku] = (int)(/*ic + (jc - 1)*ctab_size_x*/ precalc_patial_position + (kc - 1)*ctab_size_xy);
//	
//						++nlooku;
//	
//					}
//				}
//			}
//		}
//	
//		// Finished setting up the look-up table, now order the nodes such
//		// that the closest ones, according to variogram distance, are searched
//		// first. Note: the "loc" array is used because I didn't want to make
//		// special allowance for 2 byte integers in the sorting subroutine:
//		utils.sortit(0, nlooku - 1, temp1, 1, temp2, c, c, c, c, c, c);
//	
//	
//		//ijk_node = new cl_coord<int>[ncts];
//		ijk_node.reset(new cl_coord<int>[ncts]);
//	
//		for (long il = 0; il < nlooku; ++il) {
//			loc = (long)temp2[il];
//	
//			ijk_node[il].z = (int)((loc - 1) / (ctab_size_xy)) + 1;
//	
//			int	z_aux = (ijk_node[il].z - 1)*(ctab_size_xy);
//	
//			ijk_node[il].y = (int)((loc - z_aux - 1) / (nctxyz.x * 2 + 1)) + 1;
//			ijk_node[il].x = (int)loc - z_aux - (ijk_node[il].y - 1)*ctab_size_x;
//	
//	
//		}
//	
//		// Check for maximum number of nodes
//		//if (nodmax > 128 /*MAXNOD*/) {
//		//	logger << " WARNING SGSIM: Maximum close nodes = "<< std::to_strin(nodmax) << "\n ";
//		//	logger << "  Reset from your specification due to storage limitations - thats too much!. \n\n";
//		//	nodmax = 128/*MAXNOD*/;
//		//}
//	
//		// All finished:
//	
//		delete[]temp1;
//		delete[]temp2;
//	
//		return;
//	};
//	
//
//
//
//	//destructor
//	~covariance_table()
//	{
//		try
//		{
//			int ctab_x_size = (2 * nctxyz.x + 1);
//			int ctab_y_size = (2 * nctxyz.y + 1);
//			
//			for (int i = 0; i <= ctab_x_size; ++i) {
//				for (int j = 0; j <= ctab_y_size; ++j) {
//					delete[] covtab[i][j];
//				}
//				delete[] covtab[i];
//			}
//			delete[] covtab;
//			//throw 20;
//		}
//		catch (int e)
//		{
//			std::cout << "An exception occurred while trying to deallocate covtab. Ex nº " << e << '\n';
//		}
//	};
//

//	void	Ctable(GSLibGridPars& grid, misc& utils)
//	{
//		long ncts; // max number of points in the covariance table
//		double hsqd;
//
//		int ic, jc, kc;
//		//float *temp1;
//		float c[1];
//		long loc;
//
//		logger << " BUILDING COVARIANCE TABLE \n";
//
//
//		// Size of the look-up table:
//		nctxyz.x = std::min((unsigned int)((MAXCT.x - 1) / 2), (grid.get_nx() - 1));
//		nctxyz.y = std::min((unsigned int)((MAXCT.y - 1) / 2), (grid.get_ny() - 1));
//		nctxyz.z = std::min((unsigned int)((MAXCT.z - 1) / 2), (grid.get_nz() - 1));
//
//		ncts = (2 * nctxyz.x + 1)*(2 * nctxyz.y + 1)*(2 * nctxyz.z + 1);
//		//ncts = (2 * (nctxyz.x + 2))*(2 * (nctxyz.y + 2))*(2 * (nctxyz.z + 2));
//
//
//
//		//covtab = new float**[MAXCT.x +1];
//		//for (unsigned int i = 0; i <= (MAXCT.x); ++i) {
//		//	covtab[i]=  new float*[MAXCT.y +1];
//		//	for (unsigned int j = 0; j <= (MAXCT.y); ++j) {
//		//		covtab[i][j]=new float[MAXCT.z +1];
//		//		
//		//		for (unsigned int k = 0; k <= (MAXCT.z); ++k)
//		//			covtab[i][j][k] = 0;
//		//	}
//		//}
//
//
//		//covtab.reset(new boost::shared_array<boost::shared_array<float>>[(nctxyz.x * 2 + 1) + 1]);
//
//		covtab = new float**[(nctxyz.x * 2 + 1) + 1]; // +1 cause fortran
//		for (unsigned int i = 0; i <= nctxyz.x * 2 + 1; ++i)
//		{ // +1 cause fortran
//			//(covtab[i]).reset(new boost::shared_array<float>[nctxyz.y * 2 + 1 + 1]);
//			covtab[i] = new float*[nctxyz.y * 2 + 1 + 1];
//			for (unsigned int j = 0; j <= nctxyz.y * 2 + 1; ++j)
//			{ // +1 cause fortran
//				//(covtab[i])[j].reset(new float[nctxyz.z * 2 + 1 + 1]);
//				covtab[i][j] = new float[nctxyz.z * 2 + 1 + 1];
//				for (unsigned int k = 0; k <= nctxyz.z * 2 + 1; ++k)
//					covtab[i][j][k] = 0;
//			}
//		}
//
//
//		// Allocalte size for grid variables
//		//temp1= new float[ncts];
//		//temp2= new int[ncts];
//
//		boost::scoped_array<float> temp1(new float[ncts]);
//		boost::scoped_array<int> temp2(new int[ncts]);
//
//
//
//		// Debugging output:
//		logger << " Covariance Look up table and search for previously\n";
//		logger << " simulated grid nodes.  The maximum range in each\n";
//		logger << " coordinate direction for covariance look up is:\n";
//		logger << "     X direction: " << utilities::make_string(nctxyz.x*grid.get_xsiz()) << "\n";
//		logger << "     Y direction: " << utilities::make_string(nctxyz.y*grid.get_ysiz()) << "\n";
//		logger << "     Z direction: " << utilities::make_string(nctxyz.z*grid.get_zsiz()) << "\n";
//		logger << " Node values are not searched beyond this distance! \n\n";
//
//
//
//		//perhaps to save on memory, check for positions where covariances are allways zero
//		// and rescale the covtable adequatly. TO DO
//
//
//		// NOTE: If dynamically allocating memory, and if there is no shortage
//		//    it would a good idea to go at least as far as the radius and
//		//    twice that far if you wanted to be sure that all covariances
//		//    in the left hand covariance matrix are within the table look-up.
//
//		// Initialize the covariance subroutine and cbb at the same time:
//		cl_coord<double> zerocoord(0, 0, 0);
//		cbb = (float)covariance(zerocoord, zerocoord);
//
//		// Now, set up the table and keep track of the node offsets that are within the search radius:
//		nlooku = 0;
//
//
//		//float xx, yy, zz;
//		cl_coord<double> xxyyzz;
//
//		float x_siz = grid.get_xsiz();
//		float y_siz = grid.get_ysiz();
//		float z_siz = grid.get_zsiz();
//
//		for (int i = -nctxyz.x; i <= nctxyz.x; ++i)
//		{
//			xxyyzz.x = i * x_siz;
//			ic = nctxyz.x + i + 1;
//			for (int j = -nctxyz.y; j <= nctxyz.y; ++j)
//			{
//				xxyyzz.y = j * y_siz;
//				jc = nctxyz.y + j + 1;
//				for (int k = -nctxyz.z; k <= nctxyz.z; ++k)
//				{
//					xxyyzz.z = k * z_siz;
//					kc = nctxyz.z + k + 1;
//
//					//cl_coord<float> xxyyzz(xx,yy,zz);
//					covtab[ic][jc][kc] = (float)covariance(zerocoord, xxyyzz);
//
//					hsqd = search_radius->sqdist(zerocoord, xxyyzz);
//
//					if ((float)hsqd <= radsqd)
//					{
//						// We want to search by closest variogram distance (and use the
//						// anisotropic Euclidean distance to break ties:
//						temp1[nlooku] = (float)-(covtab[ic][jc][kc] - TINY * hsqd);
//
//						//temp2[nlooku] = (int) (ic + (jc-1)*MAXCT.x + (kc-1)*MAXCT.x*MAXCT.y);
//						temp2[nlooku] = (int)(ic + (jc - 1)*(nctxyz.x * 2 + 1) + (kc - 1)*(nctxyz.x * 2 + 1)*(nctxyz.y * 2 + 1)); //MAXCT.x, etc index in covtable
//
//						++nlooku;
//
//					}
//				}
//			}
//		}
//
//		// Finished setting up the look-up table, now order the nodes such
//		// that the closest ones, according to variogram distance, are searched
//		// first. Note: the "loc" array is used because I didn't want to make
//		// special allowance for 2 byte integers in the sorting subroutine:
//		utils.sortit(0, nlooku - 1, &temp1[0], 1, &temp2[0], c, c, c, c, c, c);
//
//
//		//ijk_node = new cl_coord<int>[ncts]; // ijk is the offset to current node?
//		ijk_node.reset(new cl_coord<int>[ncts]);
//
//		for (long il = 0; il < nlooku; ++il) {
//			loc = (long)temp2[il]; //index of the offset in covtable
//
//
//			//if (temp1[il] >= 0)  // perhaps if here after a certain point il will catch only zeros, reduce nlooku to that.
//			//	// but temp1 doent have the covariance directly
//			//{
//			//	nlooku = il;
//			//	break;
//			//}
//
//
//			//ijk_node[il].z = (int)((loc-1)/(MAXCT.x*MAXCT.y)) + 1;
//			//ijk_node[il].y = (int)((loc-(ijk_node[il].z -1)*(MAXCT.x*MAXCT.y)-1)/ MAXCT.x) + 1;
//			//ijk_node[il].x = (int)loc - (ijk_node[il].z -1)*MAXCT.x*MAXCT.y - (ijk_node[il].y -1)*MAXCT.x;
//
//			ijk_node[il].z = (int)((loc - 1) / ((nctxyz.x * 2 + 1)*(nctxyz.y * 2 + 1))) + 1;
//			ijk_node[il].y = (int)((loc - (ijk_node[il].z - 1)*((nctxyz.x * 2 + 1)*(nctxyz.y * 2 + 1)) - 1) / (nctxyz.x * 2 + 1)) + 1;
//			ijk_node[il].x = (int)loc - (ijk_node[il].z - 1)*(nctxyz.x * 2 + 1)*(nctxyz.y * 2 + 1) - (ijk_node[il].y - 1)*(nctxyz.x * 2 + 1);
//
//
//		}
//
//
//
//		// Check for maximum number of nodes
//		//if (nodmax > 128 /*MAXNOD*/) {
//		//	logger << " WARNING SGSIM: Maximum close nodes = "<< std::to_strin(nodmax) << "\n ";
//		//	logger << "  Reset from your specification due to storage limitations - thats too much!. \n\n";
//		//	nodmax = 128/*MAXNOD*/;
//		//}
//
//		// All finished:
//
//		//delete[]temp1;
//		//delete[]temp2;
//
//		return;
//	};

//	void	Ctable_notext(GSLibGridPars& grid, misc& utils)
//	{
//		long ncts; // max number of points in the covariance table
//		double hsqd;
//
//		int ic, jc, kc, *temp2;
//		float *temp1;
//		float c[1];
//		long loc;
//
//		//logger << " BUILDING COVARIANCE TABLE \n";
//
//
//		// Size of the look-up table:
//		nctxyz.x = std::min((unsigned int)((MAXCT.x - 1) / 2), (grid.get_nx() - 1));
//		nctxyz.y = std::min((unsigned int)((MAXCT.y - 1) / 2), (grid.get_ny() - 1));
//		nctxyz.z = std::min((unsigned int)((MAXCT.z - 1) / 2), (grid.get_nz() - 1));
//
//		ncts = (2 * nctxyz.x + 1)*(2 * nctxyz.y + 1)*(2 * nctxyz.z + 1);
//
//
//		// Allocalte size for grid variables
//		temp1 = new float[ncts];
//		temp2 = new int[ncts];
//
//
//		// NOTE: If dynamically allocating memory, and if there is no shortage
//		//    it would a good idea to go at least as far as the radius and
//		//    twice that far if you wanted to be sure that all covariances
//		//    in the left hand covariance matrix are within the table look-up.
//
//		// Initialize the covariance subroutine and cbb at the same time:
//		cl_coord<double> zerocoord(0, 0, 0);
//		cbb = (float)covariance(zerocoord, zerocoord);
//
//		// Now, set up the table and keep track of the node offsets that are within the search radius:
//		nlooku = 0;
//
//
//		//float xx, yy, zz;
//		cl_coord<double> xxyyzz;
//
//		float x_siz = grid.get_xsiz();
//		float y_siz = grid.get_ysiz();
//		float z_siz = grid.get_zsiz();
//
//		int ctab_size_x = (nctxyz.x * 2 + 1);
//		int ctab_size_xy = (nctxyz.x * 2 + 1)*(nctxyz.y * 2 + 1);
//
//		for (int i = -nctxyz.x; i <= nctxyz.x; ++i)
//		{
//			xxyyzz.x = i * x_siz;
//			ic = nctxyz.x + i + 1;
//			for (int j = -nctxyz.y; j <= nctxyz.y; ++j) {
//				xxyyzz.y = j * y_siz;
//				jc = nctxyz.y + j + 1;
//
//				int precalc_patial_position = ic + (jc - 1)*ctab_size_x;
//				for (int k = -nctxyz.z; k <= nctxyz.z; ++k) {
//					xxyyzz.z = k * z_siz;
//					kc = nctxyz.z + k + 1;
//
//					//cl_coord<float> xxyyzz(xx,yy,zz);
//					//covtab[ic][jc][kc] = (float)covariance(zerocoord, xxyyzz);
//
//
//					covtab[ic][jc][kc] = covariance_and_sqdist(zerocoord, xxyyzz, hsqd);
//
//
//
//					//hsqd = search_radius->sqdist(zerocoord, xxyyzz);
//
//					if ((float)hsqd <= radsqd)
//					{
//						// We want to search by closest variogram distance (and use the
//						// anisotropic Euclidean distance to break ties:
//						temp1[nlooku] = (float)-(covtab[ic][jc][kc] - TINY * hsqd);
//
//						//temp2[nlooku] = (int) (ic + (jc-1)*MAXCT.x + (kc-1)*MAXCT.x*MAXCT.y);
//						temp2[nlooku] = (int)(/*ic + (jc - 1)*ctab_size_x*/ precalc_patial_position + (kc - 1)*ctab_size_xy);
//
//						++nlooku;
//
//					}
//				}
//			}
//		}
//
//		// Finished setting up the look-up table, now order the nodes such
//		// that the closest ones, according to variogram distance, are searched
//		// first. Note: the "loc" array is used because I didn't want to make
//		// special allowance for 2 byte integers in the sorting subroutine:
//		utils.sortit(0, nlooku - 1, temp1, 1, temp2, c, c, c, c, c, c);
//
//
//		//ijk_node = new cl_coord<int>[ncts];
//		ijk_node.reset(new cl_coord<int>[ncts]);
//
//		for (long il = 0; il < nlooku; ++il) {
//			loc = (long)temp2[il];
//
//			ijk_node[il].z = (int)((loc - 1) / (ctab_size_xy)) + 1;
//
//			int	z_aux = (ijk_node[il].z - 1)*(ctab_size_xy);
//
//			ijk_node[il].y = (int)((loc - z_aux - 1) / (nctxyz.x * 2 + 1)) + 1;
//			ijk_node[il].x = (int)loc - z_aux - (ijk_node[il].y - 1)*ctab_size_x;
//
//
//		}
//
//		// Check for maximum number of nodes
//		//if (nodmax > 128 /*MAXNOD*/) {
//		//	logger << " WARNING SGSIM: Maximum close nodes = "<< std::to_strin(nodmax) << "\n ";
//		//	logger << "  Reset from your specification due to storage limitations - thats too much!. \n\n";
//		//	nodmax = 128/*MAXNOD*/;
//		//}
//
//		// All finished:
//
//		delete[]temp1;
//		delete[]temp2;
//
//		return;
//	};
//
//
//
//};
//


class VariogramPars 
{

public:

	unsigned int  nst; //number of structures
	float	 c0;  // nugget effect
	variogram_structure *structure;
	float	 sill; //variogram sill 

	variogram_structure *search_radius;
	float	radsqd;

	float    ***covtab;

	//boost::shared_array<boost::shared_array<boost::shared_array<float>>> covtab;
	//boost::shared_array<boost::shared_array<boost::shared_array<float>>> covtab;

	cl_coord<int> nctxyz;//lookup table size

	float	 cbb;
	//cl_coord<int> *ijk_node;
	//boost::shared_array<cl_coord<int>> ijk_node;
	std::vector<cl_coord<int>> ijk_node;


	long	 nlooku;
	//int		 nodmax;


	cl_coord<unsigned int> MAXCT;//max lookup table size

	//default constructor
	VariogramPars() {}

	//constructor from registry

	VariogramPars(registry_type *reg, std::string keyword, int sufix)
	{

		reg_key *k;
		std::string  var;
		float radius, radius1, radius2;
		sill = 0.0; // initialize
		std::string group_id;


		group_id = keyword + utilities::make_string(sufix);
		// parse variogram section

		k = get_key(reg, (char*)group_id.c_str(), (char*)("NSTRUCT"));
		if (k)
			nst = get_int(k);
		else 
			exit(-99);

		if (nst < 1)
		{
			logger << "NST must be at least 1!";
			exit(0);
		}


		k = get_key(reg, (char*)group_id.c_str(), (char*)("NUGGET"));
		if (k)
			c0 = get_float(k);
		else
			exit(-99);
		if (c0 < 0)
		{
			logger << "Nugget is smaller than 0!";
			exit(0);
		}

		logger << " Nst, C0: " << utilities::make_string(nst) << " " << utilities::make_string(c0) << "\n";



		structure = new variogram_structure[nst];

		for (unsigned int i = 0; i < nst; ++i)
		{
			structure[i].get_structure(reg, group_id, i + 1);
		}


		// end the variogram section

		// search

		// = new variogram_structure;


		k = get_key(reg, (char*)("SEARCH"), (char*)("RADIUS1"));
		if (k)
			radius = get_float(k);
		else 
			exit(-99);

		k = get_key(reg, (char*)("SEARCH"), (char*)("RADIUS2"));
		if (k)
			radius1 = get_float(k);
		else 
			exit(-99);

		k = get_key(reg, (char*)("SEARCH"), (char*)("RADIUS3"));
		if (k)
			radius2 = get_float(k);
		else 
			exit(-99);
		printf(" Search radius: %6.2f %6.2f %6.2f \n", radius, radius1, radius2);

		radsqd = radius * radius;

		search_radius = new variogram_structure();


		search_radius->max_amplitude = radius;
		search_radius->anis1 = radius1 / radius;
		search_radius->anis2 = radius2 / radius;

		k = get_key(reg, (char*)("SEARCH"), (char*)("SANG1"));
		if (k)
			search_radius->ang1 = get_float(k);
		else 
			exit(-99);

		k = get_key(reg, (char*)("SEARCH"), (char*)("SANG2"));
		if (k)
			search_radius->ang2 = get_float(k);
		else 
			exit(-99);

		k = get_key(reg, (char*)("SEARCH"), (char*)("SANG3"));
		if (k)
			search_radius->ang3 = get_float(k);
		else 
			exit(-99);
		printf(" Search anisotropy angles: %6.2f %6.2f %6.2f \n", search_radius->ang1, search_radius->ang2, search_radius->ang3);
		//

		// aand do rotation
		search_radius->rotmat = new double*[3];
		for (int j = 0; j < 3; ++j)
		{
			search_radius->rotmat[j] = new double[3];
		}

		search_radius->Setrot();


		// get the covtab parameters
		k = get_key(reg, (char*)("COVTAB"), (char*)("MAXCTX"));
		if (k)
			MAXCT.x = get_int(k);
		else 
			exit(-99);
		k = get_key(reg, (char*)("COVTAB"), (char*)("MAXCTY"));
		if (k)
			MAXCT.y = get_int(k);
		else 
			exit(-99);
		k = get_key(reg, (char*)("COVTAB"), (char*)("MAXCTZ"));
		if (k)
			MAXCT.z = get_int(k);
		else 
			exit(-99);

	};

	//copy constructor
	VariogramPars(const VariogramPars& in_vario)
	{
		this->nst = in_vario.nst; //numberof structures
		this->c0 = in_vario.c0;  // nugget effect

		this->sill = in_vario.sill; //variogram sill 
		this->radsqd = in_vario.radsqd;


		this->nctxyz = in_vario.nctxyz;//lookup table size

		this->cbb = in_vario.cbb;
		this->ijk_node = in_vario.ijk_node;

		this->nlooku = in_vario.nlooku;
		

		this->MAXCT = in_vario.MAXCT;//max lookup table size


		//float    *** covtab;
		
		//covtab.reset(new boost::shared_array<boost::shared_array<float>>[(nctxyz.x * 2 + 1) + 1]);
		allocate_ctable();
		//this->covtab = new float**[(nctxyz.x * 2 + 1) + 1]; // +1 cause fortran
		for (unsigned int i = 0; i <= nctxyz.x * 2 + 1; ++i)
		{ // +1 cause fortran
			//(covtab[i]).reset(new boost::shared_array<float>[nctxyz.y * 2 + 1 + 1]);

			//this->covtab[i] = new float*[nctxyz.y * 2 + 1 + 1];
			for (unsigned int j = 0; j <= nctxyz.y * 2 + 1; ++j)
			{ // +1 cause fortran
				//(covtab[i])[j].reset(new float[nctxyz.z * 2 + 1 + 1]);
				//this->covtab[i][j] = new float[nctxyz.z * 2 + 1 + 1];
				for (unsigned int k = 0; k <= nctxyz.z * 2 + 1; ++k)
					this->covtab[i][j][k] = in_vario.covtab[i][j][k];
			}
		}



		this->search_radius = new variogram_structure(*in_vario.search_radius);
		
		this->structure = new variogram_structure[nst];
		for (int i = 0; i < nst;++i)
			this->structure[i] = in_vario.structure[i];

	
	}


	//assignment 
	VariogramPars& operator= (const VariogramPars& in_vario)
	{
		this->nst = in_vario.nst; //numberof structures
		this->c0 = in_vario.c0;  // nugget effect

		this->sill = in_vario.sill; //variogram sill 
		this->radsqd = in_vario.radsqd;


		this->nctxyz = in_vario.nctxyz;//lookup table size

		this->cbb = in_vario.cbb;
		this->ijk_node = in_vario.ijk_node;

		this->nlooku = in_vario.nlooku;


		this->MAXCT = in_vario.MAXCT;//max lookup table size


		//covtab.reset(new boost::shared_array<boost::shared_array<float>>[(nctxyz.x * 2 + 1) + 1]);

									 //float    *** covtab;

		allocate_ctable();
		//this->covtab = new float**[(nctxyz.x * 2 + 1) + 1]; // +1 cause fortran
		for (unsigned int i = 0; i <= nctxyz.x * 2 + 1; ++i)
		{ // +1 cause fortran
			//(covtab[i]).reset(new boost::shared_array<float>[nctxyz.y * 2 + 1 + 1]);
			//this->covtab[i] = new float*[nctxyz.y * 2 + 1 + 1];
			for (unsigned int j = 0; j <= nctxyz.y * 2 + 1; ++j)
			{ // +1 cause fortran
				//(covtab[i])[j].reset(new float[nctxyz.z * 2 + 1 + 1]);
				//this->covtab[i][j] = new float[nctxyz.z * 2 + 1 + 1];
				for (unsigned int k = 0; k <= nctxyz.z * 2 + 1; ++k)
					this->covtab[i][j][k] = in_vario.covtab[i][j][k];
			}
		}



		this->search_radius = new variogram_structure(*in_vario.search_radius);

		this->structure = new variogram_structure[nst];
		for (int i = 0; i < nst; ++i)
			this->structure[i] = in_vario.structure[i];


	}



	/*
	VariogramPars(int nst, float c0, float sill, float radsqd, cl_coord<int> nctxyz, float cbb, cl_coord<int> *ijk_node, long nlooku, cl_coord<unsigned int> MAXCT)
	{
		this->nst = nst; //numberof structures
		this->c0 = c0;  // nugget effect

		this->sill = sill; //variogram sill 
		this->radsqd = radsqd;


		this->nctxyz = nctxyz;//lookup table size

		this->cbb = cbb;
		this->ijk_node = ijk_node;

		this->nlooku = nlooku;


		this->MAXCT = MAXCT;//max lookup table size


		//float    *** covtab;
		this->covtab = new float**[(nctxyz.x * 2 + 1) + 1]; // +1 cause fortran
		for (unsigned int i = 0; i <= nctxyz.x * 2 + 1; ++i)
		{ // +1 cause fortran
			this->covtab[i] = new float*[nctxyz.y * 2 + 1 + 1];
			for (unsigned int j = 0; j <= nctxyz.y * 2 + 1; ++j)
			{ // +1 cause fortran
				this->covtab[i][j] = new float[nctxyz.z * 2 + 1 + 1];
				for (unsigned int k = 0; k <= nctxyz.z * 2 + 1; ++k)
					this->covtab[i][j][k] = in_vario.covtab[i][j][k];
			}
		}



		this->search_radius = new variogram_structure(*in_vario.search_radius); //variogram_structure(model, covariance_proportion, max_amplitude, ang1, ang2, ang3, anis1, anis2)


		this->structure = new variogram_structure[nst];
		
		for (int i = 0; i < nst; ++i)
			this->structure[i] = variogram_structure(model, covariance_proportion, max_amplitude, ang1, ang2, ang3, anis1, anis2);
					
	}
	*/

	//destructor
	~VariogramPars()
	{
		delete [] structure;			
		delete search_radius;
		
		
		//delete[] ijk_node;
		try
		{
			for (int i = 0; i <= (2 * nctxyz.x + 1); ++i) {
				for (int j = 0; j <= (2 * nctxyz.y + 1); ++j) {
					delete[] covtab[i][j];
				}
				delete[] covtab[i];
			}
			delete[] covtab;
			//throw 20;
		}
		catch (int e)
		{
			std::cout << "An exception occurred while trying to deallocate covtab. Ex nº " << e << '\n';
		}


		//delete ijk_node;

	};


	
	// functions

	int		get_variogram_pars(registry_type *reg, std::string keyword, int sufix) { // sufix tells wich are we loading

		reg_key *k;

		std::string  var;

		float radius, radius1, radius2;

		sill = 0.0; // initialize

		std::string group_id;

		group_id = keyword + utilities::make_string(sufix);
		// parse variogram section

		k = get_key(reg, (char*)group_id.c_str(), (char*)("NSTRUCT"));
		if (k)
			nst = get_int(k);
		else return -1;

		if (nst < 1)
		{
			logger << "NST must be at least 1!";
			exit(0);
		}


		k = get_key(reg, (char*)group_id.c_str(), (char*)("NUGGET"));
		if (k)
			c0 = get_float(k);
		else 
			return -1;
		if (c0 < 0)
		{
			logger << "Nugget is smaller than 0!";
			exit(0);
		}

		logger << " Nst, C0: " << utilities::make_string(nst) << " " << utilities::make_string(c0) << "\n";



		structure = new variogram_structure[nst];

		for (unsigned int i = 0; i < nst; ++i)
		{
			structure[i].get_structure(reg, group_id, i+1);
		}


		// end the variogram section

		// search
		
		// = new variogram_structure;
		
		
		k = get_key(reg, (char*)("SEARCH"), (char*)("RADIUS1"));
		if (k)
			radius = get_float(k);
		else return -1;

		k = get_key(reg, (char*)("SEARCH"), (char*)("RADIUS2"));
		if (k)
			radius1 = get_float(k);
		else return -1;

		k = get_key(reg, (char*)("SEARCH"), (char*)("RADIUS3"));
		if (k)
			radius2 = get_float(k);
		else return -1;
		printf (" Search radius: %6.2f %6.2f %6.2f \n",radius,radius1,radius2);

		radsqd = radius * radius;

		search_radius = new variogram_structure();


		search_radius->max_amplitude = radius;
		search_radius->anis1 = radius1 / radius;
		search_radius->anis2 = radius2 / radius;

		k = get_key(reg, (char*)("SEARCH"), (char*)("SANG1"));
		if (k)
			search_radius->ang1 = get_float(k);
		else return -1;

		k = get_key(reg, (char*)("SEARCH"), (char*)("SANG2"));
		if (k)
			search_radius->ang2 = get_float(k);
		else return -1;

		k = get_key(reg, (char*)("SEARCH"), (char*)("SANG3"));
		if (k)
			search_radius->ang3 = get_float(k);
		else return -1;
		printf (" Search anisotropy angles: %6.2f %6.2f %6.2f \n", search_radius->ang1, search_radius->ang2, search_radius->ang3);
		//

		// aand do rotation
		search_radius->rotmat = new double*[3];
		for (int j = 0; j < 3; ++j)
		{
			search_radius->rotmat[j] = new double[3];
		}

		search_radius->Setrot();


		// get the covtab parameters
		k = get_key(reg, (char*)("COVTAB"), (char*)("MAXCTX"));
		if (k)
			MAXCT.x = get_int(k);
		else return -1;
		k = get_key(reg, (char*)("COVTAB"), (char*)("MAXCTY"));
		if (k)
			MAXCT.y = get_int(k);
		else return -1;
		k = get_key(reg, (char*)("COVTAB"), (char*)("MAXCTZ"));
		if (k)
			MAXCT.z = get_int(k);
		else return -1;


		// end covtab

		return 0;

	};

	double	covariance (cl_coord<double>& point1, cl_coord<double>& point2)
	{
	
		// Check for "zero" distance, return with cmax if so:
		float h = (float) search_radius->sqdist(point1, point2);
		
		if (h < EPSLON)
		{
			// Calculate the maximum covariance value 
			//(used for zero distances and for power model covariance):
			float cmax = c0;
			
			for (unsigned int is=0; is<nst; ++is)
			{
				cmax += structure[is].covariance_proportion;
			}
			
			return cmax;
		}
		// Loop over all the structures:
		double cova = 0.0;
		for (unsigned int is=0; is<nst; ++is)
		{
			cova += structure[is].get_covariance(point1, point2);
		}

		return cova;
	}

	double	covariance_and_sqdist(cl_coord<double>& point1, cl_coord<double>& point2, double& h )
	{

		// Check for "zero" distance, return with cmax if so:
		/*float */h = /*(float)*/search_radius->sqdist(point1, point2);

		if (h < EPSLON)
		{
			// Calculate the maximum covariance value 
			//(used for zero distances and for power model covariance):
			float cmax = c0;

			for (unsigned int is = 0; is<nst; ++is)
			{
				cmax += structure[is].covariance_proportion;
			}

			return cmax;
		}
		// Loop over all the structures:
		double cova = 0.0;
		for (unsigned int is = 0; is<nst; ++is)
		{
			cova += structure[is].get_covariance(point1, point2);
		}

		return cova;
	}


	double	covariance_no_search_range (cl_coord<double> point1, cl_coord<double> point2)
	{
	
		// Loop over all the structures:
		double cova = 0.0;
		for (unsigned int is=0; is<nst; ++is)
		{
			cova += structure[is].get_covariance(point1, point2);
		}

		return cova;
	}

	void allocate_ctable()
	{
		this->covtab = new float**[(nctxyz.x * 2 + 1) + 1]; // +1 cause fortran
		for (unsigned int i = 0; i <= nctxyz.x * 2 + 1; ++i)
		{ // +1 cause fortran
			//(covtab[i]).reset(new boost::shared_array<float>[nctxyz.y * 2 + 1 + 1]);
			this->covtab[i] = new float*[nctxyz.y * 2 + 1 + 1];
			for (unsigned int j = 0; j <= nctxyz.y * 2 + 1; ++j)
			{ // +1 cause fortran
				//(covtab[i])[j].reset(new float[nctxyz.z * 2 + 1 + 1]);
				this->covtab[i][j] = new float[nctxyz.z * 2 + 1 + 1];
				for (unsigned int k = 0; k <= nctxyz.z * 2 + 1; ++k)
					this->covtab[i][j][k] = 0;
			}
		}
	};

	void	Ctable (GSLibGridPars& grid, misc& utils)
	{
		long ncts; // max number of points in the covariance table
		double hsqd;
		
		int ic,jc,kc;
		//float *temp1;
		float c[1];
		long loc;

		logger << " BUILDING COVARIANCE TABLE \n";


		// Size of the look-up table:
		nctxyz.x = std::min((unsigned int)((MAXCT.x-1)/2),(grid.get_nx()-1));
		nctxyz.y = std::min((unsigned int)((MAXCT.y-1)/2),(grid.get_ny()-1));
		nctxyz.z = std::min((unsigned int)((MAXCT.z-1)/2),(grid.get_nz()-1));

		ncts=(2*nctxyz.x+1)*(2*nctxyz.y+1)*(2*nctxyz.z+1);
		//ncts = (2 * (nctxyz.x + 2))*(2 * (nctxyz.y + 2))*(2 * (nctxyz.z + 2));



		//covtab = new float**[MAXCT.x +1];
		//for (unsigned int i = 0; i <= (MAXCT.x); ++i) {
		//	covtab[i]=  new float*[MAXCT.y +1];
		//	for (unsigned int j = 0; j <= (MAXCT.y); ++j) {
		//		covtab[i][j]=new float[MAXCT.z +1];
		//		
		//		for (unsigned int k = 0; k <= (MAXCT.z); ++k)
		//			covtab[i][j][k] = 0;
		//	}
		//}


		//covtab.reset(new boost::shared_array<boost::shared_array<float>>[(nctxyz.x * 2 + 1) + 1]);
		allocate_ctable();
		//this->covtab = new float**[(nctxyz.x*2 + 1)+1]; // +1 cause fortran
		//for (unsigned int i = 0; i <= nctxyz.x * 2+1; ++i) 
		//{ // +1 cause fortran
		//	//(covtab[i]).reset(new boost::shared_array<float>[nctxyz.y * 2 + 1 + 1]);
		//	this->covtab[i] = new float*[nctxyz.y*2 + 1+1];
		//	for (unsigned int j = 0; j <= nctxyz.y * 2+1; ++j) 
		//	{ // +1 cause fortran
		//		//(covtab[i])[j].reset(new float[nctxyz.z * 2 + 1 + 1]);
		//		this->covtab[i][j] = new float[nctxyz.z*2 + 1+1];
		//		for (unsigned int k = 0; k <= nctxyz.z * 2+1; ++k)
		//			this->covtab[i][j][k] = 0;
		//	}
		//}


		// Allocalte size for grid variables
		//temp1= new float[ncts];
		//temp2= new int[ncts];

		boost::scoped_array<float> temp1(new float[ncts]);
		boost::scoped_array<int> temp2(new int[ncts]);



		// Debugging output:
		logger  << " Covariance Look up table and search for previously\n";
		logger  << " simulated grid nodes.  The maximum range in each\n";
		logger  << " coordinate direction for covariance look up is:\n";
		logger  << "     X direction: " << utilities::make_string(nctxyz.x*grid.get_xsiz()) << "\n";
		logger  << "     Y direction: " << utilities::make_string(nctxyz.y*grid.get_ysiz()) << "\n";
		logger  << "     Z direction: " << utilities::make_string(nctxyz.z*grid.get_zsiz()) << "\n";
		logger  << " Node values are not searched beyond this distance! \n\n";
		
		

		//perhaps to save on memory, check for positions where covariances are allways zero
		// and rescale the covtable adequatly. TO DO


		// NOTE: If dynamically allocating memory, and if there is no shortage
		//    it would a good idea to go at least as far as the radius and
		//    twice that far if you wanted to be sure that all covariances
		//    in the left hand covariance matrix are within the table look-up.

		// Initialize the covariance subroutine and cbb at the same time:
		cl_coord<double> zerocoord(0,0,0);
		cbb=(float)covariance(zerocoord,zerocoord);

		// Now, set up the table and keep track of the node offsets that are within the search radius:
		nlooku = 0;


		//float xx, yy, zz;
		cl_coord<double> xxyyzz;

		float x_siz = grid.get_xsiz();
		float y_siz = grid.get_ysiz();
		float z_siz = grid.get_zsiz();

		for ( int i= -nctxyz.x ; i <= nctxyz.x; ++i) 
		{
			xxyyzz.x = i *x_siz;
			ic = nctxyz.x  + i + 1;
			for (int j=-nctxyz.y; j<=nctxyz.y; ++j) 
			{
				xxyyzz.y = j * y_siz;
				jc = nctxyz.y + j  + 1;
				for (int k=-nctxyz.z; k<=nctxyz.z; ++k) 
				{
					xxyyzz.z = k * z_siz;
					kc = nctxyz.z + k+1;

					//cl_coord<float> xxyyzz(xx,yy,zz);
					covtab[ic][jc][kc]=(float)covariance(zerocoord,xxyyzz);
					
					hsqd = search_radius->sqdist(zerocoord,xxyyzz);
					
					if ((float)hsqd <= radsqd)
					{
						// We want to search by closest variogram distance (and use the
						// anisotropic Euclidean distance to break ties:
						temp1[nlooku] = (float) - (covtab[ic][jc][kc]-TINY*hsqd);
						
						//temp2[nlooku] = (int) (ic + (jc-1)*MAXCT.x + (kc-1)*MAXCT.x*MAXCT.y);
						temp2[nlooku] = (int)(ic + (jc - 1)*(nctxyz.x *2+1) + (kc - 1)*(nctxyz.x * 2 + 1)*(nctxyz.y * 2 + 1)); //MAXCT.x, etc index in covtable

						++nlooku ;

					}
				}
			}
		}

		// Finished setting up the look-up table, now order the nodes such
		// that the closest ones, according to variogram distance, are searched
		// first. Note: the "loc" array is used because I didn't want to make
		// special allowance for 2 byte integers in the sorting subroutine:
		utils.sortit(0,nlooku-1,&temp1[0],1,&temp2[0],c,c,c,c,c,c);


		//ijk_node = new cl_coord<int>[ncts]; // ijk is the offset to current node?
		//ijk_node.reset(new cl_coord<int>[ncts]);

		for (long il=0; il<nlooku; ++il) {
			loc = (long)temp2[il]; //index of the offset in covtable

			cl_coord<int> temp_coord;
			//if (temp1[il] >= 0)  // perhaps if here after a certain point il will catch only zeros, reduce nlooku to that.
			//	// but temp1 doent have the covariance directly
			//{
			//	nlooku = il;
			//	break;
			//}


			//ijk_node[il].z = (int)((loc-1)/(MAXCT.x*MAXCT.y)) + 1;
			//ijk_node[il].y = (int)((loc-(ijk_node[il].z -1)*(MAXCT.x*MAXCT.y)-1)/ MAXCT.x) + 1;
			//ijk_node[il].x = (int)loc - (ijk_node[il].z -1)*MAXCT.x*MAXCT.y - (ijk_node[il].y -1)*MAXCT.x;
			temp_coord.z = (int)((loc - 1) / ((nctxyz.x * 2 + 1)*(nctxyz.y * 2 + 1))) + 1;
			temp_coord.y = (int)((loc - (temp_coord.z - 1)*((nctxyz.x * 2 + 1)*(nctxyz.y * 2 + 1)) - 1) / (nctxyz.x * 2 + 1)) + 1;
			temp_coord.x = (int)loc - (temp_coord.z - 1)*(nctxyz.x * 2 + 1)*(nctxyz.y * 2 + 1) - (temp_coord.y - 1)*(nctxyz.x * 2 + 1);

			ijk_node.push_back(temp_coord);

		}



		// Check for maximum number of nodes
		//if (nodmax > 128 /*MAXNOD*/) {
		//	logger << " WARNING SGSIM: Maximum close nodes = "<< std::to_strin(nodmax) << "\n ";
		//	logger << "  Reset from your specification due to storage limitations - thats too much!. \n\n";
		//	nodmax = 128/*MAXNOD*/;
		//}

		// All finished:
		
		//delete[]temp1;
		//delete[]temp2;

		return;
	};

	void	Ctable_notext(GSLibGridPars& grid, misc& utils)
	{
		long ncts; // max number of points in the covariance table
		double hsqd;

		int ic, jc, kc, *temp2;
		float *temp1;
		float c[1];
		long loc;

		//logger << " BUILDING COVARIANCE TABLE \n";


		// Size of the look-up table:
		nctxyz.x = std::min((unsigned int)((MAXCT.x - 1) / 2), (grid.get_nx() - 1));
		nctxyz.y = std::min((unsigned int)((MAXCT.y - 1) / 2), (grid.get_ny() - 1));
		nctxyz.z = std::min((unsigned int)((MAXCT.z - 1) / 2), (grid.get_nz() - 1));

		ncts = (2 * nctxyz.x + 1)*(2 * nctxyz.y + 1)*(2 * nctxyz.z + 1);


		// Allocalte size for grid variables
		temp1 = new float[ncts];
		temp2 = new int[ncts];


		// NOTE: If dynamically allocating memory, and if there is no shortage
		//    it would a good idea to go at least as far as the radius and
		//    twice that far if you wanted to be sure that all covariances
		//    in the left hand covariance matrix are within the table look-up.

		// Initialize the covariance subroutine and cbb at the same time:
		cl_coord<double> zerocoord(0, 0, 0);
		cbb = (float)covariance(zerocoord, zerocoord);

		// Now, set up the table and keep track of the node offsets that are within the search radius:
		nlooku = 0;


		//float xx, yy, zz;
		cl_coord<double> xxyyzz;

		float x_siz = grid.get_xsiz();
		float y_siz = grid.get_ysiz();
		float z_siz = grid.get_zsiz();

		int ctab_size_x = (nctxyz.x * 2 + 1);
		int ctab_size_xy = (nctxyz.x * 2 + 1)*(nctxyz.y * 2 + 1);

		for (int i = -nctxyz.x; i <= nctxyz.x; ++i)
		{
			xxyyzz.x = i *x_siz;
			ic = nctxyz.x + i + 1;
			for (int j = -nctxyz.y; j <= nctxyz.y; ++j) {
				xxyyzz.y = j * y_siz;
				jc = nctxyz.y + j + 1;

				int precalc_patial_position = ic + (jc - 1)*ctab_size_x;
				for (int k = -nctxyz.z; k <= nctxyz.z; ++k) {
					xxyyzz.z = k * z_siz;
					kc = nctxyz.z + k + 1;

					//cl_coord<float> xxyyzz(xx,yy,zz);
					//covtab[ic][jc][kc] = (float)covariance(zerocoord, xxyyzz);


					covtab[ic][jc][kc] = covariance_and_sqdist(zerocoord, xxyyzz, hsqd);
					


					//hsqd = search_radius->sqdist(zerocoord, xxyyzz);

					if ((float)hsqd <= radsqd)
					{
						// We want to search by closest variogram distance (and use the
						// anisotropic Euclidean distance to break ties:
						temp1[nlooku] = (float)-(covtab[ic][jc][kc] - TINY*hsqd);

						//temp2[nlooku] = (int) (ic + (jc-1)*MAXCT.x + (kc-1)*MAXCT.x*MAXCT.y);
						temp2[nlooku] = (int)(/*ic + (jc - 1)*ctab_size_x*/ precalc_patial_position + (kc - 1)*ctab_size_xy);

						++nlooku;

					}
				}
			}
		}

		// Finished setting up the look-up table, now order the nodes such
		// that the closest ones, according to variogram distance, are searched
		// first. Note: the "loc" array is used because I didn't want to make
		// special allowance for 2 byte integers in the sorting subroutine:
		utils.sortit(0, nlooku - 1, temp1, 1, temp2, c, c, c, c, c, c);


		//ijk_node = new cl_coord<int>[ncts];
		//ijk_node.reset(new cl_coord<int>[ncts]);

		for (long il = 0; il<nlooku; ++il) {
			loc = (long)temp2[il];

			ijk_node[il].z = (int)((loc - 1) / (ctab_size_xy)) + 1; 

			int	z_aux = (ijk_node[il].z - 1)*(ctab_size_xy);

			ijk_node[il].y = (int)((loc - z_aux - 1) / (nctxyz.x * 2 + 1)) + 1;
			ijk_node[il].x = (int)loc - z_aux - (ijk_node[il].y - 1)*ctab_size_x;


		}

		// Check for maximum number of nodes
		//if (nodmax > 128 /*MAXNOD*/) {
		//	logger << " WARNING SGSIM: Maximum close nodes = "<< std::to_strin(nodmax) << "\n ";
		//	logger << "  Reset from your specification due to storage limitations - thats too much!. \n\n";
		//	nodmax = 128/*MAXNOD*/;
		//}

		// All finished:

		delete[]temp1;
		delete[]temp2;

		return;
	};



};


#endif