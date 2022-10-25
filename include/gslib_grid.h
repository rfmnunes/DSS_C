#ifndef GSLIB_GRID_INCLUDED
#define GSLIB_GRID_INCLUDED

#include <iostream>
#include <string>
#include <fstream>

#include "registry.h"
#include "utils.h"
#include "log.h"
#include "coordinates.h"


inline void endian_swap(unsigned int& x)
{
    x = (x>>24) | 
        ((x<<8) & 0x00FF0000) |
        ((x>>8) & 0x0000FF00) |
        (x<<24);
}

inline float float_endian_swap( const float inFloat )
{
   float retVal;
   char *floatToConvert = ( char* ) & inFloat;
   char *returnFloat = ( char* ) & retVal;

   // swap the bytes into a temporary buffer
   returnFloat[0] = floatToConvert[3];
   returnFloat[1] = floatToConvert[2];
   returnFloat[2] = floatToConvert[1];
   returnFloat[3] = floatToConvert[0];

   return retVal;
}

//
//class GSLibGridPars_generic
//{
//	unsigned int number_dimensions;
//
//	CoordinatesClass<int>* n_blocks;
//	CoordinatesClass<float>* origin_blocks;
//	CoordinatesClass<float>* size_blocks;
//	CoordinatesClass<int>* cumulative_blocks;
//
// 
//public:
//
//	GSLibGridPars_generic(registry_type *reg)
//	{//get grid parameters from registry.
//	//	reg_key *k;
//	//	/* get grid parameters */
//	//	// number
//	//	//set_nx(read_int_from_reg(reg, "GRID", "NX"));
//	//	k = get_key(reg, "GRID", "NX");
//	//	if (k)
//	//		if (set_nx(get_int(k)))
//	//		{
//	//			logger << "error loading grid parameters";
//	//			//g_log.log_string("error loading grid parameters");
//	//			return;
//	//		}
//	//	k = get_key(reg, "GRID", "NY");
//	//	if (k)
//	//		if (set_ny(get_int(k)))
//	//		{
//	//			logger << "error loading grid parameters";
//	//			return;
//	//		}
//	//	k = get_key(reg, "GRID", "NZ");
//	//	if (k)
//	//		if (set_nz(get_int(k)))
//	//		{
//	//			logger << "error loading grid parameters";
//	//			return;
//	//		}
//	//	// done n
//	//	// do orig
//	//	k = get_key(reg, "GRID", "ORIGX");
//	//	if (k)
//	//		if (set_ox(get_float(k)))
//	//		{
//	//			logger << "error loading grid parameters";
//	//			return;
//	//		}
//	//	k = get_key(reg, "GRID", "ORIGY");
//	//	if (k)
//	//		if (set_oy(get_float(k)))
//	//		{
//	//			logger << "error loading grid parameters";
//	//			return;
//	//		}
//	//	k = get_key(reg, ("GRID"), ("ORIGZ"));
//	//	if (k)
//	//		if (set_oz(get_float(k)))
//	//		{
//	//			logger << "error loading grid parameters";
//	//			return;
//	//		}
//	//	// end orig
//	//	// get size
//	//	k = get_key(reg, ("GRID"), ("SIZEX"));
//	//	if (k)
//	//		if (set_xsiz(get_float(k)))
//	//		{
//	//			logger << "error loading grid parameters";
//	//			return;
//	//		}
//	//	k = get_key(reg, ("GRID"), ("SIZEY"));
//	//	if (k)
//	//		if (set_ysiz(get_float(k)))
//	//		{
//	//			logger << "error loading grid parameters";
//	//			return;
//	//		}
//	//	k = get_key(reg, ("GRID"), ("SIZEZ"));
//	//	if (k)
//	//		if (set_zsiz(get_float(k)))
//	//		{
//	//			logger << "error loading grid parameters";
//	//			return;
//	//		}
//	//	logger << " XX: nodes, min., spacing: " << utilities::make_string(nx) << " " << utilities::make_string(ox) << " " << utilities::make_string(xsiz) << "\n";
//	//	logger << " YY: nodes, min., spacing: " << utilities::make_string(ny) << " " << utilities::make_string(oy) << " " << utilities::make_string(ysiz) << "\n";
//	//	logger << " ZZ: nodes, min., spacing: " << utilities::make_string(nz) << " " << utilities::make_string(oz) << " " << utilities::make_string(zsiz) << "\n";
//	//	//
//	//	set_nxy((long)get_nx()*(long)get_ny());
//	//	set_nxyz(get_nxy()*(long)get_nz());
//	//	// done size
//	//	logger << "Grid Parameters loaded...\n";
//	//	// end grid
//	};
//
//	//copy constructor
//	GSLibGridPars_generic(GSLibGridPars_generic& grid_pars)
//	{
//		number_dimensions = grid_pars.number_dimensions;
//
//		n_blocks = grid_pars.n_blocks;
//		origin_blocks = grid_pars.origin_blocks;
//		size_blocks = grid_pars.size_blocks;
//		cumulative_blocks = grid_pars.cumulative_blocks;
//	};
//
//	GSLibGridPars_generic(CoordinatesClass<int>& in_n_blocks, CoordinatesClass<float>& in_origin_blocks, CoordinatesClass<float>& in_size_blocks, CoordinatesClass<int>& in_cumulative_blocks)
//	{
//		number_dimensions = (unsigned int)in_n_blocks.coord.size();
//
//		n_blocks = &in_n_blocks;
//		origin_blocks = &in_origin_blocks;
//		size_blocks = &in_size_blocks;
//		cumulative_blocks = &in_cumulative_blocks;
//	};
//
//	//generic cooords
//
//	unsigned long get_index_from_ijk_c(CoordinatesClass<int> &in_struct) {
//
//		long index = 0;
//
//		//do stuff
//		for (unsigned int i = 1; i < in_struct.coord.size()/*n_dimensions*/; ++i)
//		{
//			index += (in_struct.coord[i] - 1) * cumulative_blocks->coord[i - 1];
//		}
//		//do the first coordinate
//		index += in_struct.coord[0]-1;
//
//		return index;
//	}; // get array index from the block coordinates 
//
//	unsigned long get_index_from_xyz_c(const CoordinatesClass<float> &in_struct)
//	{
//		return get_index_from_ijk_c(get_ijk_form_xyz_c(in_struct));
//	}
//
//	CoordinatesClass<int> get_ijk_form_xyz_c(const CoordinatesClass<float> &in_coord)
//	{
//		// computing indexes
//
//		CoordinatesClass<float> temp_coord = ((CoordinatesClass<float>)in_coord - *origin_blocks);
//		temp_coord = (temp_coord / *size_blocks);
//
//		CoordinatesClass<int> out_coord = (temp_coord);
//		for (unsigned int i = 0; i < in_coord.coord.size()/*n_dimensions*/; ++i)
//		{
//			out_coord.coord[i] = out_coord.coord[i] + 1.5;
//
//			if (out_coord.coord[i] < 1)
//			{
//				out_coord.coord[i] = 1;
//			}
//			else if (out_coord.coord[i] > n_blocks->coord[i])
//			{
//				out_coord.coord[i] = n_blocks->coord[i];
//			}
//		}
//
//		return out_coord;
//
//	}
//
//	CoordinatesClass<int> get_ijk_from_index_c(const unsigned long &index)
//	{
//		CoordinatesClass<int> out_struct(number_dimensions);
//
//		long rest = index;
//
//		for (unsigned int i = (unsigned int)out_struct.coord.size()/*n_dimensions*/ - 1; i > 0; --i)
//		{	
//			out_struct.coord[i] = (rest / cumulative_blocks->coord[i-1])+1;
//			rest = rest % cumulative_blocks->coord[i-1];
//
//		}
//		out_struct.coord[0] = rest + 1;
//
//		return out_struct;
//
//	};
//
//	CoordinatesClass<float> get_xyz_from_ijk_c(const CoordinatesClass<int> &in_struct)
//	{
//
//		CoordinatesClass<int> str_ones((unsigned int)in_struct.coord.size()/*n_dimensions*/);
//		for (unsigned int i = 0; i < (unsigned int)in_struct.coord.size()/*n_dimensions*/; ++i)
//		{
//			str_ones.coord[i] = 1;
//		}
//
//		CoordinatesClass<float> out_struct((unsigned int)in_struct.coord.size()/*n_dimensions*/);
//
//
//		out_struct = (*origin_blocks) + ((CoordinatesClass<int>)in_struct - str_ones)* (*size_blocks);
//
//		return out_struct;
//
//	};
//
//	CoordinatesClass<float> get_xyz_from_index_c(const unsigned long &indx)
//	{
//		return get_xyz_from_ijk_c(get_ijk_from_index_c(indx));
//	};
//
//	// fortran type
//
//	unsigned long get_index_from_ijk_f(const CoordinatesClass<int> &in_struct) {
//
//		long index = 0;
//
//		//do stuff
//		for (unsigned int i = 1; i < in_struct.coord.size(); ++i)
//		{
//			index += (in_struct.coord[i] - 1) * cumulative_blocks->coord[i - 1];
//		}
//		//do the first coordinate
//		index += in_struct.coord[0] - 1;
//
//		return index;
//	}; // get array index from the block coordinates 
//
//	unsigned long get_index_from_xyz_f(const CoordinatesClass<float> &in_struct)
//	{
//		return get_index_from_ijk_f(get_ijk_form_xyz_c(in_struct));
//	}
//
//	CoordinatesClass<int> get_ijk_form_xyz_f(const CoordinatesClass<float> &in_coord)
//	{
//		// computing indexes
//
//		CoordinatesClass<float> temp_coord = ((CoordinatesClass<float>)in_coord - *origin_blocks);
//		temp_coord = (temp_coord / *size_blocks);
//
//		CoordinatesClass<int> out_coord = (temp_coord);
//		for (unsigned int i = 0; i < in_coord.coord.size()/*n_dimensions*/; ++i)
//		{
//			out_coord.coord[i] = out_coord.coord[i] + 1.5;
//
//			if (out_coord.coord[i] < 1)
//			{
//				out_coord.coord[i] = 1;
//			}
//			else if (out_coord.coord[i] > n_blocks->coord[i])
//			{
//				out_coord.coord[i] = n_blocks->coord[i];
//			}
//		}
//
//		return out_coord;
//
//	}
//
//	CoordinatesClass<int> get_ijk_from_index_f(const unsigned long &index)
//	{
//		CoordinatesClass<int> out_struct(number_dimensions);
//
//		long rest = index-1;// frotran index is +1?????? 
//
//		for (unsigned int i = (unsigned int)out_struct.coord.size() - 1; i > 0; --i)
//		{
//			out_struct.coord[i] = (rest / cumulative_blocks->coord[i - 1]) + 1;
//			rest = rest % cumulative_blocks->coord[i - 1];
//
//		}
//		out_struct.coord[0] = rest + 1;
//
//		return out_struct;
//
//	};
//
//	CoordinatesClass<float> get_xyz_from_ijk_f(const CoordinatesClass<int> &in_struct)
//	{
//
//		CoordinatesClass<int> str_ones((unsigned int)in_struct.coord.size());
//		for (unsigned int i = 0; i < in_struct.coord.size(); ++i)
//		{
//			str_ones.coord[i] = 1;
//		}
//
//		CoordinatesClass<float> out_struct((unsigned int)in_struct.coord.size());
//
//
//		out_struct = (*origin_blocks) + ((CoordinatesClass<int>)in_struct - str_ones)* (*size_blocks);
//
//		return out_struct;
//
//	};
//
//	CoordinatesClass<float> get_xyz_from_index_f(const unsigned long &indx)
//	{
//		return get_xyz_from_ijk_f(get_ijk_from_index_f(indx));
//	};
//
//	//some utils
//	//bool is_ijk_inside(cl_coord<int> ijk)
//	//{
//	//	if (ijk.x < 1 || ijk.y < 1 || ijk.z < 1 ||
//	//		ijk.x > get_nx() || ijk.y > get_ny() ||
//	//		ijk.z > get_nz())
//	//	{
//	//		return false;
//	//	}
//	//	else
//	//	{
//	//		return true;
//	//	}
//	//};
//
//};
//

class GSLibGridPars {

	int		nx, ny, nz;
	float	ox, oy, oz,
			xsiz, ysiz, zsiz;

	long long		nxy, nxyz;


public:

	// constructor
	GSLibGridPars() { }

	GSLibGridPars(registry_type *reg)
	{//get grid parameters from registry.

		reg_key *k;

		/* get grid parameters */

		// number

		set_nx(read_int_from_reg(reg, (char*)"GRID", (char*)"NX"));

		k = get_key(reg, (char*)"GRID", (char*)"NX");
		if (k)
			if (set_nx(get_int(k)))
			{
				logger << "error loading grid parameters";
				//g_log.log_string("error loading grid parameters");
				return;
			}


		k = get_key(reg, (char*)"GRID", (char*)"NY");
		if (k)
			if (set_ny(get_int(k)))
			{
				logger << "error loading grid parameters";
				return;
			}

		k = get_key(reg, (char*)"GRID", (char*)"NZ");
		if (k)
			if (set_nz(get_int(k)))
			{
				logger << "error loading grid parameters";
				return;
			}

		// done n


		// do orig
		k = get_key(reg, (char*)"GRID", (char*)"ORIGX");
		if (k)
			if (set_ox(get_float(k)))
			{
				logger << "error loading grid parameters";
				return;
			}

		k = get_key(reg, (char*)"GRID", (char*)"ORIGY");
		if (k)
			if (set_oy(get_float(k)))
			{
				logger << "error loading grid parameters";
				return;
			}

		k = get_key(reg, (char*)("GRID"), (char*)("ORIGZ"));
		if (k)
			if (set_oz(get_float(k)))
			{
				logger << "error loading grid parameters";
				return;
			}

		// end orig

		// get size
		k = get_key(reg, (char*)("GRID"), (char*)("SIZEX"));
		if (k)
			if (set_xsiz(get_float(k)))
			{
				logger << "error loading grid parameters";
				return;
			}

		k = get_key(reg, (char*)("GRID"), (char*)("SIZEY"));
		if (k)
			if (set_ysiz(get_float(k)))
			{
				logger << "error loading grid parameters";
				return;
			}

		k = get_key(reg, (char*)("GRID"), (char*)("SIZEZ"));
		if (k)
			if (set_zsiz(get_float(k)))
			{
				logger << "error loading grid parameters";
				return;
			}


		logger << " XX: nodes, min., spacing: " << utilities::make_string(nx) << " " << utilities::make_string(ox) << " " << utilities::make_string(xsiz) << "\n";
		logger << " YY: nodes, min., spacing: " << utilities::make_string(ny) << " " << utilities::make_string(oy) << " " << utilities::make_string(ysiz) << "\n";
		logger << " ZZ: nodes, min., spacing: " << utilities::make_string(nz) << " " << utilities::make_string(oz) << " " << utilities::make_string(zsiz) << "\n";


		//
		set_nxy((long)get_nx()*(long)get_ny());
		set_nxyz(get_nxy()*(long)get_nz());

		// done size

		logger << "Grid Parameters loaded...\n";
		// end grid

	};

	GSLibGridPars(const GSLibGridPars& grid_pars)
	{
		nx = grid_pars.nx; ny = grid_pars.ny; nz = grid_pars.nz;//number
		ox = grid_pars.ox; oy = grid_pars.oy; oz = grid_pars.oz;//origin
		xsiz = grid_pars.xsiz; ysiz = grid_pars.ysiz; zsiz = grid_pars.zsiz; //size

		nxy = grid_pars.nxy; nxyz = grid_pars.nxyz; //


	};

	GSLibGridPars(int nx, int ny, int nz, float ox, float oy, float oz, float xsiz, float ysiz, float zsiz)
	{
		this->nx = nx; this->ny = ny; this->nz = nz;//number
		this->ox = ox; this->oy = oy; this->oz = oz;//origin
		this->xsiz = xsiz; this->ysiz = ysiz; this->zsiz = zsiz; //size

		this->nxy = nx*ny; this->nxyz = nx*ny*nz; //


	};

	long long get_index_from_ijk_f(const cl_coord<int>& in_struct) {
		return (in_struct.x + (in_struct.y - 1)*nx + (in_struct.z - 1)*nxy);
	}; // get array index from the block coordinates 

	long long get_index_from_xyz_f(const cl_coord<double>& in_struct) {

		cl_coord<int> temp_ijk = get_ijk_form_xyz_f(in_struct);
		long long index = get_index_from_ijk_f(temp_ijk);
		//(temp_ijk.x + (temp_ijk.y-1)*nx + (temp_ijk.z-1)*nxy);

		return index;
	};

	cl_coord<int> get_ijk_form_xyz_f(const cl_coord<double>& in_struct) //gets index of coord in ijk
	{
		cl_coord<int> out_struct;

		// Compute the index of xcord
		out_struct.x = (long)((in_struct.x - ox) / xsiz + 1.5);

		// Check to see if in or out
		if (out_struct.x < 1)
			out_struct.x = 1;
		else if (out_struct.x > nx)
			out_struct.x = nx;

		// Compute the index of ycord
		out_struct.y = (long)((in_struct.y - oy) / ysiz + 1.5);

		// Check to see if in or out
		if (out_struct.y < 1)
			out_struct.y = 1;
		else if (out_struct.y > ny)
			out_struct.y = ny;

		// Compute the index of zcord
		out_struct.z = (long)((in_struct.z - oz) / zsiz + 1.5);

		// Check to see if in or out
		if (out_struct.z < 1)
			out_struct.z = 1;
		else if (out_struct.z > nz)
			out_struct.z = nz;
		return out_struct;
	}

	cl_coord<int> get_ijk_from_index_f(long long index) {

		cl_coord<int> out_struct;

		out_struct.z = ((index - 1) / nxy) + 1;
		out_struct.y = ((index - (out_struct.z - 1)*nxy - 1) / nx) + 1;
		out_struct.x = index - (out_struct.z - 1)*nxy - (out_struct.y - 1)*nx;

		return out_struct;
	};

	cl_coord<double> get_xyz_from_ijk_f(const cl_coord<int>& in_struct) {

		cl_coord<double> out_struct;

		out_struct.x = ox + (in_struct.x - 1)*xsiz;
		out_struct.y = oy + (in_struct.y - 1)*ysiz;
		out_struct.z = oz + (in_struct.z - 1)*zsiz;

		return out_struct;
	};

	cl_coord<double> get_xyz_from_index_f(long long indx) {

		cl_coord<int> ijk = get_ijk_from_index_f(indx);
		cl_coord<double> out_struct = get_xyz_from_ijk_f(ijk);

		return out_struct;

	}; // get the block coordinates from the array index 

	//some utils
	bool is_ijk_inside(cl_coord<int>& ijk)
	{
		if (ijk.x < 1 || ijk.y < 1 || ijk.z < 1 ||
			ijk.x > get_nx() || ijk.y > get_ny() ||
			ijk.z > get_nz())
		{
			return false;
		}
		else
		{
			return true;
		}
	};

	bool is_xyz_inside(cl_coord<double> xyz)
	{
		if (xyz.x < this->get_ox() - this->get_xsiz() / 2 ||
			xyz.y < this->get_oy() - this->get_ysiz() / 2 ||
			xyz.z < this->get_oz() - this->get_zsiz() / 2 ||

			xyz.x > this->get_ox() + this->get_nx()*this->get_xsiz() ||
			xyz.y > this->get_oy() + this->get_ny()*this->get_ysiz() ||
			xyz.z > this->get_oz() + this->get_nz()*this->get_zsiz())
		{
			return false;
		}
		else
		{
			return true;
		}
	};

	//setters and getters

	int set_nx(int value)
	{
		if (value >= 1)
		{
			nx = value;
			return 0;
		}
		else
			logger << "Set number of blocks gave an error! Check the input parameters!";
		return 1;
	};
	int set_ny(int value)
	{
		if (value >= 1)
		{
			ny = value;
			return 0;
		}
		else
			logger << "Set number of blocks gave an error! Check the input parameters!";
		return 1;
	};
	int set_nz(int value)
	{
		if (value >= 1)
		{
			nz = value;
			return 0;
		}
		else
			logger << "Set number of blocks gave an error! Check the input parameters!";
		return 1;
	};

	inline unsigned int get_nx() { return nx; };
	inline unsigned int get_ny() { return ny; };
	inline unsigned int get_nz() { return nz; };

	int set_ox(float value)
	{
		ox = value;
		return 0;
	};
	int set_oy(float value)
	{
		oy = value;
		return 0;
	};
	int set_oz(float value)
	{
		oz = value;
		return 0;
	};

	inline float get_ox() { return ox; };
	inline float get_oy() { return oy; };
	inline float get_oz() { return oz; };


	int set_xsiz(float value)
	{
		if (value >= 0)
		{
			xsiz = value;
			return 0;
		}
		else
			logger << "Set size of blocks gave an error! Check the input parameters!";
		return 1;
	};
	int set_ysiz(float value)
	{
		if (value >= 0)
		{
			ysiz = value;
			return 0;
		}
		else
			logger << "Set size of blocks gave an error! Check the input parameters!";
		return 1;
	};
	int set_zsiz(float value)
	{
		if (value >= 0)
		{
			zsiz = value;
			return 0;
		}
		else
			logger << "Set size of blocks gave an error! Check the input parameters!";
		return 1;
	};

	int set_nxy(long long value)
	{

		if (value >= 0)
		{
			nxy = value;
			return 0;
		}
		else
			logger << "Something went wrong setting nxy!";
		return 1;
	}
	int set_nxyz(long long value)
	{

		if (value >= 0)
		{
			nxyz = value;
			return 0;
		}
		else
			logger << "Something went wrong setting nxyz!";
		return 1;
	}


	inline float get_xsiz() { return xsiz; };
	inline float get_ysiz() { return ysiz; };
	inline float get_zsiz() { return zsiz; };

	inline long long get_nxy() { return nxy; };
	inline long long get_nxyz() { return nxyz; };

	
};


template <class templt> class GSLibGrid : public GSLibGridPars{

public:

	//templt	*grid;

	// lets try to make a smart pointer here:
 
	boost::shared_array<templt> grid;

	// file
	float	null_data;

	long long		valid_points; // points with valid data - nxyz if no nulls are used  set with read grid, or copy from somewhere

	std::string *filename; // associated file name (if we have one)

	//header stuff
	std::string	*title; // header title
	std::string	*var_name; //variable we are using name
	
	unsigned int		n_vars; // number of vars on file

	bool	headerflag;

	double	average, variance;


	GSLibGrid(int nx, int ny, int nz, float ox, float oy, float oz, float xsiz, float ysiz, float zsiz);

	GSLibGrid(GSLibGridPars&);		// construct from parameters
	GSLibGrid(GSLibGrid<templt>&);  // construct form a grid
	
	GSLibGrid(char* &); /*from sgems*/

	~GSLibGrid();					//destructor frees memory for the grid
	

	GSLibGrid<templt>& operator= (const GSLibGrid<templt>& a);

	//void copy_grid(GSLibGrid * a); //also keep the grid

	void read_grid();

	void read_grid_binary_sgems();

	void get_average();
	
	void get_variance();



	void write_grid_ascii(); //puts *grid on file
	void write_grid_binary();

	void write_grid_binary_sgems(); 

	void zero_the_grid();

	int get_valid_points();
protected:
	
	void read_header(std::ifstream& f_in); // reads gslib header
		
	void write_header(std::ofstream& f_out);// writes gslib header

	//bool has_nulls(); // returns true if we have nulls, false otherwise

//private:	
//	
//	void put_grid_pars(GSLibGridPars& a);

	// void allocate_grid();

	//	mmap_grid(); // mmaps a file...   TO DO LATER IF USEFUL
};

//grid pars constructor

template <class templt>	 GSLibGrid<templt>::GSLibGrid(int nx, int ny, int nz, float ox, float oy, float oz, float xsiz, float ysiz, float zsiz)//constructor allocates memory for the grid (+1 for fortran) 
{

	this->set_nx(nx); this->set_ny(ny); this->set_nz(nz);//number
	this->set_ox( ox); this->set_oy (oy); this->set_oz(oz);//origin
	this->set_xsiz( xsiz); this->set_ysiz( ysiz); this->set_zsiz (zsiz); //size

	this->set_nxy(nx*ny); this->set_nxyz (nx*ny*nz); //


	
	//this->grid_pars = a;
	null_data = 0;
	n_vars = 0;
	headerflag = 0;
	valid_points = 0;
	average = 0;
	variance = 0;
	//grid = NULL;
	filename = NULL;
	title = NULL;
	var_name = NULL;

	//grid = new templt[(get_nxyz()) + 1];
	grid.reset(new templt[(get_nxyz()) + 1]);

	//grid[0] = -999.25; // just to mark the fortran type grids
	//for (unsigned int i= 0; i <= get_nxyz(); i++) grid[i] = NULL;//0; // and initialize to 0
	//allocate_grid();

};

template <class templt>	 GSLibGrid<templt>::GSLibGrid(GSLibGridPars& a)//constructor allocates memory for the grid (+1 for fortran) 
	{

	set_nx(a.get_nx()),set_ny(a.get_ny()),set_nz(a.get_nz()),
	set_ox(a.get_ox()),set_oy(a.get_oy()),set_oz(a.get_oz()),
	set_xsiz(a.get_xsiz()), set_ysiz(a.get_ysiz()),set_zsiz(a.get_zsiz()),
	set_nxy(a.get_nxy()),set_nxyz(a.get_nxyz());

	//this->grid_pars = a;
	null_data = 0;
	n_vars = 0;
	headerflag = 0;
	valid_points = 0;
	average = 0;
	variance = 0;
	//grid = NULL; 
	filename = NULL;
	title = NULL;
	var_name = NULL;

	//grid = new templt [(get_nxyz())+1]; 
	grid.reset(new templt[(get_nxyz()) + 1]);

	//grid[0] = -999.25; // just to mark the fortran type grids
	//for (unsigned int i= 0; i <= get_nxyz(); i++) grid[i] = NULL;//0; // and initialize to 0
	//allocate_grid();

	};

template <class templt>	 GSLibGrid<templt>::GSLibGrid(char* & file)//constructor from a sgems file 
	{
	filename = NULL;
	title = NULL;
	var_name = NULL;

	filename = new std::string;
	//filename(file);
	filename->append(file);
	
	std::ifstream fdata(file,std::ios::binary);

	if (!fdata.is_open()){
		logger << file << " failed to open!";
		return ;
	}

	// Write a header with a "magic number" and the grid type
	//stream << (Q_UINT32)0xB211175D;
			
	// first, load the magic number
		
	char magic[4];
	fdata.read (magic, 4);


	//then object_type

	unsigned int n_chars;
	fdata.read ((char*)&n_chars,4);
	endian_swap(n_chars);


	char *buffer = new char[n_chars];
	fdata.read (buffer,n_chars);
	std::string classname(buffer);
	delete [] buffer;

	fdata.read ((char*)&n_chars,4);
	endian_swap(n_chars);

	char *buffer2 = new char[n_chars];
	fdata.read (buffer2,n_chars);
	title = new std::string(buffer2);
	//title(buffer2);
	delete [] buffer2;

		

	//now read some stuff


	unsigned int version; 
	unsigned int nx, ny, nz;

	fdata.read((char*)&version, sizeof(unsigned int));
	//endian_swap(version);
		
	fdata.read((char*)&nx, sizeof(unsigned int));
	endian_swap(nx);
	fdata.read((char*)&ny, sizeof(unsigned int));
	endian_swap(ny);
	fdata.read((char*)&nz, sizeof(unsigned int));
	endian_swap(nz);

	set_nx(nx);
	set_ny(ny);
	set_nz(nz);

	logger << "nx: " << nx << "\n";
	logger << "ny: " << ny<< "\n";
	logger << "nz: " << nz<< "\n";

	set_nxy(nx*ny);
	set_nxyz(nx*ny*nz);

	float sizex, sizey, sizez;

	fdata.read((char*)&sizex, sizeof(float));
	sizex=float_endian_swap(sizex);
	fdata.read((char*)&sizey, sizeof(float));
	sizey=float_endian_swap(sizey);
	fdata.read((char*)&sizez, sizeof(float));
	sizez=float_endian_swap(sizez);

	set_xsiz(sizex);
	set_ysiz(sizey);
	set_zsiz(sizez);

	logger << "x_size: " << sizex<< "\n";
	logger << "y_size: " << sizey<< "\n";
	logger << "z_size: " << sizez<< "\n";


	float ox, oy, oz;

	fdata.read((char*)&ox, sizeof(float));
	ox=float_endian_swap(ox);
	fdata.read((char*)&oy, sizeof(float));
	oy=float_endian_swap(oy);
	fdata.read((char*)&oz, sizeof(float));
	oz=float_endian_swap(oz);

	set_ox(ox);
	set_oy(oy);
	set_oz(oz);

	logger << "x_origin: " <<ox<< "\n";
	logger << "y_origin: " << oy<< "\n";
	logger << "z_origin: " << oz<< "\n";

	//unsigned int n_vars;
	fdata.read((char*)&n_vars, sizeof(unsigned int));
	endian_swap(n_vars);


	var_name = new std::string[n_vars];
	for( unsigned int i = 0; i < n_vars; i++ ) {


		fdata.read ((char*)&n_chars,4);
		endian_swap(n_chars);

		char *buffer3 = new char[n_chars];
		fdata.read (buffer3,n_chars);
		var_name[i]=buffer3;
		delete [] buffer3;

	}

	//grid = NULL; 
	//grid = new templt [(get_nxyz())+1]; 
	grid.reset(new templt[(get_nxyz()) + 1]);

	float value;

	for( unsigned int i = 0; i < n_vars; ++i ) {
		for (long long j= 1;j<=get_nxyz();++j){
		
			fdata.read((char*)&value, sizeof(float));
			value=float_endian_swap(value);
		
			grid[j]=value;
		}
	}

	fdata.close();

	null_data = 0;
	//n_vars = 0;
	//headerflag = 0;
	valid_points = 0;
	average = 0;
	variance = 0;
	//grid = NULL; 
	//
	//grid = new templt [(get_nxyz())+1]; 
	
	//grid[0] = -999.25; // just to mark the fortran type grids
	//for (unsigned int i= 0; i <= get_nxyz(); i++) grid[i] = NULL;//0; // and initialize to 0
	//allocate_grid();

	};

//other grid constructor (copy) 
template <class templt>	 GSLibGrid<templt>::GSLibGrid(GSLibGrid<templt>& a)//constructor allocates memory for the grid (+1 for fortran) 
{
	filename = NULL;
	title = NULL;
	var_name = NULL;

	set_nx(a.get_nx()), set_ny(a.get_ny()), set_nz(a.get_nz()),
	set_ox(a.get_ox()),set_oy(a.get_oy()),set_oz(a.get_oz()),
	set_xsiz(a.get_xsiz()), set_ysiz(a.get_ysiz()),set_zsiz(a.get_zsiz()),
	set_nxy(a.get_nxy()),set_nxyz(a.get_nxyz());
	//grid_pars = a.grid_pars;


	null_data = a.null_data;
	n_vars = a.n_vars;
	headerflag = a.headerflag;
	valid_points = a.valid_points;
	average = a.average;
	variance = a.variance;
	grid = a.grid; 


	if (a.filename!=NULL)
	filename= new std::string(*a.filename); // associated file name (if we have one)

	if (a.title != NULL)
	title= new std::string(*a.title); // header title
	
	var_name = new std::string[n_vars];
	
	for (int i = 0; i < n_vars; ++i)
	{
		var_name[i] =a.var_name[i]; //variable we are using name
	}

	//allocate_grid();
	//grid = new templt [(get_nxyz())+1];
	grid.reset(new templt[(get_nxyz()) + 1]);
	
	//grid[0] = -999.25; // just to mark the fortran type grids
	for (long long i= 1; i <= get_nxyz(); i++)
	{
		grid[i] = a.grid[i]; 
	}
};

// also do assignment operator 
template <class templt> GSLibGrid<templt>&  GSLibGrid<templt>::operator= (const GSLibGrid<templt>& a)//constructor allocates memory for the grid (+1 for fortran) 
{
	filename = NULL;
	title = NULL;
	var_name = NULL;

	set_nx(a.get_nx()); set_ny(a.get_ny()); set_nz(a.get_nz());
	set_ox(a.get_ox()); set_oy(a.get_oy()); set_oz(a.get_oz());
	set_xsiz(a.get_xsiz()); set_ysiz(a.get_ysiz()); set_zsiz(a.get_zsiz());
	set_nxy(a.get_nxy()); set_nxyz(a.get_nxyz());
	//grid_pars = a.grid_pars;


	null_data = a.null_data;
	n_vars = a.n_vars;
	headerflag = a.headerflag;
	valid_points = a.valid_points;
	average = a.average;
	variance = a.variance;
	grid = a.grid;

	filename = new std::string(*a.filename); // associated file name (if we have one)

	title = new std::string(*a.title); // header title

	var_name = new std::string[n_vars];

	for (int i = 0; i < n_vars; ++i)
	{
		var_name[i] = a.var_name[i]; //variable we are using name
	}

	//allocate_grid();
	//grid = new templt[(get_nxyz()) + 1];
	grid.reset(new templt[(get_nxyz()) + 1]);

	//grid[0] = -999.25; // just to mark the fortran type grids
	for (long long i = 1; i <= get_nxyz(); i++)
	{
		grid[i] = a.grid[i];
	}

	return *this;
};


template <class templt> GSLibGrid<templt>::~GSLibGrid()
{//destructor frees memory for the grid
	
	//delete[] this->grid;

	//delete	filename; 
	
	if (title != NULL) 
	{ delete title; }

	if (var_name != NULL) 
	{ delete[] var_name; }
	
	if (filename != NULL)
	{ delete filename; }

	
};

//template <class templt> void GSLibGrid<templt>::allocate_grid()
//{
//		// allocating +1 for fortran 
//	this->grid = new templt [(this->grid_pars.get_nxyz())+1];
//
//};

// READS AND WRITES


//template <class templt> void GSLibGrid<templt>::put_grid_pars(GSLibGridPars& a){ // copy from a GSLibGridPars to a specific grid
//		nx=a.nx;ny=a.ny;nz=a.nz;//number
//		ox=a.ox;oy=a.oy;oz=a.oz;//origin
//		xsiz=a.xsiz;ysiz=a.ysiz;zsiz=a.zsiz; //size
//		nxy=a.nxy;
//		nxyz=a.nxyz;
//	}

//template <class templt> void GSLibGrid<templt>::copy_grid(GSLibGrid * a){// copy all properties from one grid to other (except grid itself)
//	int _getnxyz = a->get_nxyz();
//
//	for (int i= 1; i <= _getnxyz; i++) grid[i] = a->grid[i];
//}; //also keep the grid

template <class templt> void GSLibGrid<templt>::read_grid(){ // reads grid from file and puts it on *grid
		

		std::ifstream fdata(filename->c_str());

		long long nxyz = get_nxyz();

		if (!fdata.is_open()){
			logger << filename << " failed to open!";
			return ;
		}

		int icollvm =1 ; // TO DO remove this, we must have automatically TODO see par file

		if (headerflag==1)
			read_header(fdata);

		//float *var = new float[n_vars]; // if we have no header this fails - TODO see par file

		float *var = new float[n_vars]; 

		valid_points = 0;
		for (long long l = 1; l <= nxyz; l++){
			
			for (unsigned int j=0; j<n_vars; j++) 
				fdata >> var[j]; 
			grid[l] = (templt)var[icollvm-1]; //todo icollvm
			
			if(grid[l]!=null_data ) 
				++valid_points;
		
		}

		delete[] var;
		fdata.close();
	};

template <class templt> void GSLibGrid<templt>::read_grid_binary_sgems(){ // reads grid from file and puts it on *grid
		

		std::ifstream fdata(filename->c_str(),std::ios::binary);


		//unsigned int nxyz = get_nxyz();

		if (!fdata.is_open()){
			logger << filename << " failed to open!";
			return ;
		}

		//int icollvm =1 ; // TO DO remove this, we must have automatically TODO see par file
		
		// Write a header with a "magic number" and the grid type
		//stream << (Q_UINT32)0xB211175D;
			
		// first, load the magic number
		
		char magic[4];
		fdata.read (magic, 4);


		//then object_type


		unsigned int n_chars;
		fdata.read ((char*)&n_chars,4);
		endian_swap(n_chars);


		char *buffer = new char[n_chars];
		fdata.read (buffer,n_chars);
		std::string classname(buffer);
		delete [] buffer;

		fdata.read ((char*)&n_chars,4);
		endian_swap(n_chars);

		char *buffer2 = new char[n_chars];
		fdata.read (buffer2,n_chars);
		std::string title(buffer2);
		delete [] buffer2;
		

		//now read some stuff


		unsigned int version; 
		unsigned int nx, ny, nz;

		fdata.read((char*)&version, sizeof(unsigned int));
		//endian_swap(version);
		
		fdata.read((char*)&nx, sizeof(unsigned int));
		endian_swap(nx);
		fdata.read((char*)&ny, sizeof(unsigned int));
		endian_swap(ny);
		fdata.read((char*)&nz, sizeof(unsigned int));
		endian_swap(nz);

		set_nx(nx);
		set_ny(ny);
		set_nz(nz);

		set_nxy( nx*ny);
		set_nxyz( nx*ny*nz);

		float sizex, sizey, sizez;

		fdata.read((char*)&sizex, sizeof(float));
		sizex=float_endian_swap(sizex);
		fdata.read((char*)&sizey, sizeof(float));
		sizey=float_endian_swap(sizey);
		fdata.read((char*)&sizez, sizeof(float));
		sizez=float_endian_swap(sizez);

		set_xsiz(sizex);
		set_ysiz(sizey);
		set_zsiz(sizez);

		float ox, oy, oz;

		fdata.read((char*)&ox, sizeof(float));
		ox=float_endian_swap(ox);
		fdata.read((char*)&oy, sizeof(float));
		oy=float_endian_swap(oy);
		fdata.read((char*)&oz, sizeof(float));
		oz=float_endian_swap(oz);

		set_ox(ox);
		set_oy(oy);
		set_oz(oz);

		//unsigned int n_vars;
		fdata.read((char*)&n_vars, sizeof(unsigned int));
		endian_swap(n_vars);


		var_name = new std::string[n_vars];
		for( unsigned int i = 0; i < n_vars; i++ ) {
			
			fdata.read ((char*)&n_chars,4);
			endian_swap(n_chars);

			char *buffer3 = new char[n_chars];
			fdata.read (buffer3,n_chars);
			var_name[i]=buffer3;
			delete [] buffer3;
		}

		//grid = NULL; 
	
		//grid = new templt [(get_nxyz())+1]; 
		grid.reset(new templt[get_nxyz() + 1]);

		float value;

		for( unsigned int j = 0; j < n_vars; ++j ) 
		{
			valid_points = 0;

			for (long long	 i= 1;i<=get_nxyz();++i)
			{ 
		
				fdata.read((char*)&value, sizeof(float));
				value=float_endian_swap(value);
				if (value == -9966699)
					value = null_data;
				grid[i]=(templt)value;
				if (grid[i] != null_data)
					++valid_points;
			}
		}

		fdata.close();
	};

template <class templt> void GSLibGrid<templt>::write_grid_ascii()
    {
        std::ofstream fout;
		
        fout.open(filename->c_str());

		
        if (!fout.is_open())
		{
			logger << filename << " failed to open!";
			return ;
		}


        write_header(fout);

        // now write the grid data
		// fortran like, TODO PREPARE FOR C++
		long long here_nxyz = get_nxyz();
        
		for (long long i=1; i<here_nxyz+1;i++)
		{ 
			fout << std::setprecision (4) << std::fixed  << grid[i] << "\n";
		}
		fout.close();



    };//puts *grid on file

template <class templt> void GSLibGrid<templt>::write_grid_binary() // BIG ENDIAN
	{
		std::ofstream fout;


        fout.open(filename->c_str(), std::ios::out | std::ios::binary);


        if (!fout.is_open())
		{
			logger << filename << " failed to open!";
			return ;
		}


        // PPL linhas de header fixo, com tamanho 120, para escrita binaria
        fout << "VAR_PRI_Simulation_1                                                                                                    ";
        char text[120];
        int n;
        n = sprintf(text, "%4d %4d %4d %4d %100s", 1, get_nx(), get_ny(), get_nz(), " ");
        fout << text;
        fout << "value                                                                                                                   ";


        // now write the grid data
		// fortran like, TODO PREPARE FOR C++

        for (long long i=1; i<get_nxyz()+1;i++)
		{

            // PPL nao portavel
            union
            {
              float f;
              unsigned char b[4];
            } dat1, dat2;

            dat1.f = grid[i];
            dat2.b[0] = dat1.b[3];
            dat2.b[1] = dat1.b[2];
            dat2.b[2] = dat1.b[1];
            dat2.b[3] = dat1.b[0];

            fout.write((const char *)&dat2.f, sizeof(float));
        }


        fout.close();

	};

template <class templt> void GSLibGrid<templt>::write_grid_binary_sgems() 
	{
		std::ofstream fout;
		

        fout.open(filename->c_str(), std::ios::out | std::ios::binary);


        if (!fout.is_open()){
			logger << filename << " failed to open!";
			return ;
		}

		//

				// Write a header with a "magic number" and the grid type
		//stream << (Q_UINT32)0xB211175D;
			
		// first, load the magic number
		//unsigned int magic =0x;
		unsigned int magic =0x5D1711B2;
		fout.write((const char *)&magic,4);
		

		unsigned int ack=6;
		endian_swap(ack);
		fout.write((const char *)&ack,4);

		//then object_type
		std::string cgrd="Cgrid";

		fout.write(cgrd.c_str(),cgrd.size()+1);//+1stores the null


		unsigned int buffer= (unsigned int) title->size()+1;
		endian_swap(buffer);
		fout.write((const char *)&buffer,4);


		//we start the grid name

		//fout.write((const char *)&null_buffer,4);
		//fdata.read (buffer4, 4);

	//	std::string gridname="logistic_predicted";

	//	fout.write(gridname.c_str(),gridname.size()+1);
		fout.write(title->c_str(),title->size()+1);//+1stores the null
		//gridname_not var, but ok for now

		//now read some stuff


		unsigned int version=100; 
		//unsigned int nx, ny, nz;
		endian_swap(version);
		fout.write((const char *)&version,4);

		unsigned int nx= get_nx();
		unsigned int ny= get_ny();
		unsigned int nz= get_nz();

		endian_swap(nx);
		endian_swap(ny);	
		endian_swap(nz);

		fout.write((const char *)&nx,4);
		fout.write((const char *)&ny,4);
		fout.write((const char *)&nz,4);

		float sizex= get_xsiz();
		float sizey= get_ysiz();
		float sizez= get_zsiz();

		sizex=float_endian_swap(sizex);
		sizey=float_endian_swap(sizey);	
		sizez=float_endian_swap(sizez);


		fout.write((const char *)&sizex,4);
		fout.write((const char *)&sizey,4);
		fout.write((const char *)&sizez,4);


		float ox= get_ox();
		float oy= get_oy();
		float oz= get_oz();

		ox=float_endian_swap(ox);
		oy=float_endian_swap(oy);	
		oz=float_endian_swap(oz);


		fout.write((const char *)&ox,4);
		fout.write((const char *)&oy,4);
		fout.write((const char *)&oz,4);

		//grid has one property only (at least for now)
		//unsigned int n_props=1;
		endian_swap(n_vars);

		fout.write((const char *)&n_vars,4);
		//it is wrtten, lets keep it
		endian_swap(n_vars);

		for( unsigned int j = 0; j < n_vars; j++ ){
			unsigned int bs=(unsigned int)var_name[j].size()+1;
			endian_swap(bs);
			fout.write((const char *)&bs,4);
		}
		for( unsigned int j = 0; j < n_vars; j++ )
			fout.write(var_name[j].c_str(),var_name[j].size()+1);

		//std::vector< std::string > property_name( n_props );
		//for( unsigned int i = 0; i < n_props; i++ ) {
		//	fdata.read (buffer4,4); //this should be a BEL  - that is 7
		//
		//	property_name[i]="";
		//	buffer=1;
		//	while (buffer!=0){
		//		fdata.read (&buffer,1);
		//		property_name[i].push_back(buffer);
		//	}
		//}

		
		float value;

		long long nxyz=get_nx()*get_ny()*get_nz();

		for( unsigned int j = 0; j < n_vars; j++ ) { //one because we are loading only one property
			for (long long i= 1;i<=nxyz;++i){
				
				value=grid[i];

				if (value == null_data)
					value = -9966699;

				value = float_endian_swap(value);


				fout.write((const char *)&value,4);

			}
		}







  //      // PPL linhas de header fixo, com tamanho 120, para escrita binaria
  //      fout << "VAR_PRI_Simulation_1                                                                                                    ";
  //      char text[120];
  //      int n;
  //      n = sprintf(text, "%4d %4d %4d %4d %100s", 1, get_nx(), get_ny(), get_nz(), " ");
  //      fout << text;
  //      fout << "value                                                                                                                   ";


  //      // now write the grid data
		//// fortran like, TODO PREPARE FOR C++

  //      for (i=1; i<get_nxyz()+1;i++) {

  //          // PPL nao portavel
  //          union
  //          {
  //            float f;
  //            unsigned char b[4];
  //          } dat1, dat2;

  //          dat1.f = grid[i];
  //          dat2.b[0] = dat1.b[3];
  //          dat2.b[1] = dat1.b[2];
  //          dat2.b[2] = dat1.b[1];
  //          dat2.b[3] = dat1.b[0];

  //          fout.write((const char *)&dat2.f, sizeof(float));
  //      }


        fout.close();

	};

template <class templt> void GSLibGrid<templt>::read_header(std::ifstream& f_in )
{ // reads gslib header
	
		std::string buffer;

		//f_in >> title;

		title = new std::string;
		std::getline(f_in,*title);
		//f_in >> n_vars;
		//f_in >> var_name;
		std::getline(f_in,buffer);
		n_vars= atoi(buffer.c_str());

		var_name =  new std::string[n_vars];
		for (unsigned int l = 0; l< n_vars;++l){
			//f_in >> var_name[l];
			std::getline(f_in,var_name[l]);
		}

};

template <class templt> void GSLibGrid<templt>::write_header(std::ofstream& f_out){// writes gslib header

		f_out << *title << "\n";
		f_out << n_vars<< "\n";
		for (unsigned int l = 0; l< n_vars;++l)
		{
			f_out << var_name[l]<< "\n";
		}

	};

//----------------------------------------

template <class templt> void GSLibGrid<templt>::get_average()
	{
		
		double av=0;
		long long nxyz= get_nxyz();
		long long valids = get_valid_points();

		for (long long	 i=1; i<=nxyz; ++i)
		{

			if(grid[i]!=null_data ) av+=grid[i];

		}

		//average=av/std::max(nxyz,1);
		average = av / std::max(valids, (long long)1);
	};

template <class templt> void GSLibGrid<templt>::get_variance()
	{
		double ss=0;
		long long nxyz= get_nxyz();

		for (long long i=1; i<=nxyz; ++i)
		{

			if(grid[i]!=null_data ) ss+=pow((double)(grid[i]-average),2.0);;

		}

		variance=ss/std::max(nxyz-1 ,(long long)1); 

	};

template <class templt> int GSLibGrid<templt>::get_valid_points()
{
	valid_points = 0;
	long long nxyz = get_nxyz();
	for (long long i = 1; i <= nxyz; ++i)
	{
		if (grid[i] != null_data) ++valid_points;
	}

	return valid_points;
}
//template <class templt> bool GSLibGrid<templt>::has_nulls()
//	{ // returns true if we have nulls, false otherwise
//		bool a;
//		(nx*ny*nz)==valid_points ? a=0:a=1; return a;
//	};

template <class templt> void GSLibGrid<templt>::zero_the_grid()
{ // copy from a GSLibGridPars to a specific grid
	for (long long i= 0; i <= get_nxyz(); i++) grid[i] = 0;
}


#endif

																																																																																																																									