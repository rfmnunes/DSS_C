#ifndef SUPERBLOCKS_INCLUDED
#define SUPERBLOCKS_INCLUDED

#define  MAXSBX                       8
#define  MAXSBY                       8
#define  MAXSBZ                       8


#include "hardata.h"
#include "variogram.h"
#include "gslib_grid.h"


class search_results
{
public:
	unsigned int n_close_data;
	long    *inclose;
	int nodmax;
	int ndmax;
	//cl_data_point<double> *found_sample_values;
	std::vector<cl_data_point<double>> found_sample_values;

	search_results(int n_data)
	{
		n_close_data = 0;
		inclose = new long[n_data];
		//found_sample_values = new cl_data_point<double>[n_data];
	}
	~search_results() 
	{ 
		delete[] inclose; 
		//delete[] found_sample_values;
	}
};



class SuperBlocksPars
{

public:

	GSLibGridPars *superblock_grid;

	GSLibGrid< std::vector<int> > *super_blks;

	int	*n_samples_in_SB;
	cl_coord<int> *ixyzsbtosr;
 
	unsigned int nsbtosr;

	int noct;


	void setsupr(misc& utils, GSLibGridPars& grid_def,std::shared_ptr<HarddataPars> &harddata){
		long ii;
		//unsigned long i;
		float a[1];
		//float *tmp;
		int nx_1,ny_1,nz_1;
		float ox_1,oy_1,oz_1,xsiz_1,ysiz_1,zsiz_1;


		//tmp= new float[harddata.harddata_distr.n_data+1];
		boost::scoped_array<float> tmp(new float[harddata->harddata_distr.n_data + 1]);

		// Establish the number and size of the super blocks:
		nx_1   = std::min(grid_def.get_nx(),(unsigned int) (MAXSBX));
		ny_1   = std::min(grid_def.get_ny(),(unsigned int) (MAXSBY));
		nz_1   = std::min(grid_def.get_nz(),(unsigned int) (MAXSBZ));

		n_samples_in_SB = new int[(nx_1*ny_1*nz_1) + 1];

		xsiz_1 = grid_def.get_nx()*grid_def.get_xsiz()/nx_1;
		ysiz_1 = grid_def.get_ny()*grid_def.get_ysiz()/ny_1;
		zsiz_1 = grid_def.get_nz()*grid_def.get_zsiz()/nz_1;
		ox_1  = (float)( (grid_def.get_ox()-0.5*grid_def.get_xsiz()) +0.5*xsiz_1);
		oy_1  = (float)( (grid_def.get_oy()-0.5*grid_def.get_ysiz()) +0.5*ysiz_1);
		oz_1  = (float)( (grid_def.get_oz()-0.5*grid_def.get_zsiz()) +0.5*zsiz_1);

		superblock_grid= new GSLibGridPars(nx_1,ny_1,nz_1,ox_1,oy_1,oz_1,xsiz_1,ysiz_1,zsiz_1);

		// builds a grid to store the values of points belonging there
		//GSLibGrid<std::vector<int>> super_blks(*superblock_grid);

		// Initialize the extra super block array to zeros:
		for (unsigned long i=0; i<=superblock_grid->get_nxyz(); ++i)
			n_samples_in_SB[i] = 0;

		// Loop over all the data assigning the data to a super block and
		// accumulating how many data are in each super block:
		for (unsigned long i=1; i<=harddata->harddata_distr.n_data; ++i) {

			// we have two different cood systems, check if wells are outside or inside
			//if (superblock_grid->is_xyz_inside(harddata.point[i].xyz))
			//{

				ii = superblock_grid->get_index_from_xyz_f(harddata->point[i].xyz);

				tmp[i]   = (float)ii;
				++n_samples_in_SB[ii];
			//}		
			// i can make the superblocks grid a grid with a vector of the indexes of datapoints in the block.
			
				//super_blks.grid[ii].push_back(i)
				
		}







		// Sort the data by ascending super block number:
		//nsort = 4;
		//utils.sortit(1,harddata.harddata_distr.n_data,tmp,nsort, &harddata.x[0], &harddata.y[0], &harddata.z[0],&harddata.vr[0],a,a,a);
		
		utils.sortit(1,harddata->harddata_distr.n_data,&tmp[0],1, &harddata->point[0], a, a, a, a, a, a);



		// Set up array n_samples_in_SB with the starting address of the block data:
		for (unsigned long i=0; i<=(superblock_grid->get_nxyz()-1); ++i)
			n_samples_in_SB[i+1] += n_samples_in_SB[i];


		//delete[] tmp;

		return;
	};

	void new_setsupr(misc& utils, GSLibGridPars& grid_def, HarddataPars& harddata) {
		long ii;
		//unsigned long i;
		//float a[1];
		float *tmp;
		int nx_1, ny_1, nz_1;
		float ox_1, oy_1, oz_1, xsiz_1, ysiz_1, zsiz_1;


		tmp = new float[harddata.harddata_distr.n_data + 1];

		// Establish the number and size of the super blocks:
		nx_1 = std::min(grid_def.get_nx(), (unsigned int)(MAXSBX));
		ny_1 = std::min(grid_def.get_ny(), (unsigned int)(MAXSBY));
		nz_1 = std::min(grid_def.get_nz(), (unsigned int)(MAXSBZ));

		n_samples_in_SB = new int[(nx_1*ny_1*nz_1) + 1];

		xsiz_1 = grid_def.get_nx()*grid_def.get_xsiz() / nx_1;
		ysiz_1 = grid_def.get_ny()*grid_def.get_ysiz() / ny_1;
		zsiz_1 = grid_def.get_nz()*grid_def.get_zsiz() / nz_1;
		ox_1 = (float)((grid_def.get_ox() - 0.5*grid_def.get_xsiz()) + 0.5*xsiz_1);
		oy_1 = (float)((grid_def.get_oy() - 0.5*grid_def.get_ysiz()) + 0.5*ysiz_1);
		oz_1 = (float)((grid_def.get_oz() - 0.5*grid_def.get_zsiz()) + 0.5*zsiz_1);

		superblock_grid = new GSLibGridPars(nx_1, ny_1, nz_1, ox_1, oy_1, oz_1, xsiz_1, ysiz_1, zsiz_1);

		// builds a grid to store the values of points belonging there
		super_blks = new GSLibGrid< std::vector<int> >(*superblock_grid);


		for (unsigned long i = 1; i <= harddata.harddata_distr.n_data; ++i) {

			ii = superblock_grid->get_index_from_xyz_f(harddata.point[i].xyz);
			super_blks->grid[ii].push_back(i);

		}

		return;
	};

	void new_picksupr(VariogramPars& variogram) {
		double hsqd, shortest;

		float xo, yo, zo, xdis, ydis, zdis;

		ixyzsbtosr = new cl_coord<int>[8 * (MAXSBX*MAXSBY*MAXSBZ)];

		cl_coord<double> zerocoord(0, 0, 0);

		int nx = superblock_grid->get_nx();
		int ny = superblock_grid->get_ny();
		int nz = superblock_grid->get_nz();

		// Main loop over all possible super blocks:
		nsbtosr = 0;
		for (int i = -(nx - 1); i <= (nx - 1); ++i) {
			for (int j = -(ny - 1); j <= (ny - 1); ++j) {
				for (int k = -(nz - 1); k <= (nz - 1); ++k) {
					xo = i*superblock_grid->get_xsiz();
					yo = j*superblock_grid->get_ysiz();
					zo = k*superblock_grid->get_zsiz();

					// Find the closest distance between the corners of the super blocks:
					shortest = 1.0E21;

					for (int i1 = -1; i1 <= 1; ++i1) {
						for (int j1 = -1; j1 <= 1; ++j1) {
							for (int k1 = -1; k1 <= 1; ++k1) {
								for (int i2 = -1; i2 <= 1; ++i2) {
									for (int j2 = -1; j2 <= 1; ++j2) {
										for (int k2 = -1; k2 <= 1; ++k2) {
											if (i1 != 0 || j1 != 0 || k1 != 0 || i2 != 0 || j2 != 0 || k2 != 0) {
												xdis = (float)((i1 - i2)*0.5*superblock_grid->get_xsiz() + xo);
												ydis = (float)((j1 - j2)*0.5*superblock_grid->get_ysiz() + yo);
												zdis = (float)((k1 - k2)*0.5*superblock_grid->get_zsiz() + zo);
												cl_coord<double> xyzdis(xdis, ydis, zdis); // dist entre 2 blocos

												hsqd = variogram.search_radius->sqdist(zerocoord, xyzdis);//0.0,0.0,0.0,xdis,ydis,zdis,0);
												if (hsqd < shortest) shortest = hsqd;
											}
										}
									}
								}
							}
						}
					}

					// Keep this super block if it is close enoutgh:
					if (shortest <= variogram.radsqd)
					{
						ixyzsbtosr[nsbtosr].x = i;
						ixyzsbtosr[nsbtosr].y = j;
						ixyzsbtosr[nsbtosr].z = k;
						++nsbtosr;
					}
				}
			}
		}

		return;
	};


	void new_srchsupr(misc& utils, HarddataPars& harddata, VariogramPars& variogram,
		cl_coord<double> xyzloc, search_results& srch_res, variogram_structure& search_radius, int ndmax)
	{

		double hsqd;
		float  hh, c[1];


		cl_coord<int> ixyzsup;

		//inclose = new long[harddata.harddata_distr.n_data];

		//int nt;
		float *temporary = new float[harddata.harddata_distr.n_data/*+1*/];
		//tmp;

		for (unsigned int isup = 0; isup < harddata.harddata_distr.n_data/* + 1 + 1*/; ++isup)
		{
			temporary[isup] = 0;
		}
		
		// Determine the super block location of point being estimated:
		cl_coord<int> ijk = superblock_grid->get_ijk_form_xyz_f(xyzloc);


		//perhaps prepare a vector of the superblock search order
		for (unsigned int isup = 0; isup < nsbtosr; ++isup)
		{

			// Is this super block within the grid system:

			ixyzsup = ijk + ixyzsbtosr[isup];

			if (!superblock_grid->is_ijk_inside(ixyzsup)) continue;

			int curr_spr_blk = superblock_grid->get_index_from_ijk_f(ixyzsup);

			// now get data from this superblock
			for (int n = 0; n < super_blks->grid[curr_spr_blk].size(); n++)
			{

				int currsample = super_blks->grid[curr_spr_blk].at(n);
				cl_coord<double> sample_coord = harddata.point[currsample].xyz;

				// Check squared distance:

				hsqd = search_radius.sqdist(xyzloc, sample_coord);

				if (hsqd > variogram.radsqd)
					continue;

				// Accept this sample:

				srch_res.inclose[srch_res.n_close_data] = (long)currsample;
				temporary[srch_res.n_close_data] = (float)hsqd;
				++srch_res.n_close_data;
			}
			
		}

		// Sort the nearby samples by distance to point being estimated:
		utils.sortit(0, srch_res.n_close_data - 1, temporary, 1, (float *)srch_res.inclose, c, c, c, c, c, c);


		// If we aren't doing an octant search then just return:
		if (noct > 0)
		{
			//	delete[] temporary;
			//	return;
			//}
			int iq, inoct[8];
			int nt, na;
			float dx, dy, dz;

			// Partition the data into octants:
			for (int i = 0; i < 8; ++i)
				inoct[i] = 0;

			// Now pick up the closest samples in each octant:
			nt = 8 * noct;
			na = 0;

			for (unsigned int j = 0; j < srch_res.n_close_data; ++j)
			{
				int i = srch_res.inclose[j];
				hh = temporary[j];
				dx = harddata.point[i].xyz.x - xyzloc.x;
				dy = harddata.point[i].xyz.y - xyzloc.y;
				dz = harddata.point[i].xyz.z - xyzloc.z;
				if (dz < 0.)
				{
					iq = 7;
					if (dx <= 0.0 && dy < 0.0) iq = 4;
					if (dx > 0.0 && dy >= 0.0) iq = 5;
					if (dx < 0.0 && dy <= 0.0) iq = 6;
				}
				else
				{
					iq = 3;
					if (dx <= 0.0 && dy > 0.0) iq = 0;
					if (dx > 0.0 && dy >= 0.0) iq = 1;
					if (dx < 0.0 && dy <= 0.0) iq = 2;
				}


				++inoct[iq];

				// Keep this sample if the maximum has not been exceeded:
				if (inoct[iq] <= noct)
				{

					srch_res.inclose[na] = i;
					temporary[na] = hh;
					++na;
					if (na == nt) break;
				}
			}

			// End of data selection. Compute number of informed octants and return:
			srch_res.n_close_data = na;
		}
		
		delete[] temporary;
		
		return;

	};

	void picksupr(std::shared_ptr<VariogramPars> &variogram){
		double hsqd,shortest;
	
		float xo,yo,zo,xdis,ydis,zdis;

		ixyzsbtosr = new cl_coord<int>[8*(MAXSBX*MAXSBY*MAXSBZ)];

		cl_coord<double> zerocoord(0,0,0);

		int nx= superblock_grid->get_nx();
		int ny= superblock_grid->get_ny();
		int nz= superblock_grid->get_nz();

		// Main loop over all possible super blocks:
		nsbtosr = 0;
		for (int i=-(nx-1); i<=(nx-1); ++i) {
			for (int j=-(ny-1); j<=(ny-1); ++j) {
				for (int k=-(nz-1); k<=(nz-1); ++k) {
					xo = i*superblock_grid->get_xsiz();
					yo = j*superblock_grid->get_ysiz();
					zo = k*superblock_grid->get_zsiz();

					// Find the closest distance between the corners of the super blocks:
					shortest = 1.0E21;

					for (int i1=-1; i1<=1; ++i1) {
						for (int j1=-1; j1<=1; ++j1) {
							for (int k1=-1; k1<=1; ++k1) {
								for (int i2=-1; i2<=1; ++i2) {
									for (int j2=-1; j2<=1; ++j2) {
										for (int k2=-1; k2<=1; ++k2) {
											if(i1 != 0 || j1 != 0 || k1 != 0 || i2 != 0 || j2 != 0 || k2 != 0) {
												xdis = (float)((i1-i2)*0.5*superblock_grid->get_xsiz() + xo);
												ydis = (float)((j1-j2)*0.5*superblock_grid->get_ysiz() + yo);
												zdis = (float)((k1-k2)*0.5*superblock_grid->get_zsiz() + zo);
												cl_coord<double> xyzdis(xdis,ydis,zdis); // dist entre 2 blocos

												hsqd = variogram->search_radius->sqdist( zerocoord,xyzdis);//0.0,0.0,0.0,xdis,ydis,zdis,0);
												if (hsqd < shortest) shortest = hsqd;
											}
										}
									}
								}
							}
						}
					}

					// Keep this super block if it is close enoutgh:
					if (shortest <= variogram->radsqd) 
					{				
						ixyzsbtosr[nsbtosr].x =i;
						ixyzsbtosr[nsbtosr].y =j;
						ixyzsbtosr[nsbtosr].z =k;
						++nsbtosr;
					}
				}
			}
		}

		return;
	};

	void srchsupr(misc& utils, std::shared_ptr<HarddataPars> &harddata, std::shared_ptr<VariogramPars> &variogram,
		cl_coord<double>& xyzloc, search_results& srch_res, variogram_structure& search_radius, int ndmax)
	{

		double hsqd;
		float   c[1];//hh,
		int i, ii, nums;


		cl_coord<int> ixyzsup;

		//inclose = new long[harddata.harddata_distr.n_data];

		//int nt;
		float *temporary = new float[harddata->harddata_distr.n_data/*+1*/];
		//tmp;

		for (unsigned int isup = 0; isup < harddata->harddata_distr.n_data/* + 1 + 1*/; ++isup)
		{
			temporary[isup] = 0;
		}
		// Determine the super block location of point being estimated:

		cl_coord<int> ijk = superblock_grid->get_ijk_form_xyz_f(xyzloc);


		std::vector<cl_data_point<double>> temp_found_sample_values;

		// Loop over all the possible Super Blocks:
		//*nclose = 0;
		for (unsigned int isup = 0; isup<nsbtosr; ++isup)
		{

			// Is this super block within the grid system:

			ixyzsup = ijk + ixyzsbtosr[isup];

			if (!superblock_grid->is_ijk_inside(ixyzsup)) continue;
				
		
			//if (ixyzsup.x <= 0 || ixyzsup.x > superblock_grid->get_nx() ||
			//	ixyzsup.y <= 0 || ixyzsup.y > superblock_grid->get_ny() ||
			//	ixyzsup.z <= 0 || ixyzsup.z > superblock_grid->get_nz()) continue;

			// Figure out how many samples in this super block:
			ii = superblock_grid->get_index_from_ijk_f(ixyzsup);

			if (ii == 1)
			{
				nums = n_samples_in_SB[ii];
				i = 0;
			}
			else
			{
				nums = n_samples_in_SB[ii] - n_samples_in_SB[ii - 1];
				i = n_samples_in_SB[ii - 1];
			}


			// Loop over all the data in this super block:
			for (int l = 1; l <= nums; ++l)
			{
				++i;

				// Check squared distance:

				hsqd = search_radius.sqdist(xyzloc, harddata->point[i].xyz);

				if (hsqd > variogram->radsqd)
					continue;

				// Accept this sample:

				srch_res.inclose[srch_res.n_close_data] = (long)i;
				temporary[srch_res.n_close_data] = (float)hsqd;
				
				//srch_res.found_sample_values[srch_res.n_close_data] = harddata->point[i];
				temp_found_sample_values.push_back(harddata->point[i]);

				++srch_res.n_close_data;
			}
			//if (srch_res.n_close_data > ndmax){
			////	srch_res.n_close_data = ndmax;
			//	break; // leave the cycle 
			//}
		}

		// Sort the nearby samples by distance to point being estimated:
		utils.sortit(0, srch_res.n_close_data - 1, temporary, 2, (float *)srch_res.inclose, &temp_found_sample_values[0], c, c, c, c, c);


		// If we aren't doing an octant search then just return:
		if (noct > 0)
		{
		//	delete[] temporary;
		//	return;
		//}
			int iq, inoct[8];
			int nt, na ;
			double dx, dy, dz;

			// Partition the data into octants:
			for (int i = 0; i < 8; ++i)
				inoct[i] = 0;

			// Now pick up the closest samples in each octant:
			nt = 8 * noct;
			na = 0;

			for (unsigned int j = 0; j < srch_res.n_close_data; ++j)
			{
				int i = srch_res.inclose[j];
				//hh = temporary[j];
				dx = harddata->point[i].xyz.x - xyzloc.x;
				dy = harddata->point[i].xyz.y - xyzloc.y;
				dz = harddata->point[i].xyz.z - xyzloc.z;
				if (dz < 0.)
				{
					iq = 7;
					if (dx <= 0.0 && dy < 0.0) iq = 4;
					if (dx > 0.0 && dy >= 0.0) iq = 5;
					if (dx < 0.0 && dy <= 0.0) iq = 6;
				}
				else
				{
					iq = 3;
					if (dx <= 0.0 && dy > 0.0) iq = 0;
					if (dx > 0.0 && dy >= 0.0) iq = 1;
					if (dx < 0.0 && dy <= 0.0) iq = 2;
				}


				++inoct[iq];

				// Keep this sample if the maximum has not been exceeded:
				if (inoct[iq] <= noct)
				{
					
					srch_res.inclose[na] = i;

					//index = srch_res.inclose[j];
					//srch_res.found_sample_values[na] = harddata->point[i];
					srch_res.found_sample_values.push_back(harddata->point[i]);

 					//temporary[na] = hh;
					++na;
					if (na == nt) break;
				}
			}

			// End of data selection. Compute number of informed octants and return:
			srch_res.n_close_data = na;
		}
		else 
		{ // the temp_faound samples passes to official
			for (i = 0; i<ndmax;++i)
				if (temp_found_sample_values.size()>i)
					srch_res.found_sample_values.push_back(temp_found_sample_values[i]);
		}
		
		delete[] temporary;
		return;

	};
};

#endif