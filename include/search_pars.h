#ifndef SEARCH_INCLUDED
#define SEARCH_INCLUDED

#include "boost_include.h"
#include "superblocks.h"
#include "variogram.h"
//#include "SSdir.h"
//#define _REPRODUCTIBLE 	// if active does the exact path as the sequential, but is slower

//#define  TINY                    0.0001
//#define  EPSLON                 1.0E-20






class Nodes_search{

public:

	int ncnode; //number of found close nodes // can probably be replaced by a call to size one of the vectors
	std::vector<long> *icnode; // index in lookup table //can be changed by a call to a method in the covariance table that recieves the points and gives covariance ?
	std::vector<long> *close_node_index; // what it says on the tin
	std::vector<float> *close_node_value; // what it says on the tin


	//bool previous_needed;

	//short unsigned int is_needed_int;

	Nodes_search(){	
		close_node_index = new std::vector<long>;
		icnode = new std::vector<long>;
		close_node_value = new std::vector<float>;
	};


	~Nodes_search()
	{
		delete close_node_index;
		delete icnode; 
		delete close_node_value;
	};



	void  Srchnd(boost::mutex& m1, boost::condition_variable *hEvent,
		bool *hEvent_aux, long nodemax, GSLibGridPars& grid, std::shared_ptr<HarddataPars> &harddata,
		std::shared_ptr<VariogramPars> &vario, boost::shared_array<float> &sim, boost::shared_array<unsigned int> &order,
		int nodmax, search_results& srch_res , float nosvalue,
		long long  in,int numprocess, long thread_in[], misc& utils,int sstrat, int noct)
	{
		long long inde;
		
		cl_coord<int> ijk_close;
		cl_coord<int> coord_ones(1, 1, 1);

		//cl_coord<int> vario_nct(vario.nctx,vario.ncty,vario.nctz);
		//long index= order[in];

		cl_coord<int> ijk_current = grid.get_ijk_from_index_f(order[in]);

		// Consider all the nearby nodes until enough have been found
		ncnode = 0;

		int ninoct[9];
		int octants_passed = 0;
		if (noct > 0) 
		{
			
			for (unsigned int i = 1; i <= 8; ++i) 
			{
				ninoct[i] = 0; // preparing octants
				 
			}
		}


#ifdef _REPRODUCTIBLE
		int ik;

		long long min_in = grid.get_nxyz() + 1;

		for (ik = 0; ik < numprocess; ik++) {

			if ((thread_in)[ik] < min_in && (thread_in)[ik] >= 0)
				min_in = (thread_in)[ik];
		}

#endif // DEBUG
		
		cl_coord<int> recentered_ijk_close= ijk_current - vario->nctxyz - coord_ones;

		for (long il = 1; il < vario->nlooku; ++il) 
		{


			if (ncnode == nodmax || octants_passed == 8)
			{
				break;
			}

			//ijk_close = ijk_current + (vario.ijk_node[il] - vario.nctxyz - coord_ones);
			ijk_close = recentered_ijk_close + vario->ijk_node[il] ;

			if (!grid.is_ijk_inside(ijk_close)) 
			{
				continue;
			}

			inde = grid.get_index_from_ijk_f(ijk_close);

#ifdef _REPRODUCTIBLE

			for (long ik1 = in - 1; ik1 >= min_in; ik1--) //corre todos os entre o minimo e este 

				if ((ik1 - nodemax > 0) && (order[ik1] == inde))
				{ // precisa de um, vai esperar 

					{
						boost::mutex::scoped_lock lk(m1);

						while (!hEvent_aux[ik1 - nodemax - 1])
						{ // espera pelo evento
							hEvent[ik1 - nodemax - 1].wait(lk); // evento, continua
						}
					}

				}
#endif


			// Previsously simulated node is required
			if ((fabs(sim[inde] - nosvalue) > TINY))
			{ // TO DO CREATE GRIDWITH SIM VAL


				// octants
				if (noct > 0) {
					int iq;
					cl_coord<int> idxyz = ijk_current - ijk_close;
		/*			idx = ix - i;
					idy = iy - j;
					idz = iz - k;
					*/
					if (idxyz.z > 0)
					{
						iq = 4;
						if (idxyz.x <= 0 && idxyz.y > 0) iq = 1;
						if (idxyz.x > 0 && idxyz.y >= 0) iq = 2;
						if (idxyz.x < 0 && idxyz.y <= 0) iq = 3;
					}
					else 
					{
						iq = 8;
						if (idxyz.x <= 0 && idxyz.y > 0) iq = 5;
						if (idxyz.x > 0 && idxyz.y >= 0) iq = 6;
						if (idxyz.x < 0 && idxyz.y <= 0) iq = 7;
					}
					++ninoct[iq];
					if (ninoct[iq] > noct) 
					{ 
						if (ninoct[iq] == noct + 1)  
							++octants_passed;
						continue;
					}
				}

				// Check for duplicated coordinates (nodes and samples)
				int iteste = 0;

				if (sstrat==0)
				{
					cl_coord<double> xyz_coor = grid.get_xyz_from_ijk_f(ijk_close);
					//index;
					for (unsigned int kk = 0; kk < srch_res.n_close_data; ++kk)
					{
						long index = srch_res.inclose[kk];

						if ((harddata->point[index].xyz == xyz_coor)) 
						{
							iteste = 1;
							break;
						}
					}
				}

				// also test for zero covariance points
				//cl_coord<int> ijk_j = ijk_current + (vario.ijk_node[il] - vario.nctxyz - coord_ones);
				
				cl_coord<int> iijjkk = vario->nctxyz + coord_ones + (ijk_current - ijk_close);

				//if ((double)vario.covtab[iijjkk.x][iijjkk.y][iijjkk.z] == 0.0) // redundant since limited on nlookup 
				//{
				//	iteste = 1;
				//}

				//if(vario.covariance(grid.get_xyz_from_ijk_f(ijk_close), grid.get_xyz_from_ijk_f(ijk_current)) == 0.0)
				//	iteste = 1;


				
				if (iteste == 0) 
				{					
					icnode->push_back(il);
					close_node_index->push_back(inde);
					close_node_value->push_back(sim[inde]);
					++ncnode;
				}
			}
		}

		////have all nclose, lets order by covariance??
		////for that we must have covariances...  lets do R matrix. 
		////index = j - search.n_close_data;
		//float c[1];
		//float *temp_cov = new float[1000001];//[nodmax + 1];
		//for (int j = 1; j <= ncnode; ++j) {
		//int ind = icnode[j];
		//cl_coord<int> ijk_j = ijk_current + (vario.ijk_node[ind] - vario.nctxyz - coord_ones);
		//cl_coord<float> xyz_j = grid.get_xyz_from_index_f(close_node_index[j]);
		//cl_coord<int> iijjkk = vario.nctxyz + coord_ones + (ijk_current - ijk_j);
		//cl_coord<float> current_pnt_xyz = grid.get_xyz_from_index_f(order[in]);
		//if (iijjkk.x < 1 || iijjkk.x >   vario.MAXCTX ||
		//	iijjkk.y < 1 || iijjkk.y >   vario.MAXCTY ||
		//	iijjkk.z < 1 || iijjkk.z >   vario.MAXCTZ)
		//	//k_system.r[j] = (double)vario.covariance(current_pnt_xyz, xyz_j);
		//	temp_cov[j] == 1-(double)vario.covariance(current_pnt_xyz, xyz_j);
		//else
		//	//k_system.r[j] = (double)vario.covtab[iijjkk.x][iijjkk.y][iijjkk.z];
		//	temp_cov[j] = 1-(double)vario.covtab[iijjkk.x][iijjkk.y][iijjkk.z];
		//}
		//// now order  
		//utils.sortit(1, ncnode, temp_cov, 2, icnode, close_node_index, c, c, c, c, c);
		//
		//if (ncnode > nodmax) ncnode = nodmax;
		//
		//for (int j = 1; j <= nodmax; ++j) 
		//{
		//	if (temp_cov[j] == 1) {
		//		ncnode = j - 1;
		//		break;
		//	}
		//
		//}
		//delete[] temp_cov;


		// Return to calling program
		return;
	};


};


class SearchPars{
public:

	int sstrat;

	unsigned int	ndmax,
					ndmin,
					nodmax;

	// create a superblock

	SuperBlocksPars* superblocks;

	int get_search_pars(registry *reg)
	{
		reg_key *k;
		// parsing search section

		k = get_key(reg, (char*)("SEARCH"), (char*)("NDMIN"));
		if (k)
			ndmin = get_int(k);
		else return -1;

		k = get_key(reg, (char*)("SEARCH"), (char*)("NDMAX"));
		if (k)
			ndmax = get_int(k);
		else return -1;
		printf (" Number of hard samples (min/max): %d %d \n",ndmin,ndmax);

		k = get_key(reg, (char*)("SEARCH"), (char*)("NODMAX"));
		if (k)
			nodmax = get_int(k);
		else return -1;
		logger << " Max. previous simulated nodes: "<< utilities::make_string(nodmax) <<"\n";

		k = get_key(reg, (char*)("SEARCH"), (char*)("SSTRAT"));
		if (k)
			sstrat = get_int(k);
		else return -1;
		if (sstrat == 1) ndmax = 0;
		logger << " (0)two-part search/(1)data nodes flag: "<< utilities::make_string(sstrat) <<"\n";

		k = get_key(reg, (char*)("SEARCH"), (char*)("NOCT"));
		if (k)
			superblocks->noct = get_int(k);
		else return -1;
		logger << " Number of samples per octant: "<< utilities::make_string(superblocks->noct) << "\n";

		return 1; // exit ok
	
	}

	std::shared_ptr<search_results> /*search_results**/ search_data(boost::mutex& m1,boost::condition_variable *hEvent,bool *hEvent_aux,long nodemax,
						std::shared_ptr<VariogramPars> &vario, GSLibGridPars& grid, std::shared_ptr<HarddataPars> &harddata,
						Nodes_search& nodes, misc& utils, boost::shared_array<float> &sim,
						boost::shared_array<unsigned int> &order, long in, int numprocess, long thread_in[],
						variogram_structure& search_radius){
		

		cl_coord<double> xyz = grid.get_xyz_from_index_f(order[in]);
		// search stuff 

		//search_results* srch_rslts =new search_results(harddata->harddata_distr.n_data);
		std::shared_ptr<search_results> srch_rslts(new search_results(harddata->harddata_distr.n_data));
		srch_rslts->ndmax = this->ndmax;
		srch_rslts->nodmax = this->nodmax;


		//n_close_data=0;

		boost::mutex critical_sim;
		if (sstrat == 0) { // sample search

			//critical_sim.lock();
			
			
			superblocks->srchsupr (utils, harddata, vario, xyz, *srch_rslts, search_radius, ndmax);
			//superblocks->new_srchsupr(utils, harddata, vario, xyz, *srch_rslts, search_radius, ndmax);


			if (srch_rslts->n_close_data < ndmin) 
			{
				srch_rslts->n_close_data=0;
				//nodes.ncnode=0; //go to ksol also
				//return srch_rslts;
			}
			if (srch_rslts->n_close_data > ndmax)
				srch_rslts->n_close_data=ndmax;

			//critical_sim.unlock();

		}

		// we need the previous node so we must use a lock

		nodes.Srchnd( m1, hEvent,hEvent_aux,nodemax,grid, harddata, vario, sim, order,
			nodmax, *srch_rslts, utils.nosvalue, in, numprocess, thread_in, utils,sstrat, superblocks->noct ); // node search


		return srch_rslts;
	};

	search_results* search_data_krige(std::shared_ptr<VariogramPars> &vario, GSLibGridPars& grid, std::shared_ptr<HarddataPars> &harddata,
		misc& utils, long currpoint, variogram_structure& search_radius)
	{


		cl_coord<double> xyz = grid.get_xyz_from_index_f(currpoint);
		// search stuff 

		search_results* srch_rslts = new search_results(harddata->harddata_distr.n_data);

		if (sstrat == 0) { // sample search

			//superblocks->new_srchsupr(utils, harddata, vario, xyz, *srch_rslts, search_radius, ndmax);
			superblocks->srchsupr(utils, harddata, vario, xyz, *srch_rslts, search_radius, ndmax);

			if (srch_rslts->n_close_data < ndmin) {
				srch_rslts->n_close_data = 0;
				return srch_rslts;
			}
			if (srch_rslts->n_close_data > ndmax) 
				srch_rslts->n_close_data = ndmax;
		}

		return srch_rslts;

	};




};



/*
// apply A* algorithm to obtain a map of paths.

basic pseudocode is

search nearest datapoint
	while n_found < n_to_find
		move a tick in every direction 
		if datapoint found
			++n_found

trouble is getting a datastructure for this

each movement should cost something. the squared difference to the previous in the path, or something else??? 
maybe if we are modelling stuff like flows, the inverse of the permeasbility

needs a grid to store a passed non passed flag, and one to get the previous linked path 



*/





#endif