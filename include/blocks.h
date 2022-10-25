#ifndef BLOCKS_INCLUDED
#define BLOCKS_INCLUDED

#include "search_pars.h"
#include "log.h"

class BlocksPars{

public:

	//struct block_coord
	//{

	//	float *x;
	//	float *y;
	//	float *z;

	//} b_coord_type;



	unsigned int maxblocks;

	int		blocksflag;

	std::string	blocksfile;

	std::string fileid, blockid;

	unsigned int n_blocks;

	float *block_value;
	float *block_error;
	unsigned int *block_n_points;

	//block_coord 
	cl_coord<double> **blockcoords;

	double *a_b2b;
	double *a_b2p;
	double *r_b2b;
	double *rr_b2b;

	float **block_cov;
	float **block2point_covt;
	
	unsigned int ncloseblocks;
	

	//functions

	int read_blocks_parameters(registry *reg);

	void allocate_blocks_matrixes(int na);

	void deallocate_blocks_matrixes();

	//void allocate_full_matrixes(int na,int ncloseblocks, double **full_a, double **full_r, double **full_rr, double **full_s ){
	//
	//	// and allocate full matrixes
	//	int temp =(na+ncloseblocks)*(na+ncloseblocks)/2.0+1+(na+ncloseblocks)/2.0;
	//	*full_a= new double[temp];
	//	*full_r= new double[(na+1)+ncloseblocks];
	//	*full_rr=new double[(na+1)+ncloseblocks];
	//	*full_s= new double[(na+1)+ncloseblocks];
	//};

	//void deallocate_full_matrixes(double **full_a, double **full_r, double **full_rr, double **full_s ){
	//
	//	// and allocate full matrixes
	//	delete[] *full_a;
	//	delete[] *full_r;
	//	delete[] *full_rr;
	//	delete[] *full_s;
	//
	//};


	void read_block_file(std::ifstream& fbdata);

	void get_blocks_covtable(std::shared_ptr<VariogramPars> &vario);

	void get_block_to_point_covtable(std::shared_ptr<VariogramPars> &vario, GSLibGridPars& grid_def );

	void combine_matrixes(KrigeSystem& k_system, int ino, int na, KrigeSystem& full_system );


	//void put_error_term(KrigeSystem& k_system,int *nclose, long inclose[], int ncloseblocks, long *close_blocks){
	//
	//	// we now have a combined matrix a[], and need to put the error term.
	//
	//};


	int Set_matrixes_b2b( long *incloseblocks, long indpoint);


	int Set_matrixes_b2p( unsigned int na, search_results& srch_res	, long *incloseblocks, Nodes_search& nodes);

	int set_block_krig_matrix(search_results& srch_res/*SearchPars& search*/,Nodes_search& nodes,
		long * incloseblocks,long currpoint,unsigned int na,
		int ino, KrigeSystem& k_syst, KrigeSystem& full_ksyst);

private:




};

#endif