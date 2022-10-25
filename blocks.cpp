#include "blocks.h"

int BlocksPars::read_blocks_parameters(registry *reg){
		
	reg_key *k;
	// Start blocks pars

	k = get_key(reg, (char*)("BLOCKS"), (char*)("USEBLOCKS"));
	if (k)
		blocksflag = get_int(k);
	else return -1;

	if (blocksflag==1)
	{

		k = get_key(reg, (char*)("BLOCKS"), (char*)("MAXBLOCKS"));
		if (k)
			maxblocks = get_int(k);
		else return -1;


		if ((k = get_key(reg, (char*)("BLOCKS"), (char*)("BLOCKSFILE"))) != NULL)
			blocksfile= get_string(k);
		boost::algorithm::trim(blocksfile);
		logger<<" Blocks file: "<< blocksfile <<"\n";

		std::ifstream fbdata;

		fbdata.open(blocksfile.c_str());
		if ((fbdata.is_open() == 0))
		{
			logger << "blocks file: " << blocksfile << " wont open! \n";

		}


		read_block_file(fbdata);
		
		
		
	} // end blocksflag 



	return 0;
};

void BlocksPars::allocate_blocks_matrixes(int na){

	int temp=(ncloseblocks*(ncloseblocks/2.))+ncloseblocks/2.;

	a_b2b= new double[temp];

	a_b2p= new double[ncloseblocks*na];

	r_b2b= new double[ncloseblocks];

	rr_b2b= new double[ncloseblocks];

};

void BlocksPars::deallocate_blocks_matrixes(){

	delete[] a_b2b;
	delete[] a_b2p;
	delete[] r_b2b;
	delete[] rr_b2b;

};

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


void BlocksPars::read_block_file(std::ifstream& fbdata){

	fbdata >> fileid;
	fbdata >> n_blocks;

	// allocate space to keep block info

	block_value = new float[n_blocks]; // allocating values
	block_error = new float[n_blocks]; // allocating values
	block_n_points = new unsigned int[n_blocks]; // allocating values


	//cycle read_blocks

	blockcoords = new cl_coord<double>*[n_blocks];;

	for (unsigned int l = 0 ; l < n_blocks; l++)
	{

		fbdata >> blockid;
		fbdata >> block_value[l];
		fbdata >> block_error[l];
		fbdata >> block_n_points[l];

		// allocate space to put the block points

	//	blockcoords[l].x = new float[block_n_points[l]]; // allocating values

		blockcoords[l] = new cl_coord<double>[block_n_points[l]];

	//	blockcoords[l].y = new float[block_n_points[l]]; // allocating values
	//	blockcoords[l].z = new float[block_n_points[l]]; // allocating values

		for (unsigned int k = 0 ; k < block_n_points[l]; k++)
		{

			fbdata >> blockcoords[l][k].x;
			fbdata >> blockcoords[l][k].y;
			fbdata >> blockcoords[l][k].z;
		}

	}
};

void BlocksPars::get_blocks_covtable(std::shared_ptr<VariogramPars> &vario){

	double cov,covasum;

	//cria uma tab de covariancia - triangular serve, acho
	//int temp = (n_blocks*(n_blocks/2.))+n_blocks/2.;

	block_cov = new float*[n_blocks] ;

	for (unsigned int l =0; l < n_blocks; l++)
	{
		block_cov[l] = new float[n_blocks] ;
	}

	for (unsigned int l =0; l< n_blocks; l++)
	{
		for (unsigned int k = 0; k < n_blocks; k++)
		{

			covasum=0;

			for (unsigned int n=0; n < block_n_points[l];n++ )
			{ // points of first block

				//cl_coord<float> blk_1_xyz(blockcoords[l].x[n],blockcoords[l].y[n],blockcoords[l].z[n]);
				cl_coord<double> blk_1_xyz(blockcoords[l][n]);
				
				for(unsigned int m =0; m < block_n_points[k];m++ )
				{ // points of 2nd block

					//cl_coord<float> blk_2_xyz(blockcoords[k].x[m],blockcoords[k].y[m],blockcoords[k].z[m]);
					cl_coord<double> blk_2_xyz(blockcoords[k][m]);
					
					// cov entre pontos
					cov = vario->covariance_no_search_range(blk_1_xyz,blk_2_xyz);

					covasum= covasum+cov;

				}
			}

			// now average the covariance and save covariance on a table

			block_cov[l][k]=covasum/(block_n_points[l]*block_n_points[k]);

		}
	}


};

void BlocksPars::get_block_to_point_covtable(std::shared_ptr<VariogramPars> &vario, GSLibGridPars& grid_def ){


	double cov, covasum;

	long long nxyz= grid_def.get_nxyz();


	block2point_covt = new  float*[n_blocks] ;

	for (unsigned int l =0; l < n_blocks; l++)
	{
		block2point_covt[l] = new float[nxyz+1] ;
	}

	for (unsigned int l =0; l< n_blocks; l++)
	{
		for (unsigned int k = 1; k < nxyz+1; k++)
		{

			covasum=0;
			
			// now get the point

			cl_coord<double>xyz= grid_def.get_xyz_from_index_f( k);

				
			for(unsigned int m =0; m < block_n_points[l];m++ )
			{ // points of  block

				//cl_coord<float> blkxyz(blockcoords[l].x[m],blockcoords[l].y[m],blockcoords[l].z[m]);
				cl_coord<double> blkxyz(blockcoords[l][m]);
				// cov entre pontos
				cov = vario->covariance_no_search_range(xyz, blkxyz);

				covasum= covasum+cov;

			}
			
			// now average the covariance and save covariance on a table
			block2point_covt[l][k]=covasum/(block_n_points[l]);
		}
	}



};

void BlocksPars::combine_matrixes(KrigeSystem& k_system, int ino, int na, KrigeSystem& full_system ){

	//int l,i, j ;

	//	KrigeSystem *full_ksyst= new KrigeSystem(((na+1)+ncloseblocks),((na+1)+ncloseblocks));

	// first do the R and rr matrixes

	for (int l = 0; l<= na; ++l)
	{
		full_system .r[l]=k_system.r[l];
		full_system.rr[l]=k_system.rr[l];
	}
	for (unsigned int l = 0; l< ncloseblocks; ++l)
	{
		full_system.r[na+l+1]=r_b2b[l];
		full_system.rr[na+l+1]=rr_b2b[l];
	}

	// now for the A matrix
	for (int l = 0; l<= ino; ++l) // not na
	{
		full_system.a[l]=k_system.a[l];
	}

	int aux_b2b = 0;

	for (unsigned int l = 0; l< ncloseblocks; ++l)
	{
		int i;
		for (i = 0; i< na; ++i)
		{
			full_system.a[(ino+1)+(l*na+i)+aux_b2b]=a_b2p[l*na+i];
		}



		for (unsigned int j = 0; j<= l; ++j)
		{
			full_system.a[(ino+1)+(l*na+i)+aux_b2b]=a_b2b[aux_b2b];
			//full_a[(ino+1)+(l*na+i)+aux_b2b]=a_b2p[l*na+i];
			++aux_b2b;
		}

	}

	return;

};


//void put_error_term(KrigeSystem& k_system,int *nclose, long inclose[], int ncloseblocks, long *close_blocks){
//
//	// we now have a combined matrix a[], and need to put the error term.
//
//};


int BlocksPars::Set_matrixes_b2b(  long *incloseblocks, long indpoint )
{
	unsigned int j, i;//,n ; // cycle control
	int ino;

	long index1, index2;
//		float cov;

	// Set up kriging matrices:

	// get the correct one


	ino=0;
	for (j=0; j < ncloseblocks; ++j) {
		// Sort out the actual location of point "j"

		index1  = incloseblocks[j+1];


		// Sort out the actual location of point "i"
		for (i=0; i<=j; ++i) {

			index2  = incloseblocks[i+1];

			// Now, get the covariance value:

			if (index1 == index2)//TO DO CHECK THIS
				a_b2b[ino]=(1+block_error[index1]);
			else
				a_b2b[ino]=block_cov[index2-1][index1-1]; // estava -1


			////a_b2b[ino] = (double)cov- max(jerror,ierror)*cov; // we apply error here

			++ino;

		}

		// now to point

		//lets get it from covariance table


		// calculate point index 



		r_b2b[j] = block2point_covt[index1-1][indpoint];


		//cov =0 ;
		//for (n =0; n < block_n_points[index1-1];n++ ) {

		//	cov +=vario.covariance(xx,yy,zz,blockcoords[index1-1].x[n],blockcoords[index1-1].y[n],blockcoords[index1-1].z[n], utils);

		//}


		// Get the RHS value:

		//cov = vario.covariance(xx,yy,zz,x1,y1,z1, utils);
		//r_b2b[j] =(double)cov/ block_n_points[index1-1] /* - jerror*cov*/; // error applies
		rr_b2b[j]=r_b2b[j];

	}

	return ino;
};


int BlocksPars::Set_matrixes_b2p( unsigned int na, search_results& srch_res	, long *incloseblocks, Nodes_search& nodes)
{
		
	int ino;
	long index, index2;

	
	// Set up kriging matrices:

	ino=0;


	for (unsigned int j=0; j<ncloseblocks; ++j) {
		// Sort out the actual location of block "j"

		index2  = incloseblocks[j+1];

		// Sort out the actual location of point "i"
		for (unsigned int i=0; i<na; ++i) {

			if (i < srch_res.n_close_data) {

				index  = srch_res.inclose[i+1];
				//x2     = harddata.x[index];
				//y2     = harddata.y[index];
				//z2     = harddata.z[index];
			} else {
				index  = i - srch_res.n_close_data;

				//grid.get_ijk_from_index_f(nodes.close_node_index[index+1], &iix, &iiy, &iiz);

				//grid.get_xyz_from_ijk_f(iix ,iiy,iiz, &x2, &y2, &z2);

			}

			a_b2p[ino] = block2point_covt[index2-1][nodes.close_node_index->at(index)];

			//cov=0;

			//for (n =0; n < block_n_points[index2-1];n++ ) {


			//	// get the point
			//	// Now, get the covariance value:


			//	// we cant use the covariance look-up table:

			//	cov +=vario.covariance(blockcoords[index2-1].x[n],blockcoords[index2-1].y[n],blockcoords[index2-1].z[n],x2,y2,z2, utils);

			//}

			//a_b2p[ino] = (double)cov/ block_n_points[index2-1] /*- jerror*cov*/; // error term

			++ino;
		}

	}

	return ino;
};

int  BlocksPars::set_block_krig_matrix(  search_results& srch_res,
	Nodes_search& nodes, long *incloseblocks,long currpoint,
	unsigned int na,int ino, KrigeSystem& k_syst, KrigeSystem& full_ksyst  )
{
	
		if (ncloseblocks>0)	
			allocate_blocks_matrixes(na);

		int ino_b2b=Set_matrixes_b2b(  incloseblocks, currpoint);
			
		int ino_b2p=Set_matrixes_b2p( na, srch_res, incloseblocks, nodes);
			
			
		// now combine all matrixes

		combine_matrixes(k_syst, ino, na, full_ksyst );

		if (ncloseblocks>0)	
			deallocate_blocks_matrixes();

		return ino_b2b+ino_b2p;
}
