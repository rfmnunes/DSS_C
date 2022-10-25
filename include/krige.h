#ifndef KRIGE_INCLUDED
#define KRIGE_INCLUDED

#include "blocks.h"
#include <sstream>

#include "logging.h"
#include "log.h"
#include "zones.h"

class krig_results 
{
public:
	double cmean;
	double cstdev;
	double sumwts;


	krig_results() 
	{
		cmean = -999.25;
		cstdev = -999.25;
		sumwts=-999.25;
	}

	krig_results(float cbb)
	{
		cmean = 0.0;
		cstdev = cbb;
		sumwts = 0.0;
	}

	krig_results(const krig_results& inresult)
	{
		cmean = inresult.cmean;
		cstdev = inresult.cstdev;
		sumwts = inresult.sumwts;
	}

	krig_results& operator=(const krig_results& inresult)
	{
		cmean = inresult.cmean;
		cstdev = inresult.cstdev;
		sumwts = inresult.sumwts;
	}

	~krig_results()
	{

	}

};


class KrigePars{


	//   Size of the kriging system (na - number of equations)
	//
	//   Compute the estimate and kriging variance.  Recall that kriging type
	//     0 = Simple Kriging
	//     1 = Ordinary Kriging
	//     2 = Locally Varying Mean
	//     3 = External Drift
	//     4 = Collocated Cosimulation
	//     5 = Collocated Cosimulation with local corr. coef.

public:

	int		ktype;
	
	float colocorr;


	KrigeSystem *k_syst;
	KrigeSystem *full_ksyst; // for blocks


	std::string lvmfl, corrfl, secfl;

	GSLibGrid<float> *secondary_grid;
	GSLibGrid<float> *cc_grid;
	GSLibGrid<float> *LVM;

	int ising;

	int rescale;

	int icollvm,nvaril;

	KrigePars(registry* reg, GSLibGridPars& grid_def,
		misc& utils, Zones_Pars& Zones);


	int Read_krige_pars(registry *reg);


	//void Simple_krg	(  BlocksPars& blocks,
	//		HarddataPars& harddata ,SearchPars &search,
	//		VariogramPars& vario, GSLibGridPars& grid, Nodes_search& nodes,
	//		float *sim, long currpoint,krig_results& krgres, long *incloseblocks, float average);

	void Simple_krg(BlocksPars& blocks,
		HarddataPars& harddata, search_results&,
		VariogramPars& vario, GSLibGridPars& grid, Nodes_search& nodes,
		/*float *sim,*/ long currpoint, krig_results& krgres, long *incloseblocks, float current_average, boost::scoped_array<float> &rhs_average);

	void Ordinary_krg(BlocksPars& blocks, 
			HarddataPars& harddata ,search_results& srch_res,
			VariogramPars& vario, GSLibGridPars& grid, Nodes_search& nodes,
			/*float *sim,*/ long currpoint, krig_results& krgres, long *incloseblocks);

	void LVM_krg( BlocksPars& blocks 
			,HarddataPars& harddata ,search_results& srch_res,
			VariogramPars& vario, GSLibGridPars& grid, Nodes_search& nodes,
			boost::shared_array<float> &sim, long currpoint , krig_results& krgres, long *incloseblocks);

	void Ext_drift_krg( BlocksPars& blocks,
			HarddataPars& harddata ,search_results& srch_res,
			VariogramPars& vario, GSLibGridPars& grid, Nodes_search& nodes,
			boost::shared_array<float> &sim, long currpoint, boost::shared_array<float> &secondary , krig_results& krgres, long *incloseblocks);

	void ColocCosim_krg	(BlocksPars& blocks
			,HarddataPars& harddata ,search_results& srch_res,
			VariogramPars& vario, GSLibGridPars& grid, Nodes_search& nodes,
			boost::shared_array<float> &sim, long currpoint, krig_results& krgres, long *incloseblocks);

	void LocalCC_ColocCosim_krg	(BlocksPars& blocks,
			HarddataPars& harddata ,search_results& srch_res,
			VariogramPars& vario, GSLibGridPars& grid, Nodes_search& nodes,
			boost::shared_array<float> &sim, long currpoint, float localCC  , krig_results& krgres, long *incloseblocks);

	void LVM_ColocCosim_krg	(BlocksPars& blocks
		,HarddataPars& harddata ,search_results& srch_res,
		VariogramPars& vario, GSLibGridPars& grid, Nodes_search& nodes,
		boost::shared_array<float> &sim, long currpoint, float localCC  , krig_results& krgres, long *incloseblocks);

	void O_ColocCosim_krg(BlocksPars &blocks, HarddataPars &harddata, search_results& srch_res,
		VariogramPars &vario, GSLibGridPars &grid, Nodes_search &nodes, boost::shared_array<float> &sim, long currpoint,
		krig_results &krgres, long *incloseblocks, float secondary_avg);

	void O_LocalCC_ColocCosim_krg(BlocksPars &blocks, HarddataPars &harddata, search_results& srch_res,
		VariogramPars &vario, GSLibGridPars &grid, Nodes_search &nodes, boost::shared_array<float> &sim, long currpoint,
		float localCC, krig_results & krgres, long * incloseblocks, float secondary_avg);


private:
	
	//int Set_matrixes( HarddataPars& harddata , VariogramPars& vario,
	//	GSLibGridPars& grid,Nodes_search& nodes, unsigned int na, 
	//	SearchPars& search,	long currentpoint,	KrigeSystem& k_system );

	int Set_matrixes(HarddataPars& harddata, VariogramPars& vario,
		GSLibGridPars& grid, Nodes_search& nodes, unsigned int na,
		search_results&, long& currentpoint, KrigeSystem& k_system);


	// Constraints addition

	//	void Simple_constr(){} // has none

	void Ordinary_constr(int& ino, unsigned int nlines,KrigeSystem& k_system );

	//	void LVM_constr(){}// has none here

	void Ext_drift_constr(int& ino, unsigned int nb, search_results& srch_res, std::vector<long> &close_node_index,unsigned int *neq,
		long &currpoint, float *vrea, boost::shared_array<float> &secondary, float *vsec,KrigeSystem& k_system, unsigned int* );

	void ColocCosim_constr(int& ino, unsigned int nb, float &colocorr, KrigeSystem& k_system);

	//void LocalCC_ColocCosim_constr(int& ino, unsigned int nb, float &localCC,KrigeSystem& k_system);


	//void Simple_K_Local_cdf		(BlocksPars& blocks, int& idbg, GSLibGridPars& grid_def, /*HarddataPars& harddata,*/ search_results& srch_res,/*float& cbb, std::vector<long>& close_node_index,*/
	//		KrigeSystem& k_syst, krig_results& krgres,  long *incloseblocks,long& currpoint, std::vector<float>& close_node_value, unsigned int *, float& average, boost::scoped_array<float> &rhs_average
	//		, std::vector<float> &sample_values);

	void Simple_K_Local_cdf(int& idbg, GSLibGridPars& grid_def, KrigeSystem& k_syst, 
		krig_results& krgres,long& currpoint, float& average, 
		boost::scoped_array<float> &rhs_average, std::vector<float> &sample_values);

	void Ordinary_K_Local_cdf	(BlocksPars& blocks, int idbg, GSLibGridPars& grid_def, HarddataPars& harddata, search_results& srch_res,float cbb, std::vector<long> close_node_index,
			KrigeSystem& k_syst, krig_results& krgres,  long *incloseblocks,long currpoint, std::vector<float>& close_node_value, unsigned int*);

	void LVM_K_Local_cdf		(BlocksPars& blocks, int idbg, GSLibGridPars& grid_def, HarddataPars& harddata, search_results& srch_res,float cbb, std::vector<long> close_node_index,
			KrigeSystem& k_syst, krig_results& krgres,  long *incloseblocks,long currpoint, boost::shared_array<float> &sim, unsigned int*);

	void ExtrD_K_Local_cdf		(BlocksPars& blocks, int idbg, GSLibGridPars& grid_def, HarddataPars& harddata,  search_results& srch_res,float cbb, std::vector<long> close_node_index,
			KrigeSystem& k_syst, krig_results& krgres,  long *incloseblocks,long currpoint, boost::shared_array<float> &sim, unsigned int *);

	void Colloc_CoK_Local_cdf	(BlocksPars& blocks, int idbg, GSLibGridPars& grid_def, HarddataPars& harddata,  search_results& srch_res,float cbb,  std::vector<long> close_node_index,
			KrigeSystem& k_syst, krig_results& krgres,  long *incloseblocks,long currpoint, std::vector<float>& close_node_value, boost::shared_array<float> &sim, unsigned int*);

	void test_negative_variance(double *cstdev, cl_coord<int> &currpoint);

	void LVM_CoK_Local_cdf (BlocksPars& blocks, int idbg, GSLibGridPars& grid_def, HarddataPars& harddata,  search_results& srch_res,float cbb, std::vector<long> close_node_index,
			KrigeSystem& k_syst, krig_results& krgres,  long *incloseblocks,long currpoint, boost::shared_array<float> &sim, unsigned int* na );

	void O_Colloc_CoK_Local_cdf(BlocksPars & blocks, int idbg, GSLibGridPars & grid_def, HarddataPars & harddata, search_results& srch_res, float cbb, std::vector<long> close_node_index,
		KrigeSystem & k_syst, krig_results & krgres, long * incloseblocks, long currpoint, boost::shared_array<float> &sim, unsigned int * na, float secondary_avg);


	void calculate_sk_kmean_kvar(KrigeSystem& k_syst, krig_results& krgres, boost::scoped_array<float> &rhs_avg, std::vector<float> &sample_values, float& local_average);
	void calculate_ok_kmean_kvar(KrigeSystem& k_syst, krig_results& krgres, std::vector<float> &sample_values);
	void calculate_o_co_k_kmean_kvar(KrigeSystem& k_syst, krig_results& krgres, std::vector<float> &sample_values, float &secondary_value);
	void calculate_s_co_k_kmean_kvar(KrigeSystem& k_syst, krig_results& krgres, float &secondary_val, std::vector<float> &sample_values, float& local_average);

	void Get_samples_values(BlocksPars& blocks, search_results& srch_res, std::vector<float>& close_node_value,
		unsigned int* na, long *incloseblocks, std::vector<float> &sample_values);

};




#endif