#include "krige.h"
// LOCAL CDF

void KrigePars::calculate_sk_kmean_kvar(KrigeSystem& k_syst, krig_results& krgres, boost::scoped_array<float> &rhs_avg,
	std::vector<float> &sample_values, float& local_average)
{

	for (int i = 0; i < sample_values.size(); ++i)
	{
		krgres.cmean += (float)(k_syst.s[i + 1]) * (sample_values[i] - rhs_avg[i]);
		krgres.cstdev -= (float)(k_syst.s[i + 1] * k_syst.rr[i + 1]);
		krgres.sumwts += (float)(k_syst.s[i + 1]);
	}

	krgres.cmean += local_average;
	
	return;
};

void KrigePars::Get_samples_values(	BlocksPars& blocks, search_results& srch_res,std::vector<float>& close_node_value, 
									unsigned int* na, long *incloseblocks, std::vector<float> &sample_values)
{

	for (unsigned int j = 0; j < srch_res.n_close_data; ++j)
	{
		sample_values.push_back(srch_res.found_sample_values[j].value);
	}

	for (unsigned int j = srch_res.n_close_data; j < *na; ++j)
	{
		int ind = j - srch_res.n_close_data;
		sample_values.push_back(close_node_value.at(ind));
	}

	if (blocks.blocksflag == 1)
	{
		for (unsigned int k = *na; k < *na + blocks.ncloseblocks; ++k)
		{
			int ind = incloseblocks[k - *na + 1]; // inclose starts on 1
			sample_values.push_back(blocks.block_value[ind - 1]);
		}
	}

	return;
};

//void KrigePars::Simple_K_Local_cdf(/*BlocksPars& blocks,*/ int& idbg, GSLibGridPars& grid_def,
//	/*search_results& srch_res, float& cbb, std::vector<long>& close_node_index,*/
//	KrigeSystem& k_syst, krig_results& krgres, /*long *incloseblocks, */long& currpoint,
//	/*std::vector<float>& close_node_value, unsigned int* na,*/ float& local_average, boost::scoped_array<float> &rhs_avg
//	,std::vector<float> &sample_values) // float* average
//
	void KrigePars::Simple_K_Local_cdf(	int& idbg, GSLibGridPars& grid_def,
										KrigeSystem& k_syst, krig_results& krgres, long& currpoint,
										float& local_average, boost::scoped_array<float> &rhs_avg, 
										std::vector<float> &sample_values) // float* average
{
	//float vra;


	// careful with the na

	//if ( blocks.blocksflag == 1) 
	//{
	//	vra  =new float[(search.nodmax+search.ndmax+1+1)+blocks.ncloseblocks];
	//}
	//else
	//{
	//	vra  =new float[(search.nodmax+search.ndmax+1+1)];
	//}

	// Identify row, col, level and coordinates
	//long index = currpoint;
	cl_coord<int> ijk = grid_def.get_ijk_from_index_f(currpoint);


	//float gmean = average;

	if (ising == 0)
	{

		//std::vector<float> sample_values;
		//
		//Get_samples_values(blocks, srch_res, close_node_value,	na, incloseblocks, sample_values);

		//now calculate kriging mean and variance
		calculate_sk_kmean_kvar(k_syst, krgres, rhs_avg, sample_values, local_average);

		// Error message if negative variance:
		test_negative_variance(&krgres.cstdev, ijk);

	}
	else {

		// to do check if possible to ising!=0

		if (ising > 0 && idbg > 1)
			logger << " WARNING SGSIM: singular matrix for node: " << utilities::make_string(ijk.x)
			<< " " << utilities::make_string(ijk.y) << " " << utilities::make_string(ijk.z) << " " << utilities::make_string(currpoint) << "\n";

		krgres.cmean = -999.25;//gmean;//krgres.cmean = average; //gmean;
		krgres.cstdev = 1.0;


	}

	//delete[] vra;

	return;
};





void KrigePars::calculate_ok_kmean_kvar(KrigeSystem& k_syst, krig_results& krgres, std::vector<float> &sample_values) 
{

	for (int i = 0; i < sample_values.size(); ++i)
	{
		krgres.cmean += (float)(k_syst.s[i + 1]) * (sample_values[i]);
		krgres.cstdev -= (float)(k_syst.s[i + 1] * k_syst.rr[i + 1]);
		krgres.sumwts += (float)(k_syst.s[i + 1]);
	}
	//ok constraint is on the last line
	krgres.cstdev -= (double)k_syst.s[sample_values.size() + 1];

	return;

};

void KrigePars::Ordinary_K_Local_cdf(BlocksPars& blocks, int idbg, GSLibGridPars& grid_def, HarddataPars& harddata, search_results& srch_res, float cbb, std::vector<long> close_node_index,
	KrigeSystem& k_syst, krig_results& krgres, long *incloseblocks, long currpoint, std::vector<float>& close_node_value, /*float *sim,*/ unsigned int* na)
{
	long index;
	float gmean;

	float vra;

	// careful with the na

	//if ( blocks.blocksflag == 1) {
	//	vra  =new float[(search->nodmax+search->ndmax+1)+blocks.ncloseblocks];
	//}else{
	//	vra  =new float[(search->nodmax+search->ndmax+1)];
	//}
	// Identify row, col, level and coordinates
	index = currpoint;
	cl_coord<int> ijk;
	ijk = grid_def.get_ijk_from_index_f(index);

	// Medias locais

	gmean = (double)harddata.vmedexp;


	if (ising == 0)
	{
		krgres.cmean = 0.0;
		krgres.cstdev = cbb;
		krgres.sumwts = 0.0;

		std::vector<float> sample_values;
		//sample_values.push_back(srch_res.found_sample_values[j].value);
		for (unsigned int j = 0; j < srch_res.n_close_data; ++j)
		{
			index = srch_res.inclose[j];
			sample_values.push_back(harddata.point[index].value);

		}

		for (unsigned int j = srch_res.n_close_data; j < *na; ++j)
		{
			index = j - srch_res.n_close_data;
			//vra = sim[close_node_index.at(index)];//-1 fortran to C
			sample_values.push_back(close_node_value.at(index));

		}

		if (blocks.blocksflag == 1)
		{
			for (unsigned int k = *na + 1; k < *na + 1 + blocks.ncloseblocks; ++k)
			{
				index = incloseblocks[k - *na]; // inclose starts on 1
				sample_values.push_back(vra = blocks.block_value[index]);
			}
			//krgres.cstdev -= (double)k_syst.s[*na + 1 + blocks.ncloseblocks];
		}

		//else
		//	krgres.cstdev -= (double)k_syst.s[*na + 1];

		/*
	for (unsigned int j=1; j<=*na; ++j)
	{
		if (j <= srch_res.n_close_data)
		{
			index   = srch_res.inclose[j-1];
			vra  = harddata.point[index].value;
		}
		else
		{
			index   = j-srch_res.n_close_data;
			vra  = sim[close_node_index.at(index-1)];//-1 fortran to C
		}
		krgres.cmean  += (double)(k_syst.s[j]) * (vra);
		krgres.cstdev -= (double)(k_syst.s[j] * k_syst.rr[j]);
		krgres.sumwts	+= (double)(k_syst.s[j]);
		//if (k_syst.s[j]>1 || k_syst.s[j]<0)
		//	krgres->weight = k_syst.s[j];

	}



	// now for blocks, if we have them
	if (blocks.blocksflag==1)
	{
		for (unsigned int j=*na+1; j<=*na+blocks.ncloseblocks; ++j)
		{

			index   = incloseblocks[j-*na];// inclose starts on 1
			vra  = blocks.block_value[index];

			krgres.cmean  += (double)(k_syst.s[j]) * (vra);
			krgres.cstdev -= (double)(k_syst.s[j] * k_syst.rr[j]);
			krgres.sumwts	+= (double)(k_syst.s[j]);
		}

		krgres.cstdev -= (double)k_syst.s[*na+1+blocks.ncloseblocks]; // changed
	}
	else
	{
		krgres.cstdev -= (double)k_syst.s[*na+1];
	}
	*/

		calculate_ok_kmean_kvar(k_syst, krgres, sample_values);
	// Error message if negative variance:
		test_negative_variance(&krgres.cstdev, ijk);

	}
	else
	{

		// to do check if possible to ising!=0

		if (ising > 0 && idbg > 1)
			logger << " WARNING SGSIM: singular matrix for node: " << " " << utilities::make_string(ijk.x)
			<< " " << utilities::make_string(ijk.y) << " " << utilities::make_string(ijk.z) << " " << utilities::make_string(index) << "\n";
		krgres.cmean = gmean;
		krgres.cstdev = 1.0;
	}

	//delete [] vra;

	return;


};











void KrigePars::LVM_K_Local_cdf(BlocksPars& blocks, int idbg, GSLibGridPars& grid_def, HarddataPars& harddata, search_results& srch_res, float cbb, std::vector<long> close_node_index,
	KrigeSystem& k_syst, krig_results& krgres, long *incloseblocks, long currpoint, boost::shared_array<float> &sim, unsigned int* na)
{
	long index;

	float gmean;

	float vra, vrea;


	// careful with the na

	// Identify row, col, level and coordinates
	index = currpoint;

	cl_coord<int> ijk;
	ijk = grid_def.get_ijk_from_index_f(index);


	// Medias locais
	gmean = LVM->grid[index];

	if (ising == 0)
	{
		krgres.cmean = 0.0;
		krgres.cstdev = cbb;
		krgres.sumwts = 0.0;
		for (unsigned int j = 1; j <= *na; ++j)
		{
			if (j <= srch_res.n_close_data)
			{
				index = srch_res.inclose[j - 1];

				vra = harddata.point[index].value;
				vrea = harddata.vsec[index];
			}
			else
			{
				index = j - srch_res.n_close_data;
				vra = sim[close_node_index.at(index - 1)];//-1 fortran to C
				vrea = LVM->grid[close_node_index.at(index - 1)];//-1 fortran to C
			}
			krgres.cmean += (float)(k_syst.s[j]) * (vra - vrea);
			krgres.cstdev -= (float)(k_syst.s[j] * k_syst.rr[j]);
			krgres.sumwts += (float)(k_syst.s[j]);
		}
		// now for blocks, if we have them
		if (blocks.blocksflag == 1)
		{
			for (unsigned int j = *na + 1; j < *na + blocks.ncloseblocks + 1; ++j)
			{

				index = incloseblocks[j - *na];// inclose starts on 1
				vra = blocks.block_value[index];
				vrea = LVM->grid[index];

				krgres.cmean += (float)(k_syst.s[j]) * (vra - vrea);
				krgres.cstdev -= (float)(k_syst.s[j] * k_syst.rr[j]);
				krgres.sumwts += (float)(k_syst.s[j]);
			}

		}
		krgres.cmean += gmean;


		// Error message if negative variance:
		test_negative_variance(&(krgres).cstdev, ijk);

	}
	else {

		// to do check if possible to ising!=0

		if (ising > 0 && idbg > 1)
			logger << " WARNING SGSIM: singular matrix for node: " << utilities::make_string(ijk.x)
			<< " " << utilities::make_string(ijk.y) << " " << utilities::make_string(ijk.z) << " " << utilities::make_string(index) << "\n";
		krgres.cmean = gmean;
		krgres.cstdev = 1.0;
	}

	//delete[]vrea;
	//delete[]vra;

	return;
};

void KrigePars::ExtrD_K_Local_cdf(BlocksPars& blocks, int idbg, GSLibGridPars& grid_def, HarddataPars& harddata, search_results& srch_res, float cbb, std::vector<long> close_node_index,
	KrigeSystem& k_syst, krig_results& krgres, long *incloseblocks, long currpoint, boost::shared_array<float> &sim, unsigned int* na)
{
	long index;
	float gmean;

	float *vra;

	// careful with the na


	// malloc vra, vrea
	if (blocks.blocksflag == 1)
	{
		vra = new float[(srch_res.nodmax + srch_res.ndmax + 1 + 1) + blocks.ncloseblocks];
	}
	else
	{
		vra = new float[(srch_res.nodmax + srch_res.ndmax + 1 + 1)];

	}
	// Identify row, col, level and coordinates
	index = currpoint; //order[in];


	cl_coord<int> ijk = grid_def.get_ijk_from_index_f(index);


	// Medias locais
	gmean = (float)harddata.vmedexp;

	if (ising == 0)
	{
		krgres.cmean = 0.0;
		krgres.cstdev = cbb;
		krgres.sumwts = 0.0;

		for (unsigned int j = 1; j <= *na; ++j)
		{
			if (j <= srch_res.n_close_data)
			{
				index = srch_res.inclose[j - 1];
				//vra[j]  = harddata.vr[index];
				vra[j] = harddata.point[index].value;
			}
			else
			{
				index = j - srch_res.n_close_data;
				vra[j] = sim[close_node_index.at(index - 1)];//-1 fortran to C
			}

			krgres.cmean += (float)(k_syst.s[j]) * (vra[j] - harddata.vmedexp);

			krgres.cstdev -= (float)(k_syst.s[j] * k_syst.rr[j]);
			krgres.sumwts += (float)(k_syst.s[j]);
		}
		// now for blocks, if we have them
		if (blocks.blocksflag == 1)
		{
			for (unsigned int j = *na + 1; j < *na + 1 + blocks.ncloseblocks; ++j)
			{

				index = incloseblocks[j - *na];// inclose starts on 1
				vra[j] = blocks.block_value[index];

				krgres.cmean += (float)(k_syst.s[j]) * (vra[j] - harddata.vmedexp);

				krgres.cstdev -= (float)(k_syst.s[j] * k_syst.rr[j]);
				krgres.sumwts += (float)(k_syst.s[j]);
			}

		}


		// Error message if negative variance:
		test_negative_variance(&krgres.cstdev, ijk);

	}
	else {

		// to do check if possible to ising!=0

		if (ising > 0 && idbg > 1)
			logger << " WARNING SGSIM: singular matrix for node: " << utilities::make_string(ijk.x)
			<< " " << utilities::make_string(ijk.y) << " " << utilities::make_string(ijk.z) << " " << utilities::make_string(index) << "\n";
		krgres.cmean = gmean;
		krgres.cstdev = 1.0;
	}

	delete[] vra;

	return;
};



void KrigePars::calculate_s_co_k_kmean_kvar(KrigeSystem& k_syst, krig_results& krgres, float &secondary_val,
	std::vector<float> &sample_values, float& local_average)
{

	for (int i = 0; i < sample_values.size(); ++i)
	{
		krgres.cmean += (float)(k_syst.s[i + 1]) * (sample_values[i] - local_average);
		krgres.cstdev -= (float)(k_syst.s[i + 1] * k_syst.rr[i + 1]);
		krgres.sumwts += (float)(k_syst.s[i + 1]);
	}

	krgres.cmean += (float)k_syst.s[sample_values.size() + 1] * (secondary_val - local_average);
	krgres.cstdev -= (float)(k_syst.s[sample_values.size() + 1] * k_syst.rr[sample_values.size() + 1]);

	krgres.cmean += local_average;

	return;
};



void KrigePars::Colloc_CoK_Local_cdf(BlocksPars& blocks, int idbg, GSLibGridPars& grid_def,
	HarddataPars& harddata, search_results& srch_res, float cbb, std::vector<long> close_node_index,
	KrigeSystem& k_syst, krig_results& krgres, long *incloseblocks, long currpoint,
	std::vector<float>& close_node_value, boost::shared_array<float> &sim, unsigned int* na)
{
	float vra;

	// Identify row, col, level and coordinates
	long index = currpoint;

	cl_coord<int> ijk;
	ijk = grid_def.get_ijk_from_index_f(index);


	// Medias locais

	//gmean = harddata.vmedexp;

	if (ising == 0)
	{
		krgres.cmean = 0.0;
		krgres.cstdev = cbb;
		krgres.sumwts = 0.0;

		std::vector<float> sample_values;
		for (unsigned int j = 0; j < srch_res.n_close_data; ++j)
		{
			index = srch_res.inclose[j];
			//vra = 
			sample_values.push_back(harddata.point[index].value);

		}

		for (unsigned int j = srch_res.n_close_data; j < *na; ++j)
		{
			index = j - srch_res.n_close_data;
			//vra = sim[close_node_index.at(index)];//-1 fortran to C
			//vra = 
			sample_values.push_back(close_node_value.at(index));
		}

		// now for blocks, if we have them
		if (blocks.blocksflag == 1)
		{
			for (unsigned int j = *na + 1; j < *na + blocks.ncloseblocks + 1; ++j)
			{
				index = incloseblocks[j - *na];// inclose starts on 1
				//vra = 
				sample_values.push_back(blocks.block_value[index]);
			}

			index = grid_def.get_index_from_ijk_f(ijk);


		}

		////sample_values.push_back(srch_res.found_sample_values[j].value);
		//for (unsigned int j = 0; j < srch_res.n_close_data; ++j)
		//{
		//	index = srch_res.inclose[j];
		//	sample_values.push_back(harddata.point[index].value);

		float gmean = harddata.vmedexp;

		calculate_s_co_k_kmean_kvar(k_syst, krgres, secondary_grid->grid[currpoint],
			sample_values, gmean);

		//calculate_ok_kmean_kvar(k_syst, krgres, sample_values);

		// Error message if negative variance:
		test_negative_variance(&krgres.cstdev, ijk);

	}
	else
	{

		// to do check if possible to ising!=0

		if (ising > 0 && idbg > 1)
			logger << " WARNING SGSIM: singular matrix for node: " << utilities::make_string(ijk.x)
			<< " " << utilities::make_string(ijk.y) << " " << utilities::make_string(ijk.z) << " " << utilities::make_string(index) << "\n";
		krgres.cmean = harddata.vmedexp;
		krgres.cstdev = 1.0;
	}


	//delete [] vra;

	return;
};





void KrigePars::LVM_CoK_Local_cdf(BlocksPars& blocks, int idbg, GSLibGridPars& grid_def, HarddataPars& harddata, search_results& srch_res, float cbb, std::vector<long> close_node_index,
	KrigeSystem& k_syst, krig_results& krgres, long *incloseblocks, long currpoint, boost::shared_array<float> &sim, unsigned  int* na)
{
	long index;
	unsigned int  j;
	float gmean;

	float *vra, *vrea;


	// careful with the na

	if (blocks.blocksflag == 1)
	{
		vra = new float[(srch_res.nodmax + srch_res.ndmax + 1 + 1) + blocks.ncloseblocks];
		vrea = new float[(srch_res.nodmax + srch_res.ndmax + 1 + 1) + blocks.ncloseblocks];
	}
	else
	{
		vra = new float[(srch_res.nodmax + srch_res.ndmax + 1 + 1)];
		vrea = new float[(srch_res.nodmax + srch_res.ndmax + 1 + 1)];
	}

	// Identify row, col, level and coordinates
	index = currpoint; //order[in];

	//grid_def.get_ijk_from_index_f(index,&ix,&iy,&iz);
	cl_coord<int> ijk;
	ijk = grid_def.get_ijk_from_index_f(index);


	// Medias locais
	gmean = LVM->grid[index];

	if (ising == 0)
	{
		krgres.cmean = 0.0;
		krgres.cstdev = cbb;
		krgres.sumwts = 0.0;
		for (j = 1; j <= *na; ++j)
		{
			if (j <= srch_res.n_close_data)
			{
				index = srch_res.inclose[j - 1];
				//vra[j]  = harddata.vr[index];
				vra[j] = harddata.point[index].value;
				vrea[j] = harddata.vsec[index];
			}
			else
			{
				index = j - srch_res.n_close_data;
				vra[j] = sim[close_node_index.at(index - 1)];//-1 fortran to C
				vrea[j] = LVM->grid[close_node_index.at(index - 1)];//-1 fortran to C
			}
			krgres.cmean += (float)(k_syst.s[j]) * (vra[j] - vrea[j]);
			krgres.cstdev -= (float)(k_syst.s[j] * k_syst.rr[j]);
			krgres.sumwts += (float)(k_syst.s[j]);
		}
		// now for blocks, if we have them
		if (blocks.blocksflag == 1)
		{
			for (j = *na + 1; j < *na + blocks.ncloseblocks + 1; ++j)
			{

				index = incloseblocks[j - *na];// inclose starts on 1
				vra[j] = blocks.block_value[index];

				krgres.cmean += (float)(k_syst.s[j]) * (vra[j] - vrea[j]);
				krgres.cstdev -= (float)(k_syst.s[j] * k_syst.rr[j]);
				krgres.sumwts += (float)(k_syst.s[j]);
			}

			index = grid_def.get_index_from_ijk_f(ijk);

			krgres.cmean += (float)k_syst.s[*na + 1 + blocks.ncloseblocks] * (secondary_grid->grid[index] - LVM->grid[index]);//-harddata.vmedexp);
			krgres.cstdev -= (float)(k_syst.s[*na + 1 + blocks.ncloseblocks] * k_syst.rr[*na + 1 + blocks.ncloseblocks]);

		}
		else
		{
			index = grid_def.get_index_from_ijk_f(ijk);
			krgres.cmean += (float)k_syst.s[*na + 1] * (secondary_grid->grid[index] - LVM->grid[index]);//-harddata.vmedexp);
			krgres.cstdev -= (float)(k_syst.s[*na + 1] * k_syst.rr[*na + 1]);



		}
		krgres.cmean += gmean;


		// Error message if negative variance:
		test_negative_variance(&(krgres).cstdev, ijk);

	}
	else {

		// to do check if possible to ising!=0

		if (ising > 0 && idbg > 1)
			logger << " WARNING SGSIM: singular matrix for node: " << utilities::make_string(ijk.x)
			<< " " << utilities::make_string(ijk.y) << " " << utilities::make_string(ijk.z) << " " << utilities::make_string(index) << "\n";
		krgres.cmean = -999.25;//gmean;
		krgres.cstdev = 1.0;
	}

	delete[]vrea;
	delete[]vra;

	return;
};




void KrigePars::calculate_o_co_k_kmean_kvar(KrigeSystem& k_syst, krig_results& krgres, std::vector<float> &sample_values, float &secondary_value)
{

	for (int i = 0; i < sample_values.size(); ++i)
	{
		krgres.cmean += (float)(k_syst.s[i + 1]) * (sample_values[i]);
		krgres.cstdev -= (float)(k_syst.s[i + 1] * k_syst.rr[i + 1]);
		krgres.sumwts += (float)(k_syst.s[i + 1]);
	}
	
	//ok constraint is on the last line
	krgres.cstdev -= (double)k_syst.s[sample_values.size() + 1 + 1];// cause OK 
	krgres.sumwts += (float)(k_syst.s[sample_values.size() + 1 + 1]);// is this true? CHECK TO DO
	
	// this is co-kriging
	krgres.cstdev -= (float)(k_syst.s[sample_values.size() + 1] * k_syst.rr[sample_values.size() + 1]);
	
	// secondary and primary average are the same if we rescaled
	krgres.cmean += (float)k_syst.s[sample_values.size() + 1] * (secondary_value/*- secondary_avg + harddata.vmedexp*/);
	
	krgres.sumwts += (float)(k_syst.s[sample_values.size() + 1]); // CHECK IF THE SUM IS CORRECT

	return;
	
};

void KrigePars::O_Colloc_CoK_Local_cdf(BlocksPars& blocks, int idbg, GSLibGridPars& grid_def, HarddataPars& harddata, search_results& srch_res, float cbb, std::vector<long> close_node_index,
	KrigeSystem& k_syst, krig_results& krgres, long *incloseblocks, long currpoint, boost::shared_array<float> &sim, unsigned int* na, float secondary_avg)
{

	//float gmean;

	float vra;


	// careful with the na


	// Identify row, col, level and coordinates
	long index = currpoint;

	cl_coord<int> ijk;
	ijk = grid_def.get_ijk_from_index_f(index);


	// Medias locais

	//gmean = harddata.vmedexp;

	if (ising == 0)
	{
		krgres.cmean = 0.0;
		krgres.cstdev = cbb;
		krgres.sumwts = 0.0;
		for (unsigned int j = 1; j <= *na; ++j)
		{

			if (j <= srch_res.n_close_data)
			{
				index = srch_res.inclose[j - 1];
				vra = harddata.point[index].value;
			}
			else
			{
				index = j - srch_res.n_close_data;
				vra = sim[close_node_index[index - 1]];//-1 fortran to C
			}

			//krgres.cmean += (float)(k_syst.s[j]) * (vra - harddata.vmedexp);
			//krgres.cstdev -= (float)(k_syst.s[j] * k_syst.rr[j]);
			//krgres.sumwts += (float)(k_syst.s[j]);
			krgres.cmean += (float)(k_syst.s[j]) * (vra);
			krgres.cstdev -= (float)(k_syst.s[j] * k_syst.rr[j]);
			krgres.sumwts += (float)(k_syst.s[j]);

		}


		// now for blocks, if we have them
		if (blocks.blocksflag == 1)
		{
			for (unsigned int j = *na + 1; j < *na + blocks.ncloseblocks + 1; ++j)
			{
				index = incloseblocks[j - *na];// inclose starts on 1
				vra = blocks.block_value[index];

				//krgres.cmean += (float)(k_syst.s[j]) * (vra - harddata.vmedexp);

				//krgres.cstdev -= (float)(k_syst.s[j] * k_syst.rr[j]);
				//krgres.sumwts += (float)(k_syst.s[j]);

				krgres.cmean += (float)(k_syst.s[j]) * (vra);
				krgres.cstdev -= (float)(k_syst.s[j] * k_syst.rr[j]);
				krgres.sumwts += (float)(k_syst.s[j]);
			}

			/*OK*/	krgres.cstdev -= (float)k_syst.s[*na + 1 + 1 + blocks.ncloseblocks]; // changed

			index = grid_def.get_index_from_ijk_f(ijk);

			krgres.cmean += (float)k_syst.s[*na + 1 + blocks.ncloseblocks] * (secondary_grid->grid[index] - harddata.vmedexp);
			krgres.cstdev -= (float)(k_syst.s[*na + 1 + blocks.ncloseblocks] * k_syst.rr[*na + 1 + blocks.ncloseblocks]);

		}
		else
		{

			index = currpoint;//grid_def.get_index_from_ijk_f(&ijk);
							  // secondary and primary are the same if we rescaled
			krgres.cmean += (float)k_syst.s[*na + 1] * (secondary_grid->grid[index]/*- secondary_avg + harddata.vmedexp*/);
			krgres.cstdev -= (float)(k_syst.s[*na + 1] * k_syst.rr[*na + 1]);
			krgres.sumwts += (float)(k_syst.s[*na + 1]);
			/*OK*/	krgres.cstdev -= (float)k_syst.s[*na + 2];
			krgres.sumwts += (float)(k_syst.s[*na + 2]);
		}
		//krgres.cmean += harddata.vmedexp;

		// Error message if negative variance:
		test_negative_variance(&krgres.cstdev, ijk);

	}
	else
	{

		// to do check if possible to ising!=0

		if (ising > 0 && idbg > 1)
			logger << " WARNING SGSIM: singular matrix for node: " << utilities::make_string(ijk.x)
			<< " " << utilities::make_string(ijk.y) << " " << utilities::make_string(ijk.z) << " " << utilities::make_string(index) << "\n";
		krgres.cmean = harddata.vmedexp;
		krgres.cstdev = 1.0;
	}


	//delete [] vra;

	return;
};




inline void KrigePars::test_negative_variance(double *cstdev, cl_coord<int> &currpoint)
{
	// Error message if negative variance:
	if (*cstdev < 0.0) {

		std::ostringstream ix_str, iy_str, iz_str, cstdev_str;
		ix_str << currpoint.x;
		iy_str << currpoint.y;
		iz_str << currpoint.z;
		cstdev_str << *cstdev;


		logger << " ERROR: Negative Variance: " << ix_str.str()
			<< " " << iy_str.str() << " " << iz_str.str()
			<< " " << cstdev_str.str() << "\n";
		*cstdev = (float)0.0;
	}
	//o problema da variancia é aqui: agora so falta convencer o mico que assim è q esta bwm; com a imagem suave
	*cstdev = (float)sqrt(std::max(*cstdev, 0.0));
	return;

};
