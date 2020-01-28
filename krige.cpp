#include "krige.h"


int KrigePars::Read_krige_pars(registry *reg)
{
	reg_key *k;
		
// parsing kriging section

k = get_key(reg,  ("KRIGING"),  ("KTYPE"));
if (k)
	ktype = get_int(k);
else return -1;

k = get_key(reg,  ("KRIGING"),  ("COLOCORR"));
if (k)
	colocorr = get_float(k);
else return -1;
logger <<  " Kriging type (0)SK (1)OK (2)SK local means (3)KED (4)ColK-Gcc (5)ColK-Lcc (6)LVM-ColK-Gcc (7)LVM-ColK-Lcc: "<< utilities::make_string(ktype) <<"\n";
if (ktype < 0 || ktype > 7) ktype=0;
if (ktype == 4||ktype == 6) 
	logger << " Correlation coefficient: "<< utilities::make_string(colocorr)<<"\n";

if (ktype > 2){

	if ((k = get_key(reg,  ("KRIGING"),  ("SOFTFILE"))) != NULL)
		secfl= get_string(k);
	
	boost::algorithm::trim(secfl);

	logger<<" Secondary data (grid): "<< secfl <<"\n";

	k = get_key(reg,  ("KRIGING"),  ("RESCALE"));
	if (k)
		rescale = get_int(k);
	else return -1;

}

if (ktype == 2 || ktype >= 6){

	if ((k = get_key(reg,  ("KRIGING"),  ("LVMFILE"))) != NULL)
		lvmfl= get_string(k);
	
	boost::algorithm::trim(lvmfl);

	logger<<" Local Mean data (grid): "<< lvmfl <<"\n";

}


k = get_key(reg,  ("KRIGING"),  ("NVARIL"));
if (k)
	nvaril = get_int(k);
else return -1;

k = get_key(reg,  ("KRIGING"),  ("ICOLLVM"));
if (k)
	icollvm = get_int(k);
else return -1;
printf (" Number of columns/column in sec. file: %d %d \n",nvaril,icollvm);

if (ktype > 4){
	if ((k = get_key(reg,  ("KRIGING"),  ("CCFILE"))) != NULL)
		corrfl= get_string(k);
	boost::algorithm::trim(corrfl);

	logger<<" Local correlation coeficients: "<< corrfl <<"\n";

}

// end kriging section

return 0;
}


KrigePars::KrigePars(registry *reg, GSLibGridPars& grid_def,
	misc& utils, Zones_Pars& Zones)
{
	Read_krige_pars(reg);

	secondary_grid = 0;
	cc_grid = 0;
	LVM = 0;

	if (ktype == 2 || ktype >= 6) 
	{
		// prepare secondary data grid object
		LVM = new GSLibGrid<float>(grid_def);

		//	strcpy ((*secondary_grid)->filename,lvmfl); // get filename
		LVM->filename = new std::string;
		*(LVM)->filename = lvmfl;

		(LVM)->headerflag = utils.headerflag; // check if we use header

		(LVM)->null_data = utils.nosvalue; // set nulls

		logger << " Reading Local Mean datafile \n";

		switch (utils.filetype)
		{
		case GEOEAS: // case 0 is an ascii file
			(LVM)->read_grid();//reading grid, puts on object.grid
			break;
		case SGEMS: // case 1 is sgems binary
			(LVM)->read_grid_binary_sgems();
			break;
		}
		(LVM)->get_average();
		(LVM)->get_variance();

		logger << "   acceptable data  = " << (LVM)->valid_points << "\n";
		logger << "  weighted average  = " << (LVM)->average << "\n";
		logger << " weighted variance  = " << (LVM)->variance << "\n\n";


	}


	if (ktype > 2) {
		// prepare secondary data grid object
		secondary_grid = new GSLibGrid<float>(grid_def);

		//	strcpy ((*secondary_grid)->filename,lvmfl); // get filename
		secondary_grid->filename = new std::string;
		*secondary_grid->filename = secfl;

		(secondary_grid)->headerflag = utils.headerflag; // check if we use header

		(secondary_grid)->null_data = utils.nosvalue; // set nulls

		logger << " Reading secondary datafile \n";

		//(*secondary_grid)->read_grid();//reading grid, puts on object.grid

		switch (utils.filetype)
		{
		case GEOEAS: // case 0 is an ascii file
			(secondary_grid)->read_grid();//reading grid, puts on object.grid
			break;
		case SGEMS: // case 1 is sgems binary
			(secondary_grid)->read_grid_binary_sgems();
			break;
		}
		(secondary_grid)->get_average();
		(secondary_grid)->get_variance();

		logger << "   acceptable data  = " << (secondary_grid)->valid_points << "\n";
		logger << "  weighted average  = " << (secondary_grid)->average << "\n";
		logger << " weighted variance  = " << (secondary_grid)->variance << "\n\n";


	}

	// Read local correlation coeficients if necessary:
	if (ktype == 5 || ktype == 7) {
		// icollcc=1;
		// nvarilcc=1;

		cc_grid = new GSLibGrid<float>(grid_def);

		//strcpy ((*cc_grid)->filename,corrfl); // get filename
		cc_grid->filename = new std::string;
		*(cc_grid)->filename = corrfl;

		cc_grid->headerflag = utils.headerflag; // check if we use header

		cc_grid->null_data = utils.nosvalue; // set nulls

		logger << " Reading local correlation coeficients \n";

		//(*cc_grid)->read_grid();//reading grid, puts on object.grid

		switch (utils.filetype)
		{
		case GEOEAS: // case 0 is an ascii file
			(cc_grid)->read_grid();//reading grid, puts on object.grid
			break;
		case SGEMS: // case 1 is sgems binary
			(cc_grid)->read_grid_binary_sgems();
			break;
		}

		(cc_grid)->get_average();
		(cc_grid)->get_variance();

		logger << "   acceptable data  = " << (cc_grid)->valid_points << "\n";
		logger << "  weighted average  = " << (cc_grid)->average << "\n";
		logger << " weighted variance  = " << (cc_grid)->variance << "\n\n";
	}

	// Re-scale secondary variable to mean and variance of the primary variable if colocated

	//If we have secondary info and zoning, get zone dependent averages for secondary
	if (ktype >= 4) 
	{

		//Zones.n_zones = 0;
		Zones.get_zonal_avg_var(*secondary_grid);
	}



	if (ktype >= 4) {

		unsigned long   index;
		if (rescale == 1)
		{
			// do it zone to zone
			int currzone;
			for (index = 1; index <= grid_def.get_nxyz(); ++index)
			{
				currzone = Zones.get_zone_from_index(index);
							
				if ((secondary_grid->grid[index] != utils.nosvalue)&& currzone!=-1)
					secondary_grid->grid[index] =
					(float)(((secondary_grid->grid[index] - Zones.zone[currzone]->med_exp_sec)
						/
					(float)sqrt(Zones.zone[currzone]->var_exp_sec))
						*
						(float)sqrt(Zones.zone[currzone]->harddata->vvarexp)) + Zones.zone[currzone]->harddata->vmedexp;
			}
		}

	}

	// external drift

	if (ktype == 3)
	{
		for (unsigned int i = 0; i < Zones.n_zones; i++)
		{
			long   index;

			Zones.zone[i]->harddata->vsec.insert(Zones.zone[i]->harddata->vsec.end(), -999.25);

			for (unsigned int j = 1; j <= Zones.zone[i]->harddata->harddata_distr.n_data; ++j)
			{
				index = grid_def.get_index_from_xyz_f(Zones.zone[i]->harddata->point[j].xyz);
				Zones.zone[i]->harddata->vsec.insert(Zones.zone[i]->harddata->vsec.end(), secondary_grid->grid[index]);
			}
		}
	}
}


void KrigePars::Simple_krg	(	BlocksPars& blocks,
								HarddataPars& harddata ,search_results& srch_res,
								VariogramPars& vario, GSLibGridPars& grid,Nodes_search& nodes,
								long currpoint, krig_results& krgres, long *incloseblocks, float curr_loc_avg,
								boost::scoped_array<float> &rhs_avg) //average
{
	// Identify row, col, level and coordinates
	unsigned int na = srch_res.n_close_data + nodes.ncnode;
	
	int neq = na;
		
	// prepare ksol matrixes
	//int temp = na;

	KrigeSystem k_syst(na);
	
	if ( blocks.blocksflag == 1)
	{
		if (na <1)
		{
			krgres.cmean = curr_loc_avg; //average
			krgres.cstdev = 1.0;
			return;
		}


		KrigeSystem full_ksyst((na+blocks.ncloseblocks));
		
		//KrigeSystem_matrx kramt(na+blocks.ncloseblocks);
		//if(na==1) ising= full_ksyst.ksol(neq);
		//else if (na>1)
		
		// more than one, prepare system	
		int ino=Set_matrixes(harddata, vario, grid, nodes, na, srch_res, currpoint, k_syst );
		// go to ksol

		// we must also set the block to block and block to point stuff
		int block_ino=blocks.set_block_krig_matrix( srch_res, nodes, incloseblocks, currpoint, na,
													ino, k_syst, full_ksyst);
		neq += blocks.ncloseblocks;
						
		ising= full_ksyst.ksol(neq);

		std::vector<float> sample_values;


		// must be careful here when parallelizing
		Get_samples_values(blocks, srch_res, *nodes.close_node_value, &na, incloseblocks, sample_values);


		//Simple_K_Local_cdf (blocks,g_log.idbg, grid, srch_res,/* vario.cbb, *nodes.close_node_index,*/
		//					full_ksyst, krgres,  incloseblocks, currpoint, 
		//					*nodes.close_node_value , &na, curr_loc_avg, rhs_avg, sample_values); //average

		Simple_K_Local_cdf( g_log.idbg, grid, 
							full_ksyst, krgres, currpoint, 
							curr_loc_avg, rhs_avg, sample_values); //average


	}
	else
	{ // no blocks
		
		if (na <1 )
		{
			krgres.cmean= curr_loc_avg; //average
			krgres.cstdev=1.0;
			return;
		}
		
		/*else if(na==1)
		
			ising=(k_syst).ksol (neq);
		else*/

		// more than one, prepare system	
		
		int ino = Set_matrixes(harddata, vario, grid, nodes, na, srch_res, currpoint, k_syst);
		// go to ksol

		ising = k_syst.ksol (neq);

		if (ising != 0)
		{
			krgres.cmean = curr_loc_avg;//average
			krgres.cstdev = 1.0;
			return;
		}

		// simplekrig has no constraints
		std::vector<float> sample_values;

		Get_samples_values(blocks, srch_res, *nodes.close_node_value, &na, incloseblocks, sample_values);


//		Simple_K_Local_cdf(blocks, g_log.idbg, grid, srch_res,/* vario.cbb,
//			*nodes.close_node_index, */(k_syst), krgres, incloseblocks, currpoint, *nodes.close_node_value, &na, curr_loc_avg, rhs_avg, sample_values);//average
		Simple_K_Local_cdf( g_log.idbg, grid, k_syst, krgres, currpoint, 
							curr_loc_avg, rhs_avg, sample_values);//average

	}

		

};

void KrigePars::Ordinary_krg(BlocksPars& blocks, 
		HarddataPars& harddata ,search_results& srch_res,
		VariogramPars& vario, GSLibGridPars& grid, Nodes_search& nodes,
		long currpoint, krig_results& krgres, long *incloseblocks)
{
	//int *ino = new int;

	int temp;
//	int ino_b2b, ino_b2p;

	cl_coord<int> ijk	= grid.get_ijk_from_index_f(currpoint);
	cl_coord<double> xyz	= grid.get_xyz_from_ijk_f(ijk);

	unsigned int na = srch_res.n_close_data + nodes.ncnode;
	unsigned int neq = na+1;
	//if (na < 1) return;

	//alloc_matrixes(search, blocks.blocksflag);
	//temp= search.nodmax+search.ndmax/*+1+1*/;
	temp = na+1;

	KrigeSystem k_syst(temp);
		//KrigeSystem_matrx kramt(temp);
	if ( blocks.blocksflag != 1) blocks.ncloseblocks = 0;
		
	KrigeSystem full_ksyst((na+blocks.ncloseblocks));
				
	if (na <1 )
	{
		krgres.cmean=harddata.vmedexp;
		krgres.cstdev=1.0;
		return;
	}
	//else if(na==1)
	//{
	//	if ( blocks.blocksflag == 1) ising= (full_ksyst).ksol(neq+blocks.ncloseblocks);
	//	else ising=(k_syst).ksol (neq);
	//	
	//}
// more than one, prepare system	


	int ino = Set_matrixes(harddata, vario, grid,  nodes, na, srch_res,currpoint , k_syst );

	if ( blocks.blocksflag == 1) {


		int block_ino=blocks.set_block_krig_matrix( srch_res, nodes, incloseblocks, currpoint, na,
				ino, k_syst, full_ksyst);

		neq += blocks.ncloseblocks;

		int temp = ino + block_ino;
		Ordinary_constr((temp), (na+blocks.ncloseblocks),full_ksyst);
		ising=(full_ksyst).ksol (neq);

	}
	else
	{

		Ordinary_constr(ino, na,k_syst);
		ising=(k_syst).ksol (neq);

	}

	if (ising != 0)
	{
		krgres.cmean = harddata.vmedexp;
		krgres.cstdev = 1.0;
		return;
	}


	if (blocks.blocksflag==1)
	{
		Ordinary_K_Local_cdf (blocks,g_log.idbg, grid, harddata, srch_res, vario.cbb, 
				*nodes.close_node_index, (full_ksyst), krgres,  incloseblocks, currpoint, *nodes.close_node_value/*sim*/, &na);

	}
	else
		Ordinary_K_Local_cdf (blocks,g_log.idbg, grid, harddata, srch_res, vario.cbb,
				*nodes.close_node_index, (k_syst), krgres,  incloseblocks, currpoint, *nodes.close_node_value/*sim*/, &na);



}

void KrigePars::LVM_krg( BlocksPars& blocks 
		,HarddataPars& harddata ,search_results& srch_res,
		VariogramPars& vario, GSLibGridPars& grid, Nodes_search& nodes,
		boost::shared_array<float> &sim, long currpoint , krig_results& krgres, long *incloseblocks)
{
	
	int ino;

	int ino_b2b, ino_b2p;

	// Identify row, col, level and coordinates
	
	cl_coord<int> ijk	= grid.get_ijk_from_index_f(currpoint);
	cl_coord<double> xyz	= grid.get_xyz_from_ijk_f(ijk);

	unsigned int na = srch_res.n_close_data + nodes.ncnode;
	int neq = na;


	KrigeSystem k_syst(na);
		//KrigeSystem_matrx kramt(na);
	if ( blocks.blocksflag == 1)
		KrigeSystem full_ksyst(na+blocks.ncloseblocks);
	//full_ksyst= new KrigeSystem((na+blocks.ncloseblocks));

	if (na <1 )
	{
		krgres.cmean=harddata.vmedexp;
		krgres.cstdev=1.0;

	//}else if(na==1)
	//{
	//	if ( blocks.blocksflag == 1) ising= (*full_ksyst).ksol(neq);
	//	else ising=(k_syst).ksol (neq);
	//	
	}
	else
	{// more than one, prepare system	

		ino = Set_matrixes(harddata, vario, grid,  nodes,na, srch_res,currpoint, k_syst );



		if ( blocks.blocksflag == 1) 
		{
			if (blocks.ncloseblocks>0)	
				blocks.allocate_blocks_matrixes(na);

			ino_b2b=blocks.Set_matrixes_b2b( incloseblocks,currpoint );

			ino_b2p=blocks.Set_matrixes_b2p(na, srch_res,
			incloseblocks, nodes);

		// now combine all matrixes

			blocks.combine_matrixes(k_syst, ino, na, *full_ksyst);

			neq += blocks.ncloseblocks;

			ising=(*full_ksyst).ksol (neq);

			if (blocks.ncloseblocks>0)	
				blocks.deallocate_blocks_matrixes();

		}
		else
		{
		// LVM has no constraints here 
			ising=(k_syst).ksol (neq);
		}
	}


	if ( blocks.blocksflag == 1) 
		LVM_K_Local_cdf (blocks,g_log.idbg, grid, harddata,srch_res, vario.cbb, 
				*nodes.close_node_index, (*full_ksyst), krgres,  incloseblocks, currpoint,sim, &na);

	else 
		LVM_K_Local_cdf (blocks,g_log.idbg, grid, harddata,srch_res, vario.cbb, 
				*nodes.close_node_index, (k_syst), krgres,  incloseblocks, currpoint,sim, &na);

	//delete k_syst;

	//if (blocks.blocksflag==1){
	//	delete full_ksyst;
	//}
}

void KrigePars::Ext_drift_krg( BlocksPars& blocks,
		HarddataPars& harddata ,search_results& srch_res,
		VariogramPars& vario, GSLibGridPars& grid, Nodes_search& nodes,
		boost::shared_array<float> &sim, long currpoint, boost::shared_array<float> &secondary , krig_results& krgres, long *incloseblocks)
{

	int *ino = new int;
	/////////////// todo put this working, check blocks and constraint interactions

	
	float *vrea;

	int temp;
	int ino_b2b, ino_b2p;



	vrea =(float*) malloc (sizeof(float)* (srch_res.nodmax+srch_res.ndmax+1+1));


	//grid.get_ijk_from_index_f(currpoint, &ix, &iy, &iz);
	//grid.get_xyz_from_ijk_f(ix ,iy,iz, &xx, &yy, &zz);
	
	cl_coord<int> ijk	= grid.get_ijk_from_index_f(currpoint);
	cl_coord<double> xyz	= grid.get_xyz_from_ijk_f(ijk);



	unsigned int na = srch_res.n_close_data + nodes.ncnode;
	unsigned int neq = na+1;

	if (na <1)
	{
		krgres.cmean = harddata.vmedexp;
		krgres.cstdev = 1.0;
		return;
	}

	temp = na+1; 

	KrigeSystem k_syst(temp);
		//KrigeSystem_matrx kramt(temp);
	if ( blocks.blocksflag == 1)
		KrigeSystem full_ksyst((na+blocks.ncloseblocks));


	//if (na <1 )
	//{
	//	*cmean=harddata.vmedexp;
	//	*cstdev=1.0;

	//}else if(na==1)
	//{
	//	if ( blocks.blocksflag == 1) ising= (*full_ksyst).ksol(neq);
	//	else ising=(k_syst).ksol (neq);
	//	
	//}
	//else
	//{// more than one, prepare system	

	*ino = Set_matrixes(harddata,vario, grid,nodes, na, srch_res,currpoint , k_syst);

	if ( blocks.blocksflag == 1) {


		if (blocks.ncloseblocks>0)	blocks.allocate_blocks_matrixes( na);


		ino_b2b=blocks.Set_matrixes_b2b(  incloseblocks, currpoint );
		(void)ino_b2b;

		ino_b2p=blocks.Set_matrixes_b2p(na, srch_res,
				incloseblocks, nodes);
		(void)ino_b2p;


		// now combine all matrixes

		blocks.combine_matrixes(k_syst, *ino, na, *full_ksyst );

		neq += blocks.ncloseblocks;


		Ext_drift_constr(*(ino+ino_b2b+ino_b2p), (na+blocks.ncloseblocks), srch_res,
			*nodes.close_node_index,&neq,  currpoint, vrea, secondary,&harddata.vsec[0], *full_ksyst ,&na);
		//Ordinary_constr(ino+ino_b2b+ino_b2p, (na+blocks.ncloseblocks),*full_ksyst);

		ising=full_ksyst->ksol (neq);

		if (blocks.ncloseblocks>0)	blocks.deallocate_blocks_matrixes();

	}else{


		Ext_drift_constr(*ino, na, srch_res,  *nodes.close_node_index, &neq,  currpoint,
			vrea, secondary,&harddata.vsec[0], k_syst, &na );

		ising=k_syst.ksol (neq);
	}


//}

	if ( blocks.blocksflag == 1) 
		ExtrD_K_Local_cdf  (blocks,g_log.idbg, grid, harddata, srch_res, vario.cbb,
				*nodes.close_node_index, *full_ksyst, krgres,  incloseblocks, currpoint,sim, &na);
	else 
		ExtrD_K_Local_cdf (blocks,g_log.idbg, grid, harddata,srch_res, vario.cbb, 
				*nodes.close_node_index, k_syst, krgres,  incloseblocks, currpoint,sim, &na);


	//delete k_syst;
	//if (blocks.blocksflag==1){
	//	delete full_ksyst;
	//}

	free(vrea);
	//dealloc_matrixes(blocks.blocksflag);


}

void KrigePars::ColocCosim_krg	(BlocksPars& blocks
		,HarddataPars& harddata ,search_results& srch_res,
		VariogramPars& vario, GSLibGridPars& grid,  Nodes_search& nodes,
		boost::shared_array<float> &sim, long currpoint, krig_results& krgres, long *incloseblocks)
{
	unsigned int na = srch_res.n_close_data + nodes.ncnode;
	int neq = na + 1;
	
	//if (na < 1) return;
	//temp= search.nodmax+search.ndmax+1+1;
	
	int temp= na + 1;

	KrigeSystem k_syst(temp);
		
	if (blocks.blocksflag == 1) 
	{
		int ino_b2b, ino_b2p;

		KrigeSystem full_ksyst((na + blocks.ncloseblocks));
		int ino = Set_matrixes(harddata, vario, grid, nodes, na, srch_res, currpoint, k_syst);

		if (blocks.ncloseblocks>0)	blocks.allocate_blocks_matrixes(na);


		ino_b2b = blocks.Set_matrixes_b2b(incloseblocks, currpoint);
		(void)ino_b2b;

		ino_b2p = blocks.Set_matrixes_b2p(na, srch_res, incloseblocks, nodes);
		(void)ino_b2p;

		// now combine all matrixes

		blocks.combine_matrixes(k_syst, ino, na, full_ksyst);

		neq += blocks.ncloseblocks;

		int totino = ino + ino_b2b + ino_b2p;
		ColocCosim_constr(totino, (na + blocks.ncloseblocks), colocorr, full_ksyst);

		ising = (full_ksyst).ksol(neq);

		if (blocks.ncloseblocks>0)	
			blocks.deallocate_blocks_matrixes();

		Colloc_CoK_Local_cdf(blocks, g_log.idbg, grid, harddata, srch_res, vario.cbb,
			*nodes.close_node_index, (full_ksyst), krgres, incloseblocks, currpoint, *nodes.close_node_value, sim, &na);

	}
	else //no blocks
	{
		if (na < 1)
		{
			krgres.cmean = harddata.vmedexp;
			krgres.cstdev = 1.0;
			return;
		}
		int ino = Set_matrixes(harddata, vario, grid, nodes, na, srch_res, currpoint, k_syst);
		
		ColocCosim_constr(ino, na, colocorr, k_syst);

		ising = (k_syst).ksol(neq);

		Colloc_CoK_Local_cdf(blocks, g_log.idbg, grid, harddata, srch_res, vario.cbb,
			*nodes.close_node_index, (k_syst), krgres, incloseblocks, currpoint, *nodes.close_node_value, sim, &na);
	}

	//}else if(na==1)
	//{
	//	if ( blocks.blocksflag == 1) ising= (*full_ksyst).ksol(neq);
	//	else ising=(k_syst).ksol (neq);
	//	
	//}
	//else
	//{// more than one, prepare system	


	//delete k_syst;
	//if (blocks.blocksflag==1){
	//	delete full_ksyst;
	//}

};

void KrigePars::LocalCC_ColocCosim_krg	(BlocksPars& blocks,
		HarddataPars& harddata ,search_results& srch_res,
		VariogramPars& vario, GSLibGridPars& grid,  Nodes_search& nodes,
	boost::shared_array<float> &sim, long currpoint, float localCC, krig_results& krgres, long *incloseblocks)
{
	unsigned int na = srch_res.n_close_data + nodes.ncnode;
	int neq = na + 1;

	int temp = na + 1;
	KrigeSystem k_syst(temp);
	

	if (blocks.blocksflag == 1)
	{
		int ino_b2b, ino_b2p;

		KrigeSystem full_ksyst((na + blocks.ncloseblocks));
		int ino = Set_matrixes(harddata, vario, grid, nodes, na, srch_res, currpoint, k_syst );
		
		if (blocks.ncloseblocks>0)	blocks.allocate_blocks_matrixes(na);


		ino_b2b = blocks.Set_matrixes_b2b(incloseblocks, currpoint);
		(void)ino_b2b;

		ino_b2p = blocks.Set_matrixes_b2p(na, srch_res,	incloseblocks, nodes);
		(void)ino_b2p;

		// now combine all matrixes

		blocks.combine_matrixes(k_syst, ino, na, full_ksyst);

		neq += blocks.ncloseblocks;

		int totino = ino + ino_b2b + ino_b2p;
		ColocCosim_constr(totino, (na + blocks.ncloseblocks), localCC, full_ksyst);


		ising = full_ksyst.ksol(neq);

		if (blocks.ncloseblocks>0)	
			blocks.deallocate_blocks_matrixes();
		
		Colloc_CoK_Local_cdf(blocks, g_log.idbg, grid, harddata, srch_res, vario.cbb,
			*nodes.close_node_index, (full_ksyst), krgres, incloseblocks, currpoint, *nodes.close_node_value, sim, &na);

	}
	else //no blocks
	{
		if (na < 1)
		{
			krgres.cmean = harddata.vmedexp;
			krgres.cstdev = 1.0;
			return;
		}
		int ino = Set_matrixes(harddata, vario, grid, nodes, na, srch_res, currpoint, k_syst );
		
		ColocCosim_constr(ino, na, localCC, k_syst);
		
		ising = (k_syst).ksol(neq);

		Colloc_CoK_Local_cdf(blocks, g_log.idbg, grid, harddata, srch_res, vario.cbb,
			*nodes.close_node_index, (k_syst), krgres, incloseblocks, currpoint, *nodes.close_node_value, sim, &na);
	}
	//if (na <1 )
	//{
	//	*cmean=harddata.vmedexp;
	//	*cstdev=1.0;
	//if (na < 1)
	//{
	//	krgres.cmean = harddata.vmedexp;
	//	krgres.cstdev = 1.0;
	//	return;
	//}
	//	////if(neq==1){
	//	//*cmean = /*localCC**/(secondary_grid.grid[currpoint] - harddata.vmedexp);
	//	//*cstdev = 1-localCC;
	//}else
	//if(na==1)
	//{
	//	if ( blocks.blocksflag == 1) ising= (*full_ksyst).ksol(neq);
	//	else ising=(k_syst).ksol (neq);
	//	
	//}
	//else
	//{// more than one, prepare system	
	//	int ino = Set_matrixes( harddata,vario, grid, nodes,na, srch_res,currpoint , k_syst /*,kramt*/);
	//	if ( blocks.blocksflag == 1) 
	//	{
	//		if (blocks.ncloseblocks>0)	blocks.allocate_blocks_matrixes(na);
	//		ino_b2b=blocks.Set_matrixes_b2b(	incloseblocks, currpoint  );
	//		(void)ino_b2b;
	//		ino_b2p=blocks.Set_matrixes_b2p(na, srch_res,
	//			incloseblocks, nodes);
	//		(void)ino_b2p;
	//		// now combine all matrixes
	//		blocks.combine_matrixes(k_syst, ino, na,  *full_ksyst );
	//		neq += blocks.ncloseblocks;
	//		int totino = ino + ino_b2b + ino_b2p;
	//		LocalCC_ColocCosim_constr(totino,(na+blocks.ncloseblocks),  localCC,*full_ksyst);
	//		ising = full_ksyst->ksol (neq);
	//		if (blocks.ncloseblocks>0)	blocks.deallocate_blocks_matrixes();
	//	}
	//	else
	//	{
	//		LocalCC_ColocCosim_constr(ino, na,  localCC,k_syst);
	//		ising=k_syst.ksol (neq);
	//	}
	////}
	//if ( blocks.blocksflag == 1) 
	//	Colloc_CoK_Local_cdf (blocks,g_log.idbg, grid, harddata,srch_res, vario.cbb, 
	//			*nodes.close_node_index, (*full_ksyst), krgres,  incloseblocks, currpoint,sim, &na);
	//else 
	//	Colloc_CoK_Local_cdf (blocks,g_log.idbg, grid, harddata,srch_res, vario.cbb, 
	//			*nodes.close_node_index, (k_syst), krgres,  incloseblocks, currpoint,sim, &na);
	//delete k_syst;
	//if (blocks.blocksflag==1)
	//{
	//	delete full_ksyst;
	//}

}


void KrigePars::LVM_ColocCosim_krg	(BlocksPars& blocks
    ,HarddataPars& harddata ,search_results& srch_res,
		VariogramPars& vario, GSLibGridPars& grid, Nodes_search& nodes,
		boost::shared_array<float> &sim, long currpoint, float localCC  , krig_results& krgres, long *incloseblocks)
{
	int *ino = new int;

	int temp;
	int ino_b2b, ino_b2p;

	cl_coord<int> ijk	= grid.get_ijk_from_index_f(currpoint);
	cl_coord<double> xyz	= grid.get_xyz_from_ijk_f(ijk);

	unsigned int na = srch_res.n_close_data + nodes.ncnode;
	int neq = na + 1;

	temp = na + 1;
	
	KrigeSystem k_syst(temp);
	//KrigeSystem_matrx kramt(temp);


	//if (na <1 )
	//{				
	//	if ( blocks.blocksflag == 1) 
	//		ising = full_ksyst->ksol(neq);
	//	else 
	//		ising = k_syst.ksol (neq);
	//	
	//} 
	//else
	//{// more than one, prepare system	

	*ino = Set_matrixes( harddata, vario, grid, nodes,na, srch_res,currpoint , k_syst );

		if ( blocks.blocksflag == 1) 
		{		
			KrigeSystem full_ksyst((na+blocks.ncloseblocks));

			if (blocks.ncloseblocks>0)	blocks.allocate_blocks_matrixes(na);


			ino_b2b=blocks.Set_matrixes_b2b( incloseblocks, currpoint  );
			(void)ino_b2b;

			ino_b2p=blocks.Set_matrixes_b2p(na, srch_res,
				incloseblocks, nodes);
			(void)ino_b2p;

			
			// now combine all matrixes

			blocks.combine_matrixes(k_syst, *ino, na,  full_ksyst );

			neq += blocks.ncloseblocks;
			int *totino = new int;
			*totino = *ino + ino_b2b + ino_b2p;
			ColocCosim_constr(*totino,(na+blocks.ncloseblocks),  localCC,full_ksyst);


			ising = full_ksyst.ksol (neq);

			if (blocks.ncloseblocks>0)	
				blocks.deallocate_blocks_matrixes();

		}	
		else
		{
			ColocCosim_constr(*ino, na,  localCC,k_syst);
			ising=k_syst.ksol (neq);
		}

	//}


	if ( blocks.blocksflag == 1) 
		LVM_CoK_Local_cdf (blocks,g_log.idbg, grid, harddata,srch_res, vario.cbb, 
				*nodes.close_node_index, (*full_ksyst), krgres,  incloseblocks, currpoint,sim, &na);
	else 
		LVM_CoK_Local_cdf (blocks,g_log.idbg, grid, harddata,srch_res, vario.cbb, 
				*nodes.close_node_index, (k_syst), krgres,  incloseblocks, currpoint,sim, &na);

}



void KrigePars::O_ColocCosim_krg(BlocksPars& blocks,
	HarddataPars& harddata, search_results& srch_res,
	VariogramPars& vario, GSLibGridPars& grid, Nodes_search& nodes,
	boost::shared_array<float> &sim, long currpoint, krig_results& krgres, long *incloseblocks, float secondary_avg)
{
	int *ino=new int;

	int temp;
	int ino_b2b, ino_b2p;


	cl_coord<int> ijk = grid.get_ijk_from_index_f(currpoint);
	cl_coord<double> xyz = grid.get_xyz_from_ijk_f(ijk);


	unsigned int na = srch_res.n_close_data + nodes.ncnode;
	int neq = na + 1 + 1;
	//if (na < 1) return;

	//temp= search.nodmax+search.ndmax+1+1;
	temp = na + 1 + 1;

	KrigeSystem k_syst(temp);
	//KrigeSystem_matrx kramt(temp);
	if (blocks.blocksflag == 1)
		KrigeSystem full_ksyst((na + blocks.ncloseblocks));

	//if (na <1 )
	//{
	//	*cmean=harddata.vmedexp;
	//	*cstdev=1.0;
	if (na < 1)
	{
		krgres.cmean = harddata.vmedexp;
		krgres.cstdev = 1.0;
		return;
	}

	//}else if(na==1)
	//{
	//	if ( blocks.blocksflag == 1) ising= (*full_ksyst).ksol(neq);
	//	else ising=(k_syst).ksol (neq);
	//	
	//}
	//else
	//{// more than one, prepare system	



	*ino = Set_matrixes(harddata, vario, grid, nodes, na, srch_res, currpoint, k_syst);

	// CUIDADO COM A MUDANÇA PARA OK
	if (blocks.blocksflag == 1) 
	{


		if (blocks.ncloseblocks>0)	blocks.allocate_blocks_matrixes(na);


		ino_b2b = blocks.Set_matrixes_b2b(incloseblocks, currpoint);
		(void)ino_b2b;

		ino_b2p = blocks.Set_matrixes_b2p(na, srch_res,
			incloseblocks, nodes);
		(void)ino_b2p;

		// now combine all matrixes

		blocks.combine_matrixes(k_syst, *ino, na, *full_ksyst);

		neq += blocks.ncloseblocks;
		
		int *totino = new int;
		*totino = *ino + ino_b2b + ino_b2p;

		ColocCosim_constr(*totino, (na + blocks.ncloseblocks), colocorr, *full_ksyst);
		Ordinary_constr(*totino, (na + blocks.ncloseblocks) + 1, k_syst);

		ising = (*full_ksyst).ksol(neq);

		if (blocks.ncloseblocks>0)	blocks.deallocate_blocks_matrixes();


	}
	else
	{
		ColocCosim_constr(*ino, na, colocorr, k_syst);
		Ordinary_constr(*ino, na+1, k_syst);
		

		ising = (k_syst).ksol(neq);
	}

	//}
	if (blocks.blocksflag == 1)
		O_Colloc_CoK_Local_cdf(blocks, g_log.idbg, grid, harddata, srch_res, vario.cbb,
			*nodes.close_node_index, (*full_ksyst), krgres, incloseblocks, currpoint, sim, &na, secondary_avg);
	else
		O_Colloc_CoK_Local_cdf(blocks, g_log.idbg, grid, harddata, srch_res, vario.cbb,
			*nodes.close_node_index, (k_syst), krgres, incloseblocks, currpoint, sim, &na, secondary_avg);

	//dealloc_matrixes(blocks.blocksflag);

	//delete k_syst;
	//if (blocks.blocksflag==1){
	//	delete full_ksyst;
	//}

};

void KrigePars::O_LocalCC_ColocCosim_krg(BlocksPars& blocks,
	HarddataPars& harddata, search_results& srch_res,
	VariogramPars& vario, GSLibGridPars& grid, Nodes_search& nodes,
	boost::shared_array<float> &sim, long currpoint, float localCC, krig_results& krgres, long *incloseblocks, float secondary_avg)
{
	int *ino = new int;

	int temp;
	int ino_b2b, ino_b2p;

	cl_coord<int> ijk = grid.get_ijk_from_index_f(currpoint);
	cl_coord<double> xyz = grid.get_xyz_from_ijk_f(ijk);

	unsigned int na = srch_res.n_close_data + nodes.ncnode;
	int neq = na + 1+1;

	temp = na + 1+1;
	//k_syst = new KrigeSystem(temp);
	KrigeSystem k_syst(temp);
	//KrigeSystem_matrx kramt(temp);
	if (blocks.blocksflag == 1)
		KrigeSystem full_ksyst((na + blocks.ncloseblocks));
	//full_ksyst= new KrigeSystem((na+1+blocks.ncloseblocks));

	//if (na <1 )
	//{
	//	*cmean=harddata.vmedexp;
	//	*cstdev=1.0;
	if (na < 1)
	{
		krgres.cmean = harddata.vmedexp;
		krgres.cstdev = 1.0;
		return;
	}

	//	////if(neq==1){
	//	//*cmean = /*localCC**/(secondary_grid.grid[currpoint] - harddata.vmedexp);
	//	//*cstdev = 1-localCC;

	//}else
	//if(na==1)
	//{
	//	if ( blocks.blocksflag == 1) ising= (*full_ksyst).ksol(neq);
	//	else ising=(k_syst).ksol (neq);
	//	
	//}
	//else
	//{// more than one, prepare system	

	*ino = Set_matrixes(harddata, vario, grid, nodes, na, srch_res, currpoint, k_syst );

	if (blocks.blocksflag == 1)
	{

		if (blocks.ncloseblocks>0)	blocks.allocate_blocks_matrixes(na);


		ino_b2b = blocks.Set_matrixes_b2b(incloseblocks, currpoint);
		(void)ino_b2b;

		ino_b2p = blocks.Set_matrixes_b2p(na, srch_res,
			incloseblocks, nodes);
		(void)ino_b2p;

		// now combine all matrixes

		blocks.combine_matrixes(k_syst, *ino, na, *full_ksyst);

		neq += blocks.ncloseblocks;

		int *totino = new int;
		*totino = *ino + ino_b2b + ino_b2p;
		ColocCosim_constr(*totino, (na + blocks.ncloseblocks), localCC, *full_ksyst);
		Ordinary_constr(*totino, (na + blocks.ncloseblocks + 1), k_syst);
		
		delete totino;

		ising = full_ksyst->ksol(neq);

		if (blocks.ncloseblocks>0)	blocks.deallocate_blocks_matrixes();

	}
	else
	{

		ColocCosim_constr(*ino, na, localCC, k_syst);
		Ordinary_constr(*ino, na + 1, k_syst);
		
		delete ino;
		ising = k_syst.ksol(neq);
	}

	//}


	if (blocks.blocksflag == 1)
		O_Colloc_CoK_Local_cdf(blocks, g_log.idbg, grid, harddata, srch_res, vario.cbb,
			*nodes.close_node_index, (*full_ksyst), krgres, incloseblocks, currpoint, sim, &na, secondary_avg);
	else
		O_Colloc_CoK_Local_cdf(blocks, g_log.idbg, grid, harddata, srch_res, vario.cbb,
			*nodes.close_node_index, (k_syst), krgres, incloseblocks, currpoint, sim, &na, secondary_avg);

	return;
}


int KrigePars::Set_matrixes( HarddataPars& harddata , VariogramPars& vario,
							GSLibGridPars& grid, Nodes_search& nodes,unsigned int na, 
							search_results& srch_res, long& currentpoint, KrigeSystem& k_system )
{
	// i can get the values matrix here and save a cycle, i think 
	//std::vector<float> sample_values;


	cl_coord<int> current_pnt_ijk	= grid.get_ijk_from_index_f(currentpoint);
	cl_coord<double> current_pnt_xyz	= grid.get_xyz_from_ijk_f(current_pnt_ijk);
	
	cl_coord<int> ijk_j;
	cl_coord<int> ijk2;

	//cl_coord<int> vario_nct(vario.nctx,vario.ncty,vario.nctz);
	cl_coord<int> coord_ones(1,1,1);
	cl_coord<int> iijjkk;

	long index;
	long ind;

	//float cov;

	// Set up kriging matrices:
		
	int ino=0;

	cl_coord<double> xyz_j;
	cl_coord<double> xyz2;

	//separated the for loops
	
	//hard data loop
	for (unsigned int j = 0; j < srch_res.n_close_data; ++j)
	{
		// Sort out the actual location of point "j"
		// we get the coordinate from the hardata, as this is a sample
		//index = srch_res.inclose[j];
		//xyz_j = harddata.point[index].xyz;
		xyz_j =srch_res.found_sample_values[j].xyz;
	
		// Sort out the actual location of point "i"			

		for (unsigned int i = 0; i <= j; ++i)
		{
			//if (i <= srch_res.n_close_data) //useless
			//{
				//index = srch_res.inclose[i];
				//xyz2 = harddata.point[index].xyz;
				xyz2 = srch_res.found_sample_values[i].xyz;
			//}

			// Now, get the covariance value:
			++ino;

			k_system.a[ino] = (double)vario.covariance(xyz_j, xyz2);
			
		}

		// Get the RHS value (possibly with covariance look-up table):
	
		k_system.r[j+1] = (double)vario.covariance(current_pnt_xyz, xyz_j);
		k_system.rr[j+1] = k_system.r[j+1];
		
		
		//save the values
		//sample_values.push_back(srch_res.found_sample_values[j].value);
	}

	//nodes loop
	for (unsigned int j = srch_res.n_close_data; j < na; ++j)
	{
		// Sort out the actual location of point "j"
		// It is a previously simulated node (keep index for table look-up):
		// we get index from lookup table

		index = j - srch_res.n_close_data;// -1;// -1 is for c
		ind = nodes.icnode->at(index);
//gridcoord found // curr p grid coord //				// halfcovtab //fortran offset?? 
		//ijk_j = current_pnt_ijk + (vario.ijk_node[ind] - vario.nctxyz - coord_ones);// or i can just use  
		ijk_j = grid.get_ijk_from_index_f(nodes.close_node_index->at(index));

		xyz_j = grid.get_xyz_from_index_f(nodes.close_node_index->at(index));
		

		// Sort out the actual location of point "i"			

		for (unsigned int i = 0; i <= j; ++i)
		{
			if (i < srch_res.n_close_data)// this is harddata
			{
				//index = srch_res.inclose[i/* - 1*/];
				//xyz2 = harddata.point[index].xyz;
				xyz2 = srch_res.found_sample_values[i].xyz;
			}
			else // It is a previously simulated node (keep index for table look-up):
			{
				index = i - srch_res.n_close_data;// -1;// -1 is for c
				//ind = nodes.icnode->at(index);
				
				ijk2 = grid.get_ijk_from_index_f(nodes.close_node_index->at(index));
				//ijk2 = current_pnt_ijk + (vario.ijk_node[ind] - vario.nctxyz - coord_ones);

				xyz2 = grid.get_xyz_from_index_f(nodes.close_node_index->at(index));
			}

			// Now, get the covariance value:
			++ino;

			// Decide whether or not to use the covariance look-up table:
			if ( i < srch_res.n_close_data)// dealing with harddata
			{
				k_system.a[ino] = (double)vario.covariance(xyz_j, xyz2);
			}
			else // Try to use the covariance look-up (if the distance is in range): // only nodes
			{

				iijjkk = vario.nctxyz + coord_ones + (ijk_j - ijk2);

				if (iijjkk.x < 1 || iijjkk.x >   vario.nctxyz.x * 2 + 1 ||
					iijjkk.y < 1 || iijjkk.y >   vario.nctxyz.y * 2 + 1 ||
					iijjkk.z < 1 || iijjkk.z >   vario.nctxyz.z * 2 + 1)
				{
					k_system.a[ino] = (double)vario.covariance(xyz_j, xyz2);
				}
				else
				{
					k_system.a[ino] = (double)vario.covtab[iijjkk.x][iijjkk.y][iijjkk.z];
				}
			}


			// remove covarinace to do point anisotropy.
			//k_system.a[ino] = (double)vario.covariance(xyz_j,xyz2);


		}

		// Get the RHS value (possibly with covariance look-up table):

		//iijjkk = vario.nctxyz + coord_ones + (current_pnt_ijk - ijk_j);//pq é que não é o vario.ijk_node apenas?? -palpita-me que isto tá mal


		iijjkk = vario.ijk_node[ind]; //TO DO RETEST THIS

		if (iijjkk.x < 1 || iijjkk.x >   vario.nctxyz.x * 2 + 1 ||
				iijjkk.y < 1 || iijjkk.y >   vario.nctxyz.y * 2 + 1 ||
				iijjkk.z < 1 || iijjkk.z >   vario.nctxyz.z * 2 + 1)
		{
			k_system.r[j+1] = (double)vario.covariance(current_pnt_xyz, xyz_j);
		}
		else
		{
			k_system.r[j+1] = (double)vario.covtab[iijjkk.x][iijjkk.y][iijjkk.z];
		}
		k_system.rr[j+1] = k_system.r[j+1];
		

		index = j - srch_res.n_close_data;
		// sample_values.push_back(nodes.close_node_value->at(index));
	}



	//for (unsigned int j = 1; j<=na; ++j)
	//{
	//	// Sort out the actual location of point "j"
	//	if (j <= srch_res.n_close_data)
	//	{
	//		// we get the coordinate from the hardata, as this is a sample
	//		index = srch_res.inclose[j - 1];
	//		xyz_j = harddata.point[index].xyz; 
	//	}
	//	else 
	//	{
	//		// It is a previously simulated node (keep index for table look-up):
	//		// we get index from lookup table

	//		index = j - srch_res.n_close_data - 1;// -1 is for c
	//		ind    = nodes.icnode->at(index);

	//		ijk_j = current_pnt_ijk+(vario.ijk_node[ind]-vario.nctxyz-coord_ones);

	//		xyz_j = grid.get_xyz_from_index_f(nodes.close_node_index->at(index));
	//	}

	//	// Sort out the actual location of point "i"			
	//	
	//	for (unsigned int i=1; i<=j; ++i)
	//	{
	//		if (i <= srch_res.n_close_data)
	//		{
	//			index = srch_res.inclose[i - 1];
	//			xyz2 = harddata.point[index].xyz;
	//		}
	//		else // It is a previously simulated node (keep index for table look-up):
	//		{
	//			index = i - srch_res.n_close_data - 1;// -1 is for c
	//			ind    = nodes.icnode->at(index);

	//			ijk2 = current_pnt_ijk + (vario.ijk_node[ind]-vario.nctxyz-coord_ones);

	//			xyz2 = grid.get_xyz_from_index_f(nodes.close_node_index->at(index));
	//		}

	//		// Now, get the covariance value:
	//		++ino;

	//		// Decide whether or not to use the covariance look-up table:
	//		if (j <= srch_res.n_close_data || i <= srch_res.n_close_data)
	//		{
	//			k_system.a[ino] = (double)vario.covariance(xyz_j,xyz2);
	//		}
	//		else // Try to use the covariance look-up (if the distance is in range):
	//		{

	//			iijjkk= vario.nctxyz + coord_ones +(ijk_j-ijk2);
	//			
	//			if (iijjkk.x < 1 || iijjkk.x >   vario.nctxyz.x*2+1 ||
	//				iijjkk.y < 1 || iijjkk.y >   vario.nctxyz.y*2+1 ||
	//					iijjkk.z < 1 || iijjkk.z >   vario.nctxyz.z*2+1)
	//			{
	//				k_system.a[ino] =(double)vario.covariance(xyz_j,xyz2);
	//			}
	//			else
	//			{
	//				k_system.a[ino] =  (double)vario.covtab[iijjkk.x][iijjkk.y][iijjkk.z];
	//			}
	//		}


	//		// remove covarinace to do point anisotropy.
	//		//k_system.a[ino] = (double)vario.covariance(xyz_j,xyz2);


	//	}

	//	// Get the RHS value (possibly with covariance look-up table):
	//	if (j <= srch_res.n_close_data)
	//	{
	//		k_system.r[j]  = (double)vario.covariance(current_pnt_xyz,xyz_j);
	//		k_system.rr[j]=k_system.r[j];
	//	}
	//	else // Try to use the covariance look-up (if the distance is in range):
	//	{

	//		iijjkk= vario.nctxyz + coord_ones +(current_pnt_ijk-ijk_j);

	//		if (iijjkk.x < 1 || iijjkk.x >   vario.nctxyz.x*2+1 ||
	//			iijjkk.y < 1 || iijjkk.y >   vario.nctxyz.y*2+1 ||
	//			iijjkk.z < 1 || iijjkk.z >   vario.nctxyz.z*2+1)
	//		{
	//			k_system.r[j] = (double)vario.covariance(current_pnt_xyz,xyz_j);
	//		}
	//		else
	//		{
	//			k_system.r[j] = (double)vario.covtab[iijjkk.x][iijjkk.y][iijjkk.z];
	//		}
	//		k_system.rr[j] = k_system.r[j];
	//	}

	//	// dont use covariance hero too 

	//	//k_system.r[j] = (double)vario.covariance(current_pnt_xyz, xyz_j);
	//	//k_system.rr[j] = k_system.r[j];


	//}

	return ino;
}




	//void calculate_point_local_CDF(SearchPars& search, HarddataPars& harddata, 
	//	KrigeSystem& k_syst, krig_results *krgres, long *close_node_index, float *sim, int* na, float mean_multiplier)
	//{
	//	long index;
	//	float vra;

	//	for (int j = 1; j <= *na; ++j)
	//	{
	//		if (j <= search.n_close_data)
	//		{
	//			index = search.inclose[j];
	//			vra = harddata.point[index].value;
	//		}
	//		else
	//		{
	//			index = j - search.n_close_data;
	//			vra = sim[close_node_index[index]];
	//		}
	//		krgres->cmean += (float)(k_syst.s[j]) * (vra - average)mean_multiplier;
	//		krgres->cstdev -= (float)(k_syst.s[j] * k_syst.rr[j]);
	//		krgres->sumwts += (float)(k_syst.s[j]);

	//	}
	//};

