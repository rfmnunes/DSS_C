#include "krige.h"

// Constraints addition

//	void Simple_constr(){} // has none

void KrigePars::Ordinary_constr(int& ino, unsigned int nlines, KrigeSystem& k_system)	     // Addition of OK constraint:
{

	for (int i = 0; i < nlines; ++i)
	{
		++ino;
		k_system.a[ino] = 1.0;
	}

	++ino;
	k_system.a[ino] = 0.0;
	k_system.r[nlines + 1] = 1.0;
	k_system.rr[nlines + 1] = 1.0;
}

//	void LVM_constr(){}// has none here

void KrigePars::Ext_drift_constr(int& ino, unsigned int nlines, search_results& srch_res, 
	std::vector<long> &close_node_index, unsigned int *neq, long &currpoint, float *vrea, 
	boost::shared_array<float> &secondary, float *vsec, KrigeSystem& k_system, unsigned int* na)   // Addition of the External Drift Constraint:
{

	long index;

	for (unsigned int j = 0; j < nlines; ++j)
	{
		if (j + 1 <= srch_res.n_close_data)
		{
			index = srch_res.inclose[j];
			vrea[j+1] = vsec[index];
		}
		else
		{
			index = j+1 - srch_res.n_close_data;
			vrea[j+1] = secondary[close_node_index.at(index-1)];
		}
	}
	float edmin = (float) 1.0e21;
	float edmax = (float)-1.0e21;

	for (unsigned int i = 0; i < nlines; ++i)
	{
		++ino;
		k_system.a[ino] = vrea[i +1 ];
		
		if (k_system.a[ino] < edmin) 
			edmin = (float)k_system.a[ino];
		if (k_system.a[ino] > edmax) 
			edmax = (float)k_system.a[ino];
	}
	k_system.r[nlines + 1] = (double)secondary[currpoint];
	k_system.rr[nlines + 1] = k_system.r[*na + 1];
	if ((edmax - edmin) < EPSLON) 
		*neq = *neq - 1;

}

void KrigePars::ColocCosim_constr(int& ino, unsigned int nlines, float &colocorr, KrigeSystem& k_system)   // Addition of Collocated Cosimulation Constraint:
{
	int  ii;
	float edmin = (float) 1.0e21;
	float edmax = (float)-1.0e21;
	for (int i = 0; i < nlines; ++i)
	{
		++ino;
		k_system.a[ino] = (double)(colocorr*k_system.r[i + 1]);
		
		if (k_system.a[ino] < edmin) 
			edmin = (float)k_system.a[ino];
		if (k_system.a[ino] > edmax) 
			edmax = (float)k_system.a[ino];
	}
	++ino;
	k_system.a[ino] = 1.0;
	ii = nlines + 1;
	k_system.r[ii] = (double)colocorr;
	k_system.rr[ii] = k_system.r[ii];

}

//only one is needed
//void KrigePars::LocalCC_ColocCosim_constr(int& ino, unsigned int nlines, float &localCC, KrigeSystem& k_system)
//// Addition of Collocated Cosimulation Constraint: - local cc
//{
//	int ii;
//	float edmin = (float) 1.0e21;
//	float edmax = (float)-1.0e21;
//	// index  = ix + (iy-1)*grid_def.nx + (iz-1)*grid_def.nxy;
//
//	for (int i = 0; i < nlines; ++i)
//	{
//		++ino;
//		k_system.a[ino] = (double)(localCC*k_system.r[i + 1]);
//		
//		if (k_system.a[ino] < edmin) 
//			edmin = (float)k_system.a[ino];
//		if (k_system.a[ino] > edmax) 
//			edmax = (float)k_system.a[ino];
//	}
//	++ino;
//	k_system.a[ino] = 1.0;
//	ii = nlines + 1;
//	k_system.r[ii] = (double)localCC;
//	k_system.rr[ii] = k_system.r[ii];
//
//}
//
