#ifndef UTILS_INCLUDED
#define UTILS_INCLUDED


#define GEOEAS 0 // this for gslib geoeas file type
#define SGEMS 1 // this for sgems binary file type

//#include <math.h>
#include <sstream>
//
#include <iostream>     // logger
#include <fstream> 

#include "registry.h"
#include "boost_include.h"

#include "log.h"
#include "logging.h"



namespace utilities
{
	template <class templt>
	inline std::string make_string(templt in_number)
	{
		std::ostringstream out_str;
		out_str << in_number;

		return out_str.str();
	};

};


class misc{
public:

	float nosvalue;

	//use headers?
	int headerflag;


	// simulation parameters
	std::string outfl;
	int      icmean,icvar;
	unsigned int      ntry;

	// multiple grid par
	int mults, nmult;

	int nsim;

	int filetype;

	int read_miscelaneous_parameters(registry *reg)
	{
		reg_key *k;

		k = get_key(reg, (char*)("SEARCH"), (char*)("MULTS"));
		if (k)
			mults = get_int(k);
		else return -1;

		k = get_key(reg, (char*)("SEARCH"), (char*)("NMULTS"));
		if (k)
			nmult = get_int(k);
		else return -1;
		printf (" Multiple grid search (0/1)/number nodes: %d %d \n",mults,nmult);


		// parse simulation parameters
		if ((k = get_key(reg, (char*)"SIMULATION", (char*)"OUTFILE")) != NULL)
			outfl= get_string(k);

		boost::algorithm::trim(outfl);

		logger << " Output file: " << outfl << "\n";

		k = get_key(reg, (char*)("SIMULATION"), (char*)("NSIMS"));
		if (k)
			nsim = get_int(k);
		else return -1;

		std::ostringstream nsim_str;
		nsim_str << nsim;

		logger << " Realizations (#): " << nsim_str.str()/*std::to_strin(nsim)*/ << "\n";

		k = get_key(reg, (char*)("SIMULATION"), (char*)("NTRY"));
		if (k)
			ntry = get_int(k);
		else return -1;

		k = get_key(reg, (char*)("SIMULATION"), (char*)("AVGCORR"));\
		if (k)
			icmean = get_int(k);
		else return -1;

		k = get_key(reg, (char*)("SIMULATION"), (char*)("VARCORR"));
		if (k)
			icvar = get_int(k);
		else return -1;

		std::ostringstream ntry_str;
		ntry_str <<  ntry;
		std::ostringstream icmean_str;
		icmean_str <<  icmean;
		std::ostringstream icvar_str;
		icvar_str << icvar;

		logger << " Bias Correction: number of simulations = " <<  ntry_str.str() << "\n"; 
		logger << " Correction for mean = " <<  icmean_str.str() << "\n";
		logger << " Correction for variance = " <<  icvar_str.str() << "\n";


		k = get_key(reg, (char*)("GENERAL"), (char*)("USEHEADERS"));
		if (k)
			headerflag = get_int(k);
		else return -1;

		// parsing general section
		k = get_key(reg, (char*)("GENERAL"), (char*)("NULLVAL"));
		if (k)
			nosvalue= get_float(k);
		else return -1;		

		filetype = SGEMS;;
		if ((k = get_key(reg, (char*)"GENERAL", (char*)"FILETYPE")) != NULL)
		{
			std::string buffer = get_string(k);
			boost::algorithm::trim(buffer);

			if (!strcmp(buffer.c_str(), "GEOEAS")) filetype = GEOEAS;
			else if (!strcmp(buffer.c_str(), "SGEMS"))filetype = SGEMS;
			
		}

		std::ostringstream nosvalue_str;
		nosvalue_str << nosvalue;

		logger << " Not simulated nodes: "<<nosvalue_str.str()/*std::to_strin(nosvalue) */<<"\n";

		
		//k = get_key(reg,  ("DEBUG"),  ("DBGLEVEL"));
		//if (k)
		//	g_log.idbg = get_int(k);
		//else return -1;

		//if ((k = get_key(reg,  ("DEBUG"),  ("DBGFILE"))) != NULL)
		//	g_log.dbgfl= get_string(k);
	
		//boost::algorithm::trim(g_log.dbgfl);

		//logger<< " Debugging file (1/2/3): " << g_log.dbgfl << " "<< g_log.idbg << "\n";

		return 0;

	
	}

	
	float average_correction(float current_point, double &current_average, double &target_average )
	{
		float corrected = current_point - (current_average - target_average);
		return corrected;

		//COSTUMAVA SER
		//cmean -= (float)M_curr[currzone] - zones.zone[currzone].harddata.vmedexp;

	}
	float variance_correction(float in_value, float current_variance, float target_variance) 
	{
		//float curr_vari = (float)(current_variance / current_n);
		if (current_variance > 0.0)
		{
			//cpdev=(float)sqrt(corr_vari);
			float corrected = in_value* (float)(/*sqrt*/(target_variance) / /*sqrt*/(current_variance));
			return corrected;
		}
		else
			return in_value;
		
		//COSTUMAVA SER
		//float corr_vari = (float)(S_curr[currzone] / (nnodesim[currzone]));

		//if (corr_vari > 0.0)
		//{
		//	//cpdev=(float)sqrt(corr_vari);
		//	cstdev *= (float)sqrt(zones.zone[currzone].harddata.vvarexp) / sqrt(corr_vari);
		//}

	}



	
	template <class SomeType,class OtherType,class stillothertype >
	//template <class OtherType>
	void sortit (long long ib, long long ie, OtherType *a, int iperm,  SomeType *b, stillothertype *c, float d[],
			float e[], float f[], float g[], float h[])
{


	//void sortit(long ib, long ie, float a[], int iperm, float b[], float c[], float d[], 
	//		float e[], float f[], float g[], float h[])
	//{

		// The dimensions for lt and ut have to be at least log (base 2) n
		int iring,lt[64],ut[64],i,j,k,m,p,q;
		float ta,td,te,tf,tg,th,xa,xd,xe,xf,xg,xh;

		SomeType tb,xb;
		stillothertype tc, xc;

		// Initialize:
		j     = ie;
		m     = 1;
		i     = ib;
		iring = iperm+1;
		if (iperm > 7) iring=1;

		// If this segment has more than two elements  we split it
		A: if (j-i-1 < 0)
			goto B;
		else if (j-i-1 == 0)
			goto D;
		else if (j-i-1 > 0)
			goto K;

		// p is the position of an arbitrary element in the segment we choose the
		// middle element. Under certain circumstances it may be advantageous
		// to choose p at random.
		K: p    = (j+i)/2;
		ta   = a[p];
		a[p] = a[i];
		if (iring >= 2) {
			tb   = b[p];
			b[p] = b[i];
		}
		if (iring >= 3) {
			tc   = c[p];
			c[p] = c[i];
		}
		if (iring >= 4) {
			td   = d[p];
			d[p] = d[i];
		}
		if (iring >= 5) {
			te   = e[p];
			e[p] = e[i];
		}
		if (iring >= 6) {
			tf   = f[p];
			f[p] = f[i];
		}
		if (iring >= 7) {
			tg   = g[p];
			g[p] = g[i];
		}
		if (iring >= 8) {
			th   = h[p];
			h[p] = h[i];
		}

		//Start at the beginning of the segment, search for k such that a[k]>t
		q = j;
		k = i;
		H: k = k+1;
		if(k > q) goto F;
		if(a[k] <= ta) goto H;

		// Such an element has now been found now search for a q such that a[q]<t
		// starting at the end of the segment.
		I:
		if(a[q] < ta) goto J;
		q = q-1;
		if(q > k)     goto I;
		goto E;

		// a[q] has now been found. we interchange a[q] and a[k]
		J: xa   = a[k];
		a[k] = a[q];
		a[q] = xa;
		if (iring >= 2) {
			xb   = b[k];
			b[k] = b[q];
			b[q] = xb;
		}
		if (iring >= 3) {
			xc   = c[k];
			c[k] = c[q];
			c[q] = xc;
		}
		if (iring >= 4) {
			xd   = d[k];
			d[k] = d[q];
			d[q] = xd;
		}
		if (iring >= 5) {
			xe   = e[k];
			e[k] = e[q];
			e[q] = xe;
		}
		if (iring >= 6) {
			xf   = f[k];
			f[k] = f[q];
			f[q] = xf;
		}
		if (iring >= 7) {
			xg   = g[k];
			g[k] = g[q];
			g[q] = xg;
		}
		if (iring >= 8) {
			xh   = h[k];
			h[k] = h[q];
			h[q] = xh;
		}

		// Update q and search for another pair to interchange:
		q = q-1;
		goto H;
		E: q = k-1;
		F:

		// The upwards search has now met the downwards search:
		a[i]=a[q];
		a[q]=ta;
		if (iring >= 2) {
			b[i] = b[q];
			b[q] = tb;
		}
		if (iring >= 3) {
			c[i] = c[q];
			c[q] = tc;
		}
		if (iring >= 4) {
			d[i] = d[q];
			d[q] = td;
		}
		if (iring >= 5) {
			e[i] = e[q];
			e[q] = te;
		}
		if (iring >= 6) {
			f[i] = f[q];
			f[q] = tf;
		}
		if (iring >= 7) {
			g[i] = g[q];
			g[q] = tg;
		}
		if (iring >= 8) {
			h[i] = h[q];
			h[q] = th;
		}

		// The segment is now divided in three parts: (i,q-1),[q],(q+1,j)
		// store the position of the largest segment in lt and ut
		if (2*q <= i+j) goto G;
		lt[m] = i;
		ut[m] = q-1;
		i = q+1;
		goto C;
		G: lt[m] = q+1;
		ut[m] = j;
		j = q-1;

		//Update m and split the new smaller segment
		C: m = m+1;
		goto A;

		// We arrive here if the segment has  two elements we test to see if
		// the segment is properly ordered if not, we perform an interchange
		D:
		if (a[i] <= a[j]) goto B;
		xa=a[i];
		a[i]=a[j];
		a[j]=xa;
		if (iring >= 2) {
			xb   = b[i];
			b[i] = b[j];
			b[j] = xb;
		}
		if (iring >= 3) {
			xc   = c[i];
			c[i] = c[j];
			c[j] = xc;
		}
		if (iring >= 4) {
			xd   = d[i];
			d[i] = d[j];
			d[j] = xd;
		}
		if (iring >= 5) {
			xe   = e[i];
			e[i] = e[j];
			e[j] = xe;
		}
		if (iring >= 6) {
			xf   = f[i];
			f[i] = f[j];
			f[j] = xf;
		}
		if (iring >= 7) {
			xg   = g[i];
			g[i] = g[j];
			g[j] = xg;
		}
		if (iring >= 8) {
			xh   = h[i];
			h[i] = h[j];
			h[j] = xh;
		}

		// If lt and ut contain more segments to be sorted repeat process:
		B: m = m-1;
		if (m <= 0) 
		{
			return;
		}
		i = lt[m];
		j = ut[m];
		goto A;

		return;
	};

	/*

	//C quicksort , from numerical recipes
	template<class T, class U> //sort.h 
	void sort2(NRvector<T> &arr, NRvector<U> &brr) //Sort an array arr[0..n - 1] into ascending order using Quicksort, while making the corresponding rearrangement of the array brr[0..n - 1]. 
	{ 
		const Int M = 7, NSTACK = 64; 
		Int i, ir, j, k, jstack = -1, l = 0, n = arr.size(); 
		T a; 
		U b; 
		VecInt istack(NSTACK); 
		ir = n - 1; 
		for (;;) 
		{     //Insertion sort when subarray small enough. 
			if (ir - l < M) 
			{ 
				for (j = l + 1; j <= ir; j++) 
				{ 
					a = arr[j]; 
					b = brr[j]; 
					for (i = j - 1; i >= l; i--) 
					{ 
						if (arr[i] <= a) break; 
						arr[i + 1] = arr[i]; 
						brr[i + 1] = brr[i]; 
					} 
					arr[i + 1] = a; 
					brr[i + 1] = b; 
				} 
				if (jstack < 0) break; 
				ir = istack[jstack--]; //Pop stack and begin a new round of partitioning.
				l = istack[jstack--]; 
			}
			else 
			{
				k = (l + ir) >> 1; //Choose median of left, center, and right elements as partitioning element a.Also rearrange so that a[l]a[l + 1]a[ir].
				SWAP(arr[k], arr[l + 1]); 
				SWAP(brr[k], brr[l + 1]); 
				if (arr[l] > arr[ir]) 
				{ 
					SWAP(arr[l], arr[ir]); 
					SWAP(brr[l], brr[ir]); 
				} 
				if (arr[l + 1] > arr[ir]) 
				{ 
					SWAP(arr[l + 1], arr[ir]); 
					SWAP(brr[l + 1], brr[ir]); 
				} 
				if (arr[l] > arr[l + 1]) 
				{ 
					SWAP(arr[l], arr[l + 1]); 
					SWAP(brr[l], brr[l + 1]); 
				} 
				i = l + 1; //Initialize pointers for partitioning.
				j = ir; 
				a = arr[l + 1]; //Partitioning element.
				b = brr[l + 1]; 
				for (;;) 
				{ //Beginning of innermost loop. 
					do i++; while (arr[i] < a); //Scan up to ?nd element > a. 
					do j--; while (arr[j] > a); //Scan down to ?nd element < a. 
					if (j < i) break; //Pointers crossed.Partitioning complete.
					SWAP(arr[i], arr[j]); //Exchange elements of both arrays.
					SWAP(brr[i], brr[j]); 
				} //End of innermost loop.
				arr[l + 1] = arr[j]; //Insert partitioning element in both arrays.
				arr[j] = a; 
				brr[l + 1] = brr[j]; 
				brr[j] = b; jstack += 2; //Push pointers to larger subarray on stack; process smaller subarray immediately. 
				if (jstack >= NSTACK) throw("NSTACK too small in sort2."); 
				if (ir - i + 1 >= j - l) 
				{ 
					istack[jstack] = ir; 
					istack[jstack - 1] = i; ir = j - 1; 
				}
				else 
				{ 
					istack[jstack] = j - 1; 
					istack[jstack - 1] = l; l = i;
				}
			}
		}
	}
	*/


};


class KrigeSystem 
{

public:

	//double  *r;
	//boost::scoped_array<double> r;
	//boost::scoped_array<double> rr;
	//boost::scoped_array<double> s;
	//boost::scoped_array<double> a;
	double *r;

	double	*rr;
	double	*s;
	double	*a;

	KrigeSystem(int lines) 
	{

		int temp = lines + 1 + 1;

		// TODO remove fortran +1 the rest is constraints, i think

		r = new double[temp];

		rr= new double[temp];
		s = new double[temp];

		a = new double[temp*temp];

		//r.reset(new double[temp]);
		//rr.reset(new double[temp]);
		//s.reset(new double[temp]);
		//a.reset(new double[temp*temp]);

	};
	
	////copy
	//KrigeSystem(const KrigeSystem& a)
	//{	
	//	// TODO remove fortran +1 the rest is constraints, i think

	//	r = &a.r;

	//	rr = &a.rr;
	//	s = &a.s;

	//	a = &a.a;
	//}

	////assignment
	//KrigeSystem & operator=(const KrigeSystem& a)
	//{	// TODO remove fortran +1 the rest is constraints, i think

	//	r = &a.r;

	//	rr = &a.rr;
	//	s = &a.s;

	//	a = &a.a;
	//}

	~KrigeSystem() 
	{

		delete[] a;
		delete[] s;
		delete[] r;
		delete[] rr;
	};

	int ksol(int neq)
	{

		double ak, tolerance, pivot, ap;
		int  kk, km1, lp, ll, ll1,  nm1, ij = 0, in, ijm,  m1, nn, nm, nsb, nright;
		
		// If there is only one equation then set ising and return:
		nright = 1;
		nsb = 1;
		if (neq < 1)
			return -2;

		if (neq == 1) {
			s[1] = r[1] / a[1];
			return 0;
		}

		// Initialize:
		tolerance = 0.1E-06;
		nn = neq*(neq + 1) / 2;// triangular number - ino?
		nm = nsb*neq;
		m1 = neq - 1;
		kk = 0;
		
		//Start triangulation:
		for (int k = 1; k <= m1; ++k) {
			kk = kk + k;
			ak = a[kk];
			if (fabs(ak) < tolerance)
				return k; // pivot problem

			km1 = k - 1;
			int ii;
			for (int iv = 1; iv <= nright; ++iv) {
				nm1 = nm*(iv - 1);
				ii = kk + nn*(iv - 1);
				pivot = 1. / a[ii];
				lp = 0;
				for (int i = k; i <= m1; ++i) {
					ll = ii;
					ii = ii + i;
					ap = a[ii] * pivot;
					lp = lp + 1;
					ij = ii - km1;
					for (int j = i; j <= m1; ++j) {
						ij = ij + j;
						ll = ll + j;
						a[ij] = a[ij] - ap*a[ll];
					}
					for (int llb = k; llb <= nm; llb = llb + neq) {
						in = llb + lp + nm1;
						ll1 = llb + nm1;
						r[in] = r[in] - ap*r[ll1];
					}
				}
			}
		}

		// Error checking - singular matrix:
		ijm = ij - nn*(nright - 1);
		if (fabs(a[ijm]) < tolerance) 
		{
			return neq;
		}
		
		// Finished triangulation, start solving back:
		//int i;
		
		for (int iv = 1; iv <= nright; ++iv) 
		{
			nm1 = nm*(iv - 1);
			ij = ijm + nn*(iv - 1);
			pivot = 1. / a[ij];
			for (int llb = neq; llb <= nm; llb = llb + neq) 
			{
				ll1 = llb + nm1;
				s[ll1] = r[ll1] * pivot;
			}
			int i = neq;
			kk = ij;
			for (int ii = 1; ii <= m1; ++ii) 
			{
				kk = kk - i;
				pivot = 1. / a[kk];
				--i ;
				for (int llb = i; llb <= nm; llb = llb + neq) 
				{
					ll1 = llb + nm1;
					in = ll1;
					ap = r[in];
					ij = kk;
					for (int j = i; j <= m1; ++j) 
					{
						ij = ij + j;
						++in ;
						ap = ap - a[ij] * s[in];
					}
					s[ll1] = ap*pivot;
				}
			}
		}
	
		
		//
		//double **aa, **b;
		//aa = new double *[neq+1];
		//b = new double *[neq+1];
		//for (int i =0; i < neq; i++)
		//{
		//	aa[i] = new double[neq];
		//}
		//for (int i = 0; i < neq; i++)
		//{
		//	b[i] = new double[1];
		//}
		//

		//int triang=0;
		//for (int i = 0; i < neq; i++) {
		//	for (int j = 0; j < i; j++) {
		//		triang++;
		//		aa[i][j] = a[triang];
		//	}
		//}
		//	for (int i = 1; i <= neq; i++)
		//{
		//	b[0][i-1] =r[i];
		//}
	

		//int i, icol, irow, j, k, l, ll;
		//int n = neq, m = 1;
		//double big, dum, pivinv; 
		////VecInt indxc(n), indxr(n), ipiv(n);

		//int *indxc, *indxr, *ipiv;

		//indxc=new int[n];
		//indxr = new int[n];
		//ipiv = new int[n];
		//double buff;

		////These integer arrays are used for bookkeeping on the pivoting.
		//	for (j = 0; j<n; j++) 
		//		ipiv[j] = 0; 
		//	for (i = 0; i<n; i++) 
		//	{
		//	//This is the main loop over the columns to be reduced.
		//		big = 0.0; 
		//		for (j = 0; j<n; j++) 
		//			//This is the outer loop of the search for a pivot element.
		//			if (ipiv[j] != 1) 
		//				for (k = 0; k<n; k++) 
		//				{ 
		//					if (ipiv[k] == 0) 
		//					{ 
		//						if (abs(aa[j][k]) >= big) 
		//						{ 
		//							big = abs(aa[j][k]); 
		//							irow = j; 
		//							icol = k; 
		//						} 
		//					} 
		//				} 
		//		++(ipiv[icol]);
		//		if (irow != icol) 
		//		{ 
		//			for (l = 0; l < n; l++) 
		//			{
		//				buff = aa[icol][l];
		//				aa[icol][l] = aa[irow][l];
		//				aa[irow][l] = buff;
		//				//SWAP(aa[irow][l], aa[icol][l]);
		//			}
		//			for (l = 0; l<m; l++) 
		//			{
		//				buff = b[icol][l];
		//				b[icol][l] = b[irow][l];
		//				b[irow][l] = buff;
		//				//SWAP(b[irow][l], b[icol][l]);
		//			}
		//		} 
		//		indxr[i] = irow;
		//		indxc[i] = icol; 
		//		if (aa[icol][icol] == 0.0) 
		//			throw("gaussj: Singular Matrix");

		//		pivinv = 1.0 / aa[icol][icol]; 
		//		aa[icol][icol] = 1.0; 
		//		for (l = 0; l<n; l++) 
		//			aa[icol][l] *= pivinv; 
		//		for (l = 0; l<m; l++) 
		//			b[icol][l] *= pivinv; 
		//		for (ll = 0; ll<n; ll++) 
		//			//Next, we reduce the rows... 
		//			if (ll != icol) 
		//			{ 
		//				//...except for the pivot one, of course.
		//				dum = aa[ll][icol]; 
		//				aa[ll][icol] = 0.0; 
		//				for (l = 0; l<n; l++) 
		//					aa[ll][l] -= aa[icol][l] * dum; 
		//				for (l = 0; l<m; l++) 
		//					b[ll][l] -= b[icol][l] * dum; 
		//			}
		//	}
		//	for (l = n-1; l >= 0; l--)
		//	{
		//		if (indxr[l] != indxc[l])
		//			for (k = 0; k < n; k++) 
		//			{
		//				buff = aa[k][indxc[l]];
		//				aa[k][indxc[l]] = aa[k][indxr[l]];
		//				aa[k][indxr[l]] = buff;

		//				//SWAP(aa[k][indxr[l]], aa[k][indxc[l]]);

		//			}
		//	}
		//for (int i = 1; i <= neq; i++)
		//{
		//	s[i] = b[0][i - 1];
		//}
		//		
		//	//And we are done. b has the result

		//for (int i = 0; i < neq; i++)
		//{
		//	delete [] aa[i];
		//	delete [] b[i];
		//}

		//delete [] b;
		//delete [] aa;

		//delete [] indxc;
		//delete [] indxr;
		//delete [] ipiv;



		//for (int i = 1; i <= neq; i++) 
		//{
		//	float si = 0;
		//	for (int j = 1; j < i; j++) 
		//	{
		//		si = si + a[(((i - 1)*i) / 2) + j] * s[j];
		//	}
		//	s[i] = (r[i] - si) / a[(((i - 1)*i) / 2) + i];
		//}




		// Finished solving back, return:
		return 0;
	};


};

//
//class KrigeSystem_matrx
//{
//
//
//public:
//
//	Eigen::VectorXd r;
//	Eigen::VectorXd rr;
//	Eigen::VectorXd s;
//	Eigen::MatrixXd a;
//
//
//
//	KrigeSystem_matrx(int lines) 
//	{
//
//		r= Eigen::VectorXd(lines);
//		rr= Eigen::VectorXd(lines);
//		s= Eigen::VectorXd(lines);
//		a= Eigen::MatrixXd(lines,lines);
//
//	};
//
//
//	int ksol()
//	{
//
//		s = a.ldlt().solve(r);
//	
//		//s = a.triangularView<Eigen::Upper>().solve(r);
//		
//
//		return 0;
//
//	};
//
//
//
//};
//
//

#endif	

