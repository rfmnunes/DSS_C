#include "distribution.h"



	Distribution::Distribution()
	{
		this->orig_distr = new std::vector<double>;
		this->gauss_distr = new std::vector<double>;

	};

	Distribution::Distribution(std::vector<double> &distr, std::vector<double> &gauss_distr_in, unsigned int n_data_in, double z_min, double z_max) // check if scope is correct
	{
		this->orig_distr = new std::vector<double>;
		this->gauss_distr = new std::vector<double>;

		*orig_distr = distr;
		*gauss_distr = gauss_distr_in;
		n_data = n_data_in;
		zmin = z_min;
		zmax = z_max;
	};

	Distribution::Distribution(const Distribution &in_distr)
	{
		orig_distr = new std::vector<double>(*in_distr.orig_distr);
		gauss_distr = new std::vector<double>(*in_distr.gauss_distr);
		//orig_distr = in_distr.orig_distr;
		//gauss_distr = in_distr.gauss_distr;
		n_data = in_distr.n_data;
		zmin = in_distr.zmin;
		zmax = in_distr.zmax;
	};

	Distribution&  Distribution::operator=(const Distribution& in_distr)
	{
		orig_distr = in_distr.orig_distr;
		gauss_distr = in_distr.gauss_distr;
		n_data = in_distr.n_data;
		zmin = in_distr.zmin;
		zmax = in_distr.zmax;
		return *this;
	};

	Distribution::~Distribution()
	{
		//orig_distr.clear();
		//gauss_distr.clear();
		delete orig_distr;
		delete gauss_distr;

	};

	template <class SomeType>
	/*static */int Distribution::locate_in_array_c(SomeType &in_array, int size, double &value)
	{
		// find the index of the closest value smaller than value in array
		// this assumes the array starts in 0
		if (value <= in_array[0]) return 0;
		if (value > in_array[size - 1]) return size - 1;

		int low_boundary = 0;
		int high_boundary = size - 1;

		int mid_index;

		while ((high_boundary - low_boundary) > 1)
		{
			// get the middle value
			mid_index = low_boundary + (high_boundary - low_boundary) / 2;

			if (value > in_array[mid_index])// we change the lower bondary
				low_boundary = mid_index;
			else // change the high_boundary
				high_boundary = mid_index;
		}
		// when cycle ends, the lower boundary is the index we want
		return low_boundary;
	}

	template <class SomeType>
	/*static */int Distribution::locate_in_array_f(SomeType &in_array, int size, double &value)
	{
		// find the index of the closest value smaller than value in array
		// this assumes the array starts in 0
		if (value <= in_array[1]) return 1;
		if (value > in_array[size]) return size;

		int low_boundary = 1;
		int high_boundary = size;

		int mid_index;

		while ((high_boundary - low_boundary) > 1)
		{
			// get the middle value
			mid_index = low_boundary + (high_boundary - low_boundary) / 2;

			if (value > in_array[mid_index])// we change the lower bondary
				low_boundary = mid_index;
			else // change the high_boundary
				high_boundary = mid_index;
		}
		// when cycle ends, the lower boundary is the index we want
		return low_boundary;

	}


	template <class SomeType>
	/*static */int Distribution::locate_in_vector_f(std::vector<SomeType> &in_array, unsigned int size, double value)
	{
		// find the index of the closest value smaller than value in array
		// this assumes the array starts in 0
		if (value <= in_array[1]) return 1;
		if (value > in_array[size]) return size;

		int low_boundary = 1;
		int high_boundary = size;

		int mid_index;

		while ((high_boundary - low_boundary) > 1)
		{
			// get the middle value
			mid_index = low_boundary + (high_boundary - low_boundary) / 2;

			if (value > in_array[mid_index])// we change the lower bondary
				low_boundary = mid_index;
			else // change the high_boundary
				high_boundary = mid_index;
		}
		// when cycle ends, the lower boundary is the index we want
		return low_boundary;

	}


	//~Distribution()
	//{
	//	orig_distr.clear();
	//	gauss_distr.clear();
	//	//delete this;
	//} 

	double Distribution::get_limits(double &cmean)
		// return a gaussian value form the distribution, associated with cmean 
		// TODO it's in fortran index
	{
		int position = locate_in_vector_f(*orig_distr, n_data, cmean);

		if (cmean <= (*orig_distr)[1])
			return (*gauss_distr)[1];
		else if (cmean >= (*orig_distr)[n_data])
			return (*gauss_distr)[n_data];
		else
			return (*gauss_distr)[position] + (cmean - (*orig_distr)[position])*((*gauss_distr)[position + 1] - (*gauss_distr)[position]) / ((*orig_distr)[position + 1] - (*orig_distr)[position]);

		//else 
		//{
		//	for ( unsigned int j=1; j<n_data; ++j)
		//	{ //changed +1
		//		if (cmean >= orig_distr[j] && cmean < orig_distr[j+1])
		//			return gauss_distr[j]+(cmean-orig_distr[j])*(gauss_distr[j+1]-gauss_distr[j])/(orig_distr[j+1]-orig_distr[j]);  
		//	} 	 
		//}

		//return -999.99;
	};

	double Distribution::get_limits_c(double &cmean)
		// return a gaussian value form the distribution, associated with cmean 
	{ // when I have the distribution starting on 0, change to this.

		if (cmean <= (*orig_distr)[0])
			return (*gauss_distr)[0];
		else if (cmean >= (*orig_distr)[n_data - 1])
			return (*gauss_distr)[n_data - 1];
		else {
			for (unsigned int j = 0; j < n_data - 1; ++j) { //changed +1
				if (cmean >= (*orig_distr)[j] && cmean < (*orig_distr)[j + 1])
					return (*gauss_distr)[j] + (cmean - (*orig_distr)[j])*((*gauss_distr)[j + 1] - (*gauss_distr)[j]) / ((*orig_distr)[j + 1] - (*orig_distr)[j]);
			}
		}

		return -999.99;
	};

	double Distribution::backtr(double &simgval, int &ltail, double &ltpar, int &utail, double &utpar)

	{
		double backtra;
		int gsize = (*gauss_distr).size();

		// Value in the lower tail?    1=linear, 2=power, (3 and 4 are invalid)
		if (simgval <= (*gauss_distr)[1])
		{
			backtra = (*orig_distr)[1];
			double cdflo = gcum((*gauss_distr)[1]);
			double cdfbt = gcum(simgval);
			if (ltail == 1)
				backtra = powint((0.0), cdflo, zmin, (*orig_distr)[1], cdfbt, 1.0);
			else if (ltail == 2) {
				double cpow = (double)1.0 / ltpar;
				backtra = powint(0.0, cdflo, zmin, (*orig_distr)[1], cdfbt, cpow);
			}
		}

		// Value in the upper tail?     1=linear, 2=power, 4=hyperbolic:
		else if (simgval >= (*gauss_distr)[gsize - 1])
		{
			backtra = (*orig_distr)[gsize - 1];
			double cdfhi = gcum((*gauss_distr)[gsize - 1]);
			double cdfbt = gcum(simgval);
			if (utail == 1)
				backtra = powint(cdfhi, 1.0, (*orig_distr)[gsize - 1], zmax, cdfbt, 1.0);
			else if (utail == 2) {
				double cpow = (double)(1.0 / utpar);
				backtra = powint(cdfhi, 1.0, (*orig_distr)[gsize - 1], zmax, cdfbt, cpow);
			}
			else if (utail == 4) {
				double lambda = (double)((pow((*orig_distr)[gsize - 1], utpar))*(1.0 - gcum((*gauss_distr)[gsize - 1])));
				backtra = (double)(pow((lambda / (1.0 - gcum(simgval))), (1.0 / utpar)));
			}
		}

		// Value within the transformation table:
		else
		{
			//locate(&gauss_distr[0],n_data,1,n_data,simgval,&j);
			int loc_test = locate_in_array_f(gauss_distr[0], n_data, simgval);
			int j = std::max(std::min(((int)n_data - 1), loc_test), 1);
			backtra = powint((*gauss_distr)[j], (*gauss_distr)[j + 1], (*orig_distr)[j], (*orig_distr)[j + 1], simgval, 1.0);
		}

		return backtra;
	};

	double Distribution::backtr_trunc(double &simgval)
	{
		double backtra;

		// Value in the lower tail?    1=linear, 2=power, (3 and 4 are invalid)
		if (simgval <= (*gauss_distr)[1]) {

			return -999.25;
			//backtra = orig_distr[1];
		}

		// Value in the upper tail?     1=linear, 2=power, 4=hyperbolic:
		else if (simgval >= (*gauss_distr)[n_data]) {

			return -999.25;

			//backtra = orig_distr[n_data];
		}

		// Value within the transformation table:
		else {
			//locate(&gauss_distr[0],n_data,1,n_data,simgval,&j);
			int loc_test = locate_in_array_f(gauss_distr[0], n_data, simgval); //check this
			int j = std::max(std::min(((int)n_data - 1), loc_test), 1);
			backtra = powint((*gauss_distr)[j], (*gauss_distr)[j + 1], (*orig_distr)[j], (*orig_distr)[j + 1], simgval, 1.0);
		}

		return backtra;
	};

	double Distribution::backtr_c(double &simgval, int &ltail, double &ltpar, int &utail, double &utpar)

	{
		double backtra, cpow, cdfbt;

		// Value in the lower tail?    1=linear, 2=power, (3 and 4 are invalid)
		if (simgval <= (*gauss_distr)[0]) {
			backtra = (*orig_distr)[0];
			double cdflo = gcum((*gauss_distr)[0]);
			cdfbt = gcum(simgval);
			if (ltail == 1)
			{
				backtra = powint(0.0, cdflo, zmin, (*orig_distr)[0], cdfbt, 1.0);
			}
			else if (ltail == 2)
			{
				cpow = (double)1.0 / ltpar;
				backtra = powint(0.0, cdflo, zmin, (*orig_distr)[0], cdfbt, cpow);
			}
		}

		// Value in the upper tail?     1=linear, 2=power, 4=hyperbolic:
		else if (simgval >= (*gauss_distr)[n_data - 1]) {
			backtra = (*orig_distr)[n_data - 1];
			double cdfhi = gcum((*gauss_distr)[n_data - 1]);
			cdfbt = gcum(simgval);
			if (utail == 1)
				backtra = powint(cdfhi, 1.0, (*orig_distr)[n_data - 1], zmax, cdfbt, 1.0);
			else if (utail == 2) {
				cpow = (double)(1.0 / utpar);
				backtra = powint(cdfhi, 1.0, (*orig_distr)[n_data - 1], zmax, cdfbt, cpow);
			}
			else if (utail == 4) {
				double lambda = (double)((pow((*orig_distr)[n_data - 1], utpar))*(1.0 - gcum((*gauss_distr)[n_data - 1])));
				backtra = (double)(pow((lambda / (1.0 - gcum(simgval))), (1.0 / utpar)));
			}
		}

		// Value within the transformation table:
		else {
			//locate(&gauss_distr[0],n_data-1,1,n_data,simgval,&j);
			int loc_test = locate_in_array_c(gauss_distr[0], n_data, simgval);
			int j = std::max(std::min(((int)n_data - 1 - 1), loc_test), 0);
			backtra = powint((*gauss_distr)[j], (*gauss_distr)[j + 1], (*orig_distr)[j], (*orig_distr)[j + 1], simgval, (1.0));
		}

		return backtra;
	};

	double Distribution::gcum(double &x)
	{
		//Evaluate the standard normal cdf given a normal deviate x.  gcum is
		//the area under a unit normal curve to the left of x.  The results are
		//accurate only to about 5 decimal places.
		double e2, t, z, gcuma;

		z = x;
		if (z < 0.0)
			z = -z;
		t = (double)(1. / (1. + 0.2316419*z));
		gcuma = (double)(t*(0.31938153 + t * (-0.356563782 + t * (1.781477937 + t * (-1.821255978 + t * 1.330274429)))));
		e2 = (double)0.0;

		// 6 standard deviations out gets treated as infinity: */
		if (z <= (double)6.0)
			e2 = (double)(exp(-z * z / 2.)*0.3989422803);
		gcuma = (double)1.0 - e2 * gcuma;
		if (x >= (double)0.0)
			return gcuma;
		gcuma = (double)1.0 - gcuma;

		return gcuma;
	};

	double Distribution::powint(double xlow, double xhigh, double &ylow, double &yhigh, double &xval, double potencia)
	{
		static const double epslon = 1.0E-20;
		double value;

		if ((xhigh - xlow) < epslon)
			value = (yhigh + ylow) / (double)2.0;
		else
			value = ylow + (yhigh - ylow)*(double)pow(((xval - xlow) / (xhigh - xlow)), potencia);

		return value;
	};

	void Distribution::update_gaussian_distribution()
	{
		(*gauss_distr).clear();

		(*gauss_distr).push_back(-999.25); // to maintain fortran style

		for (int i = 1; i < (*orig_distr).size(); ++i)
		{
			(*gauss_distr).push_back(1.0);
		}

		double oldcp = 0.0;
		double cp = 0.0;

		double w;
		double vrg;

		for (int j = 1; j < (*orig_distr).size(); ++j)
		{
			cp = cp + (double)((*gauss_distr)[j] / (*gauss_distr).size());
			w = (cp + oldcp)*0.5;

			vrg = gauinv(w);

			oldcp = cp;

			// Now, reset the weight to the normal scores value:
			(*gauss_distr)[j] = vrg;

		}



	}


	//template <class SomeType>
	//static int locate_in_array_c(SomeType *in_array, int size, double &value )
	//{
	//	// find the index of the closest value smaller than value in array
	//	// this assumes the array starts in 0
	//	if (value <= in_array[0]) return 0;
	//	if (value > in_array[size-1]) return size-1;
	//	int low_boundary = 0;
	//	int high_boundary =  size-1;
	//	int mid_index;
	//	while((high_boundary - low_boundary) > 1)
	//	{
	//		// get the middle value
	//		mid_index= low_boundary +  (high_boundary - low_boundary)/2;
	//		if (value > in_array[mid_index])// we change the lower bondary
	//			low_boundary = mid_index;
	//		else // change the high_boundary
	//			high_boundary = mid_index;
	//	}
	//	// when cycle ends, the lower boundary is the index we want
	//	return low_boundary;
	//}

	//template <class SomeType>
	//static int locate_in_array_f(SomeType *in_array, int size, double &value )
	//{
	//	// find the index of the closest value smaller than value in array
	//	// this assumes the array starts in 0
	//	if (value <= in_array[1]) return 1;
	//	if (value > in_array[size]) return size;
	//	int low_boundary = 1;
	//	int high_boundary =  size;
	//	int mid_index;
	//	while((high_boundary - low_boundary) > 1)
	//	{
	//		// get the middle value
	//		mid_index= low_boundary +  (high_boundary - low_boundary)/2;
	//		if (value > in_array[mid_index])// we change the lower bondary
	//			low_boundary = mid_index;
	//		else // change the high_boundary
	//			high_boundary = mid_index;
	//	}
	//	// when cycle ends, the lower boundary is the index we want
	//	return low_boundary;
	//}

	//int locate_in_vector(std::vector<double> const&in_vector, double *value )
	//{
	//	// find the index of the closest value smaller than value in array
	//	// this assumes the array starts in 0
	//	if (value <= &in_vector[0]) 
	//		return 0;
	//	if (value > in_vector[in_vector.size-1]) 
	//		return in_vector.size-1;
	//	int low_boundary = 0;
	//	int high_boundary =  in_vector.size-1;
	//	int mid_index;
	//	while((high_boundary - low_boundary) > 1)
	//	{
	//		// get the middle value
	//		mid_index= (high_boundary - low_boundary)/2;
	//		if (value > &in_vector[mid_index])// we change the lower bondary
	//			low_boundary = mid_index;
	//		else // change the high_boundary
	//			high_boundary = mid_index;
	//	}
	//	// when cycle ends, the lower boundary is the index we want
	//	return low_boundary;
	//}

	//int recursive_locate (double *xx, int n, int is, int ie, double x)
	//{
	//	//
	//	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	//	//
	//	//    OBJECTIVES:
	//	//      Given an array "xx" of length "n", and given a value "x", this routine
	//	//      returns a value "j" such that "x" is between xx(j) and xx(j+1).  xx
	//	//      must be monotonic, either increasing or decreasing.  j=is-1 or j=ie is
	//	//      returned to indicate that x is out of range.
	//	//
	//	//      Bisection Concept From "Numerical Recipes", Press et. al. 1986  pp 90.
	//	//
	//	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - \																																																															- - - - - - - -
	//	//
	//	int jl,ju;
	//	//int jm;
	//	// Initialize lower and upper methods
	//	if (is <= 0) is = 1;
	//	jl = is-1;
	//	ju = ie;
	//	if (*(xx+n) <= x) {
	//		//*j = ie;
	//		//return j;
	//	}
	//			
	//	// If we are not done then compute a midpoint
	//	L: if (ju-jl > 1) {
	//		jm = (ju+jl)/2;
	//		// Replace the lower or upper limit with the midpoint
	//		if(*(xx+ie) > *(xx+is) && x > *(xx+jm))
	//			jl = jm;
	//		else
	//			ju = jm;
	//		goto L;
	//	}
	//	
	//	if (ju-jl > 1) // we must still find the point
	//	{ 
	//	// Return with the array index
	//	//*j = jl;
	//	//return j;
	//	
	//
	//};

	//int gauinv (double p, double *xp)
	//{
	//	
	//	// Check for an error situation:
	//	if (p < lim) {
	//		*xp = -1.0E10;
	//		return 0;
	//	}
	//	if (p > (1.0-lim)) {
	//		*xp =  1.0E10;
	//		return 0;
	//	}
	//	// Get k for an error situation:
	//	double pp   = p;
	//	if (p > 0.5) pp = 1 - pp;
	//	*xp   = (double)0.0;
	//	if (p == (double)0.5) return 1;
	//	// Approximate the function:
	//	double y  = sqrt(log(1.0/(pp*pp)));
	//	*xp = (double) (y + ((((y*p4+p3)*y+p2)*y+p1)*y+p0) / ((((y*q4+q3)*y+q2)*y+q1)*y+q0));
	//	if ((double)p == (double)pp) 
	//		*xp = -*xp;
	//	// Return with G^-1(p):
	//	return 1;
	//};

	/*static */double Distribution::gauinv(double &p)
	{
		static const double lim = 1.0E-10;
		static const double p0 = -0.322232431088;
		static const double p1 = -1.0;
		static const double p2 = -0.342242088547;
		static const double p3 = -0.0204231210245;
		static const double p4 = -0.0000453642210148;

		static const double q0 = 0.0993484626060;
		static const double q1 = 0.588581570495;
		static const double q2 = 0.531103462366;
		static const double q3 = 0.103537752850;
		static const double q4 = 0.0038560700634;

		double xp;

		// Check for an error situation:

		if (p < lim) {
			xp = -1.0E10;
			return xp;
		}
		if (p > (1.0 - lim)) {
			xp = 1.0E10;
			return xp;
		}

		// Get k for an error situation:
		double pp = p;
		if (p > 0.5) pp = 1 - pp;
		xp = (double)0.0;
		if (p == (double)0.5)
			return 0;

		// Approximate the function:
		double y = sqrt(log(1.0 / (pp*pp)));
		xp = (double)(y + ((((y*p4 + p3)*y + p2)*y + p1)*y + p0) / ((((y*q4 + q3)*y + q2)*y + q1)*y + q0));
		if (p == pp)
			xp = -xp;

		// Return with G^-1(p):
		return xp;
	};

	/*static	*/double Distribution::RationalApproximation(double t)
	{
		// Abramowitz and Stegun formula 26.2.23.
		// The absolute value of the error should be less than 4.5 e-4.
		double c[] = { 2.515517, 0.802853, 0.010328 };
		double d[] = { 1.432788, 0.189269, 0.001308 };
		return t - ((c[2] * t + c[1])*t + c[0]) /
			(((d[2] * t + d[1])*t + d[0])*t + 1.0);
	}

	/*static */double Distribution::NormalCDFInverse(double &p)
	{
		if (p <= 0.0 || p >= 1.0)
		{
			std::stringstream os;
			os << "Invalid input argument (" << p
				<< "); must be larger than 0 but less than 1.";
			throw std::invalid_argument(os.str());
		}

		// See article above for explanation of this section.
		if (p < 0.5)
		{
			// F^-1(p) = - G^-1(p)
			return -RationalApproximation(sqrt(-2.0*log(p)));
		}
		else
		{
			// F^-1(p) = G^-1(1-p)
			return RationalApproximation((sqrt(-2.0*log(1 - p))));
		}
	}

	void Distribution::weighted_add_to_distribution(misc utils, Distribution to_add, float weight)
	{
		// this adds a distribution to this object, assigning a weight to it.
		// weight is a fraction of 1. original distribution is weighted as 1-weight.
		float c[1];


		int dist1_npoints = n_data;
		int dist2_npoints = to_add.n_data;
		int distfinal_npoints = dist1_npoints + dist2_npoints;

		(*orig_distr).insert((*orig_distr).end(), to_add.orig_distr->begin() + 1, to_add.orig_distr->end());// +1 because fortran beginning

		int totalweight = distfinal_npoints;

		(*gauss_distr).clear();
		//gauss_distr.insert(gauss_distr.end(),-999.25);

		for (int j = 0; j < dist1_npoints; ++j)
			(*gauss_distr).insert((*gauss_distr).end(), ((1 - weight)*totalweight) / dist1_npoints);

		for (int j = dist1_npoints; j < dist1_npoints + dist2_npoints; ++j)
			(*gauss_distr).insert((*gauss_distr).end(), ((weight)*totalweight) / dist2_npoints);

		utils.sortit(0, distfinal_npoints - 1, &(*orig_distr)[0], 1, &(*gauss_distr)[0], c, c, c, c, c, c);

		double oldcp = 0.0;
		double cp = 0.0;
		double w;
		for (int j = 0; j < distfinal_npoints; ++j) {
			cp += (double)((*gauss_distr)[j] / totalweight);
			w = (cp + oldcp)*0.5;

			(*gauss_distr)[j] = gauinv(w);
			oldcp = cp;
		}
		n_data = distfinal_npoints;
		zmin = (*orig_distr)[0];
		zmax = (*orig_distr)[n_data - 1];

		return;

	}

	double Distribution::get_distribution_average()
	{
		double average = 0;

		for (unsigned int i = 1; i <= (*orig_distr).size() - 1; ++i)
		{
			average += (*orig_distr)[i];
		}


		return average / ((*orig_distr).size() - 1);
	};

	double Distribution::get_distribution_variance()
	{
		double variance = 0;
		double average = get_distribution_average();

		for (unsigned int i = 1; i <= (*orig_distr).size() - 1; ++i)
		{
			variance += pow((double)((*orig_distr)[i] - average), 2.0);
		}


		return variance / std::max((int)(((*orig_distr).size() - 1) - 1), 1);
	};

	double Distribution::get_distribution_median()
	{

		if (((*orig_distr).size()) % 2 == 0)
			return (*orig_distr)[(*orig_distr).size() / 2];
		else
			return ((*orig_distr)[((*orig_distr).size() - 1) / 2] + (*orig_distr)[((*orig_distr).size() + 1) / 2]) / 2;

		//for (unsigned int i = 1; i <= (*orig_distr).size() - 1; ++i)
		//{
		//	average += (*orig_distr)[i];
		//}

	};

	//double simulate_value_in_distribution(double cmean, double cstdev,HarddataPars harddata, Cl_mtrng RNG, int ntry)
	//{
	//	
	//	double p;
	//	double xp;
	//	float bias_corr_avg = 0.0;
	//	float bias_corr_var = 0.0;
	//	double simval;
	//	float gauss_krig_mean = (float)get_limits(cmean);
	//	int k;
	//	
	//	while ( k < ntry)
	//	{
	//		// Simulate a random gaussian value
	//		p = RNG.random_double();
	//		xp = Distribution::gauinv(p);
	//		//xp = Distribution::NormalCDFInverse(p);
	//		xp = xp * cstdev + gauss_krig_mean;
	//		simval = backtr(xp
	//			, harddata.ltail, harddata.ltpar
	//			, harddata.utail, harddata.utpar);
	//		bias_corr_avg += simval;
	//		++k;
	//	}
	//	bias_corr_avg = bias_corr_avg / ntry;
	//	//welford
	//	//bias_corr_var = S/(k);
	//	// Reter o valor simulado (SDSIM)
	//	//newsim = simval;
	//	if (ntry > 1)
	//	{
	//		simval = simval + (float)(cmean - bias_corr_avg);
	//	}
	//	return simval;
	//
	//}



