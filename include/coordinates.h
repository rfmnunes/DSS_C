#ifndef COORDS_INCLUDED
#define COORDS_INCLUDED

template <class T> class cl_coord 
{
public:
	
	//unsigned int dimensions;
	T x;
	T y;
	T z;

	cl_coord(T x_in, T y_in, T z_in)
	{
		x= x_in; y= y_in; z = z_in;
	};
		
	cl_coord(const cl_coord<T>& in)
	{
		x= in.x; y= in.y; z = in.z;
	};

	// default constructor
	cl_coord()
	{
		x= 0; y= 0; z = 0;
	};

	~cl_coord()
	{
	
	}

	cl_coord operator+(const cl_coord& b)
	{         
		cl_coord coord;
        coord.x = this->x + b.x;
        coord.y = this->y + b.y;
        coord.z = this->z + b.z;
        return coord;
	};
	
	cl_coord operator-(const cl_coord& b)
	{         
		cl_coord coord;
        coord.x = this->x - b.x;
        coord.y = this->y - b.y;
        coord.z = this->z - b.z;
        return coord;
	};
	
	bool operator==(const cl_coord& b)
	{         
		if (this->x == b.x && this->y == b.y && this->z == b.z)
			return true; 
        else return false;
	};

};



//TO DO SOMETIME

template <class T> class CoordinatesClass
{
public: //perhaps a vector is better

	std::vector<T> coord;

	// CONSTRUCTORS

	// default constructor
	CoordinatesClass(unsigned int n_dims)
	{
		for (unsigned int i = 0; i < n_dims; ++i)
			coord.push_back(0);
	};
	
	// construct from array
	CoordinatesClass(unsigned int n_dims, T *coord_in)
	{
		for (unsigned int i = 0; i < n_dims; ++i)
			coord.push_back(coord_in[i]);
	};
	
	// construct from vector
	CoordinatesClass( std::vector<T> *coord_in)
	{
		for (unsigned int i = 0; i < coord_in->size(); ++i)
			coord.push_back(coord_in->at(i));
	};


	// copy constructor
	CoordinatesClass(CoordinatesClass<T> *_in)
	{
		coord= _in.coord;
	};

	// copy with cast cnstructor
	template<typename otherT>
	CoordinatesClass(const CoordinatesClass<otherT> &_in)
	{
		for (unsigned int i = 0; i < _in.coord.size(); ++i)
			coord.push_back( (T) _in.coord.at(i));
	};


	//// construct 2d coordinate - TO GENERALIZE
	//CoordinatesClass(T _in_x,T _in_y)
	//{
	//	n_dimensions = 2;
	//	coord = new T[n_dimensions];
	//	coord[0]= _in_x;
	//	coord[1]= _in_y;
	//};

	//// construct 3d coordinate - TO GENERALIZE
	//CoordinatesClass(T _in_x,T _in_y,T _in_z)
	//{
	//	n_dimensions = 3;
	//	coord = new T[n_dimensions];
	//	coord[0]= _in_x;
	//	coord[1]= _in_y;
	//	coord[2]= _in_z;
	//};


	//// construct 4d coordinate - TO GENERALIZE
	//CoordinatesClass(T _in_x,T _in_y,T _in_z, T _in_t)
	//{
	//	n_dimensions = 4;
	//	coord = new T[n_dimensions];
	//	coord[0]= _in_x;
	//	coord[1]= _in_y;
	//	coord[2]= _in_z;
	//	coord[3]= _in_t;
	//};



	// OPERATORS

	CoordinatesClass<T> operator+(const CoordinatesClass<T>& b)
	{         
		CoordinatesClass<T> coord_out((unsigned int)b.coord.size());
		for (unsigned int i =0 ; i < (unsigned int)b.coord.size(); ++i)
			coord_out.coord[i]= this->coord[i] + b.coord[i];
		return coord_out;
	};
	
	CoordinatesClass<T> operator-(const CoordinatesClass<T>& b)
	{         
		CoordinatesClass<T> coord_out((unsigned int)b.coord.size());
		for (unsigned int i =0 ; i <(unsigned int)b.coord.size(); ++i)
			coord_out.coord[i]= this->coord[i] - b.coord[i];
		return coord_out;
	};
	
	CoordinatesClass<T> operator/(const CoordinatesClass<T>& b)
	{         
		CoordinatesClass<T> coord_out((unsigned int)b.coord.size());
		for (unsigned int i =0 ; i < b.coord.size(); ++i)
			coord_out.coord[i]= this->coord[i] / b.coord[i];
		return coord_out;
	};

	CoordinatesClass<T> operator*(const CoordinatesClass<T>& b)
	{
		CoordinatesClass<T> coord_out((unsigned int)b.coord.size());
		for (unsigned int i = 0; i < b.coord.size(); ++i)
			coord_out.coord[i] = this->coord[i] * b.coord[i];
		return coord_out;
	};

	template< typename otherT>
	CoordinatesClass<T>& operator=(const CoordinatesClass<otherT>& b) 
	{
		for (int i = 0; i < b.coord.size(); ++i)
			coord[i] = b.coord[i];
		return *this;
	};
	
	bool operator==(const CoordinatesClass& b)
	{         
		if (this.coord == b.coord)
			return true;
		else
			return false;
	};

	//template<typename FromT>
	//Pos2<T>& operator=(const Pos2<FromT>& from) {
	//	x = from.x;
	//	y = from.y;
	//	return *this;
	//}

	//template <typename Lo, typename Hi> Lo down_cast (const Hi& value)
	//{
 //   Lo result = Lo (a);
 //   if  (Hi (result) != value)
 //       throw exception ("loss of data in down cast");
 //   return result;
	//};

	//template<typename InT, typename OutT>
	//CoordinatesClass<OutT> convert( const CoordinatesClass<InT>& op )
	//{
	//	return CoordinatesClass<OutT>( op);
	//}
};


#endif
