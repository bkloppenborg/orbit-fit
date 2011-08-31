// C header file for MultiNest

#ifdef __INTEL_COMPILER 			// if the MultiNest library was compiled with ifort
       #define NESTRUN nested_mp_nestrun_
#elif defined __GNUC__ 				// if the MultiNest library was compiled with gfortran
       #define NESTRUN __nested_MOD_nestrun
#else
       #error "Don't know how to link to Fortran libraries, check symbol table for your platform (nm libnest3.a | grep nestrun) & edit example_eggbox_C/eggbox.c"
#endif

#ifndef MULTINEST_H
#define MULTINEST_H

/***************************************** C/C++ Interface to MultiNest *****************************************/

#ifdef  __cplusplus

template<typename type, int ndims> class farray_traits;

template<> class farray_traits<double, 1> { public: static const int id = 537; };
template<> class farray_traits<double, 2> { public: static const int id = 538; };
template<> class farray_traits<int, 1> { public: static const int id = 265; };
template<> class farray_traits<int, 2> { public: static const int id = 266; };

// the extra data for f90 that defines how arrays are arranged.
template<typename T, int ndim> class farray
{
	public:
		farray(T *_data, int w, int h = 0) : data(_data), offset(0), type(farray_traits<T, ndim>::id), 
		x_stride(1), x_lbound(1), x_ubound(w), y_stride(w), y_lbound(1), y_ubound(h) {};
		
		T *data;
		int offset;
		int type;
		int x_stride, x_lbound, x_ubound;
		int y_stride, y_lbound, y_ubound;
};

extern "C" {
#endif

void NESTRUN(int *, int *, int *, double *, 
    double *, int *, int *, int *, 
    int *, int *, double *,  char *, 
    int *, int *, int *, int *, 
    int *, int *, double *, 
    void (*Loglike)(double *, int *, int *, double *), 
    void (*dumper)(int &, int &, int &, double **, double **, double *, double &, double &, double &),
    int *context);
    
#ifdef  __cplusplus
}
#endif
    
#endif // MULTINEST_H
