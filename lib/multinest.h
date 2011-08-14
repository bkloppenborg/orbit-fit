// C header file for MultiNest

#ifdef __INTEL_COMPILER 			// if the MultiNest library was compiled with ifort
       #define NESTRUN nested_mp_nestrun_
#elif defined __GNUC__ 				// if the MultiNest library was compiled with gfortran
       #define NESTRUN __nested_MOD_nestrun
#else
       #error Don't know how to link to Fortran libraries, check symbol table for your platform (nm libnest3.a | grep nestrun) & edit example_eggbox_C/eggbox.c
#endif


/***************************************** C Interface to MultiNest **************************************************/

extern void NESTRUN(int *, int *, int *, double *, 
    double *, int *, int *, int *, 
    int *, int *, double *,  char *, 
    int *, int *, int *, int *, 
    int *, int *, double *, 
    void (*Loglike)(double *, int *, int *, double *), 
    void (*dumper)(int *, int *, int *, double **, double **, double *, double *, double *), 
    int *context);
