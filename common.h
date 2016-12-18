/*************************************************************
 *
 *		common.h
 *
 *		Devon Powell
 *		21 April 2015
 *
 *		Approximate Riemann solvers for d_euler_hydro
 *
 *************************************************************/

#ifndef COMMON_H_
#define COMMON_H_

// macros for named acces to state and flux arrays
#define RHO 0
#define ENERGY 1
#define MOMX 2 
#define MOMY 3 
#define MOMZ 4 

// macros for getting zone indices from the face index
#define STENCIL_SIZE 1 // num boundary ghosts
#define ZL(i) ((i)+STENCIL_SIZE-1)
#define ZR(i) ((i)+STENCIL_SIZE)

// boundary types
#define BC_WALL 0
#define BC_PERIODIC 1
#define BC_FREE 2

//#define flatind()

// real type
typedef double real;
typedef real rvec[3]; 
typedef int dvec[3];

typedef struct {
	real rho; // density
	rvec mom; // momentum
	real etot; // total energy
#ifdef RADIATION
#endif
} hydro_vector;

typedef struct {
	real rho; //density
	real etot; // total energy
	rvec v; // velocity
	real e; // specific internal energy
	real p; // pressure
	real c; // sound speed
#ifdef RADIATION
#endif
} hydro_derived_state; 

// the grid 
typedef struct {

	// current and stopping time
	char* name;
	real time;
	int step;

	// equation of state functions
	real eos_gamma;
	real (*eos_p)(real, real, real);
	real cfl_fac;

	// the grid itself
	real dx; // always cubical grid cells
	int bctype; // boundary condition type
	hydro_vector *grid; // the conserved state vector 

	// indexing information 
	int dim; // the problem dimension, 1, 2, or 3
	dvec nx; // the grid dimensions 
	dvec strides; // the grid strides in each dimension
	
} hydro_problem;








#endif // COMMON_H_
