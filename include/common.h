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
#define STENCIL_SIZE 2 // num boundary ghosts
#define ZL(i) ((i)+STENCIL_SIZE-1)
#define ZR(i) ((i)+STENCIL_SIZE)

// boundary types
#define BC_WALL 0
#define BC_PERIODIC 1
#define BC_FREE 2

#define TWO_PI 6.28318530718

//#define flatind()

// real type
typedef double real;
typedef real rvec[3]; 
typedef int dvec[3];

typedef struct {
	real rho; // density
	rvec mom; // momentum
	rvec com; // the first moment in position 
	real etot; // total energy
#ifdef RADIATION
#endif
} hydro_vector;

typedef struct {
	real rho_l, rho_c, rho_r;
	real etot; // total energy
	rvec v; // velocity
	real e_l, e_c, e_r; // specific internal energy
	real p_l, p_c, p_r; // pressure
	real c_l, c_c, c_r; // sound speed
#ifdef RADIATION
#endif
} hydro_derived_state; 



#define RAY_BASE_BITS(x) ((x)&((1<<hp->dim)-1))

#define CLIGHT 1000.0
typedef struct {
	real rmin, rmax; // the inner and outer radius for the timestep 
	real thmin, thmax; // upper and lower theta for this ray 
	dvec orcell; // the originating cell index
	real N; // number of photons in the wavepacket
} hydro_ray;

typedef struct {
	real mass; // total mass 
	rvec mom; // total momentum
	rvec com; // the first moment in position 
	real etot; // total energy
	int verts[4];
} hydro_lagrange_tet;

typedef struct {
	rvec pos, vel;
	real rho, etot;
} hydro_lagrange_vert;

// the grid 
typedef struct {

	// current and stopping time
	char* name;
	real time, dt;
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

	// the callback function pointer
	void (*output_callback)(void*);

	// experimental! Radiation field!
	hydro_ray *rays;
	int nrays, base_ray_lvl, max_ray_lvl;

	// experimental! Lagrangian mesh.
	hydro_lagrange_vert* mesh_verts;
	hydro_lagrange_tet* mesh_tets;
	int nverts, ntets;
	
	
} hydro_problem;








#endif // COMMON_H_
