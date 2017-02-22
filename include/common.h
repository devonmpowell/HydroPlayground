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
	// real x; // the ionization fraction
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
	long angle_id;
	real rmin, rmax; // the inner and outer radius for the timestep 
	dvec orcell; // the originating cell index
	real Ftot, Fcom;	
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

typedef struct {
	real E;
} rad_vector;

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
	rad_vector *rad_grid; // the conserved state vector 

	// experimental! Lagrangian mesh.
	hydro_lagrange_vert* mesh_verts;
	hydro_lagrange_tet* mesh_tets;
	int nverts, ntets;
	
	
} hydro_problem;




const static int num_base_3d = 48;
const static int base_bits_3d[48] = {0x0, 0x1000000, 0x2000000, 0x3000000, 0x4000000, 0x5000000, 0x10000000, 0x11000000, 0x12000000, 0x13000000, 0x14000000, 0x15000000, 0x20000000, 0x21000000, 0x22000000, 0x23000000, 0x24000000, 0x25000000, 0x30000000, 0x31000000, 0x32000000, 0x33000000, 0x34000000, 0x35000000, 0x40000000, 0x41000000, 0x42000000, 0x43000000, 0x44000000, 0x45000000, 0x50000000, 0x51000000, 0x52000000, 0x53000000, 0x54000000, 0x55000000, 0x60000000, 0x61000000, 0x62000000, 0x63000000, 0x64000000, 0x65000000, 0x70000000, 0x71000000, 0x72000000, 0x73000000, 0x74000000, 0x75000000};
const static real base_beams_3d[48][3][3] = {{{7.071067812e-01,7.071067812e-01,0.000000000e+00},{5.773502692e-01,5.773502692e-01,5.773502692e-01},{1.000000000e+00,0.000000000e+00,0.000000000e+00}},
{{7.071067812e-01,7.071067812e-01,0.000000000e+00},{0.000000000e+00,1.000000000e+00,0.000000000e+00},{5.773502692e-01,5.773502692e-01,5.773502692e-01}},
{{0.000000000e+00,7.071067812e-01,7.071067812e-01},{5.773502692e-01,5.773502692e-01,5.773502692e-01},{0.000000000e+00,1.000000000e+00,0.000000000e+00}},
{{0.000000000e+00,7.071067812e-01,7.071067812e-01},{0.000000000e+00,0.000000000e+00,1.000000000e+00},{5.773502692e-01,5.773502692e-01,5.773502692e-01}},
{{7.071067812e-01,0.000000000e+00,7.071067812e-01},{5.773502692e-01,5.773502692e-01,5.773502692e-01},{0.000000000e+00,0.000000000e+00,1.000000000e+00}},
{{7.071067812e-01,0.000000000e+00,7.071067812e-01},{1.000000000e+00,0.000000000e+00,0.000000000e+00},{5.773502692e-01,5.773502692e-01,5.773502692e-01}},
{{-7.071067812e-01,7.071067812e-01,0.000000000e+00},{-1.000000000e+00,-0.000000000e+00,-0.000000000e+00},{-5.773502692e-01,5.773502692e-01,5.773502692e-01}},
{{-7.071067812e-01,7.071067812e-01,0.000000000e+00},{-5.773502692e-01,5.773502692e-01,5.773502692e-01},{0.000000000e+00,1.000000000e+00,0.000000000e+00}},
{{0.000000000e+00,7.071067812e-01,7.071067812e-01},{0.000000000e+00,1.000000000e+00,0.000000000e+00},{-5.773502692e-01,5.773502692e-01,5.773502692e-01}},
{{0.000000000e+00,7.071067812e-01,7.071067812e-01},{-5.773502692e-01,5.773502692e-01,5.773502692e-01},{0.000000000e+00,0.000000000e+00,1.000000000e+00}},
{{-7.071067812e-01,0.000000000e+00,7.071067812e-01},{0.000000000e+00,0.000000000e+00,1.000000000e+00},{-5.773502692e-01,5.773502692e-01,5.773502692e-01}},
{{-7.071067812e-01,0.000000000e+00,7.071067812e-01},{-5.773502692e-01,5.773502692e-01,5.773502692e-01},{-1.000000000e+00,-0.000000000e+00,-0.000000000e+00}},
{{7.071067812e-01,-7.071067812e-01,0.000000000e+00},{1.000000000e+00,0.000000000e+00,0.000000000e+00},{5.773502692e-01,-5.773502692e-01,5.773502692e-01}},
{{7.071067812e-01,-7.071067812e-01,0.000000000e+00},{5.773502692e-01,-5.773502692e-01,5.773502692e-01},{-0.000000000e+00,-1.000000000e+00,-0.000000000e+00}},
{{0.000000000e+00,-7.071067812e-01,7.071067812e-01},{-0.000000000e+00,-1.000000000e+00,-0.000000000e+00},{5.773502692e-01,-5.773502692e-01,5.773502692e-01}},
{{0.000000000e+00,-7.071067812e-01,7.071067812e-01},{5.773502692e-01,-5.773502692e-01,5.773502692e-01},{0.000000000e+00,0.000000000e+00,1.000000000e+00}},
{{7.071067812e-01,0.000000000e+00,7.071067812e-01},{0.000000000e+00,0.000000000e+00,1.000000000e+00},{5.773502692e-01,-5.773502692e-01,5.773502692e-01}},
{{7.071067812e-01,0.000000000e+00,7.071067812e-01},{5.773502692e-01,-5.773502692e-01,5.773502692e-01},{1.000000000e+00,0.000000000e+00,0.000000000e+00}},
{{-7.071067812e-01,-7.071067812e-01,-0.000000000e+00},{-5.773502692e-01,-5.773502692e-01,5.773502692e-01},{-1.000000000e+00,-0.000000000e+00,-0.000000000e+00}},
{{-7.071067812e-01,-7.071067812e-01,-0.000000000e+00},{-0.000000000e+00,-1.000000000e+00,-0.000000000e+00},{-5.773502692e-01,-5.773502692e-01,5.773502692e-01}},
{{0.000000000e+00,-7.071067812e-01,7.071067812e-01},{-5.773502692e-01,-5.773502692e-01,5.773502692e-01},{-0.000000000e+00,-1.000000000e+00,-0.000000000e+00}},
{{0.000000000e+00,-7.071067812e-01,7.071067812e-01},{0.000000000e+00,0.000000000e+00,1.000000000e+00},{-5.773502692e-01,-5.773502692e-01,5.773502692e-01}},
{{-7.071067812e-01,0.000000000e+00,7.071067812e-01},{-5.773502692e-01,-5.773502692e-01,5.773502692e-01},{0.000000000e+00,0.000000000e+00,1.000000000e+00}},
{{-7.071067812e-01,0.000000000e+00,7.071067812e-01},{-1.000000000e+00,-0.000000000e+00,-0.000000000e+00},{-5.773502692e-01,-5.773502692e-01,5.773502692e-01}},
{{7.071067812e-01,7.071067812e-01,0.000000000e+00},{1.000000000e+00,0.000000000e+00,0.000000000e+00},{5.773502692e-01,5.773502692e-01,-5.773502692e-01}},
{{7.071067812e-01,7.071067812e-01,0.000000000e+00},{5.773502692e-01,5.773502692e-01,-5.773502692e-01},{0.000000000e+00,1.000000000e+00,0.000000000e+00}},
{{0.000000000e+00,7.071067812e-01,-7.071067812e-01},{0.000000000e+00,1.000000000e+00,0.000000000e+00},{5.773502692e-01,5.773502692e-01,-5.773502692e-01}},
{{0.000000000e+00,7.071067812e-01,-7.071067812e-01},{5.773502692e-01,5.773502692e-01,-5.773502692e-01},{-0.000000000e+00,-0.000000000e+00,-1.000000000e+00}},
{{7.071067812e-01,0.000000000e+00,-7.071067812e-01},{-0.000000000e+00,-0.000000000e+00,-1.000000000e+00},{5.773502692e-01,5.773502692e-01,-5.773502692e-01}},
{{7.071067812e-01,0.000000000e+00,-7.071067812e-01},{5.773502692e-01,5.773502692e-01,-5.773502692e-01},{1.000000000e+00,0.000000000e+00,0.000000000e+00}},
{{-7.071067812e-01,7.071067812e-01,0.000000000e+00},{-5.773502692e-01,5.773502692e-01,-5.773502692e-01},{-1.000000000e+00,-0.000000000e+00,-0.000000000e+00}},
{{-7.071067812e-01,7.071067812e-01,0.000000000e+00},{0.000000000e+00,1.000000000e+00,0.000000000e+00},{-5.773502692e-01,5.773502692e-01,-5.773502692e-01}},
{{0.000000000e+00,7.071067812e-01,-7.071067812e-01},{-5.773502692e-01,5.773502692e-01,-5.773502692e-01},{0.000000000e+00,1.000000000e+00,0.000000000e+00}},
{{0.000000000e+00,7.071067812e-01,-7.071067812e-01},{-0.000000000e+00,-0.000000000e+00,-1.000000000e+00},{-5.773502692e-01,5.773502692e-01,-5.773502692e-01}},
{{-7.071067812e-01,-0.000000000e+00,-7.071067812e-01},{-5.773502692e-01,5.773502692e-01,-5.773502692e-01},{-0.000000000e+00,-0.000000000e+00,-1.000000000e+00}},
{{-7.071067812e-01,-0.000000000e+00,-7.071067812e-01},{-1.000000000e+00,-0.000000000e+00,-0.000000000e+00},{-5.773502692e-01,5.773502692e-01,-5.773502692e-01}},
{{7.071067812e-01,-7.071067812e-01,0.000000000e+00},{5.773502692e-01,-5.773502692e-01,-5.773502692e-01},{1.000000000e+00,0.000000000e+00,0.000000000e+00}},
{{7.071067812e-01,-7.071067812e-01,0.000000000e+00},{-0.000000000e+00,-1.000000000e+00,-0.000000000e+00},{5.773502692e-01,-5.773502692e-01,-5.773502692e-01}},
{{-0.000000000e+00,-7.071067812e-01,-7.071067812e-01},{5.773502692e-01,-5.773502692e-01,-5.773502692e-01},{-0.000000000e+00,-1.000000000e+00,-0.000000000e+00}},
{{-0.000000000e+00,-7.071067812e-01,-7.071067812e-01},{-0.000000000e+00,-0.000000000e+00,-1.000000000e+00},{5.773502692e-01,-5.773502692e-01,-5.773502692e-01}},
{{7.071067812e-01,0.000000000e+00,-7.071067812e-01},{5.773502692e-01,-5.773502692e-01,-5.773502692e-01},{-0.000000000e+00,-0.000000000e+00,-1.000000000e+00}},
{{7.071067812e-01,0.000000000e+00,-7.071067812e-01},{1.000000000e+00,0.000000000e+00,0.000000000e+00},{5.773502692e-01,-5.773502692e-01,-5.773502692e-01}},
{{-7.071067812e-01,-7.071067812e-01,-0.000000000e+00},{-1.000000000e+00,-0.000000000e+00,-0.000000000e+00},{-5.773502692e-01,-5.773502692e-01,-5.773502692e-01}},
{{-7.071067812e-01,-7.071067812e-01,-0.000000000e+00},{-5.773502692e-01,-5.773502692e-01,-5.773502692e-01},{-0.000000000e+00,-1.000000000e+00,-0.000000000e+00}},
{{-0.000000000e+00,-7.071067812e-01,-7.071067812e-01},{-0.000000000e+00,-1.000000000e+00,-0.000000000e+00},{-5.773502692e-01,-5.773502692e-01,-5.773502692e-01}},
{{-0.000000000e+00,-7.071067812e-01,-7.071067812e-01},{-5.773502692e-01,-5.773502692e-01,-5.773502692e-01},{-0.000000000e+00,-0.000000000e+00,-1.000000000e+00}},
{{-7.071067812e-01,-0.000000000e+00,-7.071067812e-01},{-0.000000000e+00,-0.000000000e+00,-1.000000000e+00},{-5.773502692e-01,-5.773502692e-01,-5.773502692e-01}},
{{-7.071067812e-01,-0.000000000e+00,-7.071067812e-01},{-5.773502692e-01,-5.773502692e-01,-5.773502692e-01},{-1.000000000e+00,-0.000000000e+00,-0.000000000e+00}}};


#endif // COMMON_H_
