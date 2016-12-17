/*************************************************************
 *
 *		driver.c
 *
 *		Devon Powell
 *		21 April 2015
 *
 *		Main routines for d_euler_hydro
 *
 *************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "common.h"
#include "eos.h"

// the grid 
hydro_vector *grid; 
real dx; // always cubical grid cells
dvec nx;
int dim;
int bctype;
real time;

// equation of state functions
real eos_gamma;
real (*eos_p)(real, real, real);

// forward definitions
void init();
void release();
void evolve();
void write(char* filename);

int main() {
	init();
	evolve();
	write("out/test.dat");
	release();
	return 0;
}

// solves for local primitives from the conserved variables
hydro_derived_state get_derived_state(hydro_vector hvec) {
	int ax;
	real v2;
	hydro_derived_state ds;
	ds.rho = hvec.rho;
	ds.etot = hvec.etot;
	for(ax = 0, v2 = 0.0; ax < DIM; ++ax) {
		ds.v[ax] = hvec.mom[ax]/ds.rho;
		v2 += ds.v[ax]*ds.v[ax];	
	}
	ds.e = ds.etot/ds.rho - 0.5*v2;
	ds.p = eos_p(ds.rho, ds.e, eos_gamma);
	ds.c = sqrt(eos_gamma*ds.p/ds.rho);
	return ds;
}

void apply_boundary_conditions() {

	// copy BCs to ghost zones
	int ax, i;
	for(ax = 0; ax < DIM; ++ax) {
		if(bctype == BC_WALL) 
		for(i = 0; i < STENCIL_SIZE; ++i) {
			grid[ZL(-i)].rho = grid[ZR(i)].rho;
			grid[ZR(nx[ax]+i)].rho = grid[ZL(nx[ax]-i)].rho;
			grid[ZL(-i)].mom[ax] = -grid[ZR(i)].mom[ax];
			grid[ZR(nx[ax]+i)].mom[ax] = -grid[ZL(nx[ax]-i)].mom[ax];
			grid[ZL(-i)].etot = grid[ZR(i)].etot;
			grid[ZR(nx[ax]+i)].etot = grid[ZL(nx[ax]-i)].etot;
		}
		else if(bctype == BC_PERIODIC) 
		for(i = 0; i < STENCIL_SIZE; ++i) {
			grid[ZL(-i)].rho = grid[ZL(nx[ax]-i)].rho;
			grid[ZR(nx[ax]+i)].rho = grid[ZR(i)].rho;
			grid[ZL(-i)].mom[ax] = grid[ZL(nx[ax]-i)].mom[ax];
			grid[ZR(nx[ax]+i)].mom[ax] = grid[ZR(i)].mom[ax];
			grid[ZL(-i)].etot = grid[ZL(nx[ax]-i)].etot;
			grid[ZR(nx[ax]+i)].etot = grid[ZR(i)].etot;
		}	
	}
}

real get_max_char_speed() {
		
	int i, ax;
	real alpha_max;
	hydro_derived_state ds;
	
	// TODO: properly treat the (possibly different) axis sizes??
	for(i = 0, alpha_max = 0.0; i < nx[0]; ++i) {
		ds = get_derived_state(grid[ZR(i)]);
		for(ax = 0; ax < DIM; ++ax) {
			if(fabs(ds.v[ax] + ds.c) > alpha_max) alpha_max = fabs(ds.v[ax] + ds.c);
			if(fabs(ds.v[ax] - ds.c) > alpha_max) alpha_max = fabs(ds.v[ax] - ds.c);
		}
	}
	return alpha_max;
}

hydro_vector flux_solve_hll(hydro_derived_state left, hydro_derived_state right, int skax) {

	int ax;
	real left_alpha, right_alpha;
	hydro_vector flux;

	// get characteristic speeds for the HLL solver
	left_alpha = 0.0; // left wave speed
	if(left.c - left.v[skax] > left_alpha) left_alpha = left.c - left.v[skax];
	if(right.c - right.v[skax] > left_alpha) left_alpha = right.c - right.v[skax];
	right_alpha = 0.0; // right wave speed
	if(right.v[skax] + right.c > right_alpha) right_alpha = right.v[skax] + right.c;
	if(left.v[skax] + left.c > right_alpha) right_alpha = left.v[skax] + left.c;

	//printf("alphas = %f %f\n", left_alpha, right_alpha);

	// Compute HLL fluxes
	flux.rho = (right_alpha*left.rho*left.v[skax] + left_alpha*right.rho*right.v[skax] 
			- left_alpha*right_alpha*(right.rho - left.rho))/(right_alpha + left_alpha); 
	for(ax = 0; ax < DIM; ++ax) {
		flux.mom[ax] = (right_alpha*(left.rho*left.v[ax]*left.v[skax] + (ax==skax)*left.p) + left_alpha*(right.rho*right.v[ax]*right.v[skax] + (ax==skax)*right.p) 
			- left_alpha*right_alpha*(right.rho*right.v[ax] - left.rho*left.v[ax]))/(right_alpha + left_alpha); 
	}
	flux.etot = (right_alpha*(left.etot + left.p)*left.v[skax] + left_alpha*(right.etot + right.p)*right.v[skax]
			- left_alpha*right_alpha*(right.etot - left.etot))/(right_alpha + left_alpha); 
	return flux;
}

//long grind(dvec ind) [>{<]
	//long gr = 0;
	////gr += (ind[0]+STENCIL_SIZE)*(nx[1]+2*STENCIL_SIZE)*(nx[2]+2*STENCIL_SIZE);
	//return ind[0]+STENCIL_SIZE;
/*}*/

void update_skewer(dvec skewer, int iterax, real dt) {
	int i, ax;

	// get a local copy of the grid skewer 
	// and solve for local primitives
	// including ghost cells
	hydro_derived_state state_skewer[nx[iterax]+2*STENCIL_SIZE];
	for(skewer[iterax] = -STENCIL_SIZE; skewer[iterax] < nx[iterax]+STENCIL_SIZE; ++skewer[iterax])
		state_skewer[ZR(skewer[iterax])] = get_derived_state(grid[ZR(skewer[iterax])]);

		//printf("c:    ");
		//for(i = -STENCIL_SIZE; i < nx[0]+STENCIL_SIZE; ++i) {
			//if(i < 0 || i >= nx[0]) printf("(");
			//printf(" %.3f", state_skewer[ZR(i)].c);
			//if(i < 0 || i >= nx[0]) printf(")");
		//}
		//printf("\n");

		//printf("vx:   ");
		//for(i = -STENCIL_SIZE; i < nx[0]+STENCIL_SIZE; ++i) {
			//if(i < 0 || i >= nx[0]) printf("(");
			//printf(" %.3f", state_skewer[ZR(i)].v[0]);
			//if(i < 0 || i >= nx[0]) printf(")");
		//}
		//printf("\n");

	
	// we have skewer primitives in hand, now solve for fluxes and step forward in time 
	// forward Euler
	hydro_vector flux;
	real stfac = 0.5*dt/dx; // half-step factor
	for(i = 0; i <= nx[iterax]; ++i) {
		// solve the Riemann problem
		flux = flux_solve_hll(state_skewer[ZL(i)], state_skewer[ZR(i)], iterax);
		grid[ZL(i)].rho -= stfac*flux.rho;
		grid[ZR(i)].rho += stfac*flux.rho;
		for(ax = 0; ax < DIM; ++ax) {
			grid[ZL(i)].mom[ax] -= stfac*flux.mom[ax];
			grid[ZR(i)].mom[ax] += stfac*flux.mom[ax];
		}
		grid[ZL(i)].etot -= stfac*flux.etot;
		grid[ZR(i)].etot += stfac*flux.etot;
	}
}

void evolve() {

	int i, step, ax, iterax;
	dvec skewer;
	real dt, tmax, tsnap;
  
	// max characteristic speed
	real cfl_fac = 0.5;
	real alpha_max;

	dt = 0.01;
	time = 0.0;
	tmax = 0.2;
	step = 0;
	while(time < tmax) {

		// apply the boundary conditions to the grid
		apply_boundary_conditions();

		// get the max. characteristic speed on the grid
		// time step based on CFL condition
		alpha_max = get_max_char_speed();
		//printf("Step %d, alpha_max = %f\n", step, alpha_max);
		//printf("rho:  ");
		//for(i = -STENCIL_SIZE; i < nx[0]+STENCIL_SIZE; ++i) {
			//if(i < 0 || i >= nx[0]) printf("(");
			//printf(" %.3f", grid[ZR(i)].rho);
			//if(i < 0 || i >= nx[0]) printf(")");
		//}
		//printf("\n");


		dt = cfl_fac*dx/alpha_max;
		if(step%100 == 0) printf("Step %d, t = %f, dt = %f.\n", step, time, dt);


		// Solve the approximate Riemann problem for fluxes
		// Strang splitting
		for(iterax = 0; iterax < DIM; ++iterax) {

			// iterate over all skewers in this axis 
#if DIM > 1
			for(i = 0; i < nx[(iterax+1)%DIM]; ++i) {
				skewer[(iterax+1)%DIM] = i;
#if DIM > 2
				for(j = 0; j < nx[(iterax+2)%DIM]; ++j) {
					skewer[(iterax+2)%DIM] = j;

#endif
#endif
					//printf("My index = %d %d %d, iterax = \n", skewer.)
					// step forward in time for the 1D skewer 
					// described by skewer and iterax
					update_skewer(skewer, iterax, dt);
#if DIM > 1
				}
#if DIM > 2
			}
#endif
#endif
		}
		for(iterax = DIM-1; iterax >= 0; --iterax) {

			// iterate over all skewers in this axis 
#if DIM > 1
			for(i = 0; i < nx[(iterax+1)%DIM]; ++i) {
				skewer[(iterax+1)%DIM] = i;
#if DIM > 2
				for(j = 0; j < nx[(iterax+2)%DIM]; ++j) {
					skewer[(iterax+2)%DIM] = j;
#endif
#endif
					// step forward in time for the 1D skewer 
					// described by skewer and iterax
					update_skewer(skewer, iterax, dt);
#if DIM > 1
				}
#if DIM > 2
			}
#endif
#endif
		}

		step++;
		time += dt;
	}
	printf("Done: %d steps, t = %f.\n", step, time);
}

void init() {

	int i, ax, ncells, maxdim;

	// grid setup 
	int N = 256;
	real L = 1.0;
	dx = L/N;
	for(ax = 0, ncells = 1, maxdim = 0; ax < DIM; ++ax) {
		nx[ax] = N;
		ncells *= (nx[ax]+2*STENCIL_SIZE);
		if(nx[ax] > maxdim) maxdim = nx[ax];
	}

	grid = (hydro_vector*) malloc(ncells*sizeof(hydro_vector));

	// equation of state
	eos_gamma = EOS_SEVEN_FIFTHS;
	eos_p = eos_ideal_p;

#if 1
	// Sod shock tube
	bctype = BC_WALL;
	real rho_left = 1.0;
	real p_left = 1.0;
	real rho_right = 0.125;
	real p_right = 0.1;
	for(i = 0; i < nx[DIM-1]/2; ++i) {
		grid[ZR(i)].rho = rho_left; 
		for(ax = 0; ax < DIM; ++ax)
			grid[ZR(i)].mom[ax] = 0.0;
		grid[ZR(i)].etot = p_left/(eos_gamma-1.0);
	}
	for(i = nx[DIM-1]/2; i < nx[DIM-1]; ++i) {
		grid[ZR(i)].rho = rho_right; 
		for(ax = 0; ax < DIM; ++ax)
			grid[ZR(i)].mom[ax] = 0.0;
		grid[ZR(i)].etot = p_right/(eos_gamma-1.0);
	}
#else
	// advection test
	boundary = BC_PERIODIC;
	real pressure = 1.0;
	real v = 1.0;
	for(i = 0; i < nzones; ++i) {
		grid[RHO][ZR(i)] = 1.0; 
		if(i >= 3*nzones/8 && i < 5*nzones/8) grid[RHO][ZR(i)] *= 2.0;
		grid[MOMX][ZR(i)] = v*grid[RHO][ZR(i)];
		grid[ENERGY][ZR(i)] = pressure/(eos_gamma-1.0) + 1.0/(2.0*grid[RHO][ZR(i)])*grid[MOMX][ZR(i)]*grid[MOMX][ZR(i)];
	}
#endif

}

void release() {
	// free fields and fluxes
	if(grid) free(grid);
}

void write(char* filename) {

	// binary i/o for now
	FILE* output;
	output = fopen(filename, "w");
	if(!output) {
		printf("Failed to open %s.\n", filename);
		return;
	}

	// write nzones
	int dim = DIM;
	fwrite(&dim, sizeof(int), 1, output);
	fwrite(&nx, DIM*sizeof(int), 1, output);

	// write hydro quantities
	fwrite(&eos_gamma, sizeof(real), 1, output);
	fwrite(&time, sizeof(real), 1, output);

	// TODO:
	// how to copy skewers...?
	// this will work in 1D for now
	fwrite(grid+STENCIL_SIZE, sizeof(hydro_vector), nx[0], output);

	// write dim again as a check
	fwrite(&dim, sizeof(int), 1, output);

	fclose(output);

	// hdf5 crap, got fed up with it.
#if 0
	setbuf(stdout, NgridLL);

	// hdf5 output
	hid_t file_id; 
	hsize_t HDF_Dims;

	printf("filename = %s\n", filename);
	
	// create/overwrite the file
	file_id= H5Fcreate(filename, H5F_ACC_TRgridNC, H5P_DEFAgridLT, H5P_DEFAgridLT);


	// create dataspace with proper shape
	double test[1024];
	HDF_Dims = 1024;//nzones;
	HDF_DataspaceID = H5Screate_simple(1, &HDF_Dims, NgridLL);


	// write each hydro field
	HDF_Type = H5T_NATIVE_DOgridBLE;
	printf("create dataset\n");
	HDF_DatasetID = H5Dcreate1(HDF_FileID, "Density", HDF_Type, 
			HDF_DataspaceID, H5P_DEFAgridLT);
	H5Dwrite( HDF_DatasetID, HDF_Type, H5S_ALL, H5S_ALL, H5P_DEFAgridLT, &test[0]);
	printf("done, closing\n");
	H5Dclose(HDF_DatasetID);
	/*printf("create dataset\n");*/
	/*HDF_DatasetID = H5Dcreate(HDF_FileID, "Momentum", HDF_Type, */
			/*HDF_DataspaceID, H5P_DEFAgridLT);//, H5P_DEFAgridLT, H5P_DEFAgridLT);*/
	/*H5Dwrite( HDF_DatasetID, HDF_Type, H5S_ALL, H5S_ALL, H5P_DEFAgridLT, grid[MOMX]);*/
	/*printf("done, closing\n");*/
	/*H5Dclose(HDF_DatasetID);*/
	/*printf("create dataset\n");*/
	/*HDF_DatasetID = H5Dcreate(HDF_FileID, "Energy", HDF_Type, */
			/*HDF_DataspaceID, H5P_DEFAgridLT);//, H5P_DEFAgridLT, H5P_DEFAgridLT);*/
	/*H5Dwrite( HDF_DatasetID, HDF_Type, H5S_ALL, H5S_ALL, H5P_DEFAgridLT, grid[ENERGY]);*/
	/*printf("done, closing\n");*/
	/*H5Dclose(HDF_DatasetID);*/

	// close dataspace
	printf("close dataspace\n");
	H5Sclose(HDF_DataspaceID);
	// close the file
	printf("close file\n");
	H5Fclose(file_id);
	printf("done\n");
#endif
}
