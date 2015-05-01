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
/*#include "hdf5.h"*/

typedef float real;

// conserved hydro quantities
int nzones;
real eos_gamma;
real t;
real L;

// conserved variables and fluxes
// plus macros for access
#define RHO 0
#define MOMX 1 
#define MOMY 2 
#define MOMZ 3 
#define ENERGY 4
real* U[5];
real* F[5];

// macros for getting zone indices from the face index
#define NGHOST 1 // num boundary ghosts
#define ZL(i) ((i)+NGHOST-1)
#define ZR(i) ((i)+NGHOST)
#define LEFT (NGHOST-1)
#define RIGHT (NGHOST)

int boundary;
#define BC_WALL 0
#define BC_PERIODIC 1

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

void evolve() {

	int i, step;
	real dx = L/nzones;
	real dt, tmax, tsnap;

	// primitive variables for the Riemann problem
	real rho[2*NGHOST]; // U[RHO] 
	real p[2*NGHOST]; // pressure
	real v[2*NGHOST]; // velocity
	real E[2*NGHOST]; // U[ENERGY] U[RHO]
	real c[2*NGHOST]; // sound speed
	real alpha[2*NGHOST]; // left and right wave speeds

	// max characteristic speed
	real alpha_max;

	dt = 0.01;
	t = 0.0;
	tmax = 0.2;
	step = 0;


	while(t < tmax) {

		// max characteristic speed on the grid
		alpha_max = 0.0;

		// copy BCs to ghost zones
		if(boundary == BC_WALL) for(i = 0; i < NGHOST; ++i) {
			U[RHO][ZL(-i)] = U[RHO][ZR(i)];
			U[RHO][ZR(nzones+i)] = U[RHO][ZL(nzones-i)];
			U[MOMX][ZL(-i)] = -U[MOMX][ZR(i)];
			U[MOMX][ZR(nzones+i)] = -U[MOMX][ZL(nzones-i)];
			U[ENERGY][ZL(-i)] = U[ENERGY][ZR(i)];
			U[ENERGY][ZR(nzones+i)] = U[ENERGY][ZL(nzones-i)];
		}
		else if(boundary == BC_PERIODIC) for(i = 0; i < NGHOST; ++i) {
			U[RHO][ZL(-i)] = U[RHO][ZL(nzones-i)];
			U[RHO][ZR(nzones+i)] = U[RHO][ZR(i)];
			U[MOMX][ZL(-i)] = U[MOMX][ZL(nzones-i)];
			U[MOMX][ZR(nzones+i)] = U[MOMX][ZR(i)];
			U[ENERGY][ZL(-i)] = U[ENERGY][ZL(nzones-i)];
			U[ENERGY][ZR(nzones+i)] = U[ENERGY][ZR(i)];
		}

		// Solve the approximate Riemann problem for fluxes
		for(i = 0; i <= nzones; ++i) {

			// get local copies of the grid stencil
			rho[LEFT] = U[RHO][ZL(i)];
			rho[RIGHT] = U[RHO][ZR(i)];
			E[LEFT] = U[ENERGY][ZL(i)];
			E[RIGHT] = U[ENERGY][ZR(i)];

			// solve for local primitives
			v[LEFT] = U[MOMX][ZL(i)]/U[RHO][ZL(i)];
			v[RIGHT] = U[MOMX][ZR(i)]/U[RHO][ZR(i)];
			p[LEFT] = (eos_gamma - 1.0)*(E[LEFT] - 0.5*rho[LEFT]*v[LEFT]*v[LEFT]);
			p[RIGHT] = (eos_gamma - 1.0)*(E[RIGHT] - 0.5*rho[RIGHT]*v[RIGHT]*v[RIGHT]);
			c[LEFT] = sqrt(eos_gamma*p[LEFT]/rho[LEFT]);
			c[RIGHT] = sqrt(eos_gamma*p[LEFT]/rho[LEFT]);

			// get characteristic speeds
			alpha[LEFT] = 0.0; // left wave speed
			if(v[LEFT] - c[LEFT] < alpha[LEFT]) alpha[LEFT] = v[LEFT] - c[LEFT];
			if(v[RIGHT] - c[RIGHT] < alpha[LEFT]) alpha[LEFT] = v[RIGHT] - c[RIGHT];
			alpha[RIGHT] = 0.0; // right wave speed
			if(v[RIGHT] + c[RIGHT] > alpha[RIGHT]) alpha[RIGHT] = v[RIGHT] + c[RIGHT];
			if(v[LEFT] + c[LEFT] > alpha[RIGHT]) alpha[RIGHT] = v[LEFT] + c[LEFT];

			// save max characteristic speed for CFL
			if(-alpha[LEFT] > alpha_max) alpha_max = -alpha[LEFT];
			if(alpha[RIGHT] > alpha_max) alpha_max = alpha[RIGHT];

			// HLL fluxes
			F[RHO][i] = (alpha[RIGHT]*rho[LEFT]*v[LEFT] - alpha[LEFT]*rho[RIGHT]*v[RIGHT] 
					+ alpha[LEFT]*alpha[RIGHT]*(rho[RIGHT] - rho[LEFT]))/(alpha[RIGHT] - alpha[LEFT]); 
			F[MOMX][i] = (alpha[RIGHT]*(rho[LEFT]*v[LEFT]*v[LEFT] + p[LEFT]) - alpha[LEFT]*(rho[RIGHT]*v[RIGHT]*v[RIGHT] + p[RIGHT]) 
					+ alpha[LEFT]*alpha[RIGHT]*(rho[RIGHT]*v[RIGHT] - rho[LEFT]*v[LEFT]))/(alpha[RIGHT] - alpha[LEFT]); 
			F[ENERGY][i] = (alpha[RIGHT]*(E[LEFT] + p[LEFT])*v[LEFT] - alpha[LEFT]*(E[RIGHT] + p[RIGHT])*v[RIGHT] 
					+ alpha[LEFT]*alpha[RIGHT]*(E[RIGHT] - E[LEFT]))/(alpha[RIGHT] - alpha[LEFT]); 

		}

		// CFL condition
		dt = 0.8*dx/alpha_max;
		if(step%100 == 0) printf("Step %d, t = %f, dt = %f.\n", step, t, dt);

		// update conserved quantities
		// forward Euler
		for(i = 0; i <= nzones; ++i) {
			U[RHO][ZL(i)] -= dt/dx*F[RHO][i];
			U[RHO][ZR(i)] += dt/dx*F[RHO][i];
			U[MOMX][ZL(i)] -= dt/dx*F[MOMX][i];
			U[MOMX][ZR(i)] += dt/dx*F[MOMX][i];
			U[ENERGY][ZL(i)] -= dt/dx*F[ENERGY][i];
			U[ENERGY][ZR(i)] += dt/dx*F[ENERGY][i];
		}

		step++;
		t += dt;
	}
	printf("Done: %d steps, t = %f.\n", step, t);
}

void init() {

	nzones = 256;
	L = 1.0;
	U[RHO] = (real*) malloc((nzones+2*NGHOST)*sizeof(real));
	U[MOMX] = (real*) malloc((nzones+2*NGHOST)*sizeof(real));
	U[ENERGY] = (real*) malloc((nzones+2*NGHOST)*sizeof(real));
	F[RHO] = (real*) malloc((nzones+1)*sizeof(real));
	F[MOMX] = (real*) malloc((nzones+1)*sizeof(real));
	F[ENERGY] = (real*) malloc((nzones+1)*sizeof(real));

	int i;
#if 1
	boundary = BC_WALL;
	real rho_left = 1.0;
	real p_left = 1.0;
	real rho_right = 0.125;
	real p_right = 0.1;
	eos_gamma = 1.4;
	for(i = NGHOST; i < NGHOST + nzones/2; ++i) {
		U[RHO][i] = rho_left; 
		U[MOMX][i] = 0.0;
		U[ENERGY][i] = p_left/(eos_gamma-1.0);
	}
	for(i = NGHOST + nzones/2; i < NGHOST + nzones; ++i) {
		U[RHO][i] = rho_right; 
		U[MOMX][i] = 0.0;
		U[ENERGY][i] = p_right/(eos_gamma-1.0);
	}
#else
	// advection test
	boundary = BC_PERIODIC;
	eos_gamma = 1.4;
	real pressure = 1.0;
	real v = 1.0;
	for(i = NGHOST; i < NGHOST + nzones; ++i) {
		U[RHO][i] = 1.0; 
		if(i >= NGHOST + 3*nzones/8 && i < NGHOST + 5*nzones/8) U[RHO][i] = 2.0;
		U[MOMX][i] = v*U[RHO][i];
		U[ENERGY][i] = pressure/(eos_gamma-1.0) + 1.0/(2.0*U[RHO][i])*U[MOMX][i]*U[MOMX][i];
	}
#endif

}

void release() {
	// free fields and fluxes
	int i;
	for(i = 0; i < 5; ++i) {
		if(U[i]) free(U[i]);
		if(F[i]) free(F[i]);
	}
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
	fwrite(&nzones, sizeof(int), 1, output);

	// write hydro quantities
	fwrite(&eos_gamma, sizeof(real), 1, output);
	fwrite(&t, sizeof(real), 1, output);
	fwrite(U[RHO]+NGHOST, sizeof(real), nzones, output);
	fwrite(U[MOMX]+NGHOST, sizeof(real), nzones, output);
	fwrite(U[ENERGY]+NGHOST, sizeof(real), nzones, output);

	// write nzones again as a check
	fwrite(&nzones, sizeof(int), 1, output);

	fclose(output);

	// hdf5 crap, got fed up with it.
#if 0
	setbuf(stdout, NULL);

	// hdf5 output
	hid_t file_id; 
	hsize_t HDF_Dims;

	printf("filename = %s\n", filename);
	
	// create/overwrite the file
	file_id= H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);


	// create dataspace with proper shape
	double test[1024];
	HDF_Dims = 1024;//nzones;
	HDF_DataspaceID = H5Screate_simple(1, &HDF_Dims, NULL);


	// write each hydro field
	HDF_Type = H5T_NATIVE_DOUBLE;
	printf("create dataset\n");
	HDF_DatasetID = H5Dcreate1(HDF_FileID, "Density", HDF_Type, 
			HDF_DataspaceID, H5P_DEFAULT);
	H5Dwrite( HDF_DatasetID, HDF_Type, H5S_ALL, H5S_ALL, H5P_DEFAULT, &test[0]);
	printf("done, closing\n");
	H5Dclose(HDF_DatasetID);
	/*printf("create dataset\n");*/
	/*HDF_DatasetID = H5Dcreate(HDF_FileID, "Momentum", HDF_Type, */
			/*HDF_DataspaceID, H5P_DEFAULT);//, H5P_DEFAULT, H5P_DEFAULT);*/
	/*H5Dwrite( HDF_DatasetID, HDF_Type, H5S_ALL, H5S_ALL, H5P_DEFAULT, U[MOMX]);*/
	/*printf("done, closing\n");*/
	/*H5Dclose(HDF_DatasetID);*/
	/*printf("create dataset\n");*/
	/*HDF_DatasetID = H5Dcreate(HDF_FileID, "Energy", HDF_Type, */
			/*HDF_DataspaceID, H5P_DEFAULT);//, H5P_DEFAULT, H5P_DEFAULT);*/
	/*H5Dwrite( HDF_DatasetID, HDF_Type, H5S_ALL, H5S_ALL, H5P_DEFAULT, U[ENERGY]);*/
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
