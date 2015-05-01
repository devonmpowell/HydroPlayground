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

// global information
real L;
int nzones;
int boundary;
real time;

// conserved variables and fluxes
real* U[5];
real* U_old[5];

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

	// fluxes to compute
	real F[5];

	// max characteristic speed
	real cfl_fac = 0.5;
	real vcell;
	real ccell;
	real alpha_max;

	dt = 0.01;
	time = 0.0;
	tmax = 0.2;
	step = 0;


	while(time < tmax) {

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

		// get the max. characteristic speed on the grid
		alpha_max = 0.0;
		for(i = 0; i < nzones; ++i) {
			vcell = U[MOMX][ZR(i)]/U[RHO][ZR(i)];
			ccell = sqrt(eos_gamma*eos_p(U[RHO][ZR(i)], 
					U[ENERGY][ZR(i)]/U[RHO][ZR(i)] - 0.5*vcell*vcell, eos_gamma)/U[RHO][ZR(i)]);
			if(fabs(vcell + ccell) > alpha_max) alpha_max = fabs(vcell + ccell);
			if(fabs(vcell - ccell) > alpha_max) alpha_max = fabs(vcell - ccell);
		}

		// copy state vector U to a back-buffer
		for(i = 0; i < 5; ++i)
			if(U[i]) memcpy((void*) U_old[i], (void*) U[i], (nzones+2*NGHOST)*sizeof(real)); 
		
		// time step based on CFL condition
		dt = cfl_fac*dx/alpha_max;
		if(step%100 == 0) printf("Step %d, t = %f, dt = %f.\n", step, time, dt);

		// Solve the approximate Riemann problem for fluxes
		for(i = 0; i <= nzones; ++i) {

			// get local copies of the grid stencil
			// and solve for local primitives
			rho[LEFT] = U_old[RHO][ZL(i)];
			rho[RIGHT] = U_old[RHO][ZR(i)];
			E[LEFT] = U_old[ENERGY][ZL(i)];
			E[RIGHT] = U_old[ENERGY][ZR(i)];
			v[LEFT] = U_old[MOMX][ZL(i)]/rho[LEFT];
			v[RIGHT] = U_old[MOMX][ZR(i)]/rho[RIGHT];
			p[LEFT] = eos_p(rho[LEFT], (E[LEFT]/rho[LEFT] - 0.5*v[LEFT]*v[LEFT]), eos_gamma);
			p[RIGHT] = eos_p(rho[RIGHT], (E[RIGHT]/rho[RIGHT] - 0.5*v[RIGHT]*v[RIGHT]), eos_gamma);
			c[LEFT] = sqrt(eos_gamma*p[LEFT]/rho[LEFT]);
			c[RIGHT] = sqrt(eos_gamma*p[RIGHT]/rho[RIGHT]);

			// get characteristic speeds
			alpha[LEFT] = 0.0; // left wave speed
			if(v[LEFT] - c[LEFT] < alpha[LEFT]) alpha[LEFT] = v[LEFT] - c[LEFT];
			if(v[RIGHT] - c[RIGHT] < alpha[LEFT]) alpha[LEFT] = v[RIGHT] - c[RIGHT];
			alpha[RIGHT] = 0.0; // right wave speed
			if(v[RIGHT] + c[RIGHT] > alpha[RIGHT]) alpha[RIGHT] = v[RIGHT] + c[RIGHT];
			if(v[LEFT] + c[LEFT] > alpha[RIGHT]) alpha[RIGHT] = v[LEFT] + c[LEFT];

			// HLL fluxes
			F[RHO] = (alpha[RIGHT]*rho[LEFT]*v[LEFT] - alpha[LEFT]*rho[RIGHT]*v[RIGHT] 
					+ alpha[LEFT]*alpha[RIGHT]*(rho[RIGHT] - rho[LEFT]))/(alpha[RIGHT] - alpha[LEFT]); 
			F[MOMX] = (alpha[RIGHT]*(rho[LEFT]*v[LEFT]*v[LEFT] + p[LEFT]) - alpha[LEFT]*(rho[RIGHT]*v[RIGHT]*v[RIGHT] + p[RIGHT]) 
					+ alpha[LEFT]*alpha[RIGHT]*(rho[RIGHT]*v[RIGHT] - rho[LEFT]*v[LEFT]))/(alpha[RIGHT] - alpha[LEFT]); 
			F[ENERGY] = (alpha[RIGHT]*(E[LEFT] + p[LEFT])*v[LEFT] - alpha[LEFT]*(E[RIGHT] + p[RIGHT])*v[RIGHT] 
					+ alpha[LEFT]*alpha[RIGHT]*(E[RIGHT] - E[LEFT]))/(alpha[RIGHT] - alpha[LEFT]); 
			
			// update state vectors
			// forward Euler
			U[RHO][ZL(i)] -= dt/dx*F[RHO];
			U[RHO][ZR(i)] += dt/dx*F[RHO];
			U[MOMX][ZL(i)] -= dt/dx*F[MOMX];
			U[MOMX][ZR(i)] += dt/dx*F[MOMX];
			U[ENERGY][ZL(i)] -= dt/dx*F[ENERGY];
			U[ENERGY][ZR(i)] += dt/dx*F[ENERGY];

		}

		step++;
		time += dt;
	}
	printf("Done: %d steps, t = %f.\n", step, time);
}

void init() {

	int i;

	// grid allocation
	nzones = 256;
	L = 1.0;
	U[RHO] = (real*) malloc((nzones+2*NGHOST)*sizeof(real));
	U[MOMX] = (real*) malloc((nzones+2*NGHOST)*sizeof(real));
	U[ENERGY] = (real*) malloc((nzones+2*NGHOST)*sizeof(real));
	U_old[RHO] = (real*) malloc((nzones+2*NGHOST)*sizeof(real));
	U_old[MOMX] = (real*) malloc((nzones+2*NGHOST)*sizeof(real));
	U_old[ENERGY] = (real*) malloc((nzones+2*NGHOST)*sizeof(real));

	// equation of state
	eos_gamma = EOS_SEVEN_FIFTHS;
	eos_p = eos_ideal_p;

#if 1
	boundary = BC_WALL;
	real rho_left = 1.0;
	real p_left = 1.0;
	real rho_right = 0.125;
	real p_right = 0.1;
	for(i = 0; i < nzones/2; ++i) {
		U[RHO][ZR(i)] = rho_left; 
		U[MOMX][ZR(i)] = 0.0;
		U[ENERGY][ZR(i)] = p_left/(eos_gamma-1.0);
	}
	for(i = nzones/2; i < nzones; ++i) {
		U[RHO][ZR(i)] = rho_right; 
		U[MOMX][ZR(i)] = 0.0;
		U[ENERGY][ZR(i)] = p_right/(eos_gamma-1.0);
	}
#else
	// advection test
	boundary = BC_PERIODIC;
	real pressure = 1.0;
	real v = 1.0;
	for(i = 0; i < nzones; ++i) {
		U[RHO][ZR(i)] = 1.0; 
		if(i >= 3*nzones/8 && i < 5*nzones/8) U[RHO][ZR(i)] *= 2.0;
		U[MOMX][ZR(i)] = v*U[RHO][ZR(i)];
		U[ENERGY][ZR(i)] = pressure/(eos_gamma-1.0) + 1.0/(2.0*U[RHO][ZR(i)])*U[MOMX][ZR(i)]*U[MOMX][ZR(i)];
	}
#endif

}

void release() {
	// free fields and fluxes
	int i;
	for(i = 0; i < 5; ++i) {
		if(U[i]) free(U[i]);
		if(U_old[i]) free(U_old[i]);
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
	fwrite(&time, sizeof(real), 1, output);
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
