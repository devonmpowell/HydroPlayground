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


// forward definitions
int init();
void release();
void evolve();
void writeout(hydro_problem* hp);


long flat_index(dvec grind, hydro_problem* hp) {
	int ax;
	long flatind;
	flatind = 0; 
	for(ax = 0; ax < hp->dim; ++ax)
		flatind += hp->strides[ax]*(STENCIL_SIZE+grind[ax]);
	return flatind; 
}


// solves for local primitives from the conserved variables
hydro_derived_state get_derived_state(dvec cell, hydro_problem* hp) {
	int ax;
	real v2;
	hydro_vector hvec;
	hydro_derived_state ds;
	hvec = hp->grid[flat_index(cell, hp)];
	ds.rho = hvec.rho;
	ds.etot = hvec.etot;
	for(ax = 0, v2 = 0.0; ax < hp->dim; ++ax) {
		ds.v[ax] = hvec.mom[ax]/ds.rho;
		v2 += ds.v[ax]*ds.v[ax];	
	}
	ds.e = ds.etot/ds.rho - 0.5*v2;
	ds.p = hp->eos_p(ds.rho, ds.e, hp->eos_gamma);
	ds.c = sqrt(hp->eos_gamma*ds.p/ds.rho);
	return ds;
}

real get_max_char_speed(hydro_problem* hp) {
		
	int ax;
	dvec grind;
	real alpha_max;
	hydro_derived_state ds;

	alpha_max = 0.0;
	if(hp->dim == 1)
	for(grind[0] = 0; grind[0] < hp->nx[0]; ++grind[0]) {
		ds = get_derived_state(grind, hp);
		for(ax = 0; ax < hp->dim; ++ax) {
			if(fabs(ds.v[ax] + ds.c) > alpha_max) alpha_max = fabs(ds.v[ax] + ds.c);
			if(fabs(ds.v[ax] - ds.c) > alpha_max) alpha_max = fabs(ds.v[ax] - ds.c);
		}
	}
	else if(hp->dim == 2)
	for(grind[0] = 0; grind[0] < hp->nx[0]; ++grind[0]) {
		for(grind[1] = 0; grind[1] < hp->nx[1]; ++grind[1]) {
		ds = get_derived_state(grind, hp);
			for(ax = 0; ax < hp->dim; ++ax) {
				if(fabs(ds.v[ax] + ds.c) > alpha_max) alpha_max = fabs(ds.v[ax] + ds.c);
				if(fabs(ds.v[ax] - ds.c) > alpha_max) alpha_max = fabs(ds.v[ax] - ds.c);
			}
		}
	}
	else if(hp->dim == 3)
	for(grind[0] = 0; grind[0] < hp->nx[0]; ++grind[0]) {
		for(grind[1] = 0; grind[1] < hp->nx[1]; ++grind[1]) {
			for(grind[2] = 0; grind[2] < hp->nx[2]; ++grind[2]) {
				ds = get_derived_state(grind, hp);
				for(ax = 0; ax < hp->dim; ++ax) {
					if(fabs(ds.v[ax] + ds.c) > alpha_max) alpha_max = fabs(ds.v[ax] + ds.c);
					if(fabs(ds.v[ax] - ds.c) > alpha_max) alpha_max = fabs(ds.v[ax] - ds.c);
				}
			}
		}
	}
	return alpha_max;
}

hydro_vector flux_solve_hll(hydro_derived_state left, 
		hydro_derived_state right, int skax, hydro_problem* hp) {

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

	// Compute HLL fluxes
	flux.rho = (right_alpha*left.rho*left.v[skax] + left_alpha*right.rho*right.v[skax] 
			- left_alpha*right_alpha*(right.rho - left.rho))/(right_alpha + left_alpha); 
	for(ax = 0; ax < hp->dim; ++ax) {
		flux.mom[ax] = (right_alpha*(left.rho*left.v[ax]*left.v[skax] + (ax==skax)*left.p) + left_alpha*(right.rho*right.v[ax]*right.v[skax] + (ax==skax)*right.p) 
			- left_alpha*right_alpha*(right.rho*right.v[ax] - left.rho*left.v[ax]))/(right_alpha + left_alpha); 
	}
	flux.etot = (right_alpha*(left.etot + left.p)*left.v[skax] + left_alpha*(right.etot + right.p)*right.v[skax]
			- left_alpha*right_alpha*(right.etot - left.etot))/(right_alpha + left_alpha); 
	return flux;
}

void update_skewer(dvec skewer, int iterax, real fluxfac, hydro_problem* hp) {

	// get a local copy of the grid skewer 
	// and solve for local primitives
	int i, ax;
	long flatind;
	hydro_derived_state state_skewer[hp->nx[iterax]+2*STENCIL_SIZE];
	for(i = 0; i < hp->nx[iterax]; ++i) {
		skewer[iterax] = i;
		state_skewer[ZR(i)] = get_derived_state(skewer, hp);
	}

	// enforce boundary conditions here
	if(hp->bctype == BC_WALL) 
	for(i = 0; i < STENCIL_SIZE; ++i) {
		state_skewer[ZL(-i)] = state_skewer[ZR(i)];
		state_skewer[ZL(-i)].v[iterax] *= -1;
		state_skewer[ZR(hp->nx[iterax]+i)] = state_skewer[ZL(hp->nx[iterax]-i)];
		state_skewer[ZR(hp->nx[iterax]+i)].v[iterax] *= -1;
	}
	else if(hp->bctype == BC_PERIODIC) 
	for(i = 0; i < STENCIL_SIZE; ++i) {
		state_skewer[ZL(-i)] = state_skewer[ZL(hp->nx[iterax]-i)];
		state_skewer[ZR(hp->nx[iterax]+i)] = state_skewer[ZR(i)];
	}
	
	// we have skewer primitives in hand, now solve for fluxes and step forward in time 
	// forward Euler
	hydro_vector flux;
	for(i = 0; i <= hp->nx[iterax]; ++i) {
		// solve the Riemann problem
		flux = flux_solve_hll(state_skewer[ZL(i)], state_skewer[ZR(i)], iterax, hp);

		// update both grid cells
		skewer[iterax] = i-1;
		flatind = flat_index(skewer, hp);
		hp->grid[flatind].rho -= fluxfac*flux.rho;
		for(ax = 0; ax < hp->dim; ++ax)
			hp->grid[flatind].mom[ax] -= fluxfac*flux.mom[ax];
		hp->grid[flatind].etot -= fluxfac*flux.etot;
		skewer[iterax] = i;
		flatind = flat_index(skewer, hp);
		hp->grid[flatind].rho += fluxfac*flux.rho;
		for(ax = 0; ax < hp->dim; ++ax)
			hp->grid[flatind].mom[ax] += fluxfac*flux.mom[ax];
		hp->grid[flatind].etot += fluxfac*flux.etot;
	}
}

void evolve(real tstop, int max_steps, int output_every, hydro_problem* hp) {

	int i, j, iterax, done;
	dvec skewer;
	real dt, alpha_max;
  
	// loop until t > tmax
	tstop += hp->time;
	for(done = 0; hp->step < max_steps && !done; hp->step++, hp->time += dt) {

		// get the max. characteristic speed on the grid
		// time step based on CFL condition
		alpha_max = get_max_char_speed(hp);
		dt = hp->cfl_fac*hp->dx/alpha_max;
		if(hp->time + dt > tstop) {
			dt = tstop - hp->time;
			done = 1;
		}

		// output if called for
		if(hp->step%output_every == 0) {
			printf("Step %d, t = %f, dt = %f.\n", hp->step, hp->time, dt);
			writeout(hp);
		} 

		// Solve the approximate Riemann problem for fluxes
		// Strang splitting
		real stfac = 0.5*dt/hp->dx; // half-step factor
		for(iterax = 0; iterax < hp->dim; ++iterax) {
			// iterate over all skewers in this axis 
			if(hp->dim == 1)
				update_skewer(skewer, iterax, stfac, hp);
			else if(hp->dim == 2)
				for(i = 0; i < hp->nx[(iterax+1)%hp->dim]; ++i) {
					skewer[(iterax+1)%hp->dim] = i;
					update_skewer(skewer, iterax, stfac, hp);
				}
			else if(hp->dim == 3)
				for(i = 0; i < hp->nx[(iterax+1)%hp->dim]; ++i) {
					skewer[(iterax+1)%hp->dim] = i;
					for(j = 0; j < hp->nx[(iterax+2)%hp->dim]; ++j) {
						skewer[(iterax+2)%hp->dim] = j;
						update_skewer(skewer, iterax, stfac, hp);
					}
				}
		}
		for(iterax = hp->dim-1; iterax >= 0; --iterax) {
			if(hp->dim == 1)
				update_skewer(skewer, iterax, stfac, hp);
			else if(hp->dim == 2)
				for(i = 0; i < hp->nx[(iterax+1)%hp->dim]; ++i) {
					skewer[(iterax+1)%hp->dim] = i;
					update_skewer(skewer, iterax, stfac, hp);
				}
			else if(hp->dim == 3)
				for(i = 0; i < hp->nx[(iterax+1)%hp->dim]; ++i) {
					skewer[(iterax+1)%hp->dim] = i;
					for(j = 0; j < hp->nx[(iterax+2)%hp->dim]; ++j) {
						skewer[(iterax+2)%hp->dim] = j;
						update_skewer(skewer, iterax, stfac, hp);
					}
				}
		}
	}
	printf("Done: %d steps, t = %f.\n", hp->step, hp->time);
	writeout(hp);
}

void writeout(hydro_problem* hp) {

	char filename[128];
	sprintf(filename, "out/%s_%06d", hp->name, hp->step);
	printf(" Writing file %s...\n", filename);

	// binary i/o for now
	FILE* output;
	output = fopen(filename, "wb");
	if(!output) {
		printf("Failed to open %s.\n", filename);
		return;
	}

	// write nzones
	const int check = 0xD00D;
	fwrite(&check, sizeof(int), 1, output);
	fwrite(&hp->dim, sizeof(int), 1, output);
	fwrite(&hp->nx, 3*sizeof(int), 1, output);

	// write hydro quantities
	fwrite(&hp->eos_gamma, sizeof(real), 1, output);
	fwrite(&hp->time, sizeof(real), 1, output);

	// TODO:
	// how to copy skewers...?
	// this will work in 1D for now
	fwrite(hp->grid+STENCIL_SIZE, sizeof(hydro_vector), hp->nx[0], output);

	// check again
	fwrite(&check, sizeof(int), 1, output);
	fclose(output);
}

int init(hydro_problem* hp) {

	// TODO: this is the only thing to get into Python!!!
	hp->eos_p = eos_ideal_p;
	return 1;
}


