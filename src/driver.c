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
#include "radiation.h"


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
hydro_derived_state get_derived_state(dvec cell, int iterax, hydro_problem* hp) {
	int ax;
	real cm, v2;
	hydro_vector hvec;
	hydro_derived_state ds;
	hvec = hp->grid[flat_index(cell, hp)];

	// get the gradients as computed from the advected CM position
	ds.rho_l = hvec.rho-6*hvec.com[iterax]/(hp->dx);
	ds.rho_c = hvec.rho;
	ds.rho_r = hvec.rho+6*hvec.com[iterax]/(hp->dx);

	ds.rho_l = ds.rho_c; 
	ds.rho_r = ds.rho_c;



	ds.etot = hvec.etot;
	for(ax = 0, v2 = 0.0; ax < hp->dim; ++ax) {
		ds.v[ax] = hvec.mom[ax]/hvec.rho;
		v2 += ds.v[ax]*ds.v[ax];	
	}
	ds.e_l = ds.etot/ds.rho_l - 0.5*v2; // TODO: is this correct for the energy equation??
	ds.e_c = ds.etot/ds.rho_c - 0.5*v2; // TODO: is this correct for the energy equation??
	ds.e_r = ds.etot/ds.rho_r - 0.5*v2; // TODO: is this correct for the energy equation??
	ds.p_l = hp->eos_p(ds.rho_l, ds.e_l, hp->eos_gamma);
	ds.p_c = hp->eos_p(ds.rho_c, ds.e_c, hp->eos_gamma);
	ds.p_r = hp->eos_p(ds.rho_r, ds.e_r, hp->eos_gamma);
	ds.c_l = sqrt(hp->eos_gamma*ds.p_l/ds.rho_l);
	ds.c_c = sqrt(hp->eos_gamma*ds.p_c/ds.rho_c);
	ds.c_r = sqrt(hp->eos_gamma*ds.p_r/ds.rho_r);
	return ds;
}

real get_max_char_speed(hydro_problem* hp) {
		
	int ax;
	dvec grind;
	real alpha_max;
	hydro_derived_state ds;

// a macro for cells, to avoid repeating it in the loops
#define get_cell_char_speed(grind, hp) {	\
	for(ax = 0; ax < hp->dim; ++ax) {	\
		ds = get_derived_state(grind, ax, hp);	\
		if(fabs(ds.v[ax] + ds.c_l) > alpha_max) alpha_max = fabs(ds.v[ax] + ds.c_l); 	\
		if(fabs(ds.v[ax] - ds.c_l) > alpha_max) alpha_max = fabs(ds.v[ax] - ds.c_l);	\
		if(fabs(ds.v[ax] + ds.c_r) > alpha_max) alpha_max = fabs(ds.v[ax] + ds.c_r); 	\
		if(fabs(ds.v[ax] - ds.c_r) > alpha_max) alpha_max = fabs(ds.v[ax] - ds.c_r);	\
	}	\
}

	alpha_max = 0.0;
	if(hp->dim == 1) {
		for(grind[0] = 0; grind[0] < hp->nx[0]; ++grind[0])
			get_cell_char_speed(grind, hp);
	}
	else if(hp->dim == 2) {
		for(grind[0] = 0; grind[0] < hp->nx[0]; ++grind[0])
		for(grind[1] = 0; grind[1] < hp->nx[1]; ++grind[1])
			get_cell_char_speed(grind, hp);
	}
	else if(hp->dim == 3) {
		for(grind[0] = 0; grind[0] < hp->nx[0]; ++grind[0])
		for(grind[1] = 0; grind[1] < hp->nx[1]; ++grind[1])
		for(grind[2] = 0; grind[2] < hp->nx[2]; ++grind[2])
			get_cell_char_speed(grind, hp);
	}
	return alpha_max;
}

void update_source_terms(real dt, hydro_problem* hp) {

	int ax;
	long flatind;
	dvec grind;
	hydro_vector* cell;

// a macro for cells, to avoid repeating it in the loops
#define update_cell_source_terms(grind, hp) {	\
	flatind = flat_index(grind, hp);	\
	cell = &hp->grid[flatind];	\
	/* TODO: all other source terms (gravity?)*/ \
	/* TODO: implicit??*/ \
	cell = &hp->grid[flatind];	\
	for(ax = 0; ax < hp->dim; ++ax)	\
		cell->com[ax] += dt*cell->mom[ax];	\
}

	// go over all grid cells and add in source/sink terms for that cell
	// TODO: handle these implicitly or explicitly?
	if(hp->dim == 1) {
		for(grind[0] = 0; grind[0] < hp->nx[0]; ++grind[0])
			update_cell_source_terms(grind, hp);
	}
	else if(hp->dim == 2) {
		for(grind[0] = 0; grind[0] < hp->nx[0]; ++grind[0])
		for(grind[1] = 0; grind[1] < hp->nx[1]; ++grind[1])
			update_cell_source_terms(grind, hp);
	}
	else if(hp->dim == 3) {
		for(grind[0] = 0; grind[0] < hp->nx[0]; ++grind[0])
		for(grind[1] = 0; grind[1] < hp->nx[1]; ++grind[1])
		for(grind[2] = 0; grind[2] < hp->nx[2]; ++grind[2])
			update_cell_source_terms(grind, hp);
	}
}

hydro_vector flux_solve_hll(hydro_derived_state left, 
		hydro_derived_state right, int skax, hydro_problem* hp) {

	int ax;
	real left_alpha, right_alpha;
	hydro_vector flux;

	// get characteristic speeds for the HLL solver
	left_alpha = 0.0; // left wave speed
	if(left.c_r - left.v[skax] > left_alpha) left_alpha = left.c_r - left.v[skax];
	if(right.c_l - right.v[skax] > left_alpha) left_alpha = right.c_l - right.v[skax];
	right_alpha = 0.0; // right wave speed
	if(right.v[skax] + right.c_l > right_alpha) right_alpha = right.v[skax] + right.c_l;
	if(left.v[skax] + left.c_r > right_alpha) right_alpha = left.v[skax] + left.c_r;

	// Compute HLL fluxes
	flux.rho = (right_alpha*left.rho_r*left.v[skax] + left_alpha*right.rho_l*right.v[skax] 
			- left_alpha*right_alpha*(right.rho_l - left.rho_r))/(right_alpha + left_alpha); 
	for(ax = 0; ax < hp->dim; ++ax) {
		flux.mom[ax] = (right_alpha*(left.rho_r*left.v[ax]*left.v[skax] + (ax==skax)*left.p_r) + left_alpha*(right.rho_l*right.v[ax]*right.v[skax] + (ax==skax)*right.p_l) 
			- left_alpha*right_alpha*(right.rho_l*right.v[ax] - left.rho_r*left.v[ax]))/(right_alpha + left_alpha); 
	}
	flux.etot = (right_alpha*(left.etot + left.p_r)*left.v[skax] + left_alpha*(right.etot + right.p_l)*right.v[skax]
			- left_alpha*right_alpha*(right.etot - left.etot))/(right_alpha + left_alpha); 
	return flux;
}

void psi_flux_limiter(hydro_derived_state* dsll, hydro_derived_state* dsl, hydro_derived_state* dsr, hydro_derived_state* dsrr, hydro_problem* hp) {

	// for rho, vel, pressure...
	int ax;
	real grad_l, grad_r, grad_max; 


#if 1
	// generalized minmod limiter, with theta = 2.0 (minimal diffusivity)
	real theta = 2.9;
#define sgn(x) (((x) > 0)? 1 : -1)
#define min(x, y) (((x) < (y))? x : y)
#define min3(x, y, z) min(x, min(y, z)) 
#define minmod(x, y, z) (0.25*abs(sgn(x)+sgn(y))*(sgn(x)+sgn(z))*min3(abs(x), abs(y), abs(z)))

	real c_l, c_r;
	c_l = dsl->rho_c + 0.5*minmod(theta*(dsl->rho_c-dsll->rho_c), 0.5*(dsr->rho_c-dsll->rho_c), theta*(dsr->rho_c-dsl->rho_c));
	c_r = dsr->rho_c - 0.5*minmod(theta*(dsr->rho_c-dsl->rho_c), 0.5*(dsrr->rho_c-dsl->rho_c), theta*(dsrr->rho_c-dsr->rho_c));
	dsl->rho_r = c_l;
	dsr->rho_l = c_r;

	for(ax = 0; ax < hp->dim; ++ax) {
		c_l = dsl->v[ax] + 0.5*minmod(theta*(dsl->v[ax]-dsll->v[ax]), 0.5*(dsr->v[ax]-dsll->v[ax]), theta*(dsr->v[ax]-dsl->v[ax]));
		c_r = dsr->v[ax] - 0.5*minmod(theta*(dsr->v[ax]-dsl->v[ax]), 0.5*(dsrr->v[ax]-dsl->v[ax]), theta*(dsrr->v[ax]-dsr->v[ax]));
		dsl->v[ax] = c_l;
		dsr->v[ax] = c_r;
	}

	c_l = dsl->p_c + 0.5*minmod(theta*(dsl->p_c-dsll->p_c), 0.5*(dsr->p_c-dsll->p_c), theta*(dsr->p_c-dsl->p_c));
	c_r = dsr->p_c - 0.5*minmod(theta*(dsr->p_c-dsl->p_c), 0.5*(dsrr->p_c-dsl->p_c), theta*(dsrr->p_c-dsr->p_c));
	dsl->p_r = c_l;
	dsr->p_l = c_r;

#else

	// My special sauce!
	// monotized central difference limiter?

	// TODO: do we need to put the dx's back in?
	// density
	grad_l = 2.0*(dsl->rho_r-dsl->rho_c);
	grad_r = 2.0*(dsr->rho_c-dsr->rho_l);
	grad_max = (dsr->rho_c - dsl->rho_c);
	if(grad_l < grad_max)
		dsl->rho_r = dsl->rho_c+0.5*grad_max;
	if(grad_r > grad_max)
		dsr->rho_l = dsr->rho_c-0.5*grad_max;

	// velocity is first-order, so this does nothing for now.
	for(ax = 0; ax < hp->dim; ++ax) {
		grad_l = 2.0*(dsl->v[ax]-dsl->v[ax]);
		grad_r = 2.0*(dsr->v[ax]-dsr->v[ax]);
		grad_max = -(dsr->v[ax] - dsl->v[ax]);
		if(grad_l < grad_max)
			dsl->v[ax] = dsl->v[ax]+0.5*grad_max;
		if(grad_r > grad_max)
			dsr->v[ax] = dsr->v[ax]-0.5*grad_max;
	}

	// pressure
	grad_l = 2.0*(dsl->p_r-dsl->p_c);
	grad_r = 2.0*(dsr->p_c-dsr->p_l);
	grad_max = (dsr->p_c - dsl->p_c);
	if(grad_l < grad_max)
		dsl->p_r = dsl->p_c+0.5*grad_max;
	if(grad_r > grad_max)
		dsr->p_l = dsr->p_c-0.5*grad_max;
#endif
}

void update_skewer(dvec skewer, int iterax, real fluxfac, hydro_problem* hp) {

	// get a local copy of the grid skewer 
	// and solve for local primitives
	int i, ax;
	long flatind;
	hydro_derived_state state_skewer[hp->nx[iterax]+2*STENCIL_SIZE];
	for(i = 0; i < hp->nx[iterax]; ++i) {
		skewer[iterax] = i;
		state_skewer[ZR(i)] = get_derived_state(skewer, iterax, hp);
	}

	// enforce boundary conditions here
	if(hp->bctype == BC_WALL) 
	for(i = 0; i < STENCIL_SIZE; ++i) {
		state_skewer[ZL(-i)] = state_skewer[ZR(i)];
		state_skewer[ZL(-i)].v[iterax] *= -1;

		//double tmp;
		//tmp = state_skewer[ZL(-i)].rho_r;
		//state_skewer[ZL(-i)].rho_r = state_skewer[ZL(-i)].rho_l;
		//state_skewer[ZL(-i)].rho_l = tmp; 

		state_skewer[ZR(hp->nx[iterax]+i)] = state_skewer[ZL(hp->nx[iterax]-i)];
		state_skewer[ZR(hp->nx[iterax]+i)].v[iterax] *= -1;

		//tmp = state_skewer[ZR(hp->nx[iterax]+i)].rho_r;
		//state_skewer[ZR(hp->nx[iterax]+i)].rho_r = state_skewer[ZR(hp->nx[iterax]+i)].rho_l;
		//state_skewer[ZR(hp->nx[iterax]+i)].rho_l = tmp; 

	}
	else if(hp->bctype == BC_PERIODIC) 
	for(i = 0; i < STENCIL_SIZE; ++i) {
		state_skewer[ZL(-i)] = state_skewer[ZL(hp->nx[iterax]-i)];
		state_skewer[ZR(hp->nx[iterax]+i)] = state_skewer[ZR(i)];
	}
	else if(hp->bctype == BC_FREE) 
	for(i = 0; i < STENCIL_SIZE; ++i) {
		state_skewer[ZL(-i)] = state_skewer[ZR(i)];
		state_skewer[ZR(hp->nx[iterax]+i)] = state_skewer[ZL(hp->nx[iterax]-i)];
	}
	
	// we have skewer primitives in hand, now solve for fluxes and step forward in time 
	// forward Euler
	hydro_vector flux;
	for(i = 0; i <= hp->nx[iterax]; ++i) {
		// solve the Riemann problem
		
		// get the left and right-hand cells
		hydro_derived_state dsll = state_skewer[ZL(i-1)];
		hydro_derived_state dsl = state_skewer[ZL(i)];
		hydro_derived_state dsr = state_skewer[ZR(i)];
		hydro_derived_state dsrr = state_skewer[ZR(i+1)];

		// apply the flux limiter
		//psi_flux_limiter(&dsll, &dsl, &dsr, &dsrr, hp);

		flux = flux_solve_hll(dsl, dsr, iterax, hp);

		// update both grid cells
		skewer[iterax] = i-1;
		flatind = flat_index(skewer, hp);
		hp->grid[flatind].rho -= fluxfac*flux.rho;
		for(ax = 0; ax < hp->dim; ++ax)
			hp->grid[flatind].mom[ax] -= fluxfac*flux.mom[ax];
		hp->grid[flatind].etot -= fluxfac*flux.etot;
		hp->grid[flatind].com[iterax] -= 0.5*hp->dx*fluxfac*flux.rho;
		//hp->grid[flatind].com[iterax] -= fluxfac*flux.com[iterax];

		skewer[iterax] = i;
		flatind = flat_index(skewer, hp);
		hp->grid[flatind].rho += fluxfac*flux.rho;
		for(ax = 0; ax < hp->dim; ++ax)
			hp->grid[flatind].mom[ax] += fluxfac*flux.mom[ax];
		hp->grid[flatind].etot += fluxfac*flux.etot;
		hp->grid[flatind].com[iterax] += -0.5*hp->dx*fluxfac*flux.rho;
		//hp->grid[flatind].com[iterax] += fluxfac*flux.com[iterax];
	}
}

void evolve(real tstop, int max_steps, int output_every, hydro_problem* hp) {

	int i, j, iterax, done;
	dvec skewer;
	real dt, alpha_max;
  
	// loop until t > tmax
	tstop += hp->time;
	max_steps += hp->step;
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
		if(hp->output_callback && hp->step%output_every == 0)
			hp->output_callback(hp);

		// Solve the approximate Riemann problem for fluxes
		// Strang splitting for directional update operators
		
		for(i = 0; i < 5; ++i)
			update_radiation(0.1*dt, hp);

		// TODO: replace all of this skewer nonsense with a general grid!!
		// solve in all three axis directsions for 0.5*dt, 
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

		// source term update sandwiched in the middle
		update_source_terms(dt, hp);

		// solve in all three axis directsions for 0.5*dt, 
		// but in the reverse order
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

		for(i = 0; i < 5; ++i)
			update_radiation(0.1*dt, hp);

	}
	if(hp->output_callback)
		hp->output_callback(hp);
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


