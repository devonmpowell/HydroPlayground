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

#include "csparse.h"


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

		skewer[iterax] = i;
		flatind = flat_index(skewer, hp);
		hp->grid[flatind].rho += fluxfac*flux.rho;
		for(ax = 0; ax < hp->dim; ++ax)
			hp->grid[flatind].mom[ax] += fluxfac*flux.mom[ax];
		hp->grid[flatind].etot += fluxfac*flux.etot;
		hp->grid[flatind].com[iterax] += -0.5*hp->dx*fluxfac*flux.rho;
	}
}

void evolve(real tstop, int max_steps, int output_every, hydro_problem* hp) {

	int i, j, iterax, done;
	dvec skewer;
	real alpha_max;
  
	// loop until t > tmax
	tstop += hp->time;
	max_steps += hp->step;
	for(done = 0; hp->step < max_steps && !done; hp->step++, hp->time += hp->dt) {

		// get the max. characteristic speed on the grid
		// time step based on CFL condition
		alpha_max = get_max_char_speed(hp);
		hp->dt = hp->cfl_fac*hp->dx/alpha_max;

		//hp->dt = 0.0001;

		if(hp->time + hp->dt > tstop) {
			hp->dt = tstop - hp->time;
			done = 1;
		}

		// output if called for
		if(hp->output_callback && hp->step%output_every == 0)
			hp->output_callback(hp);

#if 1
		// Solve the approximate Riemann problem for fluxes
		// solve in all three axis directsions for 0.5*dt, 
		// TODO: replace Strang splitting with a general
		real stfac = 0.5*hp->dt/hp->dx; // half-step factor
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
		update_source_terms(hp->dt, hp);
#endif

		// update radiation, using just enough steps to avoid
		// violating the CFL condition for radiation
		real dtrad = hp->cfl_fac*hp->dx/CLIGHT;

		dtrad *= 10;


		int nrstep = ceil(hp->dt/dtrad);
		dtrad = hp->dt/nrstep;
		for(i = 0; i < nrstep; ++i)
			update_radiation(dtrad, hp);

		// update the Lagrange mesh
		//update_lagrange_mesh(dt, hp);

#if 1
		// solve in all three axis directions for 0.5*dt, 
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
#endif

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

#define FUZZ 1.0e-10
void update_lagrange_mesh(real dt, hydro_problem* hp) {

	int iter, t, i, c, v, ax;
	real vtot, vol, mtest, err, max_err;
	rvec dv[4];
	rvec tpos, tvel;
	hydro_lagrange_vert vert;
	hydro_lagrange_tet tet;

	real* tperv = malloc(hp->nverts*sizeof(real));
	real* drho = malloc(hp->nverts*sizeof(real));

	for(v = 0; v < hp->nverts; ++v) {
		if(hp->mesh_verts[v].pos[ax] != 0.0 && hp->mesh_verts[v].pos[ax] != 1.0)
			hp->mesh_verts[v].pos[ax] += dt*hp->mesh_verts[v].vel[ax];
	}



	max_err = 1000.0;
	real err_old = 2.0;
	for(iter = 0; max_err > FUZZ && iter < 1000; ++iter) {

		max_err = 0.0;

		memset(tperv, 0, hp->nverts*sizeof(real));
		memset(drho, 0, hp->nverts*sizeof(real));
		for(t = 0; t < hp->ntets; ++t) {
			
			// first, get the tet and update its volume
			tet = hp->mesh_tets[t];
			for(v = 0; v < hp->dim+1; ++v)
				for(i = 0; i < hp->dim; ++i)
					dv[v][i] = hp->mesh_verts[tet.verts[v]].pos[i];
			vol = 0.5*((dv[1][0]-dv[0][0])*(dv[2][1]-dv[0][1])-(dv[2][0]-dv[0][0])*(dv[1][1]-dv[0][1]));
			real vfac = vol/(hp->dim+1);

			// TODO: figure out what the relevant step factor for the density is!!!
			// get the mass of the tet as given by its current volume and densities
			// compute the error and relax the mesh towards the solution
			for(v = 0, mtest = 0.0; v < hp->dim+1; ++v)
				mtest += vfac*hp->mesh_verts[tet.verts[v]].rho;
			err = tet.mass - mtest;
			if(fabs(err/tet.mass) > max_err) 
				max_err = fabs(err/(tet.mass+FUZZ));
			for(v = 0; v < hp->dim+1; ++v) {
				drho[tet.verts[v]] += err/vfac; 
				tperv[tet.verts[v]] += 1; 
			}
		}
		for(v = 0; v < hp->nverts; ++v) {
			drho[v] /= tperv[v];
			hp->mesh_verts[v].rho += 0.1*drho[v];
		}


#if 0
		for(ax = 0; ax < 2; ++ax) {

			memset(tperv, 0, hp->nverts*sizeof(real));
			memset(drho, 0, hp->nverts*sizeof(real));
			for(t = 0; t < hp->ntets; ++t) {
				
				// first, get the tet and update its volume
				tet = hp->mesh_tets[t];

				
				for(v = 0; v < hp->dim+1; ++v)
					for(i = 0; i < hp->dim; ++i)
						dv[v][i] = hp->mesh_verts[tet.verts[v]].pos[i];
				vol = 0.5*((dv[1][0]-dv[0][0])*(dv[2][1]-dv[0][1])-(dv[2][0]-dv[0][0])*(dv[1][1]-dv[0][1]));
				real vfac = vol/(hp->dim+1);
	
				// TODO: figure out what the relevant step factor for the density is!!!
				// get the mass of the tet as given by its current volume and densities
				// compute the error and relax the mesh towards the solution
				mtest = 0.0;
				real sumrho = 0.0;
				real sumx = 0.0;
				for(v = 0, mtest = 0.0; v < hp->dim+1; ++v) {
					sumrho += hp->mesh_verts[tet.verts[v]].rho;
					sumx += hp->mesh_verts[tet.verts[v]].pos[ax];
					mtest += hp->mesh_verts[tet.verts[v]].rho*hp->mesh_verts[tet.verts[v]].pos[ax];
				}
				mtest = (vol/12.0)*(sumrho*sumx+mtest);
				err = tet.com[ax] - mtest;
				if(fabs(err/tet.mass) > max_err) 
					max_err = fabs(err/(tet.mass+FUZZ));

				// update the gradient
				for(v = 0; v < hp->dim+1; ++v) {
					drho[tet.verts[v]] += 12.0*err/(vol*(sumrho+hp->mesh_verts[tet.verts[v]].rho)); 
					tperv[tet.verts[v]] += 1;
				}
			}
			for(v = 0; v < hp->nverts; ++v) {
				drho[v] /= tperv[v];
				if(hp->mesh_verts[v].pos[ax] != 0.0 && hp->mesh_verts[v].pos[ax] != 1.0)
					hp->mesh_verts[v].pos[ax] += 0.0001*drho[v];
			}
		
		}

		for(ax = 0; ax < 2; ++ax) {

			memset(tperv, 0, hp->nverts*sizeof(real));
			memset(drho, 0, hp->nverts*sizeof(real));
			for(t = 0; t < hp->ntets; ++t) {
				
				// first, get the tet and update its volume
				tet = hp->mesh_tets[t];

				
				for(v = 0; v < hp->dim+1; ++v)
					for(i = 0; i < hp->dim; ++i)
						dv[v][i] = hp->mesh_verts[tet.verts[v]].vel[i];
				vol = 0.5*((dv[1][0]-dv[0][0])*(dv[2][1]-dv[0][1])-(dv[2][0]-dv[0][0])*(dv[1][1]-dv[0][1]));
				real vfac = vol/(hp->dim+1);
	
				// TODO: figure out what the relevant step factor for the density is!!!
				// get the mass of the tet as given by its current volume and densities
				// mompute the error and relax the mesh towards the solution
				mtest = 0.0;
				real sumrho = 0.0;
				real sumx = 0.0;
				for(v = 0, mtest = 0.0; v < hp->dim+1; ++v) {
					sumrho += hp->mesh_verts[tet.verts[v]].rho;
					sumx += hp->mesh_verts[tet.verts[v]].vel[ax];
					mtest += hp->mesh_verts[tet.verts[v]].rho*hp->mesh_verts[tet.verts[v]].vel[ax];
				}
				mtest = (vol/12.0)*(sumrho*sumx+mtest);
				err = tet.mom[ax] - mtest;
				if(fabs(err/tet.mass) > max_err) 
					max_err = fabs(err/(tet.mass+FUZZ));

				// update the gradient
				for(v = 0; v < hp->dim+1; ++v) {
					drho[tet.verts[v]] += 12.0*err/(vol*(sumrho+hp->mesh_verts[tet.verts[v]].rho)); 
					tperv[tet.verts[v]] += 1;
				}
			}
			for(v = 0; v < hp->nverts; ++v) {
				drho[v] /= tperv[v];
				//if(hp->mesh_verts[v].pos[ax] != 0.0 && hp->mesh_verts[v].pos[ax] != 1.0)
					//hp->mesh_verts[v].vel[ax] += 0.001*drho[v];
			}
		
		}
#endif



		if(max_err > err_old) {
			printf(" --> Stopped converging at iter %d !!\n", iter);
			break;
		} 
		err_old = max_err; 
	}

					//exit(0);

	real mtot = 0.0;

	// use the vertex-centered values to the cell-centered moments
	for(t = 0; t < hp->ntets; ++t) {
		tet = hp->mesh_tets[t];

		// TODO: make a 2D orient function here...
		for(v = 0; v < hp->dim+1; ++v)
			for(i = 0; i < hp->dim; ++i)
				dv[v][i] = hp->mesh_verts[tet.verts[v]].pos[i];
		vol = 0.5*((dv[1][0]-dv[0][0])*(dv[2][1]-dv[0][1])-(dv[2][0]-dv[0][0])*(dv[1][1]-dv[0][1]));

		// initialize mass and moments
		for(v = 0; v < hp->dim+1; ++v) 
			mtot += 0.3333333333*vol*hp->mesh_verts[tet.verts[v]].rho;

	}
	printf("Mtot = %f\n", mtot);



	free(tperv);
	free(drho);

	printf(" gradient descent-ish density solve complete, iter = %d, max err = %.5e\n", iter, max_err);


	// next, update all fluxes through the mesh faces, and update the cell-centered quantities
	
	// source terms, just momentum for now 
	for(t = 0; t < hp->ntets; ++t) {
		for(i = 0; i < 3; ++i)
			hp->mesh_tets[t].com[i] += dt*hp->mesh_tets[t].mom[i];
	}
	
}

int init(hydro_problem* hp) {

	// TODO: this is the only thing to get into Python!!!
	hp->eos_p = eos_ideal_p;

	// initialize the Lagrange mesh
	
	int i, c, t, jj, kk, v, iter, tind;
	real vol, vtot, mtest, err, sumrho, sumx;
	rvec dv[4];
	hydro_lagrange_vert vert;
	hydro_lagrange_tet tet;

	// use the vertex-centered values to the cell-centered moments
	for(t = 0; t < hp->ntets; ++t) {
		tet = hp->mesh_tets[t];

		// TODO: make a 2D orient function here...
		for(v = 0; v < hp->dim+1; ++v)
			for(i = 0; i < hp->dim; ++i)
				dv[v][i] = hp->mesh_verts[tet.verts[v]].pos[i];
		vol = 0.5*((dv[1][0]-dv[0][0])*(dv[2][1]-dv[0][1])-(dv[2][0]-dv[0][0])*(dv[1][1]-dv[0][1]));

		// initialize mass and moments
		tet.mass = 0.0;
		for(v = 0; v < hp->dim+1; ++v) 
			tet.mass += 0.3333333333*vol*hp->mesh_verts[tet.verts[v]].rho;

		for(i = 0; i < hp->dim; ++i) {

			// cell COM
			sumrho = 0.0, sumx = 0.0;
			for(v = 0; v < hp->dim+1; ++v) {
				sumrho += hp->mesh_verts[tet.verts[v]].rho;
				sumx += hp->mesh_verts[tet.verts[v]].pos[i];
				tet.com[i] += hp->mesh_verts[tet.verts[v]].rho*hp->mesh_verts[tet.verts[v]].pos[i];
			}
			tet.com[i] = (vol/12.0)*(sumrho*sumx+tet.com[i]);

			// Cell momentum
			sumrho = 0.0, sumx = 0.0;
			for(v = 0; v < hp->dim+1; ++v) {
				sumrho += hp->mesh_verts[tet.verts[v]].rho;
				sumx += hp->mesh_verts[tet.verts[v]].vel[i];
				tet.mom[i] += hp->mesh_verts[tet.verts[v]].rho*hp->mesh_verts[tet.verts[v]].vel[i];
			}
			tet.mom[i] = (vol/12.0)*(sumrho*sumx+tet.mom[i]);

		}
		// TODO: this only works with zero kinetic energy
		for(v = 0; v < hp->dim+1; ++v) 
			tet.etot += 0.3333333333*tet.mass*hp->mesh_verts[tet.verts[v]].etot;
		hp->mesh_tets[t] = tet;
	}
	return 1;
}


