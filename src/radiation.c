/*************************************************************
 *
 *		radiation.c
 *
 *		Devon Powell
 *		28 December 2016
 *
 *		Beam-traced radiation hydrodynamics	
 *
 *************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "common.h"

void update_radiation(real dt, hydro_problem* hp) {

	int i;
	int r;
	real rmin, rmax;

	//printf("Updating radiation!\n");

	if(hp->nrays > 1024) printf("Error! Overflowed ray buffer.\n");

	// STEP 1: source all of the rays
	for(r = 0; r < 32; ++r) {
		hydro_ray ray;
		ray.direction = r;
		ray.time = 0.0;
		for(i = 0; i < hp->dim; ++i)
			ray.origin[i] = hp->dx*((hp->nx[i]/2)+0.5);
		if(hp->nrays < 1024)
			hp->rays[hp->nrays++] = ray;
	}


	for(r = 0; r < hp->nrays; ++r) {
		hydro_ray ray = hp->rays[r];

		rmin = CLIGHT*ray.time; 
		rmax = CLIGHT*(ray.time+dt); 
	
		//printf("Ray %d:\n", r);
		//printf(" origin = %f %f\n", ray.origin[0], ray.origin[1]);
		//printf(" time = %f, dt = %f\n", ray.time, dt);
		//printf(" rmin = %f, rmax = %f\n", rmin, rmax);
		//printf(" energy = %f\n", ray.energy);
	
	
		ray.energy *= 0.9;
		ray.dt = dt;
		ray.time += dt;
		hp->rays[r] = ray;
	}

	// TODO: compact rays if they have been removed!

}
