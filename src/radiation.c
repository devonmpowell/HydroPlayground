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

// macros
#define ERRTOL 1.0e-6
#define wav3(va, wa, vb, wb, vr) {			\
	vr[0] = (wa*va[0] + wb*vb[0])/(wa + wb);	\
	vr[1] = (wa*va[1] + wb*vb[1])/(wa + wb);	\
	vr[2] = (wa*va[2] + wb*vb[2])/(wa + wb);	\
}
#define dot3(va, vb) (va[0]*vb[0] + va[1]*vb[1] + va[2]*vb[2])
#define cross3(v,u,w)          /* CROSS Vector Product */              \
{                                                                       \
    (v)[0] = (u)[1]*(w)[2] - (u)[2]*(w)[1];                             \
    (v)[1] = (u)[2]*(w)[0] - (u)[0]*(w)[2];                             \
    (v)[2] = (u)[0]*(w)[1] - (u)[1]*(w)[0];                             \
}

// forward declarations
typedef struct {
	struct {
		rvec pos;
		dvec pnbrs;
	} verts[64];
	int nverts;
} solid_angle_poly;
real domega3(rvec v1, rvec v2, rvec v3); 
void clip_beam_poly(solid_angle_poly* poly, rvec* clipnorms, int nclip);
void reduce_beam_poly(solid_angle_poly* poly, real* omega, rvec* centroid);
int sort_helper(const void *a, const void *b);

// functions

void update_radiation(real dt, hydro_problem* hp) {

	// all variables up top
	int v, r, f, ax, ornrays, flatind, face_id,
		cell_id, nsources, b, source_id, nstack;
	dvec grind;
	real sum, omega, flux, omega_in, omega_out, err, r_in, r_out,
		 flux_in, flux_out;
	rvec centroid, clipnorms[4], cubeverts[8];
	hydro_ray ray, rtmp;
	solid_angle_poly tpoly, curpoly;
	rad_source cs;
	const static int neighbor_face_indices[6][4] = {
		{3,7,5,1},{4,6,2,0},{6,7,3,2},{1,5,4,0},{5,7,6,4},{2,3,1,0}};
	const static int neighbor_grid_offsets[6][3] = {
		{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1}};
	struct {
		solid_angle_poly poly;
		dvec grind;
		int fid;
		real f, o;
	} beam_stack[128];

	// Add a single test source
	// TODO: hacky for now
	rad_source sources[10];
	nsources = 0;
	for(ax = 0; ax < 3; ++ax)
		sources[nsources].pos[ax] = hp->dx*(hp->nx[ax]+1)/2; //+ax+0.342)/2; 
	//sources[nsources].pos[1] -= 0.2;
	sources[nsources].dndt = 1.0;
	nsources++;


	// STEP 1: source new rays 
	for(source_id = 0; source_id < nsources; ++source_id) {

		// add new rays to the stack, one per face exiting the source cell
		// interleave bits to get the unique cell ID a la Morton order
		// TODO: only for true sources, not supersources!
		for(ax = 0; ax < 3; ++ax)
			grind[ax] = floor(sources[source_id].pos[ax]/hp->dx);
		for(cell_id = 0, b = 0; b < 24; ++b)
			cell_id |= ((grind[b%3]>>(b/3))&1)<<b;
		ray.source_id = source_id; 
		ray.face_id = ((cell_id&0xFFFFFF)<<4)|(7&0xF);
		ray.r = 0.0;
		ray.f = dt*sources[source_id].dndt;
		hp->rays[hp->nrays++] = ray;
		if(hp->nrays > 160000) printf("Error! Overflowed ray buffer.\n");
	}

	// Step 3: propagate all rays forward by c*dt
	memset(hp->rad_grid, 0, (hp->nx[0]+4)*(hp->nx[1]+4)*(hp->nx[2]+4)*sizeof(rad_vector));
	ornrays = hp->nrays;
	for(r = 0; r < ornrays; ++r) {
		ray = hp->rays[r];
		hp->rays[r].f = 0.0;
		if(ray.f <= 0.0) continue;

		// back out basic geometric information
		source_id = ray.source_id;
		face_id = ray.face_id&0xF;
		cell_id = (ray.face_id>>4);
		cs = sources[source_id];
		ray.r += CLIGHT*dt; 
		memset(grind, 0, sizeof(grind));
		for(b = 0; b < 24; ++b)
			grind[b%3] |= ((cell_id>>b)&1)<<(b/3);
		for(v = 0; v < 8; ++v) {
			cubeverts[v][0] = hp->dx*(grind[0]+((v>>0)&1)) - cs.pos[0]; 
			cubeverts[v][1] = hp->dx*(grind[1]+((v>>1)&1)) - cs.pos[1]; 
			cubeverts[v][2] = hp->dx*(grind[2]+((v>>2)&1)) - cs.pos[2]; 
		}

		//printf("Ray %d, id = %d, radius = %.5e\n", r, ray.face_id, ray.r);
		//printf("Starting grid index = %d %d %d\n", grind[0], grind[1], grind[2]);
		//printf("cell_id = %d, face_id = %d, source_id = %d\n", cell_id, face_id, source_id);

		// initialize the stack
		nstack = 0;
		if(face_id == 7) {
				
			// construct the ray poly as it exits each face of the source cell
			for(f = 0, sum = 0.0; f < 6; ++f) {
				curpoly.nverts = 4;
				for(v = 0; v < 4; ++v) {
					for(ax = 0; ax < 3; ++ax)
						curpoly.verts[v].pos[ax] = cubeverts[neighbor_face_indices[f][v]][ax];
					curpoly.verts[v].pnbrs[0] = (v-1+4)%4;
					curpoly.verts[v].pnbrs[1] = (v+1)%4;
				}
				reduce_beam_poly(&curpoly, &omega, &centroid);
				sum += omega;

				// push this poly onto the stack, tagging it as entering the corresponding
				// neighboring grid cell
				// TODO: must propagate and attenuate!
				// TODO: This is the subgrid model...
				beam_stack[nstack].poly = curpoly;
				beam_stack[nstack].fid = f;
				for(ax = 0; ax < 3; ++ax)
					beam_stack[nstack].grind[ax] = grind[ax] + neighbor_grid_offsets[f][ax];
				beam_stack[nstack].f = ray.f*omega/(2*TWO_PI);
				beam_stack[nstack].o = omega;
				nstack++;
			}

			// check that the solid angle adds to 4*pi
			err = fabs(1.0-sum/(2*TWO_PI));
			if(err > ERRTOL) {
				printf(" Solid angle err for source cell is %.5e\n", err);
				exit(0);
			}
		}
		else {
			// construct the beam poly from the upwind cell face
			// push this poly onto the stack
			curpoly.nverts = 4;
			for(v = 0; v < 4; ++v) {
				for(ax = 0; ax < 3; ++ax)
					curpoly.verts[v].pos[ax] = cubeverts[neighbor_face_indices[face_id][v]][ax];
				curpoly.verts[v].pnbrs[1] = (v-1+4)%4;
				curpoly.verts[v].pnbrs[0] = (v+1)%4;
			}
			reduce_beam_poly(&curpoly, &omega, &centroid);
			beam_stack[nstack].poly = curpoly;
			beam_stack[nstack].fid = face_id^1;
			for(ax = 0; ax < 3; ++ax)
				beam_stack[nstack].grind[ax] = grind[ax];
			beam_stack[nstack].f = ray.f;
			beam_stack[nstack].o = omega; 
			nstack++;
		}

		// now that we have initialized our beam polygon,
		// propagate it through the grid, recursing as needed
		while(nstack) {

			// pop the stack
			// skip if we have no flux or solid angle to work with
			nstack--;
			curpoly = beam_stack[nstack].poly;
			face_id = beam_stack[nstack].fid^1; // flip the face index around to outgoing status 
			flux = beam_stack[nstack].f;
			omega_in = beam_stack[nstack].o;
			memcpy(grind, beam_stack[nstack].grind, sizeof(grind));
			if(flux <= 0.0) continue;
			if(omega_in <= 0.0) continue;

			// check the grid cell vertices to ensure that we propagate
			// the entire way through during this one time-step
			// Otherwise pause the fragment as a new beam to be merged
			for(v = 0; v < 8; ++v) {
				for(ax = 0; ax < 3; ++ax)
					cubeverts[v][ax] = hp->dx*(grind[ax]+((v>>ax)&1))-cs.pos[ax]; 
				if(ray.r*ray.r < cubeverts[v][0]*cubeverts[v][0]
						+cubeverts[v][1]*cubeverts[v][1]+cubeverts[v][2]*cubeverts[v][2]) {
					rtmp = ray;
					for(cell_id = 0, b = 0; b < 24; ++b)
						cell_id |= ((grind[b%3]>>(b/3))&1)<<b;
					rtmp.face_id = ((cell_id&0xFFFFFF)<<4)|(face_id&0xF);
					rtmp.source_id = source_id; 
					rtmp.f = flux;
					hp->rays[hp->nrays++] = rtmp;
					goto next_fragment;
				}
			}

			// clip against each of the new downwind faces, tracking fractions 
			// of the solid angle through from the upwind face
			// attentuate exponentially via the mean optical depth
			for(f = 0, sum = 0.0; f < 6; ++f) {
				if(f == face_id) continue;
				for(v = 0; v < 4; ++v)
					cross3(clipnorms[v], cubeverts[neighbor_face_indices[f][v]], cubeverts[neighbor_face_indices[f][(v+1)%4]]);
				tpoly = curpoly;
				clip_beam_poly(&tpoly, clipnorms, 4);
				reduce_beam_poly(&tpoly, &omega_out, &centroid);
				if(omega_out > 0.0) {
					sum += omega_out;
	
					// get the mean optical depth
					// get the correct inner and outer radii, then attenuate exponentially 
					// TODO: cleaner logic expressions here
					// TODO: more stable numerical calculation too?
					// TODO: Use the cube vertices??
					r_in = (hp->dx*(grind[face_id/2]+(1-face_id%2))-cs.pos[face_id/2])/centroid[face_id/2];
					r_out = (hp->dx*(grind[f/2]+(1-f%2))-cs.pos[f/2])/centroid[f/2];

					// update the energy density on the grid
					real k = 100.0;
					flatind = flat_index(grind, hp);
					flux_in = flux*omega_out/omega_in;
					flux_out = flux_in*exp(-k*hp->grid[flatind].rho*(r_out-r_in)); 
					hp->rad_grid[flatind].E += 0.5*(flux_in+flux_out)*(r_out-r_in)/(CLIGHT*dt*hp->dx*hp->dx);

					// push this poly onto the stack, tagging it as entering the corresponding
					// neighboring grid cell
					beam_stack[nstack].poly = tpoly;
					beam_stack[nstack].fid = f;
					for(ax = 0; ax < 3; ++ax)
						beam_stack[nstack].grind[ax] = grind[ax] + neighbor_grid_offsets[f][ax];
					beam_stack[nstack].f = flux_out; 
					beam_stack[nstack].o = omega_out; 
					nstack++;			
				}
			}

			// error checking
			err = fabs(1.0-sum/omega_in);
			if(err > ERRTOL) {
				//printf("  Solid angle err, sum of all faces = %.5e, omega_in= %.5e, err = %.5e\n", sum, omega_in, err);
				//exit(0);
			}
			if(nstack > 64) {
				printf("  Stack overflow! %d\n", nstack);
				exit(0);
			}
			next_fragment: continue;
		}
	}
		
	// compress the ray list down, removing rays with no flux
	// then sort the rays by their unique face index, then loop through
	// merge all rays in the face index from the same source
	// TODO: do the sort and merge in one step?
	ornrays = hp->nrays;
	hp->nrays = 0;
	for(r = 0; r < ornrays; ++r) {
		ray = hp->rays[r];
		if(ray.f <= 0.0) continue;
		hp->rays[hp->nrays++] = ray;
	}
	printf("NUM RAYS = %d (pre-sort)\n", hp->nrays);
	qsort(hp->rays, hp->nrays, sizeof(hydro_ray), sort_helper);

	// finally, merge rays passing through the same face that come from the same source
	// TODO: make this a whole source-mergeing thing, with a tree traversal, etc.
	ornrays = hp->nrays;
	hp->nrays = 0;
	ray.face_id = -1; 
	for(r = 0; r < ornrays; ++r) {
		if(ray.face_id == hp->rays[r].face_id) {
			
			// flux-weight the combined radius from the source
			// TODO: this is also hacky
			ray.f += hp->rays[r].f;
			ray.r += hp->rays[r].r*hp->rays[r].f;
		}
		else {
			if(ray.face_id > 0) {
				ray.r /= ray.f;
				hp->rays[hp->nrays++] = ray;
			}
			ray.face_id = hp->rays[r].face_id;
			ray.f = hp->rays[r].f; 
			ray.r = hp->rays[r].f*hp->rays[r].r; 
		}
	}
	printf("NUM RAYS = %d (post-sort)\n", hp->nrays);
}


void reduce_beam_poly(solid_angle_poly* poly, real* omega, rvec* centroid) {

	int v;
	real len;
	rvec meanvec;

	// compute the final solid angle for this poly
	*omega = 0.0;
	if(!poly->nverts) return; 
	memset(meanvec, 0, sizeof(meanvec));
	for(v = 0; v < poly->nverts; ++v) {
		len = sqrt(poly->verts[v].pos[0]*poly->verts[v].pos[0]
			+poly->verts[v].pos[1]*poly->verts[v].pos[1]+poly->verts[v].pos[2]*poly->verts[v].pos[2]);
		poly->verts[v].pos[0] /= len;
		poly->verts[v].pos[1] /= len;
		poly->verts[v].pos[2] /= len;
		meanvec[0] += poly->verts[v].pos[0];
		meanvec[1] += poly->verts[v].pos[1];
		meanvec[2] += poly->verts[v].pos[2];
	}
	len = sqrt(meanvec[0]*meanvec[0]+meanvec[1]*meanvec[1]+meanvec[2]*meanvec[2]);
	meanvec[0] /= len;
	meanvec[1] /= len;
	meanvec[2] /= len;

	(*centroid)[0] = 0;
	(*centroid)[1] = 0;
	(*centroid)[2] = 0;
	for(v = 0; v < poly->nverts; ++v) {
		real mydo = domega3(poly->verts[v].pos, poly->verts[poly->verts[v].pnbrs[1]].pos, meanvec); 
		*omega += mydo; 
		(*centroid)[0] += mydo*(poly->verts[v].pos[0] +  poly->verts[poly->verts[v].pnbrs[1]].pos[0] +  meanvec[0]);
		(*centroid)[1] += mydo*(poly->verts[v].pos[1] +  poly->verts[poly->verts[v].pnbrs[1]].pos[1] +  meanvec[1]);
		(*centroid)[2] += mydo*(poly->verts[v].pos[2] +  poly->verts[poly->verts[v].pnbrs[1]].pos[2] +  meanvec[2]);
	}
	len = sqrt((*centroid)[0]*(*centroid)[0]+(*centroid)[1]*(*centroid)[1]+(*centroid)[2]*(*centroid)[2]);
	(*centroid)[0] /= len;
	(*centroid)[1] /= len;
	(*centroid)[2] /= len;
}


void clip_beam_poly(solid_angle_poly* poly, rvec* clipnorms, int nclip) {


	int f;
	real sdists[32];
	int clipped[32];

	// variable declarations
	int v, nclipped, np, onv, vstart, vcur, vnext, numunclipped; 

	// direct access to vertex buffer
	if(poly->nverts <= 0) return;
	int* nverts = &poly->nverts;
	for(f = 0; f < nclip; ++f) {

		// split vertices by their distance from the clip plane
		nclipped = 0;
		memset(&clipped, 0, sizeof(clipped));
		for(v = 0; v < poly->nverts; ++v) {
			sdists[v] = dot3(poly->verts[v].pos, clipnorms[f]);
			clipped[v] = (sdists[v] < 0.0);
			nclipped += clipped[v]; 
		}

		// skip this face if the poly lies entirely on one side of it 
		if(nclipped == 0) continue;
		if(nclipped == poly->nverts) {
			poly->nverts = 0;
			return;
		}

		// check all edges and insert new vertices on the bisected edges 
		onv = *nverts;
		for(vcur = 0; vcur < onv; ++vcur) {
			if(clipped[vcur]) continue;
			for(np = 0; np < 2; ++np) {
				vnext = poly->verts[vcur].pnbrs[np];
				if(!clipped[vnext]) continue;
				poly->verts[*nverts].pnbrs[1-np] = vcur;
				poly->verts[*nverts].pnbrs[np] = -1;
				poly->verts[vcur].pnbrs[np] = *nverts;
				wav3(poly->verts[vcur].pos, -sdists[vnext],
					poly->verts[vnext].pos, sdists[vcur],
					poly->verts[*nverts].pos);
				(*nverts)++;
			}
		}

		// for each new vert, search around the poly for its new neighbors
		// and doubly-link everything
		for(vstart = onv; vstart < *nverts; ++vstart) {
			if(poly->verts[vstart].pnbrs[1] >= 0) continue;
			vcur = poly->verts[vstart].pnbrs[0];
			do {
				vcur = poly->verts[vcur].pnbrs[0]; 
			} while(vcur < onv);
			poly->verts[vstart].pnbrs[1] = vcur;
			poly->verts[vcur].pnbrs[0] = vstart;
		}

		// go through and compress the vertex list, removing clipped verts
		// and re-indexing accordingly (reusing `clipped` to re-index everything)
		numunclipped = 0;
		for(v = 0; v < *nverts; ++v) {
			if(!clipped[v]) {
				poly->verts[numunclipped] = poly->verts[v];
				clipped[v] = numunclipped++;
			}
		}
		*nverts = numunclipped;
		for(v = 0; v < *nverts; ++v) {
			poly->verts[v].pnbrs[0] = clipped[poly->verts[v].pnbrs[0]];
			poly->verts[v].pnbrs[1] = clipped[poly->verts[v].pnbrs[1]];
		}	
	}
}

real domega3(rvec v1, rvec v2, rvec v3) {

	// Solid angle of a triangle by Oosterom and Strackee
	// Assumes v1, v2, v3 are unit vectors
	int i;
	real det, div;
	det = v1[0]*(v2[1]*v3[2]-v3[1]*v2[2])-v1[1]*(v2[0]*v3[2]-v3[0]*v2[2])+v1[2]*(v2[0]*v3[1]-v3[0]*v2[1]);
	div = 1.0;
	for(i = 0; i < 3; ++i) {
		div += v1[i]*v2[i];
		div += v2[i]*v3[i];
		div += v3[i]*v1[i];
	}
	return 2.0*atan2(det, div);
} 

int sort_helper(const void *a, const void *b) {
  if (((hydro_ray*)a)->face_id <  ((hydro_ray*)b)->face_id) return -1;
  if (((hydro_ray*)a)->face_id == ((hydro_ray*)b)->face_id) return 0;
  if (((hydro_ray*)a)->face_id >  ((hydro_ray*)b)->face_id) return 1;
}

