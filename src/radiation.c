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

#define unpack_id_bits(idbits, baseid, reflvl, refbits) { \
	baseid = (ray.angle_id>>24)&0xFF; \
	reflvl = (ray.angle_id>>20)&0xF; \
	refbits =  (ray.angle_id)&0xFFFFF; \
}

#define dot2(va, vb) (va[0]*vb[0] + va[1]*vb[1])
#define cross2(va, vb) (va[0]*vb[1]-va[1]*vb[0])
#define wav2(va, wa, vb, wb, vr) {			\
	vr[0] = (wa*va[0] + wb*vb[0])/(wa + wb);	\
	vr[1] = (wa*va[1] + wb*vb[1])/(wa + wb);	\
}

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

#define ERRTOL 1.0e-6

#define domega2(v0, v1) asin(cross2(v0, v1))
real domega3(rvec v1, rvec v2, rvec v3) {
	int i;
	real det, div;
	// Solid angle of a triangle by Oosterom and Strackee
	// Assumes v1, v2, v3 are unit vectors
	det = v1[0]*(v2[1]*v3[2]-v3[1]*v2[2])-v1[1]*(v2[0]*v3[2]-v3[0]*v2[2])+v1[2]*(v2[0]*v3[1]-v3[0]*v2[1]);
	div = 1.0;
	for(i = 0; i < 3; ++i) {
		div += v1[i]*v2[i];
		div += v2[i]*v3[i];
		div += v3[i]*v1[i];
	}
	return 2.0*atan2(det, div);
} 

typedef struct {
	struct {
		rvec pos;
		dvec pnbrs;
	} verts[64];
	int nverts;
} solid_angle_poly;

void clip_beam_poly(solid_angle_poly* poly, rvec* clipnorms, int nclip) {


	int f;
	real sdists[32];
	int clipped[32];

	// variable declarations
	int v, p, nclipped, np, onv, vstart, vcur, vnext, numunclipped; 
	real len;

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


void reduce_beam_poly(solid_angle_poly* poly, real* omega) {

	int v;
	real len;

	*omega = 0.0;
	if(!poly->nverts) return; 

	// compute the final solid angle for this poly
	rvec meanvec;
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

	for(v = 0; v < poly->nverts; ++v) {
		real mydo = domega3(poly->verts[v].pos, poly->verts[poly->verts[v].pnbrs[1]].pos, meanvec); 
		*omega += mydo; 
	}
#if 0
		if(poly->nverts) { //} && faces->nclip) {

			// compute the final solid angle for this poly
			rvec meanvec;
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
				*domega += mydo; 
				(*centroid)[0] += mydo*(poly->verts[v].pos[0] +  poly->verts[poly->verts[v].pnbrs[1]].pos[0] +  meanvec[0]);
				(*centroid)[1] += mydo*(poly->verts[v].pos[1] +  poly->verts[poly->verts[v].pnbrs[1]].pos[1] +  meanvec[1]);
				(*centroid)[2] += mydo*(poly->verts[v].pos[2] +  poly->verts[poly->verts[v].pnbrs[1]].pos[2] +  meanvec[2]);
			}
			len = sqrt((*centroid)[0]*(*centroid)[0]+(*centroid)[1]*(*centroid)[1]+(*centroid)[2]*(*centroid)[2]);
			(*centroid)[0] /= len;
			(*centroid)[1] /= len;
			(*centroid)[2] /= len;
		}
	}


#endif
}



void update_radiation(real dt, hydro_problem* hp) {

	int i, v, r, s, ax, rminbits, rmaxbits, fin, fout, ornrays, refine, baseid, reflvl, refbits,
		flatind, ii, jj, kk, quadrant, qid, rlvl, f, nbase, split, face_id;
	dvec lsgn, grind, ibox[2], nbox;
	real ray_mom_out, fmom_in, domega_in, domega_out, intensity_in, dobase, 
		 ray_flux_out, err, allmin, allmax, r_in, r_out, secthmid, len, ray_omega, flux, 
		 cscthmid, inmin, inmax, outmin, outmax, fmid_in, fmid_out, omega, sum, omega_in, omega_out;
	real rmin2[4], rmax2[4], flux_out[4], fmom_out[4];
	rvec x0, rmid, tmpverts[8], ntmp, v0, v1;
	hydro_ray ray, r0, r1, rtmp;
	solid_angle_poly beam_poly, inpoly, outpoly;


	rad_source sources[10];
	rad_source cs;
	int nsources;

	nsources = 0;
	for(ax = 0; ax < 3; ++ax)
		sources[nsources].pos[ax] = hp->dx*(hp->nx[ax]+1)/2; //+ax+0.342)/2; 
	sources[nsources].dndt = 1.0;
	nsources++;



	// STEP 1: source new rays 
	// TODO: subgrid model??
	for(s = 0; s < nsources; ++s) {

		// add new rays to the stack, one per face exiting the source cell
		// TODO: only for true sources, not supersources!
		for(ax = 0; ax < 3; ++ax)
			ray.grind[ax] = floor(sources[s].pos[ax]/hp->dx);
		ray.source_id = s;
		ray.face_id = 7;
		ray.r = 0.0;
		ray.f = dt*sources[s].dndt;
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
		cs = sources[ray.source_id];
		ray.r += CLIGHT*dt; 
		//printf("Ray %d: ", r);
		//printf("source_id = %d, face_id = %d, r = %.5e, f = %.5e\n", ray.source_id, ray.face_id, ray.r, ray.f);

		struct {
			solid_angle_poly poly;
			dvec grind;
			int fid;
			real f, o;
		} beam_stack[128];
		int nstack = 0;

		rvec cubeverts[8];
		for(v = 0; v < 8; ++v) {
			cubeverts[v][0] = hp->dx*(ray.grind[0]+((v>>0)&1)) - cs.pos[0]; 
			cubeverts[v][1] = hp->dx*(ray.grind[1]+((v>>1)&1)) - cs.pos[1]; 
			cubeverts[v][2] = hp->dx*(ray.grind[2]+((v>>2)&1)) - cs.pos[2]; 
		}

		rvec clipnorms[4];
		solid_angle_poly tpoly;
		solid_angle_poly curpoly;

		const static int neighbor_face_indices[6][4] = {
			{3,7,5,1},{4,6,2,0},{6,7,3,2},{1,5,4,0},{5,7,6,4},{2,3,1,0}};
		const static int neighbor_grid_offsets[6][3] = {
			{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1}};

		if(ray.face_id == 7) {
				
			// construct the ray poly as it exits each face of the source cell
			for(f = 0, sum = 0.0; f < 6; ++f) {
				curpoly.nverts = 4;
				for(v = 0; v < 4; ++v) {
					for(ax = 0; ax < 3; ++ax)
						curpoly.verts[v].pos[ax] = cubeverts[neighbor_face_indices[f][v]][ax];
					curpoly.verts[v].pos[ax] = cubeverts[neighbor_face_indices[f][v]][ax];
					curpoly.verts[v].pnbrs[0] = (v-1+4)%4;
					curpoly.verts[v].pnbrs[1] = (v+1)%4;
				}
				reduce_beam_poly(&curpoly, &omega);
				sum += omega;

				// push this poly onto the stack, tagging it as entering the corresponding
				// neighboring grid cell
				// TODO: must propagate and attenuate!
				beam_stack[nstack].poly = curpoly;
				beam_stack[nstack].fid = f;
				for(ax = 0; ax < 3; ++ax)
					beam_stack[nstack].grind[ax] = ray.grind[ax] + neighbor_grid_offsets[f][ax];
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
					curpoly.verts[v].pos[ax] = cubeverts[neighbor_face_indices[ray.face_id][v]][ax];
				curpoly.verts[v].pos[ax] = cubeverts[neighbor_face_indices[ray.face_id][v]][ax];
				curpoly.verts[v].pnbrs[1] = (v-1+4)%4;
				curpoly.verts[v].pnbrs[0] = (v+1)%4;
			}
			reduce_beam_poly(&curpoly, &omega);
			beam_stack[nstack].poly = curpoly;
			beam_stack[nstack].fid = 2*(ray.face_id/2)-(ray.face_id%2)+1; // flip the face index around to outgoing status 
			for(ax = 0; ax < 3; ++ax)
				beam_stack[nstack].grind[ax] = ray.grind[ax];
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
			face_id = 2*(beam_stack[nstack].fid/2)-(beam_stack[nstack].fid%2)+1; // flip the face index around to outgoing status 
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
					rtmp.face_id = face_id;
					memcpy(rtmp.grind, grind, sizeof(grind));
					rtmp.f = flux;
					hp->rays[hp->nrays++] = rtmp;
					goto next_fragment;
				}
			}

			// now, clip against each of the new downwind faces, tracking fractions 
			// of the solid angle through from the upwind face
			for(f = 0, sum = 0.0; f < 6; ++f) {
				if(f == face_id) continue;
				for(v = 0; v < 4; ++v)
					cross3(clipnorms[v], cubeverts[neighbor_face_indices[f][v]], cubeverts[neighbor_face_indices[f][(v+1)%4]]);
				tpoly = curpoly;
				clip_beam_poly(&tpoly, clipnorms, 4);
				reduce_beam_poly(&tpoly, &omega_out);
				if(omega_out > 0.0) {

					// process the fragment through the cell
					// get the mean optical depth and momentum transfer as well
					flatind = flat_index(grind, hp);
					hp->rad_grid[flatind].E += flux*omega_out/omega_in; 
					sum += omega_out;
	
					// push this poly onto the stack, tagging it as entering the corresponding
					// neighboring grid cell
					// TODO: must propagate and attenuate!
					beam_stack[nstack].poly = tpoly;
					beam_stack[nstack].fid = f;
					for(ax = 0; ax < 3; ++ax)
						beam_stack[nstack].grind[ax] = grind[ax] + neighbor_grid_offsets[f][ax];
					beam_stack[nstack].f = flux*omega_out/omega_in;
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
		next_ray: continue;
	}
		
	// compress the ray list down, removing rays with no flux
	ornrays = hp->nrays;
	hp->nrays = 0.0;
	for(r = 0; r < ornrays; ++r) {
		ray = hp->rays[r];
		if(ray.f <= 0.0) continue;
		hp->rays[hp->nrays++] = ray;
	}
}


