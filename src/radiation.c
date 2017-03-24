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
	} verts[32];
	int nverts;
} solid_angle_poly;

void get_beam_poly(hydro_ray ray, solid_angle_poly* poly, dvec* ibox, dvec* nbox, dvec* charloop, hydro_problem* hp) {

	real len;
	rvec rbox[2];
	int v, i, ax, baseid, reflvl, refbits, quadrant;
	unpack_id_bits(ray.angle_id, baseid, reflvl, refbits);

	// dimension-specific things, loop directions to follow characteristics,
	// and bounding box calculations for this beam
	poly->nverts = hp->dim;
	if(hp->dim == 2) {
		quadrant = baseid/2; // two base rays per quadrant
		for(v = 0; v < hp->dim; ++v) {
			for(ax = 0; ax < hp->dim; ++ax)
				poly->verts[v].pos[ax] = base_beams_2d[baseid%2][v][ax];
		}
		for(i = 0; i < reflvl; ++i) { // refine the ray
			v = (refbits>>i)&1;
			for(ax = 0; ax < hp->dim; ++ax)
				poly->verts[v].pos[ax] = 0.5*(poly->verts[0].pos[ax]+poly->verts[1].pos[ax]);
			len = sqrt(poly->verts[v].pos[0]*poly->verts[v].pos[0]+poly->verts[v].pos[1]*poly->verts[v].pos[1]); 
			for(ax = 0; ax < hp->dim; ++ax)
				poly->verts[v].pos[ax] /= len; 
		}
	}
	else if(hp->dim == 3) {
		quadrant = baseid/6; // six base rays per octant
		for(v = 0; v < hp->dim; ++v) {
			for(ax = 0; ax < hp->dim; ++ax)
				poly->verts[v].pos[ax] = base_beams_3d[baseid%6][v][ax];
			poly->verts[v].pnbrs[0] = (v-1+hp->dim)%hp->dim;
			poly->verts[v].pnbrs[1] = (v+1)%hp->dim;
		}
		for(i = 0; i < reflvl; ++i) { // refine the ray
			v = (refbits>>i)&1;
			for(ax = 0; ax < hp->dim; ++ax)
				poly->verts[v].pos[ax] = 0.5*(poly->verts[0].pos[ax]+poly->verts[1].pos[ax]);



			len = sqrt(poly->verts[v].pos[0]*poly->verts[v].pos[0]
				+poly->verts[v].pos[1]*poly->verts[v].pos[1]+poly->verts[v].pos[2]*poly->verts[v].pos[2]); 
			for(ax = 0; ax < hp->dim; ++ax)
				poly->verts[v].pos[ax] /= len; 
		}
	}
	else {
		printf("Dimension must be 2 or 3!\n");
		exit(0);
	}


	// calculate the bounding box 
	// get the loop directions such that we always move along characteristics
	// do it in all three dimensions for looping machinery later on
	for(ax = 0; ax < 3; ++ax)
		(*charloop)[ax] = 1-2*((quadrant>>ax)&1);
	for(ax = 0; ax < 3; ++ax) {
		rbox[0][ax] = 1.0e30;
		rbox[1][ax] = -1.0e30;
		for(v = 0; v < hp->dim+1; ++v) {
			if(ray.rmin*poly->verts[v].pos[ax] < rbox[0][ax]) rbox[0][ax] = ray.rmin*poly->verts[v].pos[ax];
			if(ray.rmin*poly->verts[v].pos[ax] > rbox[1][ax]) rbox[1][ax] = ray.rmin*poly->verts[v].pos[ax];
		}
		for(v = 0; v < hp->dim+1; ++v) {
			if(ray.rmax*poly->verts[v].pos[ax] < rbox[0][ax]) rbox[0][ax] = ray.rmax*poly->verts[v].pos[ax];
			if(ray.rmax*poly->verts[v].pos[ax] > rbox[1][ax]) rbox[1][ax] = ray.rmax*poly->verts[v].pos[ax];
		}
		ibox[0][ax] = floor(rbox[0][ax]/hp->dx);
		ibox[1][ax] = ceil(rbox[1][ax]/hp->dx)+1;
		(*nbox)[ax] = ibox[1][ax] - ibox[0][ax];
	}
}


	typedef struct {
		real Itot;
		rvec v0, v1;
		rvec clipnorms[8];
		int nclip;
	} face_info;

void add_clip_faces(rvec v0, rvec v1, face_info* face, hydro_problem* hp) {


	if(hp->dim == 2) {
		face->clipnorms[face->nclip][0] = -v0[1];
		face->clipnorms[face->nclip][1] = v0[0];
		face->nclip++; 
		face->clipnorms[face->nclip][0] = v1[1]; 
		face->clipnorms[face->nclip][1] = -v1[0];
		face->nclip++; 
	}
	else if(hp->dim == 3) {
		cross3(face->clipnorms[face->nclip], v0, v1)
			//printf("Adding clip face with normal %f %f %f\n", face->clipnorms[face->nclip][0], face->clipnorms[face->nclip][1], face->clipnorms[face->nclip][2]);
		face->nclip++; 
	}
}

void clip_beam_poly(solid_angle_poly* poly, face_info* faces, real *domega, rvec* centroid, hydro_problem* hp) {


	int f;
	real sdists[32];
	int clipped[32];

	// variable declarations
	int v, p, nclipped, np, onv, vstart, vcur, vnext, numunclipped; 
	real len;

	// direct access to vertex buffer
	if(poly->nverts <= 0) return;

	*domega = 0.0;

	int* nverts = &poly->nverts;

	// 2D case only!
	// TODO: clean this up
	if(hp->dim == 2) {
		for(f = 0; f < faces->nclip; ++f) {
	
			// split vertices by their distance from the clip plane
			nclipped = 0;
			memset(&clipped, 0, sizeof(clipped));
			for(v = 0; v < 2; ++v) {
				sdists[v] = dot2(poly->verts[v].pos, faces->clipnorms[f]);
				clipped[v] = (sdists[v] < 0.0);
				nclipped += clipped[v]; 
			}

	
			// skip this face if the poly lies entirely on one side of it 
			if(nclipped == 0) continue;
			if(nclipped == 2) {
				poly->nverts = 0;
				return;
			}
	
			// otherwise, check all edges and insert new vertices on the bisected edges 
			for(v = 0; v < 2; ++v) {
				if(clipped[v]) {
					wav2(poly->verts[v].pos, -sdists[1-v],
						poly->verts[1-v].pos, sdists[v],
						poly->verts[v].pos);
				} 
			}
		}
		if(poly->nverts && faces->nclip) {
			for(v = 0; v < 2; ++v) {
				len = sqrt(poly->verts[v].pos[0]*poly->verts[v].pos[0]+poly->verts[v].pos[1]*poly->verts[v].pos[1]);
				poly->verts[v].pos[0] /= len;
				poly->verts[v].pos[1] /= len;
			}
			*domega = asin(cross2(poly->verts[0].pos, poly->verts[1].pos));
		}
	}
	
	else if(hp->dim == 3) {

		//printf("Clipping poly, nclip = %d\n", faces->nclip);


		for(f = 0; f < faces->nclip; ++f) {

			//printf("clip face %d, norm = %f %f %f\n", f, faces->clipnorms[f][0], faces->clipnorms[f][1], faces->clipnorms[f][2]);
	
			// split vertices by their distance from the clip plane
			nclipped = 0;
			memset(&clipped, 0, sizeof(clipped));
			for(v = 0; v < poly->nverts; ++v) {
				sdists[v] = dot3(poly->verts[v].pos, faces->clipnorms[f]);
				clipped[v] = (sdists[v] < 0.0);
				nclipped += clipped[v]; 
			}

			//printf("nclipped = %d\n", nclipped);

	
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

					//printf("  vert %d is clipped, vert %d is not\n", vcur, vnext);


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
				//printf("start reconnect loop  with verts %d -> %d\n", vstart, vcur);
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


			//printf("clipped face %d, norm = %f %f %f\n", f, faces->clipnorms[f][0], faces->clipnorms[f][1], faces->clipnorms[f][2]);
			//for(v = 0; v < *nverts; ++v)
				//printf(" vert %d, pos = %f %f %f, pnbrs = %d %d, clipped = %d\n", v, poly->verts[v].pos[0], poly->verts[v].pos[1], poly->verts[v].pos[2], poly->verts[v].pnbrs[0], poly->verts[v].pnbrs[1], clipped[v]);


		}
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
}



void update_radiation(real dt, hydro_problem* hp) {

	int i, v, r, rminbits, rmaxbits, fin, fout, ornrays, refine, baseid, reflvl, refbits,
		flatind, ii, jj, kk, quadrant, qid, rlvl, f, nbase, split;
	dvec lsgn, grind, ibox[2], nbox;
	real ray_mom_out, fmom_in, domega_in, domega_out, intensity_in, dobase, 
		 ray_flux_out, err, allmin, allmax, r_in, r_out, secthmid, len, ray_omega, 
		 cscthmid, inmin, inmax, outmin, outmax, fmid_in, fmid_out;
	real rmin2[4], rmax2[4], flux_out[4], fmom_out[4];
	rvec x0, rmid, tmpverts[8], ntmp, v0, v1;
	hydro_ray ray, r0, r1;
	face_info ftmp, ftin, ftout, faces_in[4], faces_out[4];
	solid_angle_poly beam_poly, inpoly, outpoly;

				rvec centroid_in, centroid_out;

#define subflux(t0, t1, thmin, thmax, fc, fm) (fc*(t1-t0)/(thmax-thmin))
#define subfmom(t0, t1, thmin, thmax, fc, fm) (0.5*(th0+th1)*subflux(t0, t1, thmin, thmax, fc, fm))

	// dimension-specific constants
	if(hp->dim == 2) {
		nbase = num_base_2d;
		dobase = TWO_PI/nbase;
	}
	else if(hp->dim == 3) {
		nbase = num_base_3d;
		dobase = 2*TWO_PI/nbase;
	}

	// STEP 1: step all existing rays forward and refine if needed
	ornrays = hp->nrays;
	for(r = 0; r < ornrays; ++r) {
		ray = hp->rays[r];
		ray.rmin += CLIGHT*dt;
		ray.rmax += CLIGHT*dt;
		unpack_id_bits(ray.angle_id, baseid, reflvl, refbits);
		ray_omega = dobase/(1<<reflvl);
		split = (ray.rmax*ray_omega > 1*hp->dx);

		split = 0;
		
		if(split) {
			r0 = ray;
			r0.Ftot *= 0.5; 
			r0.angle_id = baseid<<24; // set the base ray ID 
			r0.angle_id |= (reflvl+1)<<20; // set the refinement level
			r0.angle_id |= refbits|(0<<reflvl); // set the refinement bits 
			r1 = ray;
			r1.Ftot *= 0.5;
			r1.angle_id = baseid<<24; 
			r1.angle_id |= (reflvl+1)<<20;
			r1.angle_id |= refbits|(1<<reflvl);
			hp->rays[r] = r0;
			hp->rays[hp->nrays++] = r1;
		}
		else {
			hp->rays[r] = ray;
		}
	}

	// STEP 2: source new rays 
	for(r = 0; r < nbase; ++r) {
		ray.rmin = 0.0;
		ray.rmax = CLIGHT*dt;
		ray.Ftot = dt*1.0/nbase; 
		for(i = 0; i < hp->dim; ++i) // TODO: source rays from their actual cells
			ray.orcell[i] = hp->nx[i]/2;
		ray.angle_id = ((r&0xFF)<<24); 
		hp->rays[hp->nrays++] = ray;
	}
	if(hp->nrays > 160000) printf("Error! Overflowed ray buffer.\n");

	// Step 3: propagate all rays forward by c*dt
	memset(hp->rad_grid, 0, (hp->nx[0]+4)*(hp->nx[1]+4)*(hp->nx[2]+4)*sizeof(rad_vector));
	ornrays = hp->nrays;
	hp->nrays = 0;
	for(r = 0; r < ornrays; ++r) {

		// get the ray angles and grid quadrants from the ID
		// also get the ray's bounding box, relative to the originating cell 
		// this works because we work in one single dummy quadrant
		ray = hp->rays[r];
		get_beam_poly(ray, &beam_poly, ibox, &nbox, &lsgn, hp);
		unpack_id_bits(ray.angle_id, baseid, reflvl, refbits);
		ray_omega = dobase/(1<<reflvl);

		//printf("Beam %d poly: %f %f %f, %f %f %f, %f %f %f\n", r, beam_poly.verts[0].pos[0], beam_poly.verts[0].pos[1], beam_poly.verts[0].pos[2], beam_poly.verts[1].pos[0], beam_poly.verts[1].pos[1], beam_poly.verts[1].pos[2], beam_poly.verts[2].pos[0], beam_poly.verts[2].pos[1], beam_poly.verts[2].pos[2]);
		//printf("  Total flux = %.5e, solid angle = %.5e\n", ray.Ftot, ray_omega);

		//printf("  ibox = %d %d %d to %d %d %d\n", ibox[0][0], ibox[0][1], ibox[0][2], ibox[1][0], ibox[1][1], ibox[1][2]);
		//printf("  nbox = %d %d %d\n", nbox[0], nbox[1], nbox[2]);
		//if(hp->dim == 2)
			//printf("  domega = %f\n", domega2(beam_poly.verts[0].pos, beam_poly.verts[1].pos));
		//else if(hp->dim == 3)
			//printf("  domega = %f\n", domega3(beam_poly.verts[0].pos, beam_poly.verts[1].pos, beam_poly.verts[2].pos));


		if(ray.rmin > 1.0) continue;
		// Hacky BCs. TODO: Figure this out better
		if(ibox[0][0] + ray.orcell[0] >= hp->nx[0]+1) continue; 
		if(ibox[0][1] + ray.orcell[1] >= hp->nx[1]+1)  continue; 
		if(ibox[1][0] + ray.orcell[0] < -1)  continue; 
		if(ibox[1][1] + ray.orcell[1] < -1)  continue; 

		if(hp->dim == 3) {
		
			// allocate an array of flux information
			// fill in grid indices and polar coordinates 
			// then loop over all rays, processing their geometry and fluxes
			struct {
				real flux_in[hp->dim+1];
				real fmom_in[hp->dim+1];
				real r2;
			} voxdata[nbox[0]+1][nbox[1]+1][nbox[2]+1];
			memset(voxdata, 0, sizeof(voxdata));
			for(ii = 0; ii <= nbox[0]; ++ii)
			for(jj = 0; jj <= nbox[1]; ++jj)
			for(kk = 0; kk <= nbox[2]; ++kk) {
				x0[0] = hp->dx*(ibox[0][0]+ii-0.5);
				x0[1] = hp->dx*(ibox[0][1]+jj-0.5);
				x0[2] = hp->dx*(ibox[0][2]+kk-0.5);
				voxdata[ii][jj][kk].r2 = x0[0]*x0[0] + x0[1]*x0[1] + x0[2]*x0[2];
			}
	
			ray_flux_out = 0.0;
			ray_mom_out = 0.0;
			for(ii = 0; ii < nbox[0]; ++ii)
			for(jj = 0; jj < nbox[1]; ++jj) 
			for(kk = 0; kk < nbox[2]; ++kk) {

				//if(ii >= 3)  continue; 
				//if(jj >= 3)  continue; 
				//if(kk >= 3)  continue; 
	
				if(ii >= hp->nx[0])  continue; 
				if(jj >= hp->nx[1])  continue; 
				if(kk >= hp->nx[2])  continue; 

				int px = 10, py = 1, pz = 1;
				int print = (ii == px && jj == py && kk == pz);
	
	
				grind[0] = lsgn[0]*(ibox[0][0] + ii) + ray.orcell[0];
				grind[1] = lsgn[1]*(ibox[0][1] + jj) + ray.orcell[1];
				grind[2] = lsgn[2]*(ibox[0][2] + kk) + ray.orcell[2];
				x0[0] = hp->dx*(ibox[0][0]+ii-0.5);
				x0[1] = hp->dx*(ibox[0][1]+jj-0.5);
				x0[2] = hp->dx*(ibox[0][2]+kk-0.5);
				flatind = flat_index(grind, hp);
	
				// get the bitmask for this voxel telling us the 
				// incoming and outgoing ray geometry
				rminbits = 0;
				rmaxbits = 0;
				rvec cubeverts[8];
				for(v = 0; v < 8; ++v) {
					rmin2[v] = voxdata[ii+((v>>0)&1)][jj+((v>>1)&1)][kk+((v>>2)&1)].r2 - ray.rmin*ray.rmin;
					rminbits |= (rmin2[v] > 0.0)<<v;
					rmax2[v] = ray.rmax*ray.rmax - voxdata[ii+((v>>0)&1)][jj+((v>>1)&1)][kk+((v>>2)&1)].r2;
			//printf("clip face %d, norm = %f %f %f\n", f, faces->clipnorms[f][0], faces->clipnorms[f][1], faces->clipnorms[f][2]);
			//for(v = 0; v < *nverts; ++v)
				//printf(" vert %d, pos = %f %f %f, pnbrs = %d %d, clipped = %d\n", v, poly->verts[v].pos[0], poly->verts[v].pos[1], poly->verts[v].pos[2], poly->verts[v].pnbrs[0], poly->verts[v].pnbrs[1], clipped[v]);
					rmaxbits |= (rmax2[v] > 0.0)<<v;
					cubeverts[v][0] = x0[0]+hp->dx*((v>>0)&1); 
					cubeverts[v][1] = x0[1]+hp->dx*((v>>1)&1); 
					cubeverts[v][2] = x0[2]+hp->dx*((v>>2)&1); 
				}
				if(ray.orcell[0] == grind[0] && ray.orcell[1] == grind[1] && ray.orcell[2] == grind[2]) {
					voxdata[ii][jj][kk].flux_in[hp->dim] = ray.Ftot;
					rminbits = 0xFFFF;
				} 
				if(!(rminbits && rmaxbits)) continue;

				// incoming ray cases
				memset(faces_in, 0, sizeof(faces_in));
				switch(rminbits) {

					case 0xFF:
						add_clip_faces(cubeverts[4], cubeverts[0], &faces_in[0], hp);
						add_clip_faces(cubeverts[6], cubeverts[4], &faces_in[0], hp);
						add_clip_faces(cubeverts[2], cubeverts[6], &faces_in[0], hp);
						add_clip_faces(cubeverts[0], cubeverts[2], &faces_in[0], hp);
						add_clip_faces(cubeverts[1], cubeverts[0], &faces_in[1], hp);
						add_clip_faces(cubeverts[5], cubeverts[1], &faces_in[1], hp);
						add_clip_faces(cubeverts[4], cubeverts[5], &faces_in[1], hp);
						add_clip_faces(cubeverts[0], cubeverts[4], &faces_in[1], hp);
						add_clip_faces(cubeverts[0], cubeverts[1], &faces_in[2], hp);
						add_clip_faces(cubeverts[2], cubeverts[0], &faces_in[2], hp);
						add_clip_faces(cubeverts[3], cubeverts[2], &faces_in[2], hp);
						add_clip_faces(cubeverts[1], cubeverts[3], &faces_in[2], hp);
						break;


					case 0xF0:
						printf("0xf0\n");

						add_clip_faces(cubeverts[4], cubeverts[0], &faces_in[0], hp);
						add_clip_faces(cubeverts[6], cubeverts[4], &faces_in[0], hp);
						add_clip_faces(cubeverts[2], cubeverts[6], &faces_in[0], hp);
						add_clip_faces(cubeverts[0], cubeverts[2], &faces_in[0], hp);






						add_clip_faces(cubeverts[1], cubeverts[0], &faces_in[1], hp);
						add_clip_faces(cubeverts[5], cubeverts[1], &faces_in[1], hp);
						add_clip_faces(cubeverts[4], cubeverts[5], &faces_in[1], hp);
						add_clip_faces(cubeverts[0], cubeverts[4], &faces_in[1], hp);
						add_clip_faces(cubeverts[0], cubeverts[1], &faces_in[2], hp);
						add_clip_faces(cubeverts[2], cubeverts[0], &faces_in[2], hp);
						add_clip_faces(cubeverts[3], cubeverts[2], &faces_in[2], hp);
						add_clip_faces(cubeverts[1], cubeverts[3], &faces_in[2], hp);

						break;


					case 0xFFFF: // ray origin, already set the flux
						//add_clip_faces(beam_poly.verts[0].pos, beam_poly.verts[1].pos, &faces_in[3], hp);
						//add_clip_faces(beam_poly.verts[1].pos, beam_poly.verts[2].pos, &faces_in[3], hp);
						//add_clip_faces(beam_poly.verts[2].pos, beam_poly.verts[0].pos, &faces_in[3], hp);
						break;
					default: 
						printf("Bad incoming bit mask! How can this be?? 0x%x\n", rminbits);
						//exit(0);
						return;
						break;
				}
	
				// outgoing ray cases
				memset(faces_out, 0, sizeof(faces_out));
				switch(rmaxbits) {
					case 0xFF:
						// TODO: double-check the ordering of these faces...
						add_clip_faces(cubeverts[1], cubeverts[3], &faces_out[0], hp);
						add_clip_faces(cubeverts[3], cubeverts[7], &faces_out[0], hp);
						add_clip_faces(cubeverts[7], cubeverts[5], &faces_out[0], hp);
						add_clip_faces(cubeverts[5], cubeverts[1], &faces_out[0], hp);
						add_clip_faces(cubeverts[3], cubeverts[2], &faces_out[1], hp);
						add_clip_faces(cubeverts[2], cubeverts[6], &faces_out[1], hp);
						add_clip_faces(cubeverts[6], cubeverts[7], &faces_out[1], hp);
						add_clip_faces(cubeverts[7], cubeverts[3], &faces_out[1], hp);
						add_clip_faces(cubeverts[4], cubeverts[5], &faces_out[2], hp);
						add_clip_faces(cubeverts[5], cubeverts[7], &faces_out[2], hp);
						add_clip_faces(cubeverts[7], cubeverts[6], &faces_out[2], hp);
						add_clip_faces(cubeverts[6], cubeverts[4], &faces_out[2], hp);
						break;
					default:
						//printf("outgoing 0x%x\n", rmaxbits);

						//int vcur, vnext, vcurgood, vnextgood, ax;

						//for(vcur = 0; vcur < 8; ++vcur) {

							//vcurgood = ((rmaxbits>>vcur)&1);
							//if(!vcurgood) continue;

							//printf(" vert %d is inside\n", vcur);



							//for(ax = 0; ax < hp->dim; ++ax) {
								//vnext = vcur^(1<<ax); // the next vertex along this axis
								//vnextgood = ((rmaxbits>>vnext)&1);

								//if(vnextgood) continue;

								//// vcur and vnext are on opposite sides!

								//printf("Found bisected edge between verts %d, %d\n", vcur, vnext);
							
							
							
							//}
						
						//}

						hp->rad_grid[flatind].E = 1.0; 

						continue;


						//printf("Bad outgoing bit mask! How can this be?? 0x%x\n", rmaxbits);
						//exit(0);
						break;
				}


	
				// compute pairwise intersections between in, out faces
				// to get the outgoing fluxes
				memset(flux_out, 0, sizeof(flux_out));
				memset(fmom_out, 0, sizeof(fmom_out));
				for(fin = 0; fin <= hp->dim; ++fin) {
					ftin = faces_in[fin];
					inpoly = beam_poly;


					//printf("Clipping input axis %d...\n", fin);
					clip_beam_poly(&inpoly, &ftin, &domega_in, &centroid_in, hp);
					if(domega_in <= 0.0) continue;
	
					// set the incoming flux from the beam front 
					// and get the intensity
					// TODO: numerical stability??
					if(fin == hp->dim && rminbits == 0xFFFF) 
						voxdata[ii][jj][kk].flux_in[fin] = ray.Ftot*domega_in/ray_omega;
					intensity_in = voxdata[ii][jj][kk].flux_in[fin]/domega_in;

					if(intensity_in <= 0.0) continue;

					//printf("intensity_in = %.5e\n", intensity_in);

					real out_tot = 0.0; 


					real influx = 0.0;
					real outflux = 0.0;
					influx += voxdata[ii][jj][kk].flux_in[fin];
	
					for(fout = 0; fout <= hp->dim; ++fout) {
						ftout = faces_out[fout];
						outpoly = inpoly;

						clip_beam_poly(&outpoly, &ftout, &domega_out, &centroid_out, hp);

						if(fout == hp->dim) continue;
						
						if(domega_out <= 0.0) continue;

						out_tot += domega_out;


						//if(print) printf("Centroid_out = %f %f %f\n", centroid_out[0], centroid_out[1], centroid_out[2]);

	
						// track a solid beam between the two faces
						// calculate the optical depth to first order in theta
						real k = 0.0;
						fmid_in = intensity_in*domega_out;

					err = fabs(1.0-fmid_in/(domega_out/ray_omega*ray.Ftot));
					if(err > 1.0e-12)
						printf("fmid_in = %.5e, expected flux = %.5e\n", fmid_in, domega_out/ray_omega*ray.Ftot);

	
						// get the correct inner and outer radii, then attenuate exponentially 
						if(fin == 0) r_in = x0[0]/centroid_out[0]; // x-facing faces
						else if(fin == 1) r_in = x0[1]/centroid_out[1]; // y-facing faces
						else if(fin == 2) r_in = x0[2]/centroid_out[2]; // z-facing faces
						else if(fin == hp->dim) r_in = ray.rmin; // incoming ray front 
						if(rminbits == 0xFFFF) r_in = 0.0; // TODO: why is this needed??
						if(fout == 0) r_out = (x0[0]+hp->dx)/centroid_out[0];
						else if(fout == 1) r_out = (x0[1]+hp->dx)/centroid_out[1];
						else if(fout == 2) r_out = (x0[2]+hp->dx)/centroid_out[2];
						else if(fout == hp->dim) r_out = ray.rmax; //outgoing ray front 
						//printf("Fmid_in = %.5e\n", fmid_in);

						//r_in = 0.0; r_out = 1.0;
	
						// get the correct inner and outer radii, then attenuate exponentially 
						//rmid[0] = 0.5*(outpoly.verts[0].pos[0]+outpoly.verts[1].pos[0]);
						//rmid[1] = 0.5*(outpoly.verts[0].pos[1]+outpoly.verts[1].pos[1]);
	
						// update grid hydro only if we are in a non-ghost cell

						fmid_out = fmid_in*exp(-k*hp->grid[flatind].rho*(r_out-r_in));

						if(grind[0] >= 0 && grind[0] < hp->nx[0]
								&& grind[1] >= 0 && grind[1] < hp->nx[1]) {
	
							// calculate the outgoing flux, write out the radiation energy density 
							// TODO: interpolate the mean value using the exponential
							hp->rad_grid[flatind].E += 0.5*(fmid_in+fmid_out)*(r_out-r_in)/(ray.rmax-ray.rmin)/(hp->dx*hp->dx);
							//hp->rad_grid[flatind].E += fmid_in-fmid_out;
		
							// update hydro terms by conserving the difference
							// between downwind and upwind fluxes
							//hp->grid[flatind].etot += fmid_in-fmid_out; 
							//hp->grid[flatind].mom[0] += lsgn[0]*cos(thmid)*(fmid_in-fmid_out)/CLIGHT; // TODO: review this!
							//hp->grid[flatind].mom[1] += lsgn[1]*sin(thmid)*(fmid_in-fmid_out)/CLIGHT; 
						}
	
						// update the outgoing ray fluxes and momenta
						flux_out[fout] += fmid_out; 
						

						outflux += fmid_out;
					}

					//err = fabs(1.0-domega_out/domega_in);
					//if(err > 1.0e-12)
					//printf("Clipped input axis %d, domega_in = %.5e, domega_out = %.5e, err = %.5e\n", fin, domega_in, domega_out, err);

					//err = fabs(1.0-outflux/influx);
					//if(err > 1.0e-10)
					//printf("Influx = %.5e, outflux = %.5e, frac = %.5e\n", influx, outflux, err);
				}
	
				// finally, compute the flux and moment through to the next cell
				voxdata[ii+1][jj][kk].flux_in[0] = flux_out[0];
				voxdata[ii][jj+1][kk].flux_in[1] = flux_out[1];
				voxdata[ii][jj][kk+1].flux_in[2] = flux_out[2];
				ray_flux_out += flux_out[hp->dim];

				//exit(0);
			}
			ray.Ftot = ray_flux_out;
			hp->rays[hp->nrays++] = ray;
		}

		else if(hp->dim == 2) {
		
			// allocate an array of flux information
			// fill in grid indices and polar coordinates 
			// then loop over all rays, processing their geometry and fluxes
			struct {
				real flux_in[hp->dim+1];
				real fmom_in[hp->dim+1];
				real r2;
			} voxdata[nbox[0]+1][nbox[1]+1];
			memset(voxdata, 0, sizeof(voxdata));
			for(ii = 0; ii <= nbox[0]; ++ii)
			for(jj = 0; jj <= nbox[1]; ++jj) {
				x0[0] = hp->dx*(ibox[0][0]+ii-0.5);
				x0[1] = hp->dx*(ibox[0][1]+jj-0.5);
				voxdata[ii][jj].r2 = x0[0]*x0[0] + x0[1]*x0[1];
			}
	
	
			printf("Looping...\n");
	
			ray_flux_out = 0.0;
			ray_mom_out = 0.0;
			for(ii = 0; ii < nbox[0]; ++ii)
			for(jj = 0; jj < nbox[1]; ++jj) {
	
				if(ii >= hp->nx[0])  continue; 
				if(jj >= hp->nx[1])  continue; 
	
	
				grind[0] = lsgn[0]*(ibox[0][0] + ii) + ray.orcell[0];
				grind[1] = lsgn[1]*(ibox[0][1] + jj) + ray.orcell[1];
				x0[0] = hp->dx*(ibox[0][0]+ii-0.5);
				x0[1] = hp->dx*(ibox[0][1]+jj-0.5);
				flatind = flat_index(grind, hp);
	
				// get the bitmask for this voxel telling us the 
				// incoming and outgoing ray geometry
				rminbits = 0;
				rmaxbits = 0;
				for(v = 0; v < (1<<hp->dim); ++v) {
					rmin2[v] = voxdata[ii+((v>>0)&1)][jj+((v>>1)&1)].r2 - ray.rmin*ray.rmin;
					rminbits |= (rmin2[v] > 0.0)<<v;
					rmax2[v] = ray.rmax*ray.rmax - voxdata[ii+((v>>0)&1)][jj+((v>>1)&1)].r2;
					rmaxbits |= (rmax2[v] > 0.0)<<v;
				}
				if(ray.orcell[0] == grind[0] && ray.orcell[1] == grind[1]) 
					rminbits = 0xFF; // special case for the origin cell 
				if(!(rminbits && rmaxbits)) continue;
	
				// incoming ray cases
				memset(faces_in, 0, sizeof(faces_in));
				switch(rminbits) {
					case 0x08:
						tmpverts[1][0] = x0[0] + hp->dx; 
						tmpverts[1][1] = sqrt(x0[1]*x0[1]-rmin2[1]); 
						tmpverts[2][0] = sqrt(x0[0]*x0[0]-rmin2[2]); 
						tmpverts[2][1] = x0[1] + hp->dx; 
						add_clip_faces(tmpverts[1], tmpverts[2], &faces_in[hp->dim], hp);
						break;
					case 0x0A:
						tmpverts[0][0] = x0[0] + hp->dx; 
						tmpverts[0][1] = x0[1]; 
						tmpverts[1][0] = sqrt(x0[0]*x0[0]-rmin2[0]); 
						tmpverts[1][1] = x0[1]; 
						tmpverts[2][0] = sqrt(x0[0]*x0[0]-rmin2[2]); 
						tmpverts[2][1] = x0[1] + hp->dx; 
						add_clip_faces(tmpverts[0], tmpverts[1], &faces_in[1], hp);
						add_clip_faces(tmpverts[1], tmpverts[2], &faces_in[hp->dim], hp);
						break;
					case 0x0C:
						tmpverts[1][0] = x0[0] + hp->dx; 
						tmpverts[1][1] = sqrt(x0[1]*x0[1]-rmin2[1]); 
						tmpverts[2][0] = x0[0];
						tmpverts[2][1] = sqrt(x0[1]*x0[1]-rmin2[0]); 
						tmpverts[3][0] = x0[0]; 
						tmpverts[3][1] = x0[1] + hp->dx; 
						add_clip_faces(tmpverts[1], tmpverts[2], &faces_in[hp->dim], hp);
						add_clip_faces(tmpverts[2], tmpverts[3], &faces_in[0], hp);
						break;
					case 0x0E:
						tmpverts[0][0] = x0[0] + hp->dx; 
						tmpverts[0][1] = x0[1]; 
						tmpverts[1][0] = sqrt(x0[0]*x0[0]-rmin2[0]); 
						tmpverts[1][1] = x0[1]; 
						tmpverts[2][0] = x0[0]; 
						tmpverts[2][1] = sqrt(x0[1]*x0[1]-rmin2[0]); 
						tmpverts[3][0] = x0[0]; 
						tmpverts[3][1] = x0[1] + hp->dx; 
						add_clip_faces(tmpverts[0], tmpverts[1], &faces_in[1], hp);
						add_clip_faces(tmpverts[1], tmpverts[2], &faces_in[hp->dim], hp);
						add_clip_faces(tmpverts[2], tmpverts[3], &faces_in[0], hp);
						break;
					case 0x0F:
						tmpverts[0][0] = x0[0] + hp->dx; 
						tmpverts[0][1] = x0[1];
						tmpverts[1][0] = x0[0];
						tmpverts[1][1] = x0[1];
						tmpverts[2][0] = x0[0];
						tmpverts[2][1] = x0[1] + hp->dx;
						add_clip_faces(tmpverts[0], tmpverts[1], &faces_in[1], hp);
						add_clip_faces(tmpverts[1], tmpverts[2], &faces_in[0], hp);
						break;
					case 0xFF: // special case for the ray origin
						// TODO: is this needed when clipping beam polys??
						add_clip_faces(beam_poly.verts[0].pos, beam_poly.verts[1].pos, &faces_in[hp->dim], hp);
						break;
					default: 
						printf("Bad incoming bit mask! How can this be?? 0x%x\n", rminbits);
						exit(0);
						break;
				}
	
				// outgoing ray cases
				memset(faces_out, 0, sizeof(faces_out));
				switch(rmaxbits) {
					case 0x07:
						tmpverts[0][0] = x0[0] + hp->dx; 
						tmpverts[0][1] = x0[1]; 
						tmpverts[1][0] = x0[0] + hp->dx; 
						tmpverts[1][1] = sqrt(x0[1]*x0[1]+rmax2[1]); 
						tmpverts[2][0] = sqrt(x0[0]*x0[0]+rmax2[2]); 
						tmpverts[2][1] = x0[1]+hp->dx; 
						tmpverts[3][0] = x0[0]; 
						tmpverts[3][1] = x0[1]+hp->dx; 
						add_clip_faces(tmpverts[0], tmpverts[1], &faces_out[0], hp);
						add_clip_faces(tmpverts[1], tmpverts[2], &faces_out[hp->dim], hp);
						add_clip_faces(tmpverts[2], tmpverts[3], &faces_out[1], hp);
						break;
					case 0x05:
						tmpverts[1][0] = sqrt(x0[0]*x0[0]+rmax2[0]); 
						tmpverts[1][1] = x0[1]; 
						tmpverts[2][0] = sqrt(x0[0]*x0[0]+rmax2[2]); 
						tmpverts[2][1] = x0[1] + hp->dx; 
						tmpverts[3][0] = x0[0]; 
						tmpverts[3][1] = x0[1] + hp->dx; 
						add_clip_faces(tmpverts[1], tmpverts[2], &faces_out[hp->dim], hp);
						add_clip_faces(tmpverts[2], tmpverts[3], &faces_out[1], hp);
						add_clip_faces(tmpverts[1], tmpverts[3], &faces_in[1], hp);
						break;
					case 0x03:
						tmpverts[0][0] = x0[0] + hp->dx; 
						tmpverts[0][1] = x0[1]; 
						tmpverts[1][0] = x0[0] + hp->dx; 
						tmpverts[1][1] = sqrt(x0[1]*x0[1]+rmax2[1]); 
						tmpverts[2][0] = x0[0];
						tmpverts[2][1] = sqrt(x0[1]*x0[1]+rmax2[0]); 
						add_clip_faces(tmpverts[0], tmpverts[1], &faces_out[0], hp);
						add_clip_faces(tmpverts[1], tmpverts[2], &faces_out[hp->dim], hp);
						add_clip_faces(tmpverts[0], tmpverts[2], &faces_in[0], hp);
						break;
					case 0x01:
						tmpverts[1][0] = sqrt(x0[0]*x0[0]+rmax2[0]); 
						tmpverts[1][1] = x0[1]; 
						tmpverts[2][0] = x0[0]; 
						tmpverts[2][1] = sqrt(x0[1]*x0[1]+rmax2[0]); 
						add_clip_faces(tmpverts[1], tmpverts[2], &faces_out[hp->dim], hp);
						add_clip_faces(tmpverts[1], tmpverts[2], &faces_in[0], hp);
						add_clip_faces(tmpverts[1], tmpverts[2], &faces_in[1], hp);
						break;
					case 0x0F:
						tmpverts[0][0] = x0[0] + hp->dx; 
						tmpverts[0][1] = x0[1];
						tmpverts[1][0] = x0[0] + hp->dx;
						tmpverts[1][1] = x0[1] + hp->dx;
						tmpverts[2][0] = x0[0];
						tmpverts[2][1] = x0[1] + hp->dx;
						add_clip_faces(tmpverts[0], tmpverts[1], &faces_out[0], hp);
						add_clip_faces(tmpverts[1], tmpverts[2], &faces_out[1], hp);
						break;
					default: 
						printf("Bad outgoing bit mask! How can this be?? 0x%x\n", rmaxbits);
						exit(0);
						break;
				}
	
				// compute pairwise intersections between in, out faces
				// to get the outgoing fluxes
				memset(flux_out, 0, sizeof(flux_out));
				memset(fmom_out, 0, sizeof(fmom_out));
				for(fin = 0; fin <= hp->dim; ++fin) {
					ftin = faces_in[fin];
					inpoly = beam_poly;
					clip_beam_poly(&inpoly, &ftin, &domega_in, &centroid_in, hp);
					if(domega_in <= 0.0) continue;
	
					// set the incoming flux from the beam front 
					// and get the intensity
					// TODO: numerical stability??
					if(fin == hp->dim) 
						voxdata[ii][jj].flux_in[fin] = ray.Ftot*domega_in/ray_omega;
					intensity_in = voxdata[ii][jj].flux_in[fin]/domega_in;
	
					for(fout = 0; fout <= hp->dim; ++fout) {
						ftout = faces_out[fout];
						outpoly = inpoly;
						clip_beam_poly(&outpoly, &ftout, &domega_out, &centroid_out, hp);
						if(domega_out <= 0.0) continue;
	
						// track a solid beam between the two faces
						// calculate the optical depth to first order in theta
						real k = 00.0;
						fmid_in = intensity_in*domega_out;
	
						// get the correct inner and outer radii, then attenuate exponentially 
						rmid[0] = 0.5*(outpoly.verts[0].pos[0]+outpoly.verts[1].pos[0]);
						rmid[1] = 0.5*(outpoly.verts[0].pos[1]+outpoly.verts[1].pos[1]);
						len = sqrt(rmid[0]*rmid[0]+rmid[1]*rmid[1]);
						secthmid = len/rmid[0]; 
						cscthmid = len/rmid[1];
						if(fin == 0) r_in = x0[0]*secthmid; // x-facing faces
						else if(fin == 1) r_in = x0[1]*cscthmid; // y-facing faces
						else if(fin == hp->dim) r_in = ray.rmin; // incoming ray front 
						if(rminbits == 0xFF) r_in = 0.0; // TODO: why is this needed??
						if(fout == 0) r_out = (x0[0]+hp->dx)*secthmid;
						else if(fout == 1) r_out = (x0[1]+hp->dx)*cscthmid;
						else if(fout == hp->dim) r_out = ray.rmax; //outgoing ray front 
	
						// update grid hydro only if we are in a non-ghost cell
						fmid_out = fmid_in;
						if(grind[0] >= 0 && grind[0] < hp->nx[0]
								&& grind[1] >= 0 && grind[1] < hp->nx[1]) {
	
							// calculate the outgoing flux, write out the radiation energy density 
							// TODO: interpolate the mean value using the exponential
							fmid_out = fmid_in*exp(-k*hp->grid[flatind].rho*(r_out-r_in));
							hp->rad_grid[flatind].E += 0.5*(fmid_in+fmid_out)*(r_out-r_in)/(ray.rmax-ray.rmin)/(hp->dx*hp->dx);
		
							// update hydro terms by conserving the difference
							// between downwind and upwind fluxes
							hp->grid[flatind].etot += fmid_in-fmid_out; 
							//hp->grid[flatind].mom[0] += lsgn[0]*cos(thmid)*(fmid_in-fmid_out)/CLIGHT; // TODO: review this!
							//hp->grid[flatind].mom[1] += lsgn[1]*sin(thmid)*(fmid_in-fmid_out)/CLIGHT; 
						}
	
						// update the outgoing ray fluxes and momenta
						flux_out[fout] += fmid_out; 
					}
				}
	
				// finally, compute the flux and moment through to the next cell
				voxdata[ii+1][jj].flux_in[0] = flux_out[0];
				voxdata[ii][jj+1].flux_in[1] = flux_out[1];
				ray_flux_out += flux_out[hp->dim];
			}
			ray.Ftot = ray_flux_out;
			hp->rays[hp->nrays++] = ray;

		
		}
	}
}



