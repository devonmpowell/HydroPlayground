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

typedef struct {
	struct {
		rvec pos;
		dvec pnbrs;
	} verts[32];
	int nverts;
} solid_angle_poly;

void get_beam_poly(hydro_ray ray, solid_angle_poly* poly, dvec* ibox, dvec* nbox, dvec* charloop, hydro_problem* hp) {

	rvec rbox[2];
	int v, ax, baseid, reflvl, refbits, quadrant;
	unpack_id_bits(ray.angle_id, baseid, reflvl, refbits);

	// dimension-specific things, loop directions to follow characteristics,
	// and bounding box calculations for this beam
	poly->nverts = hp->dim;
	if(hp->dim == 2) {
		quadrant = baseid/2; // two base rays per quadrant
		for(v = 0; v < hp->dim; ++v)
		for(ax = 0; ax < hp->dim; ++ax)
			poly->verts[v].pos[ax] = base_beams_2d[baseid%2][v][ax];
	}
	else if(hp->dim == 3) {
		quadrant = baseid/6; // six base rays per octant
		for(v = 0; v < hp->dim; ++v)
		for(ax = 0; ax < hp->dim; ++ax)
			poly->verts[v].pos[ax] = base_beams_3d[baseid%6][v][ax];
	}
	else {
		printf("Dimension must be 2 or 3!\n");
		exit(0);
	}

	// get the loop directions such that we always move along characteristics
	// calculate the bounding box 
	for(ax = 0; ax < hp->dim; ++ax)
		(*charloop)[ax] = 1-2*((quadrant>>ax)&1);
	for(ax = 0; ax < hp->dim; ++ax) {
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

void add_clip_faces(rvec v0, rvec v1, face_info* face) {

	rvec ntmp;

	ntmp[0] = -v0[1]; 
	ntmp[1] = v0[0]; 
	face->clipnorms[face->nclip][0] = ntmp[0];
	face->clipnorms[face->nclip][1] = ntmp[1];
	face->nclip++; 
	ntmp[0] = v1[1]; 
	ntmp[1] = -v1[0]; 
	face->clipnorms[face->nclip][0] = ntmp[0];	
	face->clipnorms[face->nclip][1] = ntmp[1];	
	face->nclip++; 

}

void clip_beam_poly(solid_angle_poly* poly, face_info* faces, real *domega, hydro_problem* hp) {


	int f;
	real sdists[32];
	int clipped[32];

	// variable declarations
	int v, p, nclipped, np, onv, vstart, vcur, vnext, numunclipped; 
	real len;

	// direct access to vertex buffer
	if(poly->nverts <= 0) return;

#define dot2(va, vb) (va[0]*vb[0] + va[1]*vb[1])
#define cross2(va, vb) (va[0]*vb[1]-va[1]*vb[0])
#define wav2(va, wa, vb, wb, vr) {			\
	vr[0] = (wa*va[0] + wb*vb[0])/(wa + wb);	\
	vr[1] = (wa*va[1] + wb*vb[1])/(wa + wb);	\
}

	// 2D case only!
	// TODO: clean this up
	if(hp->dim == 2) {
		*domega = 0.0;
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

#if 0


	// signed distances to the clipping plane
	r2d_real sdists[R2D_MAX_VERTS];
	r2d_real smin, smax;

	// for marking clipped vertices
	r2d_int clipped[R2D_MAX_VERTS];

	// loop over each clip plane
	for(p = 0; p < nplanes; ++p) {
	
		// calculate signed distances to the clip plane
		onv = *nverts;
		smin = 1.0e30;
		smax = -1.0e30;
		memset(&clipped, 0, sizeof(clipped));
		for(v = 0; v < onv; ++v) {
			sdists[v] = planes[p].d + dot(vertbuffer[v].pos, planes[p].n);
			if(sdists[v] < smin) smin = sdists[v];
			if(sdists[v] > smax) smax = sdists[v];
			if(sdists[v] < 0.0) clipped[v] = 1;
		}

		// skip this face if the poly lies entirely on one side of it 
		if(smin >= 0.0) continue;
		if(smax <= 0.0) {
			*nverts = 0;
			return;
		}

		// check all edges and insert new vertices on the bisected edges 
		for(vcur = 0; vcur < onv; ++vcur) {
			if(clipped[vcur]) continue;
			for(np = 0; np < 2; ++np) {
				vnext = vertbuffer[vcur].pnbrs[np];
				if(!clipped[vnext]) continue;
				vertbuffer[*nverts].pnbrs[1-np] = vcur;
				vertbuffer[*nverts].pnbrs[np] = -1;
				vertbuffer[vcur].pnbrs[np] = *nverts;
				wav(vertbuffer[vcur].pos, -sdists[vnext],
					vertbuffer[vnext].pos, sdists[vcur],
					vertbuffer[*nverts].pos);
				(*nverts)++;
			}
		}

		// for each new vert, search around the poly for its new neighbors
		// and doubly-link everything
		for(vstart = onv; vstart < *nverts; ++vstart) {
			if(vertbuffer[vstart].pnbrs[1] >= 0) continue;
			vcur = vertbuffer[vstart].pnbrs[0];
			do {
				vcur = vertbuffer[vcur].pnbrs[0]; 
			} while(vcur < onv);
			vertbuffer[vstart].pnbrs[1] = vcur;
			vertbuffer[vcur].pnbrs[0] = vstart;
		}

		// go through and compress the vertex list, removing clipped verts
		// and re-indexing accordingly (reusing `clipped` to re-index everything)
		numunclipped = 0;
		for(v = 0; v < *nverts; ++v) {
			if(!clipped[v]) {
				vertbuffer[numunclipped] = vertbuffer[v];
				clipped[v] = numunclipped++;
			}
		}
		*nverts = numunclipped;
		for(v = 0; v < *nverts; ++v) {
			vertbuffer[v].pnbrs[0] = clipped[vertbuffer[v].pnbrs[0]];
			vertbuffer[v].pnbrs[1] = clipped[vertbuffer[v].pnbrs[1]];
		}	
	}



#endif



}


void update_radiation(real dt, hydro_problem* hp) {

	int i, v, r, rminbits, rmaxbits, fin, fout,
		flatind, ii, jj, quadrant, qid, rlvl, f, nbase;
	dvec lsgn, grind, ibox[2], nbox;
	real ray_mom_out, fmom_in, domega_in, domega_out, intensity_in, 
		 ray_flux_out, err, allmin, allmax, r_in, r_out, secthmid, len, 
		 cscthmid, inmin, inmax, outmin, outmax, fmid_in, fmid_out;
	real rmin2[4], rmax2[4], flux_out[4], fmom_out[4];
	rvec x0, rmid, tmpverts[4], ntmp, v0, v1;
	hydro_ray ray, r0, r1;
	face_info ftmp, ftin, ftout, faces_in[4], faces_out[4];
	solid_angle_poly beam_poly, inpoly, outpoly;


//#define subflux(t0, t1, thmin, thmax, fc, fm)  \
	//((t0 - t1)*(-6*fm*(t0 + t1 - thmax - thmin)  + \
			   //fc*(3*(thmax + thmin)*(t0+t1) - \
				//4*(thmax*thmax + thmax*thmin + thmin*thmin)) )/((thmax-thmin)*(thmax-thmin)*(thmax-thmin)))
//#define subfmom(t0, t1, thmin, thmax, fc, fm) \
			//(fm*(-4*t0*t0*t0 + 3*t0*t0*(thmax + thmin) + \
					   //t1*t1*(4*t1 - 3*(thmax + thmin))) + \
			//2*fc*(t0*t0*t0*(thmax + thmin) - t0*t0*(thmax*thmax + thmax*thmin + thmin*thmin) + \
						//t1*t1*(thmax*thmax + thmax*thmin + thmin*thmin - t1*(thmax + thmin)))/((thmax-thmin)*(thmax-thmin)*(thmax-thmin)))
#define subflux(t0, t1, thmin, thmax, fc, fm) (fc*(t1-t0)/(thmax-thmin))
#define subfmom(t0, t1, thmin, thmax, fc, fm) (0.5*(th0+th1)*subflux(t0, t1, thmin, thmax, fc, fm))


	if(hp->dim != 2) {
		printf("2D radiation only for now!\n");
		exit(0);
	}

	// STEP 1: step all existing rays forward and refine if needed
	int ornrays = hp->nrays;
	for(r = 0; r < ornrays; ++r) {
		ray = hp->rays[r];
		ray.rmin += CLIGHT*dt;
		ray.rmax += CLIGHT*dt;
		hp->rays[r] = ray;

		// no refinement for now...

#if 0
		unpack_id_bits(idbits, baseid, reflvl, refbits);
		quadrant = ;
		
		// TODO: can test the criterion by using the reflvl,
		// no need to compute the thetas!
		thmin = (TWO_PI*qid)/(1<<(hp->dim+rlvl));
		thmax = (TWO_PI*(qid+1))/(1<<(hp->dim+rlvl));
		if(ray.rmax*(thmax-thmin) > 20*hp->dx) {
			r0 = ray; r1 = ray;
			r0.Ftot = subflux(thmin, 0.5*(thmin+thmax), thmin, thmax, ray.Ftot, ray.Fcom); 
			r1.Ftot = subflux(0.5*(thmin+thmax), thmax, thmin, thmax, ray.Ftot, ray.Fcom); 
			r0.Fcom = subfmom(thmin, 0.5*(thmin+thmax), thmin, thmax, ray.Ftot, ray.Fcom); 
			r1.Fcom = subfmom(0.5*(thmin+thmax), thmax, thmin, thmax, ray.Ftot, ray.Fcom); 
			r0.angle_id = quadrant<<24; // set the quadrant (highest byte)
			r0.angle_id |= (rlvl+1)<<16; // set the refinement level
			r0.angle_id |= (qid<<1); // set the ray id (lowest two bytes) 
			r1.angle_id = quadrant<<24; 
			r1.angle_id |= (rlvl+1)<<16;
			r1.angle_id |= (qid<<1)+1; 
			hp->rays[r] = r0;
			hp->rays[hp->nrays++] = r1;
		}
		else {
			hp->rays[r] = ray;
		}
#endif
	}

	// STEP 2: source new rays 
	if(hp->dim == 2) nbase = num_base_2d;
	else if(hp->dim == 3) nbase = num_base_3d;
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
	memset(hp->rad_grid, 0, (hp->nx[0]+4)*(hp->nx[1]+4)*sizeof(rad_vector));
	ornrays = hp->nrays;
	hp->nrays = 0;
	for(r = 0; r < ornrays; ++r) {

		// get the ray angles and grid quadrants from the ID
		// also get the ray's bounding box, relative to the originating cell 
		// this works because we work in one single dummy quadrant
		ray = hp->rays[r];
		get_beam_poly(ray, &beam_poly, ibox, &nbox, &lsgn, hp);

		// Hacky BCs. TODO: Figure this out better
		if(ibox[0][0] + ray.orcell[0] >= hp->nx[0]+1) continue; 
		if(ibox[0][1] + ray.orcell[1] >= hp->nx[1]+1)  continue; 
		if(ibox[1][0] + ray.orcell[0] < -1)  continue; 
		if(ibox[1][1] + ray.orcell[1] < -1)  continue; 

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
		ray_flux_out = 0.0;
		ray_mom_out = 0.0;
		for(ii = 0; ii < nbox[0]; ++ii)
		for(jj = 0; jj < nbox[1]; ++jj) {
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
					add_clip_faces(tmpverts[1], tmpverts[2], &faces_in[hp->dim]);
					break;
				case 0x0A:
					tmpverts[0][0] = x0[0] + hp->dx; 
					tmpverts[0][1] = x0[1]; 
					tmpverts[1][0] = sqrt(x0[0]*x0[0]-rmin2[0]); 
					tmpverts[1][1] = x0[1]; 
					tmpverts[2][0] = sqrt(x0[0]*x0[0]-rmin2[2]); 
					tmpverts[2][1] = x0[1] + hp->dx; 
					add_clip_faces(tmpverts[0], tmpverts[1], &faces_in[1]);
					add_clip_faces(tmpverts[1], tmpverts[2], &faces_in[hp->dim]);
					break;
				case 0x0C:
					tmpverts[1][0] = x0[0] + hp->dx; 
					tmpverts[1][1] = sqrt(x0[1]*x0[1]-rmin2[1]); 
					tmpverts[2][0] = x0[0];
					tmpverts[2][1] = sqrt(x0[1]*x0[1]-rmin2[0]); 
					tmpverts[3][0] = x0[0]; 
					tmpverts[3][1] = x0[1] + hp->dx; 
					add_clip_faces(tmpverts[1], tmpverts[2], &faces_in[hp->dim]);
					add_clip_faces(tmpverts[2], tmpverts[3], &faces_in[0]);
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
					add_clip_faces(tmpverts[0], tmpverts[1], &faces_in[1]);
					add_clip_faces(tmpverts[1], tmpverts[2], &faces_in[hp->dim]);
					add_clip_faces(tmpverts[2], tmpverts[3], &faces_in[0]);
					break;
				case 0x0F:
					tmpverts[0][0] = x0[0] + hp->dx; 
					tmpverts[0][1] = x0[1];
					tmpverts[1][0] = x0[0];
					tmpverts[1][1] = x0[1];
					tmpverts[2][0] = x0[0];
					tmpverts[2][1] = x0[1] + hp->dx;
					add_clip_faces(tmpverts[0], tmpverts[1], &faces_in[1]);
					add_clip_faces(tmpverts[1], tmpverts[2], &faces_in[0]);
					break;
				case 0xFF: // special case for the ray origin
					// TODO: is this needed when clipping beam polys??
					add_clip_faces(beam_poly.verts[0].pos, beam_poly.verts[1].pos, &faces_in[hp->dim]);
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
					add_clip_faces(tmpverts[0], tmpverts[1], &faces_out[0]);
					add_clip_faces(tmpverts[1], tmpverts[2], &faces_out[hp->dim]);
					add_clip_faces(tmpverts[2], tmpverts[3], &faces_out[1]);
					break;
				case 0x05:
					tmpverts[1][0] = sqrt(x0[0]*x0[0]+rmax2[0]); 
					tmpverts[1][1] = x0[1]; 
					tmpverts[2][0] = sqrt(x0[0]*x0[0]+rmax2[2]); 
					tmpverts[2][1] = x0[1] + hp->dx; 
					tmpverts[3][0] = x0[0]; 
					tmpverts[3][1] = x0[1] + hp->dx; 
					add_clip_faces(tmpverts[1], tmpverts[2], &faces_out[hp->dim]);
					add_clip_faces(tmpverts[2], tmpverts[3], &faces_out[1]);
					add_clip_faces(tmpverts[1], tmpverts[3], &faces_in[1]);
					break;
				case 0x03:
					tmpverts[0][0] = x0[0] + hp->dx; 
					tmpverts[0][1] = x0[1]; 
					tmpverts[1][0] = x0[0] + hp->dx; 
					tmpverts[1][1] = sqrt(x0[1]*x0[1]+rmax2[1]); 
					tmpverts[2][0] = x0[0];
					tmpverts[2][1] = sqrt(x0[1]*x0[1]+rmax2[0]); 
					add_clip_faces(tmpverts[0], tmpverts[1], &faces_out[0]);
					add_clip_faces(tmpverts[1], tmpverts[2], &faces_out[hp->dim]);
					add_clip_faces(tmpverts[0], tmpverts[2], &faces_in[0]);
					break;
				case 0x01:
					tmpverts[1][0] = sqrt(x0[0]*x0[0]+rmax2[0]); 
					tmpverts[1][1] = x0[1]; 
					tmpverts[2][0] = x0[0]; 
					tmpverts[2][1] = sqrt(x0[1]*x0[1]+rmax2[0]); 
					add_clip_faces(tmpverts[1], tmpverts[2], &faces_out[hp->dim]);
					add_clip_faces(tmpverts[1], tmpverts[2], &faces_in[0]);
					add_clip_faces(tmpverts[1], tmpverts[2], &faces_in[1]);
					break;
				case 0x0F:
					tmpverts[0][0] = x0[0] + hp->dx; 
					tmpverts[0][1] = x0[1];
					tmpverts[1][0] = x0[0] + hp->dx;
					tmpverts[1][1] = x0[1] + hp->dx;
					tmpverts[2][0] = x0[0];
					tmpverts[2][1] = x0[1] + hp->dx;
					add_clip_faces(tmpverts[0], tmpverts[1], &faces_out[0]);
					add_clip_faces(tmpverts[1], tmpverts[2], &faces_out[1]);
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
				clip_beam_poly(&inpoly, &ftin, &domega_in, hp);
				if(domega_in <= 0.0) continue;

				// set the incoming flux from the beam front 
				// and get the intensity
				// TODO: numerical stability??
				if(fin == hp->dim) 
					voxdata[ii][jj].flux_in[fin] = ray.Ftot*domega_in/(TWO_PI*0.125);
				intensity_in = voxdata[ii][jj].flux_in[fin]/domega_in;

				for(fout = 0; fout <= hp->dim; ++fout) {
					ftout = faces_out[fout];
					outpoly = inpoly;
					clip_beam_poly(&outpoly, &ftout, &domega_out, hp);
					if(domega_out <= 0.0) continue;

					// track a solid beam between the two faces
					// calculate the optical depth to first order in theta
					real k = 0.0;
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
					//fmom_out[fout] += fmid_out*thmid; 
				}
			}

			// finally, compute the flux and moment through to the next cell
			voxdata[ii+1][jj].flux_in[0] = flux_out[0];
			voxdata[ii+1][jj].fmom_in[0] = fmom_out[0];
			voxdata[ii][jj+1].flux_in[1] = flux_out[1];
			voxdata[ii][jj+1].fmom_in[1] = fmom_out[1];
			ray_flux_out += flux_out[hp->dim];
			ray_mom_out += fmom_out[hp->dim];
		}
		ray.Ftot = ray_flux_out;
		ray.Fcom = ray_mom_out; 
		hp->rays[hp->nrays++] = ray;
	}
}



