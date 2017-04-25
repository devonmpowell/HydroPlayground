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
#include "geometry.h"

#define unpack_id_bits(idbits, baseid, reflvl, refbits) { \
	baseid = (ray.angle_id>>24)&0xFF; \
	reflvl = (ray.angle_id>>20)&0xF; \
	refbits =  (ray.angle_id)&0xFFFFF; \
}

#define dot2(va, vb) (va[0]*vb[0] + va.y*vb.y)
#define cross2(va, vb) (va.x*vb.y-va.y*vb.x)
#define wav2(va, wa, vb, wb, vr) {			\
	vr.x = (wa*va.x + wb*vb.x)/(wa + wb);	\
	vr.y = (wa*va.y + wb*vb.y)/(wa + wb);	\
}

#define wav3(va, wa, vb, wb, vr) {			\
	vr.x = (wa*va.x + wb*vb.x)/(wa + wb);	\
	vr.y = (wa*va.y + wb*vb.y)/(wa + wb);	\
	vr.z = (wa*va.z + wb*vb.z)/(wa + wb);	\
}

#define dot3(va, vb) (va.x*vb.x + va.y*vb.y + va.z*vb.z)
#define cross3(v,u,w)          /* CROSS Vector Product */              \
{                                                                       \
    (v).x = (u).y*(w).z - (u).z*(w).y;                             \
    (v).y = (u).z*(w).x - (u).x*(w).z;                             \
    (v).z = (u).x*(w).y - (u).y*(w).x;                             \
}


#define domega2(v0, v1) asin(cross2(v0, v1))
real domega3(rvec v1, rvec v2, rvec v3) {
	int i;
	real det, div;

	// Assumes v1, v2, v3 are unit vectors
	div = sqrt(v1.x*v1.x + v1.y*v1.y + v1.z*v1.z);
	for(i = 0; i < 3; ++i) v1.xyz[i] /= div;
	div = sqrt(v2.x*v2.x + v2.y*v2.y + v2.z*v2.z);
	for(i = 0; i < 3; ++i) v2.xyz[i] /= div;
	div = sqrt(v3.x*v3.x + v3.y*v3.y + v3.z*v3.z);
	for(i = 0; i < 3; ++i) v3.xyz[i] /= div;

	// Solid angle of a triangle by Oosterom and Strackee
	det = v1.x*(v2.y*v3.z-v3.y*v2.z)-v1.y*(v2.x*v3.z-v3.x*v2.z)+v1.z*(v2.x*v3.y-v3.x*v2.y);
	div = 1.0;
	for(i = 0; i < 3; ++i) {
		div += v1.xyz[i]*v2.xyz[i];
		div += v2.xyz[i]*v3.xyz[i];
		div += v3.xyz[i]*v1.xyz[i];
	}
	return 2.0*atan2(det, div);
} 


void clip_beam_poly(psi_poly* poly, psi_plane* faces, int nclip, hydro_problem* hp);

void get_beam_poly(hydro_ray ray, psi_poly* poly, dvec* nbox, dvec* charloop, hydro_problem* hp) {

	real len;
	rvec rbox[2];
	int v, i, ax, baseid, reflvl, refbits, quadrant;
	unpack_id_bits(ray.angle_id, baseid, reflvl, refbits);

	rvec refined_vertex_pos[4];

	// get the quandrant and the refined beam vertex positions
	quadrant = baseid/6; 
	for(v = 0; v < hp->dim; ++v) {
		for(ax = 0; ax < hp->dim; ++ax)
			refined_vertex_pos[v].xyz[ax] = base_beams_3d[baseid][v][ax];
	}
	for(i = 0; i < reflvl; ++i) { // refine the ray
		v = (refbits>>i)&1;
		for(ax = 0; ax < hp->dim; ++ax)
			refined_vertex_pos[v].xyz[ax] = 0.5*(refined_vertex_pos[0].xyz[ax]+refined_vertex_pos[1].xyz[ax]);
		len = sqrt(refined_vertex_pos[v].x*refined_vertex_pos[v].x
				+ refined_vertex_pos[v].y*refined_vertex_pos[v].y
				+ refined_vertex_pos[v].z*refined_vertex_pos[v].z);
		for(ax = 0; ax < hp->dim; ++ax)
			refined_vertex_pos[v].xyz[ax] /= len; 
	}

	// fill in the beam vertices 
	for(ax = 0; ax < hp->dim; ++ax) { 
		poly->verts[0].pos.xyz[ax] = ray.origin.xyz[ax] + ray.rmin*refined_vertex_pos[0].xyz[ax]; 
		poly->verts[1].pos.xyz[ax] = ray.origin.xyz[ax] + ray.rmin*refined_vertex_pos[1].xyz[ax]; 
		poly->verts[2].pos.xyz[ax] = ray.origin.xyz[ax] + ray.rmin*refined_vertex_pos[2].xyz[ax]; 
		poly->verts[3].pos.xyz[ax] = ray.origin.xyz[ax] + ray.rmax*refined_vertex_pos[0].xyz[ax]; 
		poly->verts[4].pos.xyz[ax] = ray.origin.xyz[ax] + ray.rmax*refined_vertex_pos[1].xyz[ax]; 
		poly->verts[5].pos.xyz[ax] = ray.origin.xyz[ax] + ray.rmax*refined_vertex_pos[2].xyz[ax]; 
	}

	// initialize the connectivity and local basis
	poly->verts[0].q.xyz[0] = 0.0; 
	poly->verts[0].q.xyz[1] = 0.0; 
	poly->verts[0].q.xyz[2] = 0.0; 
	poly->verts[1].q.xyz[0] = 1.0; 
	poly->verts[1].q.xyz[1] = 0.0; 
	poly->verts[1].q.xyz[2] = 0.0; 
	poly->verts[2].q.xyz[0] = 0.0; 
	poly->verts[2].q.xyz[1] = 1.0; 
	poly->verts[2].q.xyz[2] = 0.0; 
	poly->verts[3].q.xyz[0] = 0.0; 
	poly->verts[3].q.xyz[1] = 0.0; 
	poly->verts[3].q.xyz[2] = 1.0; 
	poly->verts[4].q.xyz[0] = 1.0; 
	poly->verts[4].q.xyz[1] = 0.0; 
	poly->verts[4].q.xyz[2] = 1.0; 
	poly->verts[5].q.xyz[0] = 0.0; 
	poly->verts[5].q.xyz[1] = 1.0; 
	poly->verts[5].q.xyz[2] = 1.0; 

	poly->verts[0].flags = 1; 
	poly->verts[1].flags = 1; 
	poly->verts[2].flags = 1; 

	poly->verts[0].fnbrs[0] = -1; 
	poly->verts[0].fnbrs[1] = BEAM_IN; 
	poly->verts[0].fnbrs[2] = -1; 
	poly->verts[1].fnbrs[0] = -1; 
	poly->verts[1].fnbrs[1] = BEAM_IN; 
	poly->verts[1].fnbrs[2] = -1; 
	poly->verts[2].fnbrs[0] = -1; 
	poly->verts[2].fnbrs[1] = BEAM_IN; 
	poly->verts[2].fnbrs[2] = -1; 
	poly->verts[3].fnbrs[0] = -1; 
	poly->verts[3].fnbrs[1] = BEAM_OUT; 
	poly->verts[3].fnbrs[2] = -1; 
	poly->verts[4].fnbrs[0] = -1; 
	poly->verts[4].fnbrs[1] = BEAM_OUT; 
	poly->verts[4].fnbrs[2] = -1; 
	poly->verts[5].fnbrs[0] = -1; 
	poly->verts[5].fnbrs[1] = BEAM_OUT; 
	poly->verts[5].fnbrs[2] = -1; 

	poly->verts[0].pnbrs[0] = 2; 
	poly->verts[0].pnbrs[1] = 1; 
	poly->verts[0].pnbrs[2] = 3; 
	poly->verts[1].pnbrs[0] = 0; 
	poly->verts[1].pnbrs[1] = 2; 
	poly->verts[1].pnbrs[2] = 4; 
	poly->verts[2].pnbrs[0] = 1; 
	poly->verts[2].pnbrs[1] = 0; 
	poly->verts[2].pnbrs[2] = 5; 
	poly->verts[3].pnbrs[0] = 4; 
	poly->verts[3].pnbrs[1] = 5; 
	poly->verts[3].pnbrs[2] = 0; 
	poly->verts[4].pnbrs[0] = 5; 
	poly->verts[4].pnbrs[1] = 3; 
	poly->verts[4].pnbrs[2] = 1; 
	poly->verts[5].pnbrs[0] = 3; 
	poly->verts[5].pnbrs[1] = 4; 
	poly->verts[5].pnbrs[2] = 2; 

	// calculate the bounding box 
	// get the loop directions such that we always move along characteristics
	// do it in all three dimensions for looping machinery later on
	for(ax = 0; ax < 3; ++ax)
		(*charloop).ijk[ax] = 1-2*((quadrant>>ax)&1);
	poly->nverts = 2*hp->dim;
	for(ax = 0; ax < hp->dim; ++ax) {
		rbox[0].xyz[ax] = 1.0e30;
		rbox[1].xyz[ax] = -1.0e30;
		for(v = 0; v < poly->nverts; ++v) {
			if(poly->verts[v].pos.xyz[ax] < rbox[0].xyz[ax]) rbox[0].xyz[ax] = poly->verts[v].pos.xyz[ax];
			if(poly->verts[v].pos.xyz[ax] > rbox[1].xyz[ax]) rbox[1].xyz[ax] = poly->verts[v].pos.xyz[ax];
		}
		poly->ibox[0].ijk[ax] = floor(rbox[0].xyz[ax]/hp->dx);
		poly->ibox[1].ijk[ax] = ceil(rbox[1].xyz[ax]/hp->dx)+1;
		(*nbox).ijk[ax] = poly->ibox[1].ijk[ax] - poly->ibox[0].ijk[ax];
	}
}

void update_radiation(real dt, hydro_problem* hp) {

	int i, v, r, ax, rminbits, rmaxbits, fin, fout, ornrays, refine, baseid, reflvl, refbits,
		flatind, quadrant, qid, rlvl, f, nbase, split;
	dvec lsgn, grind, nbox, loopind;
	real ray_mom_out, fmom_in, domega_in, domega_out, intensity_in, dobase, flux_in, 
		 ray_flux_out, err, allmin, allmax, r_in, r_out, secthmid, len, ray_omega, 
		 cscthmid, inmin, inmax, outmin, outmax, fmid_in, fmid_out;
	real rmin2[4], rmax2[4], flux_out[4], fmom_out[4], moments[10];
	rvec rmid, tmpverts[8], ntmp, v0, v1;
	hydro_ray ray, r0, r1;
	psi_poly beam_poly, inpoly, outpoly, cell_poly, beamtmp;

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
			ray.origin.xyz[i] = 0.5*hp->dx*(hp->nx.ijk[i]+1);
		//ray.origin.x = 0.510283476082;
		//ray.origin.y = 0.430283476082;
		//ray.origin.z = 0.490283476082;
		ray.angle_id = ((r&0xFF)<<24); 
		hp->rays[hp->nrays++] = ray;
	}
	if(hp->nrays > 160000) printf("Error! Overflowed ray buffer.\n");

	// Step 3: propagate all rays forward by c*dt
	memset(hp->rad_grid, 0, (hp->nx.i+4)*(hp->nx.j+4)*(hp->nx.k+4)*sizeof(rad_vector));
	ornrays = hp->nrays;
	hp->nrays = 0;
	for(r = 0; r < ornrays; ++r) {

		// get the ray angles and grid quadrants from the ID
		// also get the ray's bounding box, relative to the originating cell 
		// this works because we work in one single dummy quadrant
		ray = hp->rays[r];
		get_beam_poly(ray, &beam_poly, &nbox, &lsgn, hp);
		unpack_id_bits(ray.angle_id, baseid, reflvl, refbits);
		ray_omega = dobase/(1<<reflvl);               

		if(baseid > 5) continue;

		int nfaces;
		psi_face_buffer faces[16], face_in, face_out;
		psi_voxels vox;
		psi_poly curpoly;

		printf("Ray %d\n", r);

		// translate the poly to an integer offset
		dvec offset = beam_poly.ibox[0]; 
		rvec ortmp = ray.origin; 
		for(ax = 0; ax < 3; ++ax) {
			//offset.ijk[ax] = beam_poly.ibox[(1+lsgn.ijk[ax])/2].ijk[ax];
			beam_poly.ibox[0].ijk[ax] -= offset.ijk[ax];
			beam_poly.ibox[1].ijk[ax] -= offset.ijk[ax];
			for(v = 0; v < beam_poly.nverts; ++v)
				beam_poly.verts[v].pos.xyz[ax] -= hp->dx*offset.ijk[ax];
			ortmp.xyz[ax] -= hp->dx*offset.ijk[ax];
		}

		struct {
			real flux_in[8]; // TODO: Don't need this many indices
		} transfer_info[nbox.i+1][nbox.j+1][nbox.k+1];

		memset(transfer_info, 0, sizeof(transfer_info));

		// voxelize the beam
		ray_flux_out = 0.0;
		psi_voxels_init(&vox, &beam_poly, &lsgn, hp);
		while(psi_voxels_next(&vox, &curpoly)) {

			// make sure we're inside the grid bounds before processing
			dvec locind = curpoly.ibox[0];
			for(ax = 0; ax < 3; ++ax)
				grind.ijk[ax] = locind.ijk[ax] + offset.ijk[ax];
			for(ax = 0; ax < hp->dim; ++ax)
				if(grind.ijk[ax] < 0 || grind.ijk[ax] >= hp->nx.ijk[ax])
					goto next_cell;
			flatind = flat_index(grind, hp);

			// process the solid angle intersections of the faces
			// first, translate the voxel into coordinates centered on the ray source
			int nflag = 0;
			for(v = 0; v < curpoly.nverts; ++v) {
				nflag += curpoly.verts[v].flags;
				for(ax = 0; ax < hp->dim; ++ax)
					curpoly.verts[v].pos.xyz[ax] -= ortmp.xyz[ax];
			}
			if(nflag == 3) {
				transfer_info[locind.i+1][locind.j][locind.k].flux_in[0] = 0.33333333*ray.Ftot;
				transfer_info[locind.i][locind.j+1][locind.k].flux_in[2] = 0.33333333*ray.Ftot;
				transfer_info[locind.i][locind.j][locind.k+1].flux_in[4] = 0.33333333*ray.Ftot;
				continue;
			}

			psi_poly poly_in, poly_out;
			real omega_in, omega_in_tot, omega_out_tot;

			omega_in_tot = 0.0;
			omega_out_tot = 0.0;

			psi_extract_faces(&curpoly, faces, &nfaces);

			int ornf = nfaces;
			nfaces = 0;
			for(f = 0; f < ornf; ++f) {

				// make a poly out of it and filter tiny pieces
				face_in = faces[f];
				poly_in.nverts = face_in.nverts;
				for(v = 0; v < face_in.nverts; ++v) {
					poly_in.verts[v].pos = face_in.verts[v];
					poly_in.verts[v].pnbrs[0] = (v + face_in.nverts - 1)%face_in.nverts; 
					poly_in.verts[v].pnbrs[1] = (v + 1)%face_in.nverts; 
				}
				reduce_beam_poly(&poly_in, &omega_in);

				if(face_in.face_id%2)
					omega_out_tot += -omega_in;
				else	
					omega_in_tot += omega_in;


				printf("omega = %.5e\n", omega_in);
				//if(fabs(omega_in) > 1.0e-12*ray_omega)
					faces[nfaces++] = face_in;
			}

			err = fabs(1.0 - omega_out_tot/omega_in_tot);
			if(err > 1.0e-6) {
				printf(" First solid angle check.  Omega_in_tot = %.5e, omega_out_tot = %.5e, err = %.5e\n", omega_in_tot, omega_out_tot, err);
			}



			real domega;
			rvec centroid;

			real omega_x = 0.0;

			real flux_in_tot = 0.0;
			real flux_out_tot = 0.0;

			//printf(" voxel %d %d %d\n", curpoly.ibox[0].i, curpoly.ibox[0].j, curpoly.ibox[0].k);

			for(fin = 0; fin < nfaces; ++fin) {
				face_in = faces[fin];
				if(face_in.face_id%2 != 0) continue;

				// copy the input poly into a new poly struct for clipping
				// get its solid angle as viewed from the source
				poly_in.nverts = face_in.nverts;
				for(v = 0; v < face_in.nverts; ++v) {
					poly_in.verts[v].pos = face_in.verts[v];
					poly_in.verts[v].pnbrs[0] = (v + face_in.nverts - 1)%face_in.nverts; 
					poly_in.verts[v].pnbrs[1] = (v + 1)%face_in.nverts; 

					//printf("In vert = %.5e, %.5e, %.5e\n", face_in.verts[v].x, face_in.verts[v].y, face_in.verts[v].z);
				}
				reduce_beam_poly(&poly_in, &omega_in);
				if(omega_in <= 0.0) continue;


				psi_plane clip_in[16];
				memset(clip_in, 0, sizeof(clip_in));
				int nclip = 0; 
				for(v = 0; v < face_in.nverts; ++v)
					cross3(clip_in[nclip++].n, face_in.verts[(v+1)%face_in.nverts], face_in.verts[v]);
						



				//printf("in %d, domega = %.5e\n", face_in.face_id, omega_in);

				real omega_out;
				real omega_out_tot = 0.0;
				rvec centroid_out;
				rvec cliptmp;
				for(fout = 0; fout < nfaces; ++fout) {
					face_out = faces[fout];
					if(face_out.face_id%2 == 0) continue;

					poly_out.nverts = face_out.nverts;
					for(v = 0; v < face_out.nverts; ++v) {
						poly_out.verts[v].pos = face_out.verts[v];
						poly_out.verts[v].pnbrs[0] = (v + face_out.nverts - 1)%face_out.nverts; 
						poly_out.verts[v].pnbrs[1] = (v + 1)%face_out.nverts; 
	
						//printf("In vert = %.5e, %.5e, %.5e\n", face_out.verts[v].x, face_out.verts[v].y, face_out.verts[v].z);
					}
					reduce_beam_poly(&poly_out, &omega_out);
					//if(omega_out <= 0.0) continue;



					// TODO: make these clip planes ahead of time
					psi_plane clip_out[16];
					memset(clip_out, 0, sizeof(clip_out));
					int nclip = 0; 
					for(v = 0; v < face_out.nverts; ++v) {
						cross3(cliptmp, face_out.verts[(v+1)%face_out.nverts], face_out.verts[v]);
						//printf("  Out vert = %.5e, %.5e, %.5e\n", face_out.verts[v].x, face_out.verts[v].y, face_out.verts[v].z);
						//printf("   face norm = %.5e, %.5e, %.5e\n", cliptmp.x, cliptmp.y, cliptmp.z);
						len = sqrt(cliptmp.x*cliptmp.x + cliptmp.y*cliptmp.y + cliptmp.z*cliptmp.z);

						//if(len > 1.0e-6*hp->dx)
							clip_out[nclip++].n = cliptmp;
					}
					poly_out = poly_in;
					clip_beam_poly(&poly_out, clip_out, nclip, hp); 
					reduce_beam_poly(&poly_out, &domega);

					//printf("  out %d, domega = %.5e\n", face_out.face_id, domega);
					for(v = 0; v < poly_out.nverts; ++v)
						//printf("  Clipped out vert = %.5e, %.5e, %.5e\n", poly_out.verts[v].pos.x, poly_out.verts[v].pos.y, poly_out.verts[v].pos.z);
					//if(domega <= 0.0) continue;

					omega_out += domega;
					//centroid_out[face_out.face_id] = centroid;
					omega_out_tot += domega;
				}

				err = fabs(1.0 - omega_out_tot/omega_in);
				if(err > 1.0e-6) {
					//printf("Incorrect solid angle! Omega_in = %.5e, omega_out_tot = %.5e, err = %.5e\n", omega_in, omega_out_tot, err);
					//exit(0);
				}

				//hp->rad_grid[flatind].E += err; 
				//hp->rad_grid[flatind].E += omega_in;
				hp->rad_grid[flatind].E += omega_out_tot;


#if 0

				if(face_in.face_id == BEAM_IN)
					flux_in = omega_out_tot/ray_omega*ray.Ftot; 
				else
					flux_in = transfer_info[locind.i][locind.j][locind.k].flux_in[face_in.face_id];
				if(flux_in <= 0.0) continue;


				flux_in_tot += flux_in; 

				if(omega_out_tot <= 0.0) {
				
					//hp->rad_grid[flatind].E += 1.0; 
					//
					if(flux_in > 1.0e-10*ray.Ftot)
						printf(" Omega_out = %.5e, flux_in = %.5e for incoming face %d\n", omega_out_tot, flux_in, face_in.face_id);
					continue;
				} 


				// get the fraction of solid angle belonging to each outgoing face
				for(fout = 0; fout < 8; ++fout) 
					omega_out[fout] /= omega_out_tot;

				real flux_out = 0.0;

				// fout now refers to face id 
				for(fout = 1; fout < 6; fout += 2) {

					//if(omega_out[fout] < 1.0e-6)
						//continue;

					fmid_in = flux_in*omega_out[fout];

					centroid = centroid_out[fout];

					fmid_out = fmid_in;

						//if(print) printf("Centroid_out = %f %f %f\n", centroid_out[0], centroid_out[1], centroid_out[2]);

	
						// track a solid beam between the two faces
						// calculate the optical depth to first order in theta
						//real k = 0.0;
						//fmid_in = intensity_in*domega_out;

					////err = fabs(1.0-fmid_in/(domega_out/ray_omega*ray.Ftot));
					////if(err > 1.0e-12)
						////printf("fmid_in = %.5e, expected flux = %.5e\n", fmid_in, domega_out/ray_omega*ray.Ftot);

	
						//// get the correct inner and outer radlocind.i, then attenuate exponentially 
						if(face_in.face_id == 0) r_in = (hp->dx*locind.i-ortmp.x)/centroid.x; // x-facing faces
						else if(face_in.face_id == 2) r_in = (hp->dx*locind.j-ortmp.y)/centroid.y; // y-facing faces
						else if(face_in.face_id == 4) r_in = (hp->dx*locind.k-ortmp.z)/centroid.z; // z-facing faces
						else if(face_in.face_id == hp->dim) r_in = ray.rmin; // incoming ray front 
						//if(rminbits == 0xFFFF) r_in = 0.0; // TODO: why is this needed??

						if(fout == 1) r_out = (hp->dx*(locind.i+1)-ortmp.x)/centroid.x;
						else if(fout == 3) r_out = (hp->dx*(locind.j+1)-ortmp.y)/centroid.y;
						else if(fout == 5) r_out = (hp->dx*(locind.k+1)-ortmp.z)/centroid.z;
						else if(fout == hp->dim) r_out = ray.rmax; //outgoing ray front 
						//printf("Fmid_in = %.5e\n", fmid_in);


					//hp->rad_grid[flatind].E += omega_out[fout]*omega_out_tot; 
					//hp->rad_grid[flatind].E += fmid_in;// 0.5*(fmid_in+fmid_out)*(r_out-r_in)/(ray.rmax-ray.rmin)/(hp->dx*hp->dx); 
					//hp->rad_grid[flatind].E += 0.5*(fmid_in+fmid_out)*(r_out-r_in)/(ray.rmax-ray.rmin)/(hp->dx*hp->dx); 


					transfer_info[locind.i+(fout==1)][locind.j+(fout==3)][locind.k+(fout==5)].flux_in[fout-1]
						+= fmid_in;
					flux_out_tot += fmid_in;
					flux_out += fmid_in;
				
				}

				//hp->rad_grid[flatind].E += flux_in*omega_out[BEAM_OUT];
				ray_flux_out += flux_in*omega_out[BEAM_OUT];
				flux_out_tot += flux_in*omega_out[BEAM_OUT];

				flux_out += flux_in*omega_out[BEAM_OUT];



				err = fabs(1.0 - flux_out/flux_in);
				if(err > 1.0e-15) {
					printf(" Flux error: voxel %d %d %d, face in %d\n", grind.i, grind.j, grind.k, face_in.face_id);
				//printf("flux_in = %.5e, flux_out = %.5e, err = %.5e\n", flux_in, flux_out, err);
					printf("   flux_in = %.5e, flux_out = %.5e, err = %.5e\n", flux_in, flux_out, err);

				
				}
#endif



			
			}

			//err = fabs(1.0 - flux_out_tot/flux_in_tot);
			//if(err > 1.0e-12) {
			//if(1) {
				//printf("voxel %d %d %d\n", grind.i, grind.j, grind.k);
				//printf(" flux_in = %.5e, flux_out = %.5e, err = %.5e\n", flux_in_tot, flux_out_tot, err);
			
			//}

			next_cell: continue;
		}


		//err = fabs(1.0 - ray_flux_out/ray.Ftot);
		//if(err > 1.0e-12) {
			//printf(" ray %d: flux_in = %.5e, flux_out = %.5e, err = %.5e\n", r, ray.Ftot, ray_flux_out, err);
		//}
		

		// set the outgoing ray flux
		//ray.Ftot = 1.0; 
		hp->rays[hp->nrays++] = ray;
	}

}

void reduce_beam_poly(psi_poly* poly, real *domega) {


	// variable declarations
	int f;
	real sdists[32];
	int clipped[32];
	int v, p, nclipped, np, onv, vstart, vcur, vnext, numunclipped; 
	real len;
	rvec meanvec;

	// direct access to vertex buffer
	*domega = 0.0;
	if(poly->nverts <= 0)
		return;

	// compute the final solid angle for this poly
	memset(&meanvec, 0, sizeof(meanvec));
	for(v = 0; v < poly->nverts; ++v) {
		len = sqrt(poly->verts[v].pos.xyz[0]*poly->verts[v].pos.xyz[0]
			+poly->verts[v].pos.xyz[1]*poly->verts[v].pos.xyz[1]+poly->verts[v].pos.xyz[2]*poly->verts[v].pos.xyz[2]);
		poly->verts[v].pos.xyz[0] /= len;
		poly->verts[v].pos.xyz[1] /= len;
		poly->verts[v].pos.xyz[2] /= len;
		meanvec.xyz[0] += poly->verts[v].pos.xyz[0];
		meanvec.xyz[1] += poly->verts[v].pos.xyz[1];
		meanvec.xyz[2] += poly->verts[v].pos.xyz[2];
	}
	len = sqrt(meanvec.x*meanvec.x+meanvec.y*meanvec.y+meanvec.z*meanvec.z);
	meanvec.xyz[0] /= len;
	meanvec.xyz[1] /= len;
	meanvec.xyz[2] /= len;

	for(v = 0; v < poly->nverts; ++v) {
		real mydo = domega3(poly->verts[v].pos, poly->verts[poly->verts[v].pnbrs[1]].pos, meanvec); 
		*domega += mydo; 
	}
}


void clip_beam_poly(psi_poly* poly, psi_plane* faces, int nclip, hydro_problem* hp) {


	// variable declarations
	int f;
	real sdists[32];
	int clipped[32];
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
			sdists[v] = dot3(poly->verts[v].pos, faces[f].n);
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

#if 0
	*domega = 0.0;
	if(poly->nverts) { 

		// compute the final solid angle for this poly
		rvec meanvec;
		memset(&meanvec, 0, sizeof(meanvec));
		for(v = 0; v < poly->nverts; ++v) {
			len = sqrt(poly->verts[v].pos.xyz[0]*poly->verts[v].pos.xyz[0]
				+poly->verts[v].pos.xyz[1]*poly->verts[v].pos.xyz[1]+poly->verts[v].pos.xyz[2]*poly->verts[v].pos.xyz[2]);
			poly->verts[v].pos.xyz[0] /= len;
			poly->verts[v].pos.xyz[1] /= len;
			poly->verts[v].pos.xyz[2] /= len;
			meanvec.xyz[0] += poly->verts[v].pos.xyz[0];
			meanvec.xyz[1] += poly->verts[v].pos.xyz[1];
			meanvec.xyz[2] += poly->verts[v].pos.xyz[2];
		}

		len = sqrt(meanvec.xyz[0]*meanvec.xyz[0]+meanvec.xyz[1]*meanvec.xyz[1]+meanvec.xyz[2]*meanvec.xyz[2]);
		meanvec.xyz[0] /= len;
		meanvec.xyz[1] /= len;
		meanvec.xyz[2] /= len;

		(*centroid).xyz[0] = 0;
		(*centroid).xyz[1] = 0;
		(*centroid).xyz[2] = 0;
		for(v = 0; v < poly->nverts; ++v) {
			real mydo = domega3(poly->verts[v].pos, poly->verts[poly->verts[v].pnbrs[1]].pos, meanvec); 
			*domega += mydo; 
			(*centroid).xyz[0] += mydo*(poly->verts[v].pos.xyz[0] +  poly->verts[poly->verts[v].pnbrs[1]].pos.xyz[0] +  meanvec.xyz[0]);
			(*centroid).xyz[1] += mydo*(poly->verts[v].pos.xyz[1] +  poly->verts[poly->verts[v].pnbrs[1]].pos.xyz[1] +  meanvec.xyz[1]);
			(*centroid).xyz[2] += mydo*(poly->verts[v].pos.xyz[2] +  poly->verts[poly->verts[v].pnbrs[1]].pos.xyz[2] +  meanvec.xyz[2]);
		}
		len = sqrt((*centroid).xyz[0]*(*centroid).xyz[0]+(*centroid).xyz[1]*(*centroid).xyz[1]+(*centroid).xyz[2]*(*centroid).xyz[2]);
		(*centroid).xyz[0] /= len;
		(*centroid).xyz[1] /= len;
		(*centroid).xyz[2] /= len;
	}
#endif
}

	
#if 0
	// 2D case only!
	// TODO: clean this up
	//else if(hp->dim == 2) {
		//for(f = 0; f < faces->nclip; ++f) {
	
			//// split vertices by their distance from the clip plane
			//nclipped = 0;
			//memset(&clipped, 0, sizeof(clipped));
			//for(v = 0; v < 2; ++v) {
				//sdists[v] = dot2(poly->verts[v].pos, faces->clipnorms[f]);
				//clipped[v] = (sdists[v] < 0.0);
				//nclipped += clipped[v]; 
			//}

	
			//// skip this face if the poly lies entirely on one side of it 
			//if(nclipped == 0) continue;
			//if(nclipped == 2) {
				//poly->nverts = 0;
				//return;
			//}
	
			//// otherwise, check all edges and insert new vertices on the bisected edges 
			//for(v = 0; v < 2; ++v) {
				//if(clipped[v]) {
					//wav2(poly->verts[v].pos, -sdists[1-v],
						//poly->verts[1-v].pos, sdists[v],
						//poly->verts[v].pos);
				//} 
			//}
		//}
		//if(poly->nverts && faces->nclip) {
			//for(v = 0; v < 2; ++v) {
				//len = sqrt(poly->verts[v].pos.xyz[0]*poly->verts[v].pos.xyz[0]+poly->verts[v].pos.xyz[1]*poly->verts[v].pos.xyz[1]);
				//poly->verts[v].pos.xyz[0] /= len;
				//poly->verts[v].pos.xyz[1] /= len;
			//}
			//*domega = asin(cross2(poly->verts[0].pos, poly->verts[1].pos));
		//}
	//}
	
	
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
						if(fin == 0) r_in = x0.x/centroid_out[0]; // x-facing faces
						else if(fin == 1) r_in = x0.y/centroid_out[1]; // y-facing faces
						else if(fin == 2) r_in = x0.z/centroid_out[2]; // z-facing faces
						else if(fin == hp->dim) r_in = ray.rmin; // incoming ray front 
						if(rminbits == 0xFFFF) r_in = 0.0; // TODO: why is this needed??
						if(fout == 0) r_out = (x0.x+hp->dx)/centroid_out[0];
						else if(fout == 1) r_out = (x0.y+hp->dx)/centroid_out[1];
						else if(fout == 2) r_out = (x0.z+hp->dx)/centroid_out[2];
						else if(fout == hp->dim) r_out = ray.rmax; //outgoing ray front 
						//printf("Fmid_in = %.5e\n", fmid_in);

						//r_in = 0.0; r_out = 1.0;
	
						// get the correct inner and outer radii, then attenuate exponentially 
						//rmid[0] = 0.5*(outpoly.verts[0].pos[0]+outpoly.verts[1].pos[0]);
						//rmid[1] = 0.5*(outpoly.verts[0].pos[1]+outpoly.verts[1].pos[1]);
	
						// update grid hydro only if we are in a non-ghost cell

						fmid_out = fmid_in*exp(-k*hp->grid[flatind].rho*(r_out-r_in));

						if(grind.i >= 0 && grind.i < hp->nx.i
								&& grind.j >= 0 && grind.j < hp->nx.j) {
	
							// calculate the outgoing flux, write out the radiation energy density 
							// TODO: interpolate the mean value using the exponential
							hp->rad_grid[flatind].E += 0.5*(fmid_in+fmid_out)*(r_out-r_in)/(ray.rmax-ray.rmin)/(hp->dx*hp->dx);
							//hp->rad_grid[flatind].E += fmid_in-fmid_out;
		
							// update hydro terms by conserving the difference
							// between downwind and upwind fluxes
							//hp->grid[flatind].etot += fmid_in-fmid_out; 
		psi_reduce(&beam_poly, moments, 0, 0);
		printf("Original beam, mom[0] = %.5e:\n", moments[0]*hp->nrays);
							//hp->grid[flatind].mom[0] += lsgn.i*cos(thmid)*(fmid_in-fmid_out)/CLIGHT; // TODO: review this!
							//hp->grid[flatind].mom[1] += lsgn.j*sin(thmid)*(fmid_in-fmid_out)/CLIGHT; 
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
				next_cell: continue;
			}
		}
#endif




