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



void update_radiation(real dt, hydro_problem* hp) {

	int i, v, r, ax, rminbits, rmaxbits, fin, fout, ornrays, refine, baseid, reflvl, refbits,
		flatind, quadrant, qid, rlvl, f, nbase, split, cell_good;
	dvec lsgn, grind, nbox, loopind, locind;
						dvec nextcell;
	real ray_mom_out, fmom_in, domega_in, domega_out, intensity_in, dobase, flux_in, omega_out, 
		 ray_flux_out, err, allmin, allmax, r_in, r_out, secthmid, len, ray_omega, r_old, r_new, 
		 cscthmid, inmin, inmax, outmin, outmax, fmid_in, fmid_out;
	real rmin2[4], rmax2[4], flux_out, fmom_out[4], moments[10], dr;
	rvec rmid, tmpverts[8], ntmp, v0, v1;
	hydro_ray ray, r0, r1;
	psi_poly beam_poly, inpoly, outpoly, cell_poly, beamtmp;

	rvec centroid_in, centroid_out;

	// dimension-specific constants
	if(hp->dim == 3) {
		nbase = num_base_3d;
		dobase = 2*TWO_PI/nbase;
	}
	else{
		printf("3D only!\n");
		exit(0);
	}

	real source_power = 1.577e2; // units of 10^60 photons per Myr
	real cross_section = 0;// 6.62e-1; // units of 10^-60 kpc^2
	real alpha_b = 2.78e-4; // recombination rate

	// STEP 1: source new rays 
	for(r = 0; r < nbase; ++r) {
		ray.radius = 0.0;
		ray.flux = dt*source_power/nbase; 
		for(i = 0; i < hp->dim; ++i) // TODO: source rays from their actual cells
			//ray.origin.xyz[i] = ERRTOL*hp->dx;// 0.5*hp->dx*(hp->nx.ijk[i]+1);
			ray.origin.xyz[i] = 0.5*hp->dx*(hp->nx.ijk[i]+1);
		ray.angle_id = ((r&0xFF)<<24); 
		hp->rays[hp->nrays++] = ray;
	}
	if(hp->nrays > 160000) printf("Error! Overflowed ray buffer.\n");

	// Step 3: propagate all rays forward by c*dt
	memset(hp->rad_grid, 0, (hp->nx.i+4)*(hp->nx.j+4)*(hp->nx.k+4)*sizeof(rad_vector));
	for(i = 0; i < (hp->nx.i+4)*(hp->nx.j+4)*(hp->nx.k+4); ++i)
		hp->grid[i].dN = 0.0;
	ornrays = hp->nrays;
	hp->nrays = 0;
	for(r = 0; r < ornrays; ++r) {

		// get the ray angles and grid quadrants from the ID
		// also get the ray's bounding box, relative to the originating cell 
		// this works because we work in one single dummy quadrant
		
		// get the ray and step it forward in time, saving the old and new radii from the source
		ray = hp->rays[r];
		r_old = ray.radius;
		ray.radius += dt*CLIGHT;
		r_new = ray.radius;

		// get the beam polygon and the bounding box
		get_beam_poly(ray, r_old, r_new, &beam_poly, &nbox, &lsgn, hp);
		unpack_id_bits(ray.angle_id, baseid, reflvl, refbits);
		ray_omega = dobase/(1<<reflvl);               
		if(beam_poly.ibox[0].i >= hp->nx.i) continue;
		if(beam_poly.ibox[0].j >= hp->nx.j) continue;
		if(beam_poly.ibox[0].k >= hp->nx.k) continue;
		if(beam_poly.ibox[1].i < 0) continue;
		if(beam_poly.ibox[1].j < 0) continue;
		if(beam_poly.ibox[1].k < 0) continue;

		int nfaces;
		psi_face_buffer faces[16], face_in, face_out;
		psi_voxels vox;
		psi_poly curpoly;

		// translate the poly to an integer offset 
		// based on its "upper-left" corner
		//dvec offset;
		   //for(ax = 0; ax < 3; ++ax) // black magic logic
			//offset.ijk[ax] = beam_poly.ibox[lsgn.ijk[ax]<0].ijk[ax]-2*(lsgn.ijk[ax]<0); 
		//rvec ortmp = ray.origin; 
		//for(ax = 0; ax < 3; ++ax) {
			//beam_poly.ibox[0].ijk[ax] -= offset.ijk[ax];
			//beam_poly.ibox[1].ijk[ax] -= offset.ijk[ax];
			//for(v = 0; v < beam_poly.nverts; ++v)
				//beam_poly.verts[v].pos.xyz[ax] -= hp->dx*offset.ijk[ax];
			//ortmp.xyz[ax] -= hp->dx*offset.ijk[ax];
		//}


		//// make planes for the beam poly and the voxel
		//// these will be useful for finding optical depths
		//rvec tmp0, tmp1;
		psi_plane cell_planes[8];
		//for(ax = 0; ax < 3; ++ax) {
			//tmp0.xyz[ax] = beam_poly.verts[4].pos.xyz[ax] - beam_poly.verts[3].pos.xyz[ax];
			//tmp1.xyz[ax] = beam_poly.verts[5].pos.xyz[ax] - beam_poly.verts[3].pos.xyz[ax];
		//}
		//cross3(cell_planes[BEAM_OUT].n, tmp0, tmp1);
		//cell_planes[BEAM_OUT].d = dot3(cell_planes[BEAM_OUT].n, beam_poly.verts[3].pos);

		//for(ax = 0; ax < 3; ++ax) {
			//tmp0.xyz[ax] = beam_poly.verts[1].pos.xyz[ax] - beam_poly.verts[0].pos.xyz[ax];
			//tmp1.xyz[ax] = beam_poly.verts[2].pos.xyz[ax] - beam_poly.verts[0].pos.xyz[ax];
		//}
		//cross3(cell_planes[BEAM_IN].n, tmp0, tmp1);
		//cell_planes[BEAM_IN].d = dot3(cell_planes[BEAM_IN].n, beam_poly.verts[0].pos);

#if 1

		// voxelize the beam
		ray_flux_out = 0.0;
		struct {
			real flux_in[8]; // TODO: Don't need this many indices
		} transfer_info[nbox.i+1][nbox.j+1][nbox.k+1];
		memset(transfer_info, 0, sizeof(transfer_info));
		for(locind.i = 0; locind.i < nbox.i; ++locind.i)
		for(locind.j = 0; locind.j < nbox.j; ++locind.j)
		for(locind.k = 0; locind.k < nbox.k; ++locind.k) {

			// get the relative and absolute grid indices
			// make sure we're inside the grid bounds before processing
			for(ax = 0; ax < 3; ++ax)
				grind.ijk[ax] = lsgn.ijk[ax]*locind.ijk[ax] + 32;// offset.ijk[ax];
			for(cell_good = 1, ax = 0; ax < hp->dim; ++ax)
				if(grind.ijk[ax] < 0 || grind.ijk[ax] >= hp->nx.ijk[ax])
					cell_good = 0;
			flatind = flat_index(grind, hp);

			if(cell_good)
				hp->rad_grid[flatind].E += 1; 

			printf(" ray %d, grid index %d %d %d\n", r, grind.i, grind.j, grind.k);

			continue;

			// process the solid angle intersections of the faces
			// first, translate the voxel into coordinates centered on the ray source
			int nflag = 0;
			for(v = 0; v < curpoly.nverts; ++v) {
				nflag += curpoly.verts[v].flags;
				for(ax = 0; ax < hp->dim; ++ax)
					curpoly.verts[v].pos.xyz[ax] -= ray.origin.xyz[ax];
			}
			if(r_old == 0.0 && nflag == 3) {
				// TODO: better subgrid model!!
				transfer_info[locind.i+1][locind.j][locind.k].flux_in[0] = ray.flux;
				transfer_info[locind.i][locind.j+1][locind.k].flux_in[2] = ray.flux;
				transfer_info[locind.i][locind.j][locind.k+1].flux_in[4] = ray.flux;
				continue;
			}

			psi_poly poly_in, poly_out;
			real omega_in, omega_in_tot, omega_out_tot;

			omega_in_tot = 0.0;
			omega_out_tot = 0.0;

			// Extract all upwind/downwind faces from the poly
			// Filter out all degenerate edges

			// normalize the face vertices to be like ray normals
			rvec cvert;
			for(v = 0; v < curpoly.nverts; ++v) {
				cvert = curpoly.verts[v].pos;
				len = sqrt(cvert.x*cvert.x+cvert.y*cvert.y+cvert.z*cvert.z);
				for(ax = 0; ax < hp->dim; ++ax)
					curpoly.verts[v].pos.xyz[ax] /= len;
			}
			psi_extract_faces(&curpoly, faces, &nfaces);
			rvec cliptmp;
			int ornf = nfaces;
			nfaces = 0;
			for(fin = 0; fin < ornf; ++fin) {
				face_in = faces[fin];
				int ornv = face_in.nverts;
				face_in.nverts = 0;
				rvec pvert = face_in.verts[ornv-1];
				for(v = 0; v < ornv; ++v) {
					cvert = face_in.verts[v];
					cross3(cliptmp, cvert, pvert); 
					len = sqrt(cliptmp.x*cliptmp.x + cliptmp.y*cliptmp.y + cliptmp.z*cliptmp.z);
					if(len*r_new > ERRTOL*hp->dx)  { // TODO: better tolerance
						face_in.verts[face_in.nverts++] = cvert;
						pvert = cvert;
					}
				}
				if(face_in.nverts >= hp->dim) 
					faces[nfaces++] = face_in;
			}


			real domega;
			rvec centroid;

			real omega_x = 0.0;

			real flux_in_tot = 0.0;
			real flux_out_tot = 0.0;


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
				}
				reduce_beam_poly(&poly_in, &omega_in);
				if(omega_in < ERRTOL*ray_omega) continue;

				// get the incoming beam front (TODO: do this ahead of time??)
				if(face_in.face_id == BEAM_IN)
				   flux_in = ray.flux*omega_in/ray_omega;
				else
					flux_in = transfer_info[locind.i][locind.j][locind.k].flux_in[face_in.face_id];
				if(flux_in <= 0.0) continue;
				intensity_in = flux_in/omega_in;

				real omega_out_tot = 0.0;
				for(fout = 0; fout < nfaces; ++fout) {
					face_out = faces[fout];
					if(face_out.face_id%2 == 0) continue;



					// TODO: make these clip planes ahead of time
					psi_plane clip_out[16];
					memset(clip_out, 0, sizeof(clip_out));
					int nclip = 0; 
					rvec pvert = face_out.verts[0];
					rvec cvert;
					for(v = 1; v < face_out.nverts+1; ++v) {
						cvert = face_out.verts[v%face_out.nverts];
						cross3(cliptmp, cvert, pvert); 
						clip_out[nclip++].n = cliptmp;
						pvert = cvert;
					}

					poly_out = poly_in;
					clip_beam_poly(&poly_out, clip_out, nclip, hp); 
					reduce_beam_poly_deluxe(&poly_out, locind, ray.origin, face_in.face_id, face_out.face_id, cell_planes[BEAM_IN], cell_planes[BEAM_OUT], &omega_out, &dr, hp);
					if(omega_out < ERRTOL*omega_in) continue;

					flux_in = intensity_in*omega_out;
					flux_out = flux_in;
					if(cell_good)
						flux_out *= exp(-cross_section*(1.0-hp->grid[flatind].x)*hp->grid[flatind].rho*dr);

					if(face_out.face_id == BEAM_OUT)
						ray_flux_out += flux_out; 
					else {
						memset(&nextcell, 0, sizeof(nextcell));
						nextcell.ijk[face_out.face_id/2] = 1;
						transfer_info[locind.i+nextcell.i][locind.j+nextcell.j][locind.k+nextcell.k].flux_in[face_out.face_id^1] += flux_out; 
					}

					omega_out_tot += omega_out;
					if(cell_good) {
						hp->rad_grid[flatind].E += 0.5*(flux_in+flux_out)*dr/((r_new-r_old)*(hp->dx*hp->dx)); 
						hp->grid[flatind].dN += flux_in-flux_out; 
					}
				}

				// TODO: fix this check!
				//err = fabs(1.0 - omega_out_tot/omega_in);
				//if(err > ERRTOL) {
					//printf("Incorrect solid angle!  face_in = %d, Omega_in = %.5e, omega_out_tot = %.5e, err = %.5e\n", face_in.face_id, omega_in, omega_out_tot, err);
					//exit(0);
				//}
			}
		}


#else


		// voxelize the beam
		ray_flux_out = 0.0;
		struct {
			real flux_in[8]; // TODO: Don't need this many indices
			int touched_by_ax[8];
		} transfer_info[nbox.i+1][nbox.j+1][nbox.k+1];
		memset(transfer_info, 0, sizeof(transfer_info));
		psi_voxels_init(&vox, &beam_poly, lsgn, hp);
		while(psi_voxels_next(&vox, &curpoly)) {

			// get the relative and absolute grid indices
			// make sure we're inside the grid bounds before processing
			locind = curpoly.ibox[0];
			for(ax = 0; ax < 3; ++ax)
				grind.ijk[ax] = locind.ijk[ax] + offset.ijk[ax];
			int cell_good = 1;
			for(ax = 0; ax < hp->dim; ++ax)
				if(grind.ijk[ax] < 0 || grind.ijk[ax] >= hp->nx.ijk[ax])
					cell_good = 0;
			flatind = flat_index(grind, hp);

			// process the solid angle intersections of the faces
			// first, translate the voxel into coordinates centered on the ray source
			int nflag = 0;
			for(v = 0; v < curpoly.nverts; ++v) {
				nflag += curpoly.verts[v].flags;
				for(ax = 0; ax < hp->dim; ++ax)
					curpoly.verts[v].pos.xyz[ax] -= ortmp.xyz[ax];
			}
			if(r_old == 0.0 && nflag == 3) {
				// TODO: better subgrid model!!
				transfer_info[locind.i+1][locind.j][locind.k].flux_in[0] = ray.flux;
				transfer_info[locind.i][locind.j+1][locind.k].flux_in[2] = ray.flux;
				transfer_info[locind.i][locind.j][locind.k+1].flux_in[4] = ray.flux;
				continue;
			}

			psi_poly poly_in, poly_out;
			real omega_in, omega_in_tot, omega_out_tot;

			omega_in_tot = 0.0;
			omega_out_tot = 0.0;

			// Extract all upwind/downwind faces from the poly
			// Filter out all degenerate edges

			// normalize the face vertices to be like ray normals
			rvec cvert;
			for(v = 0; v < curpoly.nverts; ++v) {
				cvert = curpoly.verts[v].pos;
				len = sqrt(cvert.x*cvert.x+cvert.y*cvert.y+cvert.z*cvert.z);
				for(ax = 0; ax < hp->dim; ++ax)
					curpoly.verts[v].pos.xyz[ax] /= len;
			}
			psi_extract_faces(&curpoly, faces, &nfaces);
			rvec cliptmp;
			int ornf = nfaces;
			nfaces = 0;
			for(fin = 0; fin < ornf; ++fin) {
				face_in = faces[fin];
				int ornv = face_in.nverts;
				face_in.nverts = 0;
				rvec pvert = face_in.verts[ornv-1];
				for(v = 0; v < ornv; ++v) {
					cvert = face_in.verts[v];
					cross3(cliptmp, cvert, pvert); 
					len = sqrt(cliptmp.x*cliptmp.x + cliptmp.y*cliptmp.y + cliptmp.z*cliptmp.z);
					if(len*r_new > ERRTOL*hp->dx)  { // TODO: better tolerance
						face_in.verts[face_in.nverts++] = cvert;
						pvert = cvert;
					}
				}
				if(face_in.nverts >= hp->dim) 
					faces[nfaces++] = face_in;
			}


			real domega;
			rvec centroid;

			real omega_x = 0.0;

			real flux_in_tot = 0.0;
			real flux_out_tot = 0.0;


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
				}
				reduce_beam_poly(&poly_in, &omega_in);
				if(omega_in < ERRTOL*ray_omega) continue;

				// get the incoming beam front (TODO: do this ahead of time??)
				if(face_in.face_id == BEAM_IN)
				   flux_in = ray.flux*omega_in/ray_omega;
				else
					flux_in = transfer_info[locind.i][locind.j][locind.k].flux_in[face_in.face_id];
				if(flux_in <= 0.0) continue;
				intensity_in = flux_in/omega_in;

				real omega_out_tot = 0.0;
				for(fout = 0; fout < nfaces; ++fout) {
					face_out = faces[fout];
					if(face_out.face_id%2 == 0) continue;



					// TODO: make these clip planes ahead of time
					psi_plane clip_out[16];
					memset(clip_out, 0, sizeof(clip_out));
					int nclip = 0; 
					rvec pvert = face_out.verts[0];
					rvec cvert;
					for(v = 1; v < face_out.nverts+1; ++v) {
						cvert = face_out.verts[v%face_out.nverts];
						cross3(cliptmp, cvert, pvert); 
						clip_out[nclip++].n = cliptmp;
						pvert = cvert;
					}

					poly_out = poly_in;
					clip_beam_poly(&poly_out, clip_out, nclip, hp); 
					reduce_beam_poly_deluxe(&poly_out, locind, ortmp, face_in.face_id, face_out.face_id, cell_planes[BEAM_IN], cell_planes[BEAM_OUT], &omega_out, &dr, hp);
					if(omega_out < ERRTOL*omega_in) continue;

					flux_in = intensity_in*omega_out;
					flux_out = flux_in;
					if(cell_good)
						flux_out *= exp(-cross_section*(1.0-hp->grid[flatind].x)*hp->grid[flatind].rho*dr);

					if(face_out.face_id == BEAM_OUT)
						ray_flux_out += flux_out; 
					else {
						memset(&nextcell, 0, sizeof(nextcell));
						nextcell.ijk[face_out.face_id/2] = 1;
						transfer_info[locind.i+nextcell.i][locind.j+nextcell.j][locind.k+nextcell.k].flux_in[face_out.face_id^1] += flux_out; 
					}

					omega_out_tot += omega_out;
					if(cell_good) {
						hp->rad_grid[flatind].E += 0.5*(flux_in+flux_out)*dr/((r_new-r_old)*(hp->dx*hp->dx)); 
						hp->grid[flatind].dN += flux_in-flux_out; 
					}
				}

				// TODO: fix this check!
				//err = fabs(1.0 - omega_out_tot/omega_in);
				//if(err > ERRTOL) {
					//printf("Incorrect solid angle!  face_in = %d, Omega_in = %.5e, omega_out_tot = %.5e, err = %.5e\n", face_in.face_id, omega_in, omega_out_tot, err);
					//exit(0);
				//}
			}
		}
#endif

		//err = fabs(1.0 - ray_flux_out/ray.flux);
		//if(err > 1.0e-12) {
			//printf(" ray %d: flux_in = %.5e, flux_out = %.5e, err = %.5e\n", r, ray.flux, ray_flux_out, err);
		//}

		// set the outgoing ray flux
		ray.flux = ray_flux_out;
		hp->rays[hp->nrays++] = ray;
	}

	// update the ionization fractions on the grid
	for(i = 0; i < (hp->nx.i+4)*(hp->nx.j+4)*(hp->nx.k+4); ++i) {
	
		real rhov = hp->grid[i].rho*hp->dx*hp->dx*hp->dx;
		hp->grid[i].x += hp->grid[i].dN/rhov - dt*alpha_b*hp->grid[i].rho*hp->grid[i].x*hp->grid[i].x;
		if(hp->grid[i].x > 1.0+ERRTOL) {
			printf("Oh no! Ionization fraction is too high! %.10f\n", hp->grid[i].x);
			hp->grid[i].x = 1.0;
		}
	}


#if 0
	// TODO: re-instate refinement and merging
	// STEP 1: step all existing rays forward and refine if needed
	ornrays = hp->nrays;
	for(r = 0; r < ornrays; ++r) {
		ray = hp->rays[r];        
		unpack_id_bits(ray.angle_id, baseid, reflvl, refbits);
		ray_omega = dobase/(1<<reflvl);
		split = (ray.radius*ray_omega > 1*hp->dx*hp->dx);
		//split = 0;
		if(split) {
			r0 = ray;
			r0.flux *= 0.5; 
			r0.angle_id = baseid<<24; // set the base ray ID 
			r0.angle_id |= (reflvl+1)<<20; // set the refinement level
			r0.angle_id |= refbits|(0<<reflvl); // set the refinement bits 
			r1 = ray;
			r1.flux *= 0.5;
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
#endif




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

void get_beam_poly(hydro_ray ray, real rmin, real rmax, psi_poly* poly, dvec* nbox, dvec* charloop, hydro_problem* hp) {

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
		poly->verts[0].pos.xyz[ax] = ray.origin.xyz[ax] + rmin*refined_vertex_pos[0].xyz[ax]; 
		poly->verts[1].pos.xyz[ax] = ray.origin.xyz[ax] + rmin*refined_vertex_pos[1].xyz[ax]; 
		poly->verts[2].pos.xyz[ax] = ray.origin.xyz[ax] + rmin*refined_vertex_pos[2].xyz[ax]; 
		poly->verts[3].pos.xyz[ax] = ray.origin.xyz[ax] + rmax*refined_vertex_pos[0].xyz[ax]; 
		poly->verts[4].pos.xyz[ax] = ray.origin.xyz[ax] + rmax*refined_vertex_pos[1].xyz[ax]; 
		poly->verts[5].pos.xyz[ax] = ray.origin.xyz[ax] + rmax*refined_vertex_pos[2].xyz[ax]; 
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


void reduce_beam_poly_deluxe(psi_poly* poly, dvec index, rvec origin, int fin, int fout, psi_plane beam_in, psi_plane beam_out, real *domega, real *dr, hydro_problem* hp) {


	// variable declarations
	int f, ax;
	real sdists[32];
	real r_in, r_out;
	int clipped[32];
	int v, p, nclipped, np, onv, vstart, vcur, vnext, numunclipped; 
	real len;
	rvec meanvec;

	// direct access to vertex buffer
	*domega = 0.0;
	*dr = 0.0;
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

		// solid angle
		real mydo = domega3(poly->verts[v].pos, poly->verts[poly->verts[v].pnbrs[1]].pos, meanvec); 
		*domega += mydo; 

		// center of this piece
		rvec centroid;
		memset(&centroid, 0, sizeof(centroid));
		for(ax = 0; ax < 3; ++ax)
			centroid.xyz[ax] += 0.33333333333333*(poly->verts[v].pos.xyz[ax] + poly->verts[poly->verts[v].pnbrs[1]].pos.xyz[ax] + meanvec.xyz[ax]); 
	
		//// get the correct inner and outer radindex.i, then attenuate exponentially 
		if(fin == 0) r_in = (hp->dx*index.i-origin.x)/centroid.x; // x-facing faces
		else if(fin == 2) r_in = (hp->dx*index.j-origin.y)/centroid.y; // y-facing faces
		else if(fin == 4) r_in = (hp->dx*index.k-origin.z)/centroid.z; // z-facing faces
		else if(fin == BEAM_IN)
			r_in = (beam_in.d-dot3(beam_in.n, origin))/(dot3(beam_in.n, centroid)); 
		//if(rminbits == 0xFFFF) r_in = 0.0; // TODO: why is this needed??

		if(fout == 1) r_out = (hp->dx*(index.i+1)-origin.x)/centroid.x;
		else if(fout == 3) r_out = (hp->dx*(index.j+1)-origin.y)/centroid.y;
		else if(fout == 5) r_out = (hp->dx*(index.k+1)-origin.z)/centroid.z;
		else if(fout == BEAM_OUT) 
			r_out = (beam_out.d-dot3(beam_out.n, origin))/(dot3(beam_out.n, centroid)); 

		*dr += (r_out-r_in)*mydo;
	}

	// properly weight the dr
	*dr /= *domega;
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




