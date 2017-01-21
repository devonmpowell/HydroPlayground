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

void hypy_get_face_flux_transfers(real* matx_angle, real* matx_depth, int rminbits, int rmaxbits, real* rmin2, real* rmax2, rvec x0, real clipmin, real clipmax, hydro_ray ray, hydro_problem* hp) {

	int fin, fout, nfin, nfout;
	real th0, th1;
	rvec tmpverts[4]; 
	typedef struct {
		real thmin, thmax;
		int dim_index;
	} face_info;
	face_info ftmp, ftin, ftout;
	face_info faces_in[4];
	face_info faces_out[4];

	// TODO: do this without the arctangents!
#define clip_face(v0, v1, ax, faces, nf) {	\
	ftmp.thmin = atan2(v0[1], v0[0]);	\
	ftmp.thmax = atan2(v1[1], v1[0]);	\
	if(ftmp.thmin < clipmax && ftmp.thmax > clipmin) {	\
		ftmp.dim_index = ax;	\
		ftmp.thmax = fmin(ftmp.thmax, clipmax);	\
		ftmp.thmin = fmax(ftmp.thmin, clipmin);	\
		if(ftmp.thmax > ftmp.thmin)	\
			faces[nf++] = ftmp;	\
	}	\
}

	// skip if the cell is outside of the beam
	if(!(rminbits && rmaxbits))
		return;

	// incoming ray cases
	nfin = 0;
	switch(rminbits) {
		case 0x08:
			tmpverts[1][0] = x0[0] + hp->dx; 
			tmpverts[1][1] = sqrt(x0[1]*x0[1]-rmin2[1]); 
			tmpverts[2][0] = sqrt(x0[0]*x0[0]-rmin2[2]); 
			tmpverts[2][1] = x0[1]+hp->dx; 
			clip_face(tmpverts[1], tmpverts[2], hp->dim, faces_in, nfin);
			break;
		case 0x0A:
			tmpverts[0][0] = x0[0] + hp->dx; 
			tmpverts[0][1] = x0[1]; 
			tmpverts[1][0] = sqrt(x0[0]*x0[0]-rmin2[0]); 
			tmpverts[1][1] = x0[1]; 
			tmpverts[2][0] = sqrt(x0[0]*x0[0]-rmin2[2]); 
			tmpverts[2][1] = x0[1] + hp->dx; 
			clip_face(tmpverts[0], tmpverts[1], 1, faces_in, nfin);
			clip_face(tmpverts[1], tmpverts[2], hp->dim, faces_in, nfin);
			break;
		case 0x0C:
			tmpverts[1][0] = x0[0] + hp->dx; 
			tmpverts[1][1] = sqrt(x0[1]*x0[1]-rmin2[1]); 
			tmpverts[2][0] = x0[0];
			tmpverts[2][1] = sqrt(x0[1]*x0[1]-rmin2[0]); 
			tmpverts[3][0] = x0[0]; 
			tmpverts[3][1] = x0[1] + hp->dx; 
			clip_face(tmpverts[1], tmpverts[2], hp->dim, faces_in, nfin);
			clip_face(tmpverts[2], tmpverts[3], 0, faces_in, nfin);
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
			clip_face(tmpverts[0], tmpverts[1], 1, faces_in, nfin);
			clip_face(tmpverts[1], tmpverts[2], hp->dim, faces_in, nfin);
			clip_face(tmpverts[2], tmpverts[3], 0, faces_in, nfin);
			break;
		case 0x0F:
			tmpverts[0][0] = x0[0] + hp->dx; 
			tmpverts[0][1] = x0[1];
			tmpverts[1][0] = x0[0];
			tmpverts[1][1] = x0[1];
			tmpverts[2][0] = x0[0];
			tmpverts[2][1] = x0[1] + hp->dx;
			clip_face(tmpverts[0], tmpverts[1], 1, faces_in, nfin);
			clip_face(tmpverts[1], tmpverts[2], 0, faces_in, nfin);
			break;
		default: 
			break;
	}

	// outgoing ray cases
	nfout = 0;
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
			clip_face(tmpverts[0], tmpverts[1], 0, faces_out, nfout);
			clip_face(tmpverts[1], tmpverts[2], hp->dim, faces_out, nfout);
			clip_face(tmpverts[2], tmpverts[3], 1, faces_out, nfout);
			break;
		case 0x05:
			tmpverts[1][0] = sqrt(x0[0]*x0[0]+rmax2[0]); 
			tmpverts[1][1] = x0[1]; 
			tmpverts[2][0] = sqrt(x0[0]*x0[0]+rmax2[2]); 
			tmpverts[2][1] = x0[1] + hp->dx; 
			tmpverts[3][0] = x0[0]; 
			tmpverts[3][1] = x0[1] + hp->dx; 
			clip_face(tmpverts[1], tmpverts[2], hp->dim, faces_out, nfout);
			clip_face(tmpverts[2], tmpverts[3], 1, faces_out, nfout);
			break;
		case 0x03:
			tmpverts[0][0] = x0[0] + hp->dx; 
			tmpverts[0][1] = x0[1]; 
			tmpverts[1][0] = x0[0] + hp->dx; 
			tmpverts[1][1] = sqrt(x0[1]*x0[1]+rmax2[1]); 
			tmpverts[2][0] = x0[0];
			tmpverts[2][1] = sqrt(x0[1]*x0[1]+rmax2[0]); 
			clip_face(tmpverts[0], tmpverts[1], 0, faces_out, nfout);
			clip_face(tmpverts[1], tmpverts[2], hp->dim, faces_out, nfout);
			break;
		case 0x01:
			tmpverts[1][0] = sqrt(x0[0]*x0[0]+rmax2[0]); 
			tmpverts[1][1] = x0[1]; 
			tmpverts[2][0] = x0[0]; 
			tmpverts[2][1] = sqrt(x0[1]*x0[1]+rmax2[0]); 
			clip_face(tmpverts[1], tmpverts[2], hp->dim, faces_out, nfout);
			break;
		case 0x0F:
			tmpverts[0][0] = x0[0] + hp->dx; 
			tmpverts[0][1] = x0[1];
			tmpverts[1][0] = x0[0] + hp->dx;
			tmpverts[1][1] = x0[1] + hp->dx;
			tmpverts[2][0] = x0[0];
			tmpverts[2][1] = x0[1] + hp->dx;
			clip_face(tmpverts[0], tmpverts[1], 0, faces_out, nfout);
			clip_face(tmpverts[1], tmpverts[2], 1, faces_out, nfout);
			break;
		default: 
			break;
	}

	// now compute pairwise intersections between in, out faces
	for(fin = 0; fin < nfin; ++fin) {
		ftin = faces_in[fin];
		for(fout = 0; fout < nfout; ++fout) {
			ftout = faces_out[fout];
			if(ftout.thmin < ftin.thmax && ftout.thmax > ftin.thmin) {
				th0 = fmax(ftin.thmin, ftout.thmin);
				th1 = fmin(ftin.thmax, ftout.thmax);
				if(th1-th0 > 0.0) { // TODO: is this "if" needed??
					real r_in, r_out;
					switch(ftin.dim_index) {
						case 0: // x-facing faces
							r_in = -x0[0]/(th1-th0)*log(((1+tan(0.5*th0))/(1-tan(0.5*th0)))/((1+tan(0.5*th1))/(1-tan(0.5*th1))));
							break;
						case 1: // y-facing faces
							r_in = -x0[1]/(th1-th0)*log(tan(0.5*th0)/tan(0.5*th1));
							break;
						case 2://hp->dim: // inner radius
							r_in = ray.rmin; 
							break;
						default: 
							r_in = 0.0;
							break;	
					}
					switch(ftout.dim_index) {
						case 0: // x-facing faces
							r_out = -(x0[0]+hp->dx)/(th1-th0)*log(((1+tan(0.5*th0))/(1-tan(0.5*th0)))/((1+tan(0.5*th1))/(1-tan(0.5*th1))));
							break;
						case 1: // y-facing faces
							r_out = -(x0[1]+hp->dx)/(th1-th0)*log(tan(0.5*th0)/tan(0.5*th1));
							break;
						case 2:// hp->dim: // inner radius
							r_out = ray.rmax; 
							break;
						default: 
							r_out = 0.0;
							break;	
					}
					matx_depth[(hp->dim+1)*ftin.dim_index+ftout.dim_index] = r_out-r_in;
					matx_angle[(hp->dim+1)*ftin.dim_index+ftout.dim_index] = th1-th0;
				}
			}
		}
	}
	//memset(matx_depth, 0, 9*sizeof(real));
}

void update_radiation(real dt, hydro_problem* hp) {

	int i, j, v, r, itmp, axin, rminbits, rmaxbits,
		flatind, nvox, locind, ii, jj, axout, quadrant, qid, rlvl;
	dvec lsgn, ibox[2], nbox, grind;
	real rmin, rmax, depth, flux_in, fangles[2],
		thmin, thmax, th0, th1, dotot, ray_out, domega, invdomega;
	real rmin2[4], rmax2[4], domegas[4];
	real face_transfer_matrix[16];
	real face_depth_matrix[16];
	rvec v0, v1, x0;
	hydro_ray ray, r0, r1;

	// STEP 1: step all existing rays forward and refine if needed
	int ornrays = hp->nrays;
	for(r = 0; r < ornrays; ++r) {
		ray = hp->rays[r];
		ray.rmin += CLIGHT*dt;
		ray.rmax += CLIGHT*dt;
		quadrant = (ray.angle_id>>24)&0xFF;
		rlvl = (ray.angle_id>>16)&0xFF;
		qid = (ray.angle_id>>0)&0xFFFF;
		thmin = (TWO_PI*qid)/(1<<(hp->dim+rlvl));
		thmax = (TWO_PI*(qid+1))/(1<<(hp->dim+rlvl));
		if(ray.rmax*(thmax-thmin) > hp->dx*sqrt(hp->dim)) {
			r0 = ray; r1 = ray;
			r0.N *= 0.5; r1.N *= 0.5;
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
	}

	// STEP 2: source new rays 
	#define BASELVL 3
	for(r = 0; r < (1<<(hp->dim+BASELVL)); ++r) {
		ray.rmin = 0.0;
		ray.rmax = CLIGHT*dt;
		ray.angle_id = (r&((1<<hp->dim)-1))<<24; // set the quadrant (highest byte)
		ray.angle_id |= BASELVL<<16; // set the refinement level
		ray.angle_id |= (r>>hp->dim)&0xFFFF; // set the ray id (lowest two bytes) 
		ray.N = 1.0e4*dt;
		for(i = 0; i < hp->dim; ++i)
			ray.orcell[i] = hp->nx[i]/2;
		hp->rays[hp->nrays++] = ray;
	}
	if(hp->nrays > 160000) printf("Error! Overflowed ray buffer.\n");

	// TODO: test only -- make a separate radiation energy field!
	memset(hp->rad_grid, 0, (hp->nx[0]+4)*(hp->nx[1]+4)*sizeof(rad_vector));
	ornrays = hp->nrays;
	hp->nrays = 0;
	for(r = 0; r < ornrays; ++r) {
		ray = hp->rays[r];

		// 2D only for now...
		if(hp->dim == 2) {

			// get the ray angles and grid quadrants from the ID
			// also get the ray's bounding box, relative to the originating cell 
			// this works because we work in one single dummy quadrant
			quadrant = (ray.angle_id>>24)&0xFF;
			rlvl = (ray.angle_id>>16)&0xFF;
			qid = (ray.angle_id>>0)&0xFFFF;
			thmin = (TWO_PI*qid)/(1<<(hp->dim+rlvl));
			thmax = (TWO_PI*(qid+1))/(1<<(hp->dim+rlvl));
			lsgn[0] = 1-2*((quadrant>>0)&1);
			lsgn[1] = 1-2*((quadrant>>1)&1);
			ibox[0][0] = floor((0.5*hp->dx + ray.rmin*cos(thmax))/hp->dx);
			ibox[0][1] = floor((0.5*hp->dx + ray.rmin*sin(thmin))/hp->dx);
			ibox[1][0] = ceil((0.5*hp->dx + ray.rmax*cos(thmin))/hp->dx);
			ibox[1][1] = ceil((0.5*hp->dx + ray.rmax*sin(thmax))/hp->dx);
			nbox[0] = ibox[1][0] - ibox[0][0];
			nbox[1] = ibox[1][1] - ibox[0][1];

			// Hacky BCs. TODO: Figure this out better
			if(ibox[0][0] + ray.orcell[0] < 0) continue; 
			if(ibox[0][1] + ray.orcell[1] < 0)  continue; 
			if(ibox[1][0] + ray.orcell[0] >= hp->nx[0])  continue; 
			if(ibox[1][1] + ray.orcell[1] >= hp->nx[1])  continue; 

			// allocate an array of flux information
			// fill in grid indices and polar coordinates 
			// then loop over all rays, processing their geometry and fluxes
			struct {
				real flux_in[4];
				real r2;
			} voxdata[nbox[0]+1][nbox[1]+1];
			memset(voxdata, 0, sizeof(voxdata));
			for(ii = 0; ii <= nbox[0]; ++ii)
			for(jj = 0; jj <= nbox[1]; ++jj) {
				x0[0] = hp->dx*(ibox[0][0]+ii-0.5);
				x0[1] = hp->dx*(ibox[0][1]+jj-0.5);
				voxdata[ii][jj].r2 = x0[0]*x0[0] + x0[1]*x0[1];
			}
			ray_out = 0.0;
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

				// get the flux transfer matrix
				// special case for the origin cell
				// TODO: make this not a special case
				memset(face_transfer_matrix, 0, sizeof(face_transfer_matrix));
				memset(face_depth_matrix, 0, sizeof(face_depth_matrix));
				if(ray.orcell[0] == grind[0] && ray.orcell[1] == grind[1]) {
					if(thmax < TWO_PI/8)
						face_transfer_matrix[(hp->dim+1)*hp->dim + 0] = thmax-thmin; 
					else if(thmin > TWO_PI/8)
						face_transfer_matrix[(hp->dim+1)*hp->dim + 1] = thmax-thmin; 
					else {
						face_transfer_matrix[(hp->dim+1)*hp->dim + 0] = TWO_PI/8-thmin;
						face_transfer_matrix[(hp->dim+1)*hp->dim + 1] = thmax-TWO_PI/8;
					}
				}
				else {
					// TODO: simplify this function call...
					hypy_get_face_flux_transfers(face_transfer_matrix, face_depth_matrix, rminbits, rmaxbits, 
						rmin2, rmax2, x0, thmin, thmax, ray, hp);
				}

				// get the incoming delta-omega to get the total
				// flux from the incoming beam front 
				// TODO: move this into its own function?
				for(i = 0, domega = 0.0; i < 3; ++i)
					domega += face_transfer_matrix[(hp->dim+1)*hp->dim + i];
				voxdata[ii][jj].flux_in[hp->dim] += ray.N*domega/(thmax-thmin);

				// calculate the flux transfers between each in, out pair of cell faces
				// or to the exiting beam front 
				for(axin = 0; axin <= hp->dim; ++axin) {
					flux_in = voxdata[ii][jj].flux_in[axin];
					if(flux_in <= 0.0) continue;
					for(i = 0, domega = 0.0; i < 3; ++i)
						domega += face_transfer_matrix[(hp->dim+1)*axin + i];
					invdomega = 1.0/domega;


					real k = 5.0e2;

					voxdata[ii+1][jj].flux_in[0] += flux_in*exp(-k*face_depth_matrix[(hp->dim+1)*axin+0]
							*hp->grid[flatind].rho)*invdomega*face_transfer_matrix[(hp->dim+1)*axin + 0];
					voxdata[ii][jj+1].flux_in[1] += flux_in*exp(-k*face_depth_matrix[(hp->dim+1)*axin+1]
							*hp->grid[flatind].rho)*invdomega*face_transfer_matrix[(hp->dim+1)*axin + 1];


					ray_out += flux_in*exp(-k*face_depth_matrix[(hp->dim+1)*axin+0]
							*hp->grid[flatind].rho)*invdomega*face_transfer_matrix[(hp->dim+1)*axin + hp->dim];

					// TODO: this is just debugging for now, but keep track of the total
					// radiation energy density!
					hp->rad_grid[flatind].E += flux_in*(1.0 - invdomega*face_transfer_matrix[(hp->dim+1)*axin + hp->dim]);///[>invdomega;

					// keep track of the outgoing flux;
					real flux_out = 0.0; 
					flux_out += flux_in*exp(-k*face_depth_matrix[(hp->dim+1)*axin+0]
							*hp->grid[flatind].rho)*invdomega*face_transfer_matrix[(hp->dim+1)*axin + 0];
					flux_out += flux_in*exp(-k*face_depth_matrix[(hp->dim+1)*axin+1]
							*hp->grid[flatind].rho)*invdomega*face_transfer_matrix[(hp->dim+1)*axin + 1];
					flux_out += flux_in*exp(-k*face_depth_matrix[(hp->dim+1)*axin+hp->dim]
							*hp->grid[flatind].rho)*invdomega*face_transfer_matrix[(hp->dim+1)*axin + hp->dim];



					// update hydro terms by conserving the difference
					// between downwind and upwind fluxes
					// TODO: better volume averaging here...
					hp->grid[flatind].etot += flux_in-flux_out; 
					hp->grid[flatind].mom[0] += cos(0.5*(thmax+thmin))*(flux_in-flux_out)/CLIGHT; 
					hp->grid[flatind].mom[1] += sin(0.5*(thmax+thmin))*(flux_in-flux_out)/CLIGHT; 
				}
			}
			ray.N = ray_out;
			hp->rays[hp->nrays++] = ray;
		}
		else {
			printf("2D only for now!\n");
			exit(0);
		}
	}

	// TODO: compact rays if they have been removed!
}
