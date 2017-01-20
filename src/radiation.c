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


int mypopc(int i) {
	 i = i - ((i >> 1) & 0x55555555);
	 i = (i & 0x33333333) + ((i >> 2) & 0x33333333);
	 return (((i + (i >> 4)) & 0x0F0F0F0F) * 0x01010101) >> 24;
}

#include "r2d.h"
#include "r3d.h"

typedef struct {
	dvec grind;
	real depth;
	real flux_in;
} vd;

int cmpfunc (const void* a, const void* b) {
	if((*(vd*)a).depth < (*(vd*)b).depth) return -1;
	if((*(vd*)a).depth == (*(vd*)b).depth) return 0;
	if((*(vd*)a).depth > (*(vd*)b).depth) return 1;
}

void hypy_get_face_flux_transfers(real* matx, int rminbits, int rmaxbits, real* rmin2, real* rmax2, rvec x0, real clipmin, real clipmax, hydro_problem* hp) {

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



#define clip_face(v0, v1, ax, faces, nf) {	\
	ftmp.thmin = atan2(v0[1], v0[0]);	\
	ftmp.thmax = atan2(v1[1], v1[0]);	\
	/*if(ftmp.thmin - ftmp.thmax >= TWO_PI) ftmp.thmax += TWO_PI;*/	\
	/*if(ftmp.thmax - ftmp.thmin >= TWO_PI) ftmp.thmin += TWO_PI;*/	\
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

	// debugging, make sure angles all sum up
	if(nfin > 0 && nfout > 0) {
		real doin, doout;
		if(doin > 1.0e-10 && fabs(1.0-doout/doin) > 1.0e-10) {
			printf("-- theta_in = %.5e, theta_out = %.5e\n", doin, doout);
			printf("Angle error! %.5e\n", fabs(1.0-doout/doin));
			exit(0);
		}
	}

	// now compute pairwise intersections between in, out faces
	for(fin = 0; fin < nfin; ++fin) {
		ftin = faces_in[fin];
		for(fout = 0; fout < nfout; ++fout) {
			ftout = faces_out[fout];
			if(ftout.thmin < ftin.thmax && ftout.thmax > ftin.thmin) {
				th0 = fmax(ftin.thmin, ftout.thmin);
				th1 = fmin(ftin.thmax, ftout.thmax);
				if(th1-th0 > 0.0)
					matx[(hp->dim+1)*ftin.dim_index+ftout.dim_index] = th1-th0;
			}
		}
	}
}


#define BASERAYS 128 //64 

void update_radiation(real dt, hydro_problem* hp) {

	int i, j, v, r, itmp, axin, rminbits, rmaxbits,
		flatind, nvox, locind, ii, jj, axout;
	dvec loopdir, ibox[2], nbox, grind;
	real rmin, rmax, depth, flux_in, fangles[2],
		thmin, thmax, th0, th1, dotot, ray_out, domega, invdomega;
	real rmin2[4], rmax2[4], domegas[4];
	real face_transfer_matrix[(hp->dim+1)*(hp->dim+1)];
	rvec verts[4], origin, v0, v1, x0;
	hydro_ray ray;

	// STEP 1: step all existing rays forward
	for(r = 0; r < hp->nrays; ++r) {
		hydro_ray ray = hp->rays[r];
		ray.rmin += CLIGHT*dt;
		ray.rmax += CLIGHT*dt;
		hp->rays[r] = ray;
	}

	// STEP 2: source new rays 
	for(r = 0; r < BASERAYS; ++r) {
		hydro_ray ray;
		ray.rmin = 0.0;
		ray.rmax = CLIGHT*dt;
		ray.angle_id = r;
		ray.N = 1.0;
		for(i = 0; i < hp->dim; ++i)
			ray.orcell[i] = hp->nx[i]/2;
		//if(r < BASERAYS/4)
			hp->rays[hp->nrays++] = ray;
	}
	if(hp->nrays > 8192) printf("Error! Overflowed ray buffer.\n");

	// TODO: test only -- make a separate radiation energy field!
	memset(hp->rad_grid, 0, (hp->nx[0]+4)*(hp->nx[1]+4)*sizeof(rad_vector));
	int ornrays = hp->nrays;
	hp->nrays = 0;
	for(r = 0; r < ornrays; ++r) {

		// get the ray and find its origin
		ray = hp->rays[r];
		origin[0] = hp->dx*(0.5+ray.orcell[0]);
		origin[1] = hp->dx*(0.5+ray.orcell[1]);

		// 2D only for now...
		if(hp->dim == 2) {


			int quadrant = ray.angle_id&((1<<(hp->dim+1))-1);
			int qid = ray.angle_id>>hp->dim;

			printf("Ray id = %d, quadrant = %d, qid = %d\n", r, quadrant, qid);

			// get the ray angles from the ID
			thmin = (TWO_PI*qid)/BASERAYS;
			thmax = (TWO_PI*(qid+1))/BASERAYS;

			// TODO: figure out the loop symmetry!!!
			loopdir[0] = (quadrant>>0)&1; 
			loopdir[1] = (quadrant>>1)&1; 

			// calculate the ray vertex positions and get the bounding box
			verts[0][0] = origin[0] + ray.rmin*cos(thmax);
			verts[0][1] = origin[1] + ray.rmin*sin(thmax);
			verts[1][0] = origin[0] + ray.rmax*cos(thmax);
			verts[1][1] = origin[1] + ray.rmax*sin(thmax);
			verts[2][0] = origin[0] + ray.rmax*cos(thmin);
			verts[2][1] = origin[1] + ray.rmax*sin(thmin);
			verts[3][0] = origin[0] + ray.rmin*cos(thmin); 
			verts[3][1] = origin[1] + ray.rmin*sin(thmin); 
			for(i = 0; i < hp->dim; ++i) {
				ibox[0][i] = 100000000;
				ibox[1][i] = -100000000;
				for(v = 0; v < 4; ++v) {
					itmp = floor(verts[v][i]/hp->dx);
					if(itmp < ibox[0][i]) ibox[0][i] = itmp;
					if(itmp+1 > ibox[1][i]) ibox[1][i] = itmp+1;
				
				}
				nbox[i] = ibox[1][i] - ibox[0][i];
			}

			// Hacky BCs. Figure this out better
			if(ibox[0][0] < 0) continue; 
			if(ibox[0][1] < 0)  continue; 
			if(ibox[1][0] >= hp->nx[0])  continue; 
			if(ibox[1][1] >= hp->nx[1])  continue; 

			// allocate an array of flux information
			// fill in grid indices and polar coordinates 
			struct {
				real flux_in[4];
				real r2;
			} voxdata[nbox[0]+1][nbox[1]+1];
			memset(voxdata, 0, sizeof(voxdata));
			for(ii = 0; ii <= nbox[0]; ++ii)
			for(jj = 0; jj <= nbox[1]; ++jj) {
				x0[0] = hp->dx*(ibox[0][0]+ii)-origin[0];
				x0[1] = hp->dx*(ibox[0][1]+jj)-origin[1];
				voxdata[ii][jj].r2 = x0[0]*x0[0] + x0[1]*x0[1];
			}

			// loop over all rays, processing their geometry now
			ray_out = 0.0;
			for(ii = 0; ii < nbox[0]; ++ii)
			for(jj = 0; jj < nbox[1]; ++jj) {
				grind[0] = ibox[0][0] + ii;
				grind[1] = ibox[0][1] + jj;
				x0[0] = hp->dx*grind[0]-origin[0];
				x0[1] = hp->dx*grind[1]-origin[1];
				//grind[0] = (grind[0]+hp->nx[0])%hp->nx[0];
				//grind[1] = (grind[1]+hp->nx[1])%hp->nx[1];
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
				memset(face_transfer_matrix, 0, sizeof(face_transfer_matrix));
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
					hypy_get_face_flux_transfers(face_transfer_matrix, rminbits, rmaxbits, 
						rmin2, rmax2, x0, thmin, thmax, hp);
				}

				// get the incoming delta-omega to get the total
				// flux from the incoming beam front 
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


					real transmit = exp(-hp->grid[flatind].rho);

					voxdata[ii+1][jj].flux_in[0] += flux_in*transmit*invdomega*face_transfer_matrix[(hp->dim+1)*axin + 0];
					voxdata[ii][jj+1].flux_in[1] += flux_in*transmit*invdomega*face_transfer_matrix[(hp->dim+1)*axin + 1];
					ray_out += flux_in*transmit*invdomega*face_transfer_matrix[(hp->dim+1)*axin + hp->dim];

					// TODO: this is just debugging for now, but keep track of the total
					// radiation energy density!
					hp->rad_grid[flatind].E += flux_in*(1.0 - invdomega*face_transfer_matrix[(hp->dim+1)*axin + hp->dim]);///[>invdomega;

					// keep track of the outgoing flux;
					real flux_out = 0.0; 
					flux_out += flux_in*transmit*invdomega*face_transfer_matrix[(hp->dim+1)*axin + 0];
					flux_out += flux_in*transmit*invdomega*face_transfer_matrix[(hp->dim+1)*axin + 1];
					flux_out += flux_in*transmit*invdomega*face_transfer_matrix[(hp->dim+1)*axin + hp->dim];



					// update hydro terms by conserving the difference
					// between downwind and upwind fluxes
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
