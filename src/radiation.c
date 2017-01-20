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

void hypy_get_face_flux_transfers(real* matx, int rminbits, int rmaxbits, real* rmin2, real* rmax2, rvec x0, hydro_ray ray, hydro_problem* hp) {

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
	if(ftmp.thmin < ray.thmax && ftmp.thmax > ray.thmin) {	\
		ftmp.dim_index = ax;	\
		ftmp.thmax = fmin(ftmp.thmax, ray.thmax);	\
		ftmp.thmin = fmax(ftmp.thmin, ray.thmin);	\
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


#define BASERAYS 28 //20 

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
		ray.rmin = ray.rmax;
		ray.rmax += CLIGHT*dt;
		hp->rays[r] = ray;
	}

	// STEP 2: source new rays 
	for(r = 0; r < BASERAYS && hp->nrays < BASERAYS; ++r) {
		hydro_ray ray;
		ray.rmin = 0.0;
		ray.rmax = CLIGHT*dt;
		ray.thmin = r*TWO_PI/BASERAYS;
		ray.thmax = (r+1)*TWO_PI/BASERAYS;
		ray.N = 1.0;//(TWO_PI/8)/(BASERAYS/4);
		for(i = 0; i < hp->dim; ++i)
			ray.orcell[i] = hp->nx[i]/2;
		if((r == 2 || r == 3) && hp->nrays < BASERAYS/4) //BASERAYS)
			hp->rays[hp->nrays++] = ray;
	}
	if(hp->nrays > 1024) printf("Error! Overflowed ray buffer.\n");

	// TODO: test only -- make a separate radiation energy field!
	memset(hp->grid, 0, (hp->nx[0]+4)*(hp->nx[1]+4)*sizeof(hydro_vector));
	for(r = 0; r < hp->nrays; ++r) {

		// get the ray and find its origin
		ray = hp->rays[r];
		origin[0] = hp->dx*(0.5+ray.orcell[0]);
		origin[1] = hp->dx*(0.5+ray.orcell[1]);

		// 2D only for now...
		if(hp->dim == 2) {

			// TODO: figure out the loop symmetry!!!
//#define sgn(x) (((x) < 0)? 1 : 0)
			//loopdir[0] = sgn(cos(0.5*(ray.thmin+ray.thmax)));
			//loopdir[1] = sgn(sin(0.5*(ray.thmin+ray.thmax)));

			// calculate the ray vertex positions and get the bounding box
			verts[0][0] = origin[0] + ray.rmin*cos(ray.thmax);
			verts[0][1] = origin[1] + ray.rmin*sin(ray.thmax);
			verts[1][0] = origin[0] + ray.rmax*cos(ray.thmax);
			verts[1][1] = origin[1] + ray.rmax*sin(ray.thmax);
			verts[2][0] = origin[0] + ray.rmax*cos(ray.thmin);
			verts[2][1] = origin[1] + ray.rmax*sin(ray.thmin);
			verts[3][0] = origin[0] + ray.rmin*cos(ray.thmin); 
			verts[3][1] = origin[1] + ray.rmin*sin(ray.thmin); 
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
				flatind = flat_index(grind, hp);
				x0[0] = hp->dx*grind[0]-origin[0];
				x0[1] = hp->dx*grind[1]-origin[1];


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
					if(ray.thmax < TWO_PI/8)
						face_transfer_matrix[(hp->dim+1)*hp->dim + 0] = ray.thmax-ray.thmin; 
					else if(ray.thmin > TWO_PI/8)
						face_transfer_matrix[(hp->dim+1)*hp->dim + 1] = ray.thmax-ray.thmin; 
					else {
						face_transfer_matrix[(hp->dim+1)*hp->dim + 0] = TWO_PI/8-ray.thmin;
						face_transfer_matrix[(hp->dim+1)*hp->dim + 1] = ray.thmax-TWO_PI/8;
					}
				}
				else {
					// TODO: simplify this function call...
					hypy_get_face_flux_transfers(face_transfer_matrix, rminbits, rmaxbits, 
						rmin2, rmax2, x0, ray, hp);
				}

				// get the incoming delta-omega to get the total
				// flux from the incoming beam front 
				for(i = 0, domega = 0.0; i < 3; ++i)
					domega += face_transfer_matrix[(hp->dim+1)*hp->dim + i];
				voxdata[ii][jj].flux_in[hp->dim] += ray.N*domega/(ray.thmax-ray.thmin);

				// calculate the flux transfers between each in, out pair of cell faces
				// or to the exiting beam front 
				for(axin = 0; axin <= hp->dim; ++axin) {
					flux_in = voxdata[ii][jj].flux_in[axin];
					if(flux_in <= 0.0) continue;
					for(i = 0, domega = 0.0; i < 3; ++i)
						domega += face_transfer_matrix[(hp->dim+1)*axin + i];
					invdomega = 1.0/domega;
					voxdata[ii+1][jj].flux_in[0] += flux_in*invdomega*face_transfer_matrix[(hp->dim+1)*axin + 0];
					voxdata[ii][jj+1].flux_in[1] += flux_in*invdomega*face_transfer_matrix[(hp->dim+1)*axin + 1];
					ray_out += flux_in*invdomega*face_transfer_matrix[(hp->dim+1)*axin + hp->dim];

					// TODO: this is just debugging for now, but keep track of the total
					// radiation energy density!
					hp->grid[flatind].rho += flux_in*(1.0 - invdomega*face_transfer_matrix[(hp->dim+1)*axin + hp->dim]);//*invdomega;
				}
			}

			// lastly, update the ray's energy using the total outgoing flux
			printf("ray_in = %f, ray_out = %f, err = %.5e\n", hp->rays[r].N, ray_out, fabs(1.0-ray_out/hp->rays[r].N));
			ray.N = ray_out;
			hp->rays[r] = ray;
		}
		else {
			printf("2D only for now!\n");
			exit(0);
		}
	}

	// TODO: compact rays if they have been removed!
}
