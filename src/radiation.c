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

void ray_info_from_id(int idbits, int* rlvl, int* quadrant, real* verts_out, hydro_problem* hp) {


}

void update_radiation(real dt, hydro_problem* hp) {

	int i, v, r, rminbits, rmaxbits, fin, fout,
		flatind, ii, jj, quadrant, qid, rlvl, f;
	dvec lsgn, nbox, grind, ibox[2];
	real ray_mom_out, fmom_in, thmin, thmax, th0, th1, 
		 ray_flux_out, err, allmin, allmax, r_in, r_out, secthmid, 
		 cscthmid, inmin, inmax, outmin, outmax, thmid, fmid_in, fmid_out;
	real rmin2[4], rmax2[4], flux_out[4], fmom_out[4];
	rvec x0, tmpverts[4];
	hydro_ray ray, r0, r1;
	typedef struct {
		real thmin, thmax;
		real Imin, Imax;
	} face_info;
	face_info ftmp, ftin, ftout, faces_in[4], faces_out[4];

	// TODO: do this without the arctangents!
	#define clip_face(v0, v1, ax, faces) {	\
		ftmp.thmin = atan2(v0[1], v0[0]);	\
		ftmp.thmax = atan2(v1[1], v1[0]);	\
		faces[ax] = ftmp;	\
	}

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
		quadrant = (ray.angle_id>>24)&0xFF;
		rlvl = (ray.angle_id>>16)&0xFF;
		qid = (ray.angle_id>>0)&0xFFFF;
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
	}

	// STEP 2: source new rays 
	#define BASELVL 3
	for(r = 0; r < (1<<(hp->dim+BASELVL)); ++r) {
		ray.rmin = 0.0;
		ray.rmax = CLIGHT*dt;
		ray.Ftot = dt*1.0/(1<<(hp->dim+BASELVL));
		ray.Fcom = ray.Ftot*TWO_PI*(0.5+((r>>hp->dim)&0xFFFF))/(1<<(hp->dim+BASELVL)); // intensity-weighted deviation from the beam center
		for(i = 0; i < hp->dim; ++i)
			ray.orcell[i] = hp->nx[i]/2;
		// TODO: clean this up with macros! Or just use a char[]... 
		ray.angle_id = (r&((1<<hp->dim)-1))<<24; // set the quadrant (highest byte)
		ray.angle_id |= BASELVL<<16; // set the refinement level
		ray.angle_id |= (r>>hp->dim)&0xFFFF; // set the ray id (lowest two bytes) 
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
		if(ibox[0][0] + ray.orcell[0] >= hp->nx[0]+1) continue; 
		if(ibox[0][1] + ray.orcell[1] >= hp->nx[1]+1)  continue; 
		if(ibox[1][0] + ray.orcell[0] < -1)  continue; 
		if(ibox[1][1] + ray.orcell[1] < -1)  continue; 

		// allocate an array of flux information
		// fill in grid indices and polar coordinates 
		// then loop over all rays, processing their geometry and fluxes
		struct {
			real flux_in[4];
			real fmom_in[4];
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
					clip_face(tmpverts[1], tmpverts[2], hp->dim, faces_in);
					inmin = faces_in[hp->dim].thmin;
					inmax = faces_in[hp->dim].thmax;
					break;
				case 0x0A:
					tmpverts[0][0] = x0[0] + hp->dx; 
					tmpverts[0][1] = x0[1]; 
					tmpverts[1][0] = sqrt(x0[0]*x0[0]-rmin2[0]); 
					tmpverts[1][1] = x0[1]; 
					tmpverts[2][0] = sqrt(x0[0]*x0[0]-rmin2[2]); 
					tmpverts[2][1] = x0[1] + hp->dx; 
					clip_face(tmpverts[0], tmpverts[1], 1, faces_in);
					clip_face(tmpverts[1], tmpverts[2], hp->dim, faces_in);
					inmin = faces_in[1].thmin;
					inmax = faces_in[hp->dim].thmax;
					break;
				case 0x0C:
					tmpverts[1][0] = x0[0] + hp->dx; 
					tmpverts[1][1] = sqrt(x0[1]*x0[1]-rmin2[1]); 
					tmpverts[2][0] = x0[0];
					tmpverts[2][1] = sqrt(x0[1]*x0[1]-rmin2[0]); 
					tmpverts[3][0] = x0[0]; 
					tmpverts[3][1] = x0[1] + hp->dx; 
					clip_face(tmpverts[1], tmpverts[2], hp->dim, faces_in);
					clip_face(tmpverts[2], tmpverts[3], 0, faces_in);
					inmin = faces_in[hp->dim].thmin;
					inmax = faces_in[0].thmax;
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
					clip_face(tmpverts[0], tmpverts[1], 1, faces_in);
					clip_face(tmpverts[1], tmpverts[2], hp->dim, faces_in);
					clip_face(tmpverts[2], tmpverts[3], 0, faces_in);
					inmin = faces_in[1].thmin;
					inmax = faces_in[0].thmax;
					break;
				case 0x0F:
					tmpverts[0][0] = x0[0] + hp->dx; 
					tmpverts[0][1] = x0[1];
					tmpverts[1][0] = x0[0];
					tmpverts[1][1] = x0[1];
					tmpverts[2][0] = x0[0];
					tmpverts[2][1] = x0[1] + hp->dx;
					clip_face(tmpverts[0], tmpverts[1], 1, faces_in);
					clip_face(tmpverts[1], tmpverts[2], 0, faces_in);
					inmin = faces_in[1].thmin;
					inmax = faces_in[0].thmax;
					break;
				case 0xFF: // special case for the ray origin
					ftmp.thmin = thmin; 
					ftmp.thmax = thmax; 
					faces_in[hp->dim] = ftmp;
					inmin = faces_in[hp->dim].thmin;
					inmax = faces_in[hp->dim].thmax;
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
					clip_face(tmpverts[0], tmpverts[1], 0, faces_out);
					clip_face(tmpverts[1], tmpverts[2], hp->dim, faces_out);
					clip_face(tmpverts[2], tmpverts[3], 1, faces_out);
					outmin = faces_out[0].thmin;
					outmax = faces_out[1].thmax;
					break;
				case 0x05:
					tmpverts[1][0] = sqrt(x0[0]*x0[0]+rmax2[0]); 
					tmpverts[1][1] = x0[1]; 
					tmpverts[2][0] = sqrt(x0[0]*x0[0]+rmax2[2]); 
					tmpverts[2][1] = x0[1] + hp->dx; 
					tmpverts[3][0] = x0[0]; 
					tmpverts[3][1] = x0[1] + hp->dx; 
					clip_face(tmpverts[1], tmpverts[2], hp->dim, faces_out);
					clip_face(tmpverts[2], tmpverts[3], 1, faces_out);
					outmin = faces_out[hp->dim].thmin;
					outmax = faces_out[1].thmax;
					break;
				case 0x03:
					tmpverts[0][0] = x0[0] + hp->dx; 
					tmpverts[0][1] = x0[1]; 
					tmpverts[1][0] = x0[0] + hp->dx; 
					tmpverts[1][1] = sqrt(x0[1]*x0[1]+rmax2[1]); 
					tmpverts[2][0] = x0[0];
					tmpverts[2][1] = sqrt(x0[1]*x0[1]+rmax2[0]); 
					clip_face(tmpverts[0], tmpverts[1], 0, faces_out);
					clip_face(tmpverts[1], tmpverts[2], hp->dim, faces_out);
					outmin = faces_out[0].thmin;
					outmax = faces_out[hp->dim].thmax;
					break;
				case 0x01:
					tmpverts[1][0] = sqrt(x0[0]*x0[0]+rmax2[0]); 
					tmpverts[1][1] = x0[1]; 
					tmpverts[2][0] = x0[0]; 
					tmpverts[2][1] = sqrt(x0[1]*x0[1]+rmax2[0]); 
					clip_face(tmpverts[1], tmpverts[2], hp->dim, faces_out);
					outmin = faces_out[hp->dim].thmin;
					outmax = faces_out[hp->dim].thmax;
					break;
				case 0x0F:
					tmpverts[0][0] = x0[0] + hp->dx; 
					tmpverts[0][1] = x0[1];
					tmpverts[1][0] = x0[0] + hp->dx;
					tmpverts[1][1] = x0[1] + hp->dx;
					tmpverts[2][0] = x0[0];
					tmpverts[2][1] = x0[1] + hp->dx;
					clip_face(tmpverts[0], tmpverts[1], 0, faces_out);
					clip_face(tmpverts[1], tmpverts[2], 1, faces_out);
					outmin = faces_out[0].thmin;
					outmax = faces_out[1].thmax;
					break;
				default: 
					printf("Bad outgoing bit mask! How can this be?? 0x%x\n", rmaxbits);
					exit(0);
					break;
			}

			// clip the min/max voxel bounds 
			allmin = fmax(inmin, outmin);
			allmax = fmin(inmax, outmax);
			for(f = 0; f <= hp->dim; ++f) {
				faces_in[f].thmin = fmax(fmax(faces_in[f].thmin, allmin), thmin);
				faces_in[f].thmax = fmin(fmin(faces_in[f].thmax, allmax), thmax);
				faces_out[f].thmin = fmax(fmax(faces_out[f].thmin, allmin), thmin);
				faces_out[f].thmax = fmin(fmin(faces_out[f].thmax, allmax), thmax);
				if(faces_in[f].thmax <= faces_in[f].thmin) {
					faces_in[f].thmin = 0.0; 
					faces_in[f].thmax = 0.0;
				}
				if(faces_out[f].thmax <= faces_out[f].thmin) {
					faces_out[f].thmin = 0.0; 
					faces_out[f].thmax = 0.0;
				}
			}

			// set the incoming flux from the beam front for this cell
			// then compute pairwise intersections between in, out faces
			// to get the outgoing fluxes
			voxdata[ii][jj].flux_in[hp->dim] = subflux(faces_in[hp->dim].thmin, faces_in[hp->dim].thmax, thmin, thmax, ray.Ftot, ray.Fcom);
			voxdata[ii][jj].fmom_in[hp->dim] = subfmom(faces_in[hp->dim].thmin, faces_in[hp->dim].thmax, thmin, thmax, ray.Ftot, ray.Fcom);
			memset(flux_out, 0, sizeof(flux_out));
			memset(fmom_out, 0, sizeof(fmom_out));
			for(fin = 0; fin <= hp->dim; ++fin) {
				ftin = faces_in[fin];
				for(fout = 0; fout <= hp->dim; ++fout) {
					ftout = faces_out[fout];
					th0 = fmax(ftin.thmin, ftout.thmin);
					th1 = fmin(ftin.thmax, ftout.thmax);
					if(th0 >= th1) continue;

					// track a solid beam between the two faces
					// calculate the optical depth to first order in theta
					// TODO: 1st order
					real k = 00.0;//1;//0.1;
					fmid_in = subflux(th0, th1, faces_in[fin].thmin, faces_in[fin].thmax, voxdata[ii][jj].flux_in[fin], voxdata[ii][jj].fmom_in[fin]);
					fmom_in = subfmom(th0, th1, faces_in[fin].thmin, faces_in[fin].thmax, voxdata[ii][jj].flux_in[fin], voxdata[ii][jj].fmom_in[fin]);
					thmid = 0.5*(th0+th1);
					secthmid = 1.0/cos(thmid);
					cscthmid = 1.0/sin(thmid);

					// get the correct inner and outer radii, then attenuate exponentially 
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

						fmid_out = fmid_in*exp(-k*hp->grid[flatind].rho*(r_out-r_in));
	
						// Write out the radiation energy density 
						// TODO: interpolate the mean value using the exponential
						hp->rad_grid[flatind].E += 0.5*(fmid_in+fmid_out)*(r_out-r_in)/(ray.rmax-ray.rmin)/(hp->dx*hp->dx);
						//printf("detot = %.5e, dr = %.5e, fmid_in = %.5e\n", 0.5*(fmid_in+fmid_out)*(r_out-r_in)/(ray.rmax-ray.rmin)/(hp->dx*hp->dx), r_out-r_in, fmid_in);
	
						// update hydro terms by conserving the difference
						// between downwind and upwind fluxes
						hp->grid[flatind].etot += fmid_in-fmid_out; 
						hp->grid[flatind].mom[0] += lsgn[0]*cos(thmid)*(fmid_in-fmid_out)/CLIGHT; // TODO: review this!
						hp->grid[flatind].mom[1] += lsgn[1]*sin(thmid)*(fmid_in-fmid_out)/CLIGHT; 
					}

					// update the outgoing ray fluxes and momenta
					flux_out[fout] += fmid_out; 
					fmom_out[fout] += fmid_out*thmid; 
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



