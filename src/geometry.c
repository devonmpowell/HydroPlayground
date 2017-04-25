#include "geometry.h"

// geometry-related constants
#define ONE_THIRD 0.3333333333333333333333333333333333
#define ONE_SIXTH 0.1666666666666666666666666666666666
#define ONE_TWLTH 0.0833333333333333333333333333333333

#define PSI_NDIM 3

#define min(x, y) (((x) < (y))? (x) : (y))
#define max(x, y) (((x) > (y))? (x) : (y))




// internal declarations for very low-level voxelization routines
void psi_split_coord(psi_poly* inpoly, psi_poly* outpolys, real coord, int ax, int splitdir);
void psi_reduce(psi_poly* poly, real* moments, int polyorder, int weight);

void psi_init_tet(psi_poly* poly, rvec* verts);

real psi_orient_tet(rvec* pos, rvec* vel);


void psi_split_coord(psi_poly* inpoly, psi_poly* outpolys, real coord, int ax, int splitdir) {

	int v, np, i, npnxt, onv, vcur, vnext, vstart, pnext, nright, fcur;
	real invw;
	rvec newpos, newq;
	int* nverts = &inpoly->nverts;
	psi_vertex* vertbuffer = inpoly->verts; 
	int side[POLYSZ];
	real sdists[POLYSZ];

	// calculate signed distances to the clip plane
	nright = 0;
	memset(&side, 0, sizeof(side));
	for(v = 0; v < *nverts; ++v) {
		sdists[v] = coord - vertbuffer[v].pos.xyz[ax];
		if(sdists[v] < 0.0) {
			side[v] = 1;
			nright++;
		}
	}

	// return if the poly lies entirely on one side of it 
	if(nright == 0) {
		outpolys[0] = *inpoly;
		outpolys[1].nverts = 0;
		return;
	}
	if(nright == *nverts) {
		outpolys[1] = *inpoly;
		outpolys[0].nverts = 0;
		return;
	}

	// check all edges and insert new vertices on the bisected edges 
	onv = inpoly->nverts;
	for(vcur = 0; vcur < onv; ++vcur) {
		if(side[vcur]) continue;
#if PSI_NDIM == 3
		for(np = 0; np < 3; ++np) {
			vnext = vertbuffer[vcur].pnbrs[np];
			if(!side[vnext]) continue;
			invw = 1.0/(sdists[vcur]-sdists[vnext]); // guaranteed to be nonzero
			for(i = 0; i < 3; ++i) {
				newpos.xyz[i] = invw*(vertbuffer[vnext].pos.xyz[i]*sdists[vcur] - vertbuffer[vcur].pos.xyz[i]*sdists[vnext]);
				newq.xyz[i] = invw*(vertbuffer[vnext].q.xyz[i]*sdists[vcur] - vertbuffer[vcur].q.xyz[i]*sdists[vnext]);
			}
			for(npnxt = 0; npnxt < 3; ++npnxt) 
				if(vertbuffer[vnext].pnbrs[npnxt] == vcur) break;
			vertbuffer[*nverts].pos = newpos;
			vertbuffer[*nverts].q = newq;
			vertbuffer[*nverts].pnbrs[0] = vcur;
			vertbuffer[*nverts].fnbrs[0] = vertbuffer[vnext].fnbrs[npnxt];
			vertbuffer[vcur].pnbrs[np] = *nverts;
			(*nverts)++;
			vertbuffer[*nverts].pos = newpos;
			vertbuffer[*nverts].q = newq;
			side[*nverts] = 1;
			vertbuffer[*nverts].pnbrs[0] = vnext;
			vertbuffer[*nverts].fnbrs[0] = vertbuffer[vcur].fnbrs[np];
			vertbuffer[vnext].pnbrs[npnxt] = *nverts;
			(*nverts)++;
		}
#elif PSI_NDIM == 2
		//for(np = 0; np < 2; ++np) {
			//vnext = vertbuffer[vcur].pnbrs[np];
			//if(!side[vnext]) continue;
			//invw = 1.0/(sdists[vcur]-sdists[vnext]); // guaranteed to be nonzero
			//for(i = 0; i < 2; ++i) {
				//newpos.xyz[i] = invw*(vertbuffer[vnext].pos.xyz[i]*sdists[vcur] - vertbuffer[vcur].pos.xyz[i]*sdists[vnext]);
				//newq.xyz[i] = invw*(vertbuffer[vnext].q.xyz[i]*sdists[vcur] - vertbuffer[vcur].q.xyz[i]*sdists[vnext]);
			//}
			//vertbuffer[*nverts].pos = newpos;
			//vertbuffer[*nverts].q = newq;
			//vertbuffer[*nverts].pnbrs[1-np] = vcur;
			//vertbuffer[*nverts].pnbrs[np] = -1;
			//vertbuffer[vcur].pnbrs[np] = *nverts;
			//(*nverts)++;
			//vertbuffer[*nverts].pos = newpos;
			//vertbuffer[*nverts].q = newq;
			//side[*nverts] = 1;
			//vertbuffer[*nverts].pnbrs[np] = vnext;
			//vertbuffer[*nverts].pnbrs[1-np] = -1;
			//vertbuffer[vnext].pnbrs[1-np] = *nverts;
			//(*nverts)++;
		//}
#endif
	}

	// for each new vert, search around the faces for its new neighbors
	// and doubly-link everything
	for(vstart = onv; vstart < *nverts; ++vstart) {
		vcur = vstart;
		vnext = vertbuffer[vcur].pnbrs[0];
		fcur = vertbuffer[vcur].fnbrs[0];
#if PSI_NDIM == 3
		do {
			for(np = 0; np < 3; ++np) if(vertbuffer[vnext].pnbrs[np] == vcur) break;
			vcur = vnext;
			pnext = (np+1)%3;
			vnext = vertbuffer[vcur].pnbrs[pnext];
			vertbuffer[vcur].fnbrs[pnext] = fcur;
		} while(vcur < onv);
		vertbuffer[vstart].pnbrs[2] = vcur;
		vertbuffer[vcur].pnbrs[1] = vstart;

		// label the new face according to the clip plane index
		vertbuffer[vstart].fnbrs[2] = 2*ax+(1-side[vstart]);

#elif PSI_NDIM == 2
		if(vertbuffer[vstart].pnbrs[1] >= 0) continue;
		do {
			vcur = vnext;
			vnext = vertbuffer[vcur].pnbrs[0];
		} while(vcur < onv);
		vertbuffer[vstart].pnbrs[1] = vcur;
		vertbuffer[vcur].pnbrs[0] = vstart;
#endif
	}

	if(*nverts >= POLYSZ) printf("WARNING: Overflowed vertex buffer. nverts = %d >= %d\n", *nverts, POLYSZ);

	// copy and compress vertices into their new buffers
	// TODO: do this beforehand, since we have two buffers ready
	// reusing side[] for reindexing
	onv = *nverts;
	outpolys[0].nverts = 0;
	outpolys[1].nverts = 0;
	for(v = 0; v < onv; ++v) {
		outpolys[side[v]].verts[outpolys[side[v]].nverts] = vertbuffer[v];
		side[v] = outpolys[side[v]].nverts++;
	}
	for(v = 0; v < outpolys[0].nverts; ++v) 
		for(np = 0; np < PSI_NDIM; ++np)
			outpolys[0].verts[v].pnbrs[np] = side[outpolys[0].verts[v].pnbrs[np]];
	for(v = 0; v < outpolys[1].nverts; ++v) 
		for(np = 0; np < PSI_NDIM; ++np)
			outpolys[1].verts[v].pnbrs[np] = side[outpolys[1].verts[v].pnbrs[np]];
}


void psi_extract_faces(psi_poly* poly, psi_face_buffer* faces, int* nfaces) {

	// reduce according to Powell and Abel (2015) for simplicity
	real mass;
	int np, vstart, pstart, vcur, vnext, pnext, fcur;
	rvec q0, q1, q2; 
	psi_vertex* vertbuffer = poly->verts; 
	int* nverts = &poly->nverts; 
	int emarks[*nverts][PSI_NDIM];
	memset(&emarks, 0, sizeof(emarks));

	*nfaces = 0;
	psi_face_buffer curface;

	for(vstart = 0; vstart < *nverts; ++vstart)
	for(pstart = 0; pstart < 3; ++pstart) {
		if(emarks[vstart][pstart]) continue;
		pnext = pstart; 
		vcur = vstart; 
		emarks[vcur][pnext] = 1;
		fcur = vertbuffer[vcur].fnbrs[pnext]; 
		vnext = vertbuffer[vcur].pnbrs[pnext];

		// skip if the face is a beam boundary (zero flux)
		if(fcur < 0) continue;
		curface.face_id = fcur;
		curface.nverts = 0;

		//printf( "New face, verts(faces) = ");

		do {
			curface.verts[curface.nverts++] = vertbuffer[vcur].pos;
			//printf( "%d(%d) ", vcur, fcur);
			// go to the next face
			for(np = 0; np < 3; ++np) if(vertbuffer[vnext].pnbrs[np] == vcur) break;
			vcur = vnext;
			pnext = (np+1)%3;
			emarks[vcur][pnext] = 1;
			vnext = vertbuffer[vcur].pnbrs[pnext];
		} while(vcur != vstart);
		faces[(*nfaces)++] = curface;
	}	

}

	
void psi_reduce(psi_poly* poly, real* moments, int polyorder, int weight) {
	// reduce according to Powell and Abel (2015) for simplicity
#if PSI_NDIM == 3
	const static int nmom[3] = {1, 4, 10};
#elif PSI_NDIM == 2
	const static int nmom[3] = {1, 3, 6};
#endif
	real mass;
	int np, vstart, pstart, vcur, vnext, pnext, fcur;
	rvec q0, q1, q2; 
	psi_vertex* vertbuffer = poly->verts; 
	int* nverts = &poly->nverts; 
	int emarks[*nverts][PSI_NDIM];
	memset(moments, 0, nmom[polyorder]*sizeof(real)); 
	memset(&emarks, 0, sizeof(emarks));
#if PSI_NDIM == 3
	for(vstart = 0; vstart < *nverts; ++vstart)
	for(pstart = 0; pstart < 3; ++pstart) {
		if(emarks[vstart][pstart]) continue;
		pnext = pstart; 
		vcur = vstart; 
		emarks[vcur][pnext] = 1;
		fcur = vertbuffer[vcur].fnbrs[pnext]; 
		vnext = vertbuffer[vcur].pnbrs[pnext];
		//q0 = vertbuffer[vcur].q;
		q0 = vertbuffer[vcur].pos;
		for(np = 0; np < 3; ++np) if(vertbuffer[vnext].pnbrs[np] == vcur) break;
		vcur = vnext;
		pnext = (np+1)%3;
		emarks[vcur][pnext] = 1;
		fcur = vertbuffer[vcur].fnbrs[pnext]; 
		vnext = vertbuffer[vcur].pnbrs[pnext];
		while(vnext != vstart) {
			//q1 = vertbuffer[vnext].q;
			//q2 = vertbuffer[vcur].q;
			q1 = vertbuffer[vnext].pos;
			q2 = vertbuffer[vcur].pos;

			mass = (-q2.x*q1.y*q0.z + q1.x*q2.y*q0.z + q2.x*q0.y*q1.z
			   	- q0.x*q2.y*q1.z - q1.x*q0.y*q2.z + q0.x*q1.y*q2.z); 
			moments[0] += mass;

			for(np = 0; np < 3; ++np) if(vertbuffer[vnext].pnbrs[np] == vcur) break;
			vcur = vnext;
			pnext = (np+1)%3;
			emarks[vcur][pnext] = 1;
			fcur = vertbuffer[vcur].fnbrs[pnext]; 
			vnext = vertbuffer[vcur].pnbrs[pnext];
		}
	}	
#elif PSI_NDIM == 2
	//for(vcur = 0; vcur < *nverts; ++vcur) {
		//vnext = vertbuffer[vcur].pnbrs[0];
		//q0 = vertbuffer[vcur].q;
		//q1 = vertbuffer[vnext].q;
		//mass = (q0.x*q1.y - q0.y*q1.x); 
		//moments[0] += mass;
		//if(polyorder > 0) {
			//moments[1] += ONE_THIRD*mass*(q0.x + q1.x);
			//moments[2] += ONE_THIRD*mass*(q0.y + q1.y);
		//}
		//if(polyorder > 1) {
			//moments[3] += ONE_SIXTH*mass*(q0.x*q0.x + q1.x*q1.x + q0.x*q1.x);
			//moments[4] += ONE_TWLTH*mass*(q0.x*(2*q0.y + q1.y) + q1.x*(q0.y + 2*q1.y));
			//moments[5] += ONE_SIXTH*mass*(q0.y*q0.y + q1.y*q1.y + q0.y*q1.y);
		//}
	//}	
#endif

}

void psi_init_tet(psi_poly* poly, rvec* pos) {
	// initialize the tet as an edge-vertex graph 
	// fill in mass coordinates too
	int v;
	memset(poly, 0, sizeof(psi_poly));
#if PSI_NDIM == 3
	poly->nverts = 4;
	poly->verts[0].pnbrs[0] = 1;	
	poly->verts[0].pnbrs[1] = 3;	
	poly->verts[0].pnbrs[2] = 2;	
	poly->verts[1].pnbrs[0] = 2;	
	poly->verts[1].pnbrs[1] = 3;	
	poly->verts[1].pnbrs[2] = 0;	
	poly->verts[2].pnbrs[0] = 0;	
	poly->verts[2].pnbrs[1] = 3;	
	poly->verts[2].pnbrs[2] = 1;	
	poly->verts[3].pnbrs[0] = 1;	
	poly->verts[3].pnbrs[1] = 2;	
	poly->verts[3].pnbrs[2] = 0;	
	// TODO: all four BC coordinates for symmetry?
	poly->verts[1].q.x = 1.0;
	poly->verts[2].q.y = 1.0;
	poly->verts[3].q.z = 1.0;
#elif PSI_NDIM == 2
	poly->nverts = 3;
	poly->verts[0].pnbrs[0] = 1;	
	poly->verts[0].pnbrs[1] = 2;	
	poly->verts[1].pnbrs[0] = 2;	
	poly->verts[1].pnbrs[1] = 0;	
	poly->verts[2].pnbrs[0] = 0;	
	poly->verts[2].pnbrs[1] = 1;	
	poly->verts[1].q.x = 1.0;
	poly->verts[2].q.y = 1.0;
#endif
	for(v = 0; v < PSI_NDIM+1; ++v) poly->verts[v].pos = pos[v];
}

void psi_voxels_init(psi_voxels* vox, psi_poly* poly, dvec* ldir, hydro_problem* grid) {

	int i, cliplo, cliphi, nplanes;
	real cmin, cmax;
	psi_plane planes[2*PSI_NDIM];

#if 0
	// first see if the poly needs to be clipped against the window boundary 
	// and compute the index bounds for the polyhedron
	memset(planes, 0, sizeof(planes));
	nplanes = 0;
	for(i = 0; i < PSI_NDIM; ++i) {
		cliplo = 0; cliphi = 0;
		cmin = rbox[0].xyz[i];
		cmax = rbox[1].xyz[i];
		if(grid->window[0][i] > cmin) {
			cmin = grid->window[0][i];
			cliplo = 1;
		}	
		if(grid->periodic && grid->box[0][i] > cmin) {
			cmin = grid->box[0][i];
			cliplo = 1;
		}
		if(grid->window[1][i] < cmax) {
			cmax = grid->window[1][i];
			cliphi = 1;
		}	
		if(grid->periodic && grid->box[1][i] < cmax) {
			cmax = grid->box[1][i];
			cliphi = 1;
		}
		if(cliplo) {
			planes[nplanes].n.xyz[i] = 1.0;
			planes[nplanes].d = -cmin;
			nplanes++;
		}
		if(cliphi) {
			planes[nplanes].n.xyz[i] = -1.0;
			planes[nplanes].d = cmax;
			nplanes++;
		}
		poly->ibox[0].ijk[i] = floor((cmin-grid->window[0][i])/grid->d[i]);
		poly->ibox[1].ijk[i] = ceil((cmax-grid->window[0][i])/grid->d[i]);
		if(poly->ibox[0].ijk[i] < 0) poly->ibox[0].ijk[i] = 0;
		if(poly->ibox[1].ijk[i] > grid->n[i]) poly->ibox[1].ijk[i] = grid->n[i];
	}
	psi_clip(poly, planes, nplanes);
#endif

	// initialize the stack
	vox->grid = grid;
	vox->splitdir = *ldir;
	vox->stack[0] = *poly;
	vox->nstack = 1;
}

int psi_voxels_next(psi_voxels* vox, psi_poly* poly) {

	int i, dmax, spax, siz;
	psi_poly curpoly;

	// pointers for easier acces
	psi_poly* stack = vox->stack;
	int* nstack = &vox->nstack;
	hydro_problem* grid = vox->grid;

	while(*nstack > 0) {
		curpoly = stack[--*nstack];
		if(curpoly.nverts <= 0) continue;
		dmax = 0; spax = 0;
		for(i = 0; i < PSI_NDIM; ++i) {
			siz = curpoly.ibox[1].ijk[i]-curpoly.ibox[0].ijk[i];
			if(siz > dmax) {
				dmax = siz; 
				spax = i;
			}	
		}
		if(dmax == 1) {
			*poly = curpoly;
			return 1;
		}

		// split the poly and push children to the stack
		psi_split_coord(&curpoly, &stack[*nstack], grid->dx*(curpoly.ibox[0].ijk[spax]+dmax/2), spax, vox->splitdir.ijk[spax]);
		memcpy(stack[*nstack].ibox, curpoly.ibox, 2*sizeof(dvec));
		stack[*nstack].ibox[1].ijk[spax] -= dmax-dmax/2; 
		memcpy(stack[*nstack+1].ibox, curpoly.ibox, 2*sizeof(dvec));
		stack[*nstack+1].ibox[0].ijk[spax] += dmax/2;
		*nstack += 2;

		// TODO: This is a hack!
		if(vox->splitdir.ijk[spax]) {
			psi_poly tmp = stack[*nstack-1];
			stack[*nstack-1] = stack[*nstack-2];
			stack[*nstack-2] = tmp; 
		}

	}
	return 0;
}

