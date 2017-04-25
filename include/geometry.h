#include "common.h"

#ifndef _GEOMETRY_H_
#define _GEOMETRY_H_


#define BEAM_IN 6
#define BEAM_OUT 7

// structure to hold a polyhedron,
// including mass coordinates
#define PSI_NDIM 3
#define POLYSZ 64 // TODO: fix this number better. Should be 16, accd to Euler Characteristic with 10 faces
typedef struct {
	rvec pos, q;
	int pnbrs[PSI_NDIM];
	int fnbrs[PSI_NDIM];
	int flags;
} psi_vertex;
typedef struct {
	int nverts;
	dvec ibox[2];
	psi_vertex verts[POLYSZ];
} psi_poly; 
typedef struct {
	rvec n;
	real d;
} psi_plane;

// struct to voxelize a polyhedron
// contains a stack and some grid information
typedef struct {
	psi_poly stack[32]; // TODO: how big to make this??
	int nstack;
	dvec splitdir;
	hydro_problem* grid;
} psi_voxels;


typedef struct {
	rvec verts[16];
	int nverts;
	int face_id;
} psi_face_buffer;


void psi_voxels_init(psi_voxels* vox, psi_poly* poly, dvec* splitdir, hydro_problem* grid);
int psi_voxels_next(psi_voxels* vox, psi_poly* poly);

//void psi_voxelize_tet(rvec* pos, rvec* vel, real mass, rvec* rbox, psi_dest_grid* dest_grid);

//void psi_point_sample_tet(rvec* pos, rvec* vel, real mass, rvec* rbox, psi_dest_grid* dest_grid);

//real psi_barycentric(rvec* pos, rvec samppos, real* bcoords);

//void psi_voxelize_annihilation(rvec* pos0, rvec* vel0, real mass0, rvec* rbox0, rvec* pos1, rvec* vel1, real mass1, rvec* rbox1, psi_dest_grid* dest_grid);



#endif // _GEOMETRY_H_
