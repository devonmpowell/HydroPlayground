/*************************************************************
 *
 *		riemann.c
 *
 *		Devon Powell
 *		21 April 2015
 *
 *		Approximate Riemann solvers for d_euler_hydro
 *
 *************************************************************/

/**
 *
 * riemann_hll
 *
 * rho, E, p are 2-vectors containing the left- and right- states.
 * ax is the axis normal to this face
 * dim is the dimensionality of the problem to be solved
 * v has length 2*dim, e.g. { vx_left, vy_left, vx_right, vy_right }
 *
 * Result goes in the vector F, of length 5
 *
 */
void riemann_hll(real* rho, real* v, real* E, real* p, int ax, int dim, real* F) {

	real F_l[5];
	real F_r[5];

	for


}
