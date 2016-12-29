/*************************************************************
 *
 *		eos.c
 *
 *		Devon Powell
 *		21 April 2015
 *
 *		Equations of state
 *
 *************************************************************/

#include "eos.h"

// ideal gas
real eos_ideal_p(real rho, real e, real gamma) {
	return (gamma - 1.0)*e*rho;
}
