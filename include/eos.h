/*************************************************************
 *
 *		eos.h
 *
 *		Devon Powell
 *		21 April 2015
 *
 *		Equations of state
 *
 *************************************************************/

#ifndef EOS_H_
#define EOS_H_

#include "common.h"

#define EOS_SEVEN_FIFTHS 1.4
#define EOS_FIVE_THIRDS 1.6666666666666666666667

// ideal gas
real eos_ideal_p(real rho, real e, real gamma);


#endif // EOS_H_
