/*
 * simplecheck.h
 *
 * This provides a basic implementation of a consistency checking class.
 *
 */

#ifndef SIMPLECHECK_H_
#define SIMPLECHECK_H_

#include "consistency.h"

#include <cmath>
// stringstreams
#include <iostream>
#include <string>
#include <sstream>

class SimpleCheck : public Consistency {
	public:
		// Inherited functions
		void checkstate (const double*, const double, IntParams&, Output&, const double*);
		void checkfinal (const double*, const double, IntParams&, Output&, const double*);

		// Constructor
		// Initialize variables
		SimpleCheck() {
			friederror = false;
			phantom = false;
			superluminal = false;
			laplacian = false;
			ghost = false;
			tsuperluminal = false;
			tlaplacian = false;
			tghost = false;
			negenergy = false;
		}

	private:
		// Store booleans on which warning have been issued

		// Error in Friedmann equation sufficiently small
		bool friederror;
		// Equation of state is not phantom
		bool phantom;
		// Speed of sound is not superluminal
		bool superluminal;
		// Speed of sound is not imaginary
		bool laplacian;
		// Perturbations are not ghost-like
		bool ghost;
		// Tensor speed is not superluminal
		bool tsuperluminal;
		// Tensor speed is not imaginary
		bool tlaplacian;
		// Tensors are not ghost-like
		bool tghost;
		// Energy density is not negative
		bool negenergy;

};

#endif /* SIMPLECHECK_H_ */
