/*
 * integrate.h
 *
 * This class provides a wrapper around the GSL ode integration routines.
 *
 * The crux of the class is the dointstep function, which takes in a function to call to obtain derivatives,
 * a set of parameters to pass to that function, the present data, the present time, and the max time to run to.
 *
 */

#ifndef INTEGRATE_H_
#define INTEGRATE_H_

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include "intparams.h"

// This defines the integrator class
class Integrator {

	public:
		Integrator(); // Constructor
		~Integrator(); // Destructor

		// Routine that actually does the integration
		int dointstep(int (*func)(double, const double*, double*, void*), IntParams&, double*, double&, const double);

		// Just a setter and getter for stepsize
		inline double getstepsize() const { return stepsize; }

		void setstepsize(const double newsize) {
			// Do some simple error checking
			if (newsize <= 0)
				return;
			stepsize = newsize;
			if (stepsize < minstepsize())
				stepsize = minstepsize();
			if (stepsize > maxstepsize())
				stepsize = maxstepsize();
		}

		// For information
		inline double maxstepsize() { return 0.005; }
		inline double minstepsize() { return 1e-7; }

	private:
		// These are initialized in the constructor, and released in the destructor
		// Step, control and evolution objects
		gsl_odeiv2_step *step;
		gsl_odeiv2_control *control;
		gsl_odeiv2_evolve *evolve;

		// This is the stepsize that has been recommended by the integrator
		double stepsize;

		// This is the number of elements that are being integrated
		size_t numelements;
};

#endif /* INTEGRATE_H_ */
