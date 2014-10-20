#include "integrate.h"

Integrator::Integrator() {
	// Constructor

	// Set the number of elements to be integrated
	numelements = 4;

	// Integration method as supplied by GSL
	const gsl_odeiv2_step_type *Type = gsl_odeiv2_step_rk8pd;
    //const gsl_odeiv2_step_type *Type = gsl_odeiv2_step_rkf45;

	// Initialize GSL
	step = gsl_odeiv2_step_alloc(Type, numelements);
	control = gsl_odeiv2_control_yp_new(0, 1.0e-8); // absolute error, relative error
	evolve = gsl_odeiv2_evolve_alloc(numelements);

    // Set the initial stepsize to be quite small, in order to get good data on derivatives at the beginning
	stepsize = minstepsize();

}

Integrator::~Integrator() {
	// Destructor

	// Release the GSL objects
	gsl_odeiv2_evolve_free(evolve);
	gsl_odeiv2_control_free(control);
	gsl_odeiv2_step_free(step);

}

int Integrator::dointstep(int (*func)(double, const double*, double*, void*), IntParams &params, double data[], double &time, const double endtime) {
	// Integrates a single step

	// Sets up the system to integrate, including the function that specifies the derivatives, NULL for the Jacobian, the number of elements
	// being integrated, and the parameters to pass through to the function
	gsl_odeiv2_system sys = { func, NULL, numelements, &params };

	// Make a copy of the stepsize - the evolution function wants to update this quantity, but a reference to the class' private information
	// makes things unhappy
	double stepsizecopy = stepsize;

	// Actually take the integration step
	int status = gsl_odeiv2_evolve_apply(evolve, control, step, &sys, &time, endtime, &stepsizecopy, data);

	// Update the stepsize according to what GSL thinks will work well, but make sure it doesn't get too big. We want some resolution on our functions!
	setstepsize(stepsizecopy);

	// Return the status: should be GSL_SUCCESS if everything worked
	return status;

}
