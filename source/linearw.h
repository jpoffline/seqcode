/*
 * linearw.h
 *
 * This implements a dark energy model with equation of state w = w0 + wa (1 - a).
 * Because such a situation can exist in a variety of models, we do not implement the
 * speed of sound or ghost checks.
 *
 */

#ifndef LINEARW_H_
#define LINEARW_H_

#include "model.h"

class LinearW : public Model {

	public:
		// Here are the functions that are overridden by the LinearW class
		int derivatives(const double data[], double derivs[], Parameters &params);
		int init(double data[], const double time, Parameters &params, IniReader &init, Output &output);
        bool isvalidconfig(const double data[]) {return true;}

	private:
		// Here are some overridden internal functions. They're pretty self-explanatory.
		double pressure(const double data[], const double hdot);
		double energydensity(const double data[]);

		// Omega_Lambda h^2 today
		double OmegaLambdah2;
		// Equation of state w(a) = w0 + wa (1 − a)
		double w0;
		double wa;

        // Each time the stuff is calculated, store both the data and it, so as not to waste computation time
		// Raising something to a power is expensive! (so is the exponential)
        double storeddata[4];
        double menergydensity;

};

#endif /* LINEARW_H_ */
