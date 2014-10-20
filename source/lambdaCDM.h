/*
 * lambdaCDM.h
 *
 * This is a very basic model - it ignores the scalar field altogether, and implements a cosmological constant with
 * Omega_DE = 0.7 (or as specified).
 *
 */

#ifndef LAMBDACDM_H_
#define LAMBDACDM_H_

#include "model.h"

class LambdaCDM : public Model {

	public:
		// Here are the base functions that are overridden by this class
		int derivatives(const double data[], double derivs[], Parameters &params);
		int init(double data[], const double time, Parameters &params, IniReader &init, Output &output);
        bool isvalidconfig(const double data[]) {return true;}

	private:
		// Here are some overridden internal functions. They're pretty self-explanatory.
		double energydensity(const double data[]);
		double pressure(const double data[], const double hdot);
		// Variable to store Omega_Lambda
		double OmegaLambdah2;
};

#endif /* LAMBDACDM_H_ */
