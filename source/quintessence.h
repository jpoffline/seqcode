/*
 * quintessence2.h
 *
 * This class is an implementation of the Model class; it implements a very basic quintessence model,
 * mostly to make sure that everything is working in the program.
 *
 * The difference between this quintessence module and quintessence.h is that this module integrates the acceleration equation
 * instead of just using the Friedmann equation to obtain the evolution of a.
 *
 */

#ifndef QUINTESSENCE_H_
#define QUINTESSENCE_H_

#include "model.h"

class Quintessence : public Model {

	public:
		// Here are the functions that are overridden by the quintessence class
		int derivatives(const double data[], double derivs[], Parameters &params);
		int init(double data[], const double time, Parameters &params, IniReader &init, Output &output);
        bool isvalidconfig(const double data[]) {return true;}

		// The speed of sound/tensors and ghost checks are all trivial for quintessence
		double speedofsound2(const double data[]) {return 1.0;} // The speed of sound in quintessence is always 1.
		bool implementsSOS() {return true;}
		bool isghost(const double data[]) {return false;} // Quintessence is never ghost-like
		bool implementsghost() {return true;}

		double speedoftensor2(const double data[]) {return 1.0;} // The speed of tensor perturbations in quintessence is always 1.
		bool implementsSOT() {return true;}
		bool istensorghost(const double data[]) {return false;} // Quintessence is never ghost-like
		bool implementstensorghost() {return true;}


	private:
		// Here are some overridden internal functions. They're pretty self-explanatory.
		double pressure(const double data[], const double hdot);
		double energydensity(const double data[]);
		// These internal functions are specific to quintessence
		double potential(const double phi);
		double potentialprime(const double phi);

		// Some parameters that are read in from params.ini
		int modeltype;
		double mass;
		double lambda;
		double alpha;
		double beta;
};

#endif /* QUINTESSENCE_H_ */
