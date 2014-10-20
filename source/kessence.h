/*
 * kessence.h
 *
 * This implements a general k-essence class. The user will need to define their own k-essence function.
 * The default implementation is a (\nabla \phi)^4 term.
 *
 */

#ifndef KESSENCE_H_
#define KESSENCE_H_

#include "model.h"

class Kessence : public Model {

	public:
		// Here are the functions that are overridden by the quintessence class
		int derivatives(const double data[], double derivs[], Parameters &params);
		int init(double data[], const double time, Parameters &params, IniReader &init, Output &output);
        bool isvalidconfig(const double data[]) {return true;}

		// The speed of sound and scalar ghost conditions need to be calculated for k-essence
		double speedofsound2(const double data[]);
		bool implementsSOS() {return true;}
		bool isghost(const double data[]);
		bool implementsghost() {return true;}

		// The speed of tensor perturbations are unchanged from GR however
		double speedoftensor2(const double data[]) {return 1.0;} // The speed of tensor perturbations in k-essence is always 1.
		bool implementsSOT() {return true;}
		bool istensorghost(const double data[]) {return false;} // k-essence is never ghost-like
		bool implementstensorghost() {return true;}

	private:
		// Here are some overridden internal functions. They're pretty self-explanatory.
		double pressure(const double data[], const double hdot);
		double energydensity(const double data[]);

		// Just some parameters
		double lambda;
		double alpha;
		double beta;
		double n;

		// Function to calculate the Lagrangian and all its appropriate derivatives:
		// U, Up, Upp, UX, UXX, UXP
		void computelagrangian(const double data[]);
		// Each time the stuff is calculated, store both the data and it, so as not to waste computation time
		double storeddata[4];
		double X;
		double L;
		double LX;
		double LXX;
		double Lp;
		double LXp;

};

#endif /* KESSENCE_H_ */
