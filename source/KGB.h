/*
 * KGB.h
 *
 * This implements a general KGB class. The user will need to define their own KGB function.
 *
 */

#ifndef KGB_H_
#define KGB_H_

#include "model.h"

class KGB : public Model {

	public:
		// Here are the functions that are overridden by the KGB class
		int derivatives(const double data[], double derivs[], Parameters &params);
		int init(double data[], const double time, Parameters &params, IniReader &init, Output &output);
        bool isvalidconfig(const double data[]) {return true;}

		// The speed of sound and scalar ghost conditions need to be calculated for KGB
		double speedofsound2(const double data[]);
		bool implementsSOS() {return true;}
		bool isghost(const double data[]);
		bool implementsghost() {return true;}
		// The speed of tensor perturbations are unchanged from GR however
		double speedoftensor2(const double data[]) {return 1.0;} 
		bool implementsSOT() {return true;}
		bool istensorghost(const double data[]) {return false;} 
		bool implementstensorghost() {return true;}

	private:
		// Here are some overridden internal functions. They're pretty self-explanatory.
		double pressure(const double data[], const double hdot);
		double energydensity(const double data[]);
		int commonterms(const double data[]);
		// Just some parameters
		double lambda;
		double alpha;
		double beta;
		double n;

		// Function to calculate the Lagrangian and all its appropriate derivatives:
		// The results array should be of length 6
		int computelagrangian(const double data[]);
		// Each time the stuff is calculated, store both the data and it, so as not to waste computation time
		double storeddata[4];
		double a2,a4,a8, hubble;
		double phidot, phidot2, phidot3, phidot4;
		double X;
		double L3;
		double V;
		double Vp;
		double VX;
		double VXX;
		double VXp;
		double L3p;
		double L3X;
		double L3pp;
		double L3Xp;
		double L3XX;
		double P1, P2, D, NV, NL3_1, NL3_2;
		double cs2denom;
};

#endif /* KGB_H_ */
