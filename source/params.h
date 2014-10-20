/*
 * params.h
 *
 * This is a storage class for cosmological parameters. It is given a variety of inputs, computes some results from that,
 * and then just holds them for future reference.
 *
 */

#ifndef PARAMS_H_
#define PARAMS_H_

#include <cmath>
#include "inireader.h"

class Parameters {
	public:

		// Constructor. Sets the parameters of the model, which cannot be changed afterwards
		Parameters (IniReader&);

		// In this code, we normalise all densities by dividing by \rho_0 = 3 H_0^2 / 8 \pi G
		// where H0 = 100 km/s/Mpc. This means that the dimensionless energy densities are Omega_x * h^2 in the usual parlance
		// Getters for the dimensionless parameters of the model
		inline double rhoM () {return mOmegaMh2;}
		inline double rhoB () {return mOmegaBh2;}
		inline double rhoR () {return mOmegaRh2;}
		inline double rhoGamma () {return mOmegaGammah2;}
		inline double rhoK () {return mOmegaKh2;}  // = -k/H_0^2
		inline double getconH0 () {return mconh0;}  // Returns c / H0 in Mpc (where H0 = 100 km/s/Mpc)

		// Getters for some basic parameters
		inline double Tgamma () {return mT;}       // Temperature of photons today
		inline double Neff () {return mNeff;}      // Effective number of relativistic species (3.046 for standard neutrinos)
		inline double z0 () {return mz0;}          // Redshift to start the evolution at

		// Routine to store (and get) the value of h that has been calculated by an evolution
		void seth(double h) {mh = h;}
		double geth() {return mh;}

		// Routine to obtain the real Omegas, after h has been obtained
		inline double getOmegaK() {return mOmegaKh2 / mh / mh;}
        inline double getOmegaM() {return mOmegaMh2 / mh / mh;}
        inline double getOmegaB() {return mOmegaBh2 / mh / mh;}
        inline double getOmegaR() {return mOmegaRh2 / mh / mh;}
        inline double getOmegaGamma() {return mOmegaGammah2 / mh / mh;}
        // Routine to return DH, after h as been obtained
        inline double getDH() {return mconh0 / mh;}

        // The bitmask used to determine which combination of data sets are used
        inline unsigned int getbitmask() {return bitmask;}


	private:
		// Internal storage for values
		double mOmegaMh2;
		double mOmegaBh2;
		double mOmegaRh2;
		double mOmegaGammah2;
		double mOmegaKh2;
		double mNeff;
		double mT;
		double mz0;
		double mh;
        double mconh0;
        unsigned int bitmask;
};

#endif /* PARAMS_H_ */
