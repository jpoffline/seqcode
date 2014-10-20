#include "KGB.h"

// The derivatives routine is given the state of the system (a, phi and \dot{\phi}) as well as the parameters of the system,
// and returns the derivatives (\dot{a}, \dot{\phi}, \ddot{\phi}, \dot{H})
int KGB::derivatives(const double data[], double derivs[], Parameters &params) {

	// The first section here is model independent

	// Extract data for easier reading of the code
	double a = data[0];
	double phi = data[1];
	double phidot = data[2];
	double phidot2 = phidot * phidot;
	double phidot3 = phidot * phidot2;
	double phidot4 = phidot2 * phidot2;
	double hubble = data[3];

	// Computing \dot{a}, knowing the hubble rate
	derivs[0] = hubble * a;
	// Computing \dot{\phi}
	derivs[1] = phidot;

	// Model dependent code below

	// Computing \ddot{\phi}. This requires solving the scalar equation of motion for \ddot{\phi}.
	// (1) Compute KGB Lagrangian quantities
	int result = computelagrangian(data);
	
	// Compute \ddot{\phi} = N/D,
	// and split the numerator N and denominator D
	// into contributions from V & L3 separately.
	// N = NV + NL3, D = DV + DL3.
	// These are computed in the "computelagrangian" call.
	

	// Computing \dot{H}. This requires solving the acceleration equation for \dot{H}.
	// Note that pressure does depend on \dot{H} in this model,
	
	// Write the pressure as P = P1 + P2\ddot(phi)
	// P1 & P2 are computed in the "computelagrangian" call

	// Split up the acceleration equation: H1 is all the "standard" contributions
	double H1 = - 0.5 * hubble * hubble - 0.5 * params.rhoR() / a2 + 0.5 * params.rhoK();
	// put it all together and compute \dot{H}
	derivs[3]=(H1-1.5*a2*(P1+P2*(NV+NL3_1)/D)) / (1.0+1.5*a2*P2*NL3_2/D);

	// Now that we've computed hubbledot, can finish off the scalar field equation of motion:
	derivs[2] = ( NV + NL3_1 + NL3_2 * derivs[3] ) / D;
	
	// To compute the pressure from the function, need to pass it \dot{H}
	double press = pressure(data, derivs[3]);
	
	// GSL_SUCCESS indicates that the computation was successful.
	// If it failed, look up the appropriate error code in the same enum as GSL_SUCCESS
	return GSL_SUCCESS;
}

// Function to calculate the Lagrangian and all its appropriate derivatives
// The results array should be of length 6
int KGB::computelagrangian(const double data[]) {
	// First, check the data against the previous data
	// This prevents the results being computed multiple times on the same data
	if (data[0] == storeddata[0] &&
			data[1] == storeddata[1] &&
			data[2] == storeddata[2] &&
			data[3] == storeddata[3]) {
		// It's the same as before, so don't recompute it
		// Success!
		return 0;
	}

	// Extract data for easier reading of the code
	double phi = data[1];
	// Compute the X quantity
	double temp = data[2] / data[0];
	X = temp * temp / 2.0;
	
	// Compute powers of a:
	a2 = data[0] * data[0];
	a4 = a2 * a2;
	a8 = a4 * a4;

    // Here we construct the Lagrangian & important derivatives thereof.
	// L = V(phi, X) - L3(phi, X) * boxphi;
	
    // Construct the two functions V(phi, X) & L3(phi, X).
    V = 1.0;
	L3 = 1.0;
	
	// Now construct their derivatives.
	
	VX = 1.0;
	VXX = 1.0;

	L3p = 1.0;
	L3pp = 1.0;
	L3X = 1.0;
	L3XX = 1.0;
	L3Xp = 1.0;
	
	//  Useful parameterizations:
	// V = c_2X^n,
	// L3 = c_3 X^m
	// .. "generalized galileon"
	// Could also easily add in potential terms.
	
	
	// Compute some commonly occuring terms in the equation of motion
	// and acceleration equation.
	
	int res = commonterms(data);
	
	
	// Store the data for which these results are correct
	for (int i = 0; i < 4; i++)
		storeddata[i] = data[i];

	// Success!
	return 0;
}

int KGB::commonterms(const double data[]){
	hubble = data[3];
	phidot = data[2];
	phidot2 = phidot * phidot;
	phidot3 = phidot * phidot2;
	phidot4 = phidot2 * phidot2;
	
	// SCALAR FIELD EQUATION OF MOTION
	
	// \ddot\phi = N / D.
	// NUMERATORS, N = NV + NL3_1 + NL3_2*hubbledot.
	// Numerator for the "V"-part
	NV = a8*(a2*Vp+hubble*(VXX*phidot3/a2-2.0*VX*phidot)-VXp*phidot2); 
	NL3_1 = a4*(4.0*a4*hubble*phidot*L3p
		+3.0*hubble*hubble*phidot4*L3XX*a4*phidot2*L3pp-4.0*a2*hubble*phidot3*L3X);
	NL3_2 = -a4*3.0*a2*phidot2*L3X;
	
	// DENOMINATOR
	
	// Denominator for "V" part
	double DV = a8*(VXX*phidot2/a2 + VX);
	// denominator for "L3" part
	double DL3 = a4*(-2.0*a4*L3p+6.0*a2*hubble*phidot*L3X-a2*phidot2*L3X+3.0*L3XX*hubble*phidot3);
	// Full denominator
	D = DV + DL3;

	// Write the pressure as P = P1 + P2\ddot(phi)
	P1 = (V-2.0*X/a2*(a2*L3p-L3X*hubble*phidot))/3.0;
	P2 = -2.0*X/a2*L3X/3.0;
	
	// denominator in the sound speed
	// Also computed to check stability.
	cs2denom=a2*(-2.0*L3p+VX+2.0*X*(VXX-L3Xp)+6.0*X*X*L3X*L3X);
	cs2denom+=6.0*hubble*(L3X+X*L3XX)*phidot;
	
	return 0;
	
} // END commonterms()

// Returns the ratio rho_Q/rho_c
double KGB::energydensity(const double data[]){
	// Extract data for easier reading of the code
	double a = data[0];
	double phi = data[1];
	double phidot = data[2];
	double hubble = data[3];
	double a2 = a * a;
	
	// Compute quantities
	int result = computelagrangian(data);

	// Now construct energy/3 for KGB
	return (2.0*X*VX-V+2.0*X/a2*(-a2*L3p+3.0*L3X*hubble*phidot))/3.0;
}

// Returns the ratio P_Q/rho_c
double KGB::pressure(const double data[], const double hdot){
	// Extract data for easier reading of the code
	double a = data[0];
	double phi = data[1];
	double phidot = data[2];
	double phidot2 = phidot * phidot;
	double phidot3 = phidot * phidot2;
	double phidot4 = phidot2 * phidot2;
	double hubble = data[3];
	
	// Compute quantities
	int result = computelagrangian(data);

	// P = P1 + P2 * phidotdot.
	// phidotdot = N / D,
	// with N = NV + NL3_1 + NL3_2*hdot,
	// and D = DV + DL3.
	// These were calculated inside "commonterms", 
	// called from "computelagrangian".
	
	
	// Now put the above bit in for \ddot{\phi}.
	return P1 + P2 * (NV + NL3_1 + NL3_2 * hdot) / D;
}

/* This function does a few things:
 * - Sets the name of the class
 * - Initializes the value of H using the Friedmann equation
 * - Returns a log output
 */
int KGB::init(double data[], double time, Parameters &params, IniReader &init, Output &output) {

	// Set the name of the class
	section = "KGB";

	// Go and get model parameters
	// lambda is just a parameter in the ini file. It can be used to do a single parameter scan, for example.
	lambda = init.getiniDouble("lambda", 1.0, section);
	alpha = init.getiniDouble("alpha", 1.0, section);
	beta = init.getiniDouble("beta", 1.0, section);
	n = init.getiniDouble("n", 2.0, section);

	// Construct H
	// Temporary variable
	double temp;
	// Scale factor
	double a = data[0];
	double a2 = a * a;

	// Calculate H^2
	
	/// NEED to modify this for computing the initial Hubble for KGB
	temp = params.rhoM() / a + params.rhoR() / a2 + params.rhoK() + a2 * energydensity(data);

	// Calculate H
	data[3] = sqrt(temp);

    // Print stuff to the log
    output.printlog("Running KGB model.");
    output.printvalue("phi0", data[1]);
    output.printvalue("phidot0", data[2]);

    // Success!
    return 0;

}

// The speedofsound2 returns the speed of sound squared, given the state of the system
double KGB::speedofsound2(const double data[]) {
	// The speed of sound in KGB can vary from 1.
	// Extract quantities for easier reading of the code
	
	double phidot = data[2];
	double hubble = data[3];
	
	
	// Compute quantities
	int result = computelagrangian(data);

	// Compute numerator of cs2
	double cs2numer = a2*(-2.0*L3p+VX+2.0*L3Xp*X-2.0*X*X*L3X*L3X);
	cs2numer+= 2.0*hubble*(L3X-X*L3XX)*phidot;
	
	// NOTE: this cs2 needs \ddot{\phi}
	//numer+=2.0*(L3X+X*L3XX)*phidotdot;
	
	// Compute denominator of cs2:
	// done in "commonterms" in "computelagrangian" call.
	
	// Return result
	return cs2numer / cs2denom;
}

// The isghost function is given the state of the system and returns whether or not the theory has become ghostlike
bool KGB::isghost(const double data[]) {
	
	// KGB becomes ghost-like when the denominator in the speed of sound squared becomes negative.
	double phidot = data[2];
	double hubble = data[3];
	
	// Compute quantities
	int result = computelagrangian(data);

	// Check if positive, and return the result
	if (cs2denom < 0)
		return true;
	else
		return false;

}
