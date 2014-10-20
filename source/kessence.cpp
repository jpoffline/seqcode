#include "kessence.h"

// The derivatives routine is given the state of the system (a, phi and \dot{\phi}) as well as the parameters of the system,
// and returns the derivatives (\dot{a}, \dot{\phi}, \ddot{\phi}, \dot{H})
int Kessence::derivatives(const double data[], double derivs[], Parameters &params) {

	// The first section here is model independent

	// Extract data for easier reading of the code
	double a = data[0];
	double a2 = a * a;
	double phidot = data[2];
	double phidot2 = phidot * phidot;
	double phidot3 = phidot * phidot2;
	double hubble = data[3];

	// Computing \dot{a}, knowing the hubble rate
	derivs[0] = hubble * a;
	// Computing \dot{\phi}
	derivs[1] = phidot;

	// Model dependent code below

	// Computing \ddot{\phi}. This requires solving the scalar equation of motion for \ddot{\phi}.
	// Compute k-essence quantities
	computelagrangian(data);

	derivs[2] = (a2 * Lp
	                   - 2 * hubble * phidot * LX
	                   - LXp * phidot2
	                   + LXX * hubble * phidot3 / a2)
	                   / (LX + LXX * phidot2 / a2);
	// This is equivalent to
	// ( 1.0 - 3.0 * SoS ) * phidot * hubble + SoS / LX * ( a2 * Lp - pow( phidot , 2.0 ) * LpX );
	// where SoS = c^2. However, the full expression doesn't have a div/0 when LX = 0.

	// Computing \dot{H}. This requires solving the acceleration equation for \dot{H}.
	// Note that pressure does not depend on \dot{H} in this model,
	// so we pass in 0 for \dot{H} when calculating pressure
	double press = pressure(data, 0.0);
    derivs[3] = 0.5 * (- params.rhoR() / a2 - 3.0 * a2 * press - hubble * hubble + params.rhoK());

	// GSL_SUCCESS indicates that the computation was successful.
	// If it failed, look up the appropriate error code in the same enum as GSL_SUCCESS
	return GSL_SUCCESS;
}

// Function to calculate the Lagrangian and all its appropriate derivatives:
// U, Up, Upp, UX, UXX, UXP
void Kessence::computelagrangian(const double data[]) {
	// First, check the data against the previous data
	// This prevents the results being computed multiple times on the same data
	if (data[0] == storeddata[0] &&
			data[1] == storeddata[1] &&
			data[2] == storeddata[2] &&
			data[3] == storeddata[3]) {
		// It's the same as before, so don't recompute it
		return;
	}

	// Extract data for easier reading of the code
	double phi = data[1];
	// Compute the X quantity
	double temp = data[2] / data[0];
	X = temp * temp / 2.0;

	// Calculate everything
	// Below is a different choice of Lagrangian
    // This Lagrangian is L = X + alpha X^n - beta e^{-lambda phi}
    // (all suitably dimensionless)
	// Quintessence has alpha = 0

    // Here is the potential
    double pot = beta * exp(- lambda * phi);
    double dpot = - lambda * pot;

    // Here we construct the Lagrangian & important derivatives thereof.

    // (1) Lagrangian
    L = X + alpha * pow( X, n ) - pot;

    // (2) dL/dX
    LX = 1.0 + alpha * n * pow ( X, n - 1.0 );

    // (3) d^2L/dX^2
    if (n == 1.0)
    	LXX = 0;
    else
        LXX = alpha * n * (n - 1.0) * pow( X, n - 2.0);

    // (4) dL/dphi
    Lp = - dpot;

    // (5) d^2L/dXdphi
    LXp = 0.0;

	// Store the data for which these results are correct
	for (int i = 0; i < 4; i++)
		storeddata[i] = data[i];

}

// Returns the ratio rho_Q/rho_c
double Kessence::energydensity(const double data[]){
	// Compute quantities
	computelagrangian(data);

	return (2.0 * X * LX - L) / 3.0;
}
// Returns the ratio P_Q/rho_c
double Kessence::pressure(const double data[], const double hdot){
	// Compute quantities
	computelagrangian(data);

	return L / 3.0;
}

/* This function does a few things:
 * - Sets the name of the class
 * - Initializes the value of H using the Friedmann equation
 * - Returns a log output
 */
int Kessence::init(double data[], double time, Parameters &params, IniReader &init, Output &output) {

	// Set the name of the class
	section = "Kessence";

	// Go and get model parameters
	// These are just parameters in the ini file.
	lambda = init.getiniDouble("lambda", 1.0, section);
	alpha = init.getiniDouble("alpha", 1.0, section);
	beta = init.getiniDouble("beta", 1.0, section);
	n = init.getiniDouble("n", 2.0, section);

    // Set the computelagrangian data to uninitialized
    for (int i = 0; i < 4; i++) storeddata[i] = -1;

	// Construct H
	// Temporary variable
	double temp;
	// Scale factor
	double a = data[0];
	double a2 = a * a;

	// Calculate H^2
	temp = params.rhoM() / a + params.rhoR() / a2 + params.rhoK() + a2 * energydensity(data);

	// Calculate H
	data[3] = sqrt(temp);

	// Print stuff to the log
	output.printlog("Running Kessence model.");
	output.printvalue("n", n);
	output.printvalue("lambda", lambda);
	output.printvalue("alpha", alpha);
	output.printvalue("beta", beta);
	output.printvalue("phi0", data[1]);
	output.printvalue("phidot0", data[2]);

	// Success!
	return 0;
}

// The speedofsound2 returns the speed of sound squared, given the state of the system
double Kessence::speedofsound2(const double data[]) {
	// The speed of sound in k-essence can vary from 1.

	// Compute quantities
	computelagrangian(data);

	// Return result
	return LX / (LX + 2.0 * X * LXX);
}

// The isghost function is given the state of the system and returns whether or not the theory has become ghostlike
bool Kessence::isghost(const double data[]) {
	// K-essence becomes ghost-like when the denominator in the speed of sound squared becomes negative.

	// Compute quantities
	computelagrangian(data);

	// Calculate and return the result
	if (LX + 2.0 * X * LXX < 0)
		return true;
	else
		return false;

}
