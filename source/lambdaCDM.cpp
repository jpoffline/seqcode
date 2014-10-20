#include "lambdaCDM.h"

// The derivatives routine is given the state of the system (a, phi and \dot{\phi}) as well as the parameters of the system,
// and returns the derivatives (\dot{a}, \dot{\phi}, \ddot{\phi}, \dot{H})
int LambdaCDM::derivatives(const double data[], double derivs[], Parameters &params) {

	// The first section here is model independent

	// Extract data for easier reading of the code
	double a = data[0];
	double a2 = a * a;   // a^2
	double phidot = data[2];
	double hubble = data[3];

	// Computing \dot{a}, knowing the hubble rate
	derivs[0] = hubble * a;
	// Computing \dot{\phi}
	derivs[1] = phidot;

	// Model dependent code below

	// Computing \ddot{\phi}. This requires solving the scalar equation of motion for \ddot{\phi}.
	derivs[2] = 0;

	// Computing \dot{H}. This requires solving the acceleration equation for \dot{H}.
	// Note that pressure does not depend on \dot{H} in this model,
	// so we pass in 0 for \dot{H} when calculating pressure
	double press = pressure(data, 0.0);
    derivs[3] = 0.5 * (- params.rhoR() / a2 - 3.0 * a2 * press - hubble * hubble + params.rhoK());

	return GSL_SUCCESS;
}

// Returns the ratio rho_Q/rho_0
double LambdaCDM::energydensity(const double data[]){
	return OmegaLambdah2;
}
// Returns the ratio P_Q/rho_0
double LambdaCDM::pressure(const double data[], const double hdot){
	return -OmegaLambdah2;
}

/* This function does four things:
 * - Sets the name of the class
 * - Reads in the value of Omega_Lambda from the ini file
 * - Initializes the value of H using the Friedmann equation
 * - Returns a log output
 */
int LambdaCDM::init(double data[], double time, Parameters &params, IniReader &init, Output &output) {

	// Set the name of the class
	section = "LambdaCDM";

	// We want to construct Omega_Lambda h^2 = Lambda 8 pi G / 3 / H0^2 (where H0 = 100 km/s/Mpc)
    // where the action S = \int d^4x \sqrt{-g} ( m_P^2/2 R - Lambda ) defines Lambda

	// There are two ways of doing this. The first is to read in Omega_Lambda, and use knowledge of all the other density fractions
	// The second is to read in the desired value of h, and use knowledge of all the other density fractions

    if (init.getiniBool("precise", false, section) == true) {
        // Use the value of h
        // We need to know the value of h to calculate this
        double h = init.getiniDouble("desiredh", 0.7, "Cosmology");
        OmegaLambdah2 = h * h - (params.rhoK() + params.rhoM() + params.rhoR());
    } else {
        // Use the given value of Omega_Lambda
        double temp = init.getiniDouble("OmegaLambda", 0.7, section);
        OmegaLambdah2 = temp / (1 - temp) * (params.rhoK() + params.rhoM() + params.rhoR());
    }


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

	// Discard phi, \dot{\phi}
	data[1] = 0;
	data[2] = 0;

	// Print stuff to the log
	output.printlog("Running LambdaCDM model.");
	output.printvalue("OmegaLambdah2", OmegaLambdah2);

    // We have success!
    return 0;

}
