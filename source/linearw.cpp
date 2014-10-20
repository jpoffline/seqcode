#include "linearw.h"

// The derivatives routine is given the state of the system (a, phi and \dot{\phi}) as well as the parameters of the system,
// and returns the derivatives (\dot{a}, \dot{\phi}, \ddot{\phi}, \dot{H})
int LinearW::derivatives(const double data[], double derivs[], Parameters &params) {

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
	// There is no scalar field in this model, so phidot = 0.
	derivs[2] = 0;

	// Computing \dot{H}. This requires solving the acceleration equation for \dot{H}.
	// Note that pressure does not depend on \dot{H} in this model,
	// so we pass in 0 for \dot{H} when calculating pressure
	double press = pressure(data, 0.0);
    derivs[3] = 0.5 * (- params.rhoR() / a2 - 3.0 * a2 * press - hubble * hubble + params.rhoK());

	// GSL_SUCCESS indicates that the computation was successful. If it failed, look up the appropriate error code in the same enum as GSL_SUCCESS
	return GSL_SUCCESS;
}

// Returns the ratio rho_Q/rho_0
double LinearW::energydensity(const double data[]){
    // First, check the data against the previous data
    // This prevents the results being computed multiple times on the same data
    if (data[0] == storeddata[0] &&
            data[1] == storeddata[1] &&
            data[2] == storeddata[2] &&
            data[3] == storeddata[3]) {
        // It's the same as before, so don't recompute it
        return menergydensity;
    }

	// Extract data for easier reading of the code
	double a = data[0];

	// The formula for energy density is rho = rho_0 a^(-3(1 + w0 + wa)) e^(3 wa (a - 1))
	// menergydensity = pow(a, - 3.0 * (1.0 + w0 + wa)) * exp(3 * wa * (a - 1)) * OmegaLambdah2;
    menergydensity = exp(3 * wa * (a - 1) - 3.0 * (1.0 + w0 + wa) * log(a)) * OmegaLambdah2;

    // Store the data for which these results are correct
    for (int i = 0; i < 4; i++)
        storeddata[i] = data[i];

    // Return the value
    return menergydensity;
}
// Returns the ratio P_Q/rho_0
double LinearW::pressure(const double data[], const double hdot){
	// Extract data for easier reading of the code
	double a = data[0];

	// Just take advantage of the equation of state
	return (w0 + wa * (1 - a)) * energydensity(data);
}

/* This function does four things:
 * - Sets the name of the class
 * - Reads in the value of Omega_Lambda from the ini file
 * - Initializes the value of H using the Friedmann equation
 * - Returns a log output
 */
int LinearW::init(double data[], double time, Parameters &params, IniReader &init, Output &output) {

	// Set the name of the class
	section = "LinearW";

    // We want to construct Omega_Lambda h^2 = Lambda 8 pi G / 3 / H0^2 (where H0 = 100 km/s/Mpc)
    // where Lambda is the energy density today

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

	// Also extract w0 and wa from the ini file
	w0 = init.getiniDouble("wnaught", -1.0, section);
	wa = init.getiniDouble("wa", 0.0, section);

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

	// Discard phi, \dot{\phi}
	data[1] = 0;
	data[2] = 0;

    // Print stuff to the log
    output.printlog("Running LinearW model.");
    output.printvalue("OmegaLambdah2", OmegaLambdah2);
    output.printvalue("wnaught", w0);
    output.printvalue("wa", wa);

    // We have success!
    return 0;

}
