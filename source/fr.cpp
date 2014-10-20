#include "fr.h"

// The derivatives routine is given the state of the system (a, phi and \dot{\phi}) as well as the parameters of the system,
// and returns the derivatives (\dot{a}, \dot{\phi}, \ddot{\phi}, \dot{H})
int Fr::derivatives(const double data[], double derivs[], Parameters &params) {

    // The first section here is model independent

    // Extract data for easier reading of the code
    double a = data[0];
    double a2 = a * a;
    double phi = data[1];
    double phidot = data[2];
    double hubble = data[3];
    double X = phidot * phidot / 2.0 / a2;

    // Computing \dot{a}, knowing the hubble rate
    derivs[0] = hubble * a;
    // Computing \dot{\phi}
    derivs[1] = phidot;

    // Model dependent code below

    // First, update the Lagrangian values
    computelagrangian(data);

    // Computing \dot{H}. This requires solving the acceleration equation for \dot{H}.
    // This turns out to be surprisingly straightforward from the phi equation of motion
    derivs[3] = - hubble * hubble + params.rhoK() + a2 * phi / 6.0;

    // Computing \ddot{\phi}. This requires solving the scalar equation of motion for \ddot{\phi}.
    // This one is a mess, however.
    double temp1 = a2 * (valF + 4.0 * valFppp * X - phi * valFp);
    double temp2 = 2.0 * ((1 + valFp) * (hubble * hubble + 2.0 * derivs[3] - params.rhoK()) + valFpp * hubble * phidot);
    double temp3 = 2.0 * params.rhoR() / a2;
    double denom = 2.0 * valFpp;
    derivs[2] = - (temp1 + temp2 + temp3) / denom;

    // GSL_SUCCESS indicates that the computation was successful. If it failed, look up the appropriate error code in the same enum as GSL_SUCCESS
    return GSL_SUCCESS;
}

// Returns the ratio rho_Q/rho_c
double Fr::energydensity(const double data[]){
    // Extract data for easier reading of the code
//    double a = data[0];
//    double phi = data[1];
//    double phidot = data[2];
//    double hubble = data[3];

    // Your code here
    return 0;
}
// Returns the ratio P_Q/rho_c
double Fr::pressure(const double data[], const double hdot){
    // Extract data for easier reading of the code
//    double a = data[0];
//    double phi = data[1];
//    double phidot = data[2];
//    double hubble = data[3];

    // Your code here
    return 0;
}

/* This function does a few things:
 * - Sets the name of the class
 * - Initializes the value of H using the Friedmann equation
 * - Returns a log output
 */
int Fr::init(double data[], double time, Parameters &params, IniReader &init, Output &output) {

    // Set the name of the class
    section = "Fr";

    // These are just parameters in the ini file.
    alpha = init.getiniDouble("alpha", 1.0, section);
    n = init.getiniDouble("n", 2.0, section);

    // Set the computelagrangian data to uninitialized
    for (int i = 0; i < 4; i++) storeddata[i] = -1;

//    // Construct H
//    // Temporary variable
//    double temp;
//    // Scale factor
//    double a = data[0];
//    double a2 = a * a;
//    // Scalar field
//    double phi = data[1];
//    double phidot = data[2];

    // There are two values of H; construct both of them and choose the positive one
    // If there are two positive roots, use the ini file to determine which to take
    // If there are no real roots, then return an error

    // First, calculate the Lagrangian values
    computelagrangian(data);

    data[3] = 0;

    // Print stuff to the log
    output.printlog("Running F(R) model.");
    output.printvalue("n", n);
    output.printvalue("alpha", alpha);
    output.printvalue("phi0", data[1]);
    output.printvalue("phidot0", data[2]);

    // Success!
    return 0;

}

// Function to calculate the Lagrangian and all its appropriate derivatives:
// F, F', F'', F'''
int Fr::computelagrangian(const double data[]) {
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

    // Our function F has the following form:
    // F(phi) = \alpha phi^n
    // Go calculating!

    valF = alpha * pow(phi, n);
    valFp = n * alpha * pow(phi, n - 1.0);
    valFpp = n * (n - 1.0) * alpha * pow(phi, n - 2.0);
    valFppp = n * (n - 1.0) * (n - 2.0) * alpha * pow(phi, n - 3.0);

    // Store the data for which these results are correct
    for (int i = 0; i < 4; i++)
        storeddata[i] = data[i];

    // Success!
    return 0;
}

// Function to determine if we are in a valid configuration
bool Fr::isvalidconfig(const double data[]) {
    // First, update the Lagrangian values
    computelagrangian(data);

    // The requirement for a valid configuration is that F'' != 0
    if (valFpp == 0.0) return false;

    return true;
}

// Function to determine if the model is ghostlike or not
bool Fr::isghost(const double data[]) {
    // The model is not ghostlike so long as 1 + F' > 0

    // First, update the Lagrangian values
    computelagrangian(data);

    // Perform the check
    if (1.0 + valFp < 0) return true;

    return false;
}
