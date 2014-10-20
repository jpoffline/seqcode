#include "quintessence.h"

// To modify this module to a given quintessence potential, the following two functions must be specified
// Returns the potential for a given phi (used in calculating energy density and pressure)
double Quintessence::potential(const double phi){


    double phi2 = phi * phi;

	// The potential is selected by a parameter that is specified in params.ini
	switch(modeltype){
		case 1 :
			// lambda phi^4
			// This is just the potential V = mass * m_P^2 H_0^2 phi^2 / 2 + lambda * m_P^2 H_0^2 phi^4
			// (for dimensionless phi)
			// In terms of the dimensionless potential U, this is U = mass * phi^2 / 2 + lambda * phi^4
			return mass * phi2 / 2 + lambda * phi2 * phi2;

		case 2 :
			// Exponential
			// This is an exponential function: V = alpha m_P^2 H_0^2 exp(- beta phi)
			// (for dimensionless phi)
			// In terms of the dimensionless potential U, this is U = alpha exp(- beta phi)
			return alpha * exp(- beta * phi);

		case 3 :
			// User defined potential
			// (Presently not defined)
			return mass * phi;

		case 0 :
		default : // Catch all
			// mass term only
			// This is just the potential V = mass * m_P^2 H_0^2 phi^2 / 2 (for dimensionless phi)
			// In terms of the dimensionless potential U, this is U = mass * phi^2 / 2
			return mass * phi2 / 2;
	}

	// It's an error if we get to here, but return something sensible
	return 0;

}
// Returns the derivative of the potential for a given phi (used in calculating the scalar equation of motion)
double Quintessence::potentialprime(const double phi){

	// The potential is selected by a parameter that is specified in params.ini
	switch(modeltype){
		case 1 :
			// lambda phi^4
			// This is just the potential V = mass * m_P^2 H_0^2 phi^2 / 2 + lambda * m_P^2 H_0^2 phi^4
			// (for dimensionless phi)
			// In terms of the dimensionless potential U, this is U = mass * phi^2 / 2 + lambda * phi^4
			// The derivative is U' = mass * phi + 4 * lambda * phi^3
			return mass * phi + 4.0 * lambda * phi * phi * phi;

		case 2 :
			// Exponential
			// This is an exponential function: V = alpha m_P^2 H_0^2 exp(- beta phi)
			// (for dimensionless phi)
			// In terms of the dimensionless potential U, this is U = alpha exp(- beta phi)
			// The derivative is U' = - alpha beta exp(- beta phi)
			return - beta * alpha * exp(- beta * phi);

		case 3 :
			// User defined potential
			// (Presently not defined)
			return mass ;

		case 0 :
		default : // Catch all
			// mass term only
			// This is just the potential V = mass * m_P^2 H_0^2 phi^2 / 2 (for dimensionless phi)
			// In terms of the dimensionless potential U, this is U = mass * phi^2 / 2
			// The derivative is U' = mass * phi
			return mass * phi;
	}

	// It's an error if we get to here, but return something sensible
	return 0;

}

/* This function does three things:
 * - Sets the name of the class
 * - Initializes the value of H using the Friedmann equation
 * - Returns a log output
 */
int Quintessence::init(double data[], double time, Parameters &params, IniReader &init, Output &output) {

	// Set class name
	section = "Quintessence";

	// Go and get model parameters
	modeltype = init.getiniInt("PotentialType", 0, section);
	if (modeltype < 0 || modeltype > 3)  // do a quick bit of error checking
		modeltype = 0;
	mass = init.getiniDouble("mass", 1.0, section);
	lambda = init.getiniDouble("lambda", 1.0, section);
	alpha = init.getiniDouble("alpha", 1.0, section);
	beta = init.getiniDouble("beta", 1.0, section);

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

	// Print stuff to the logs
	switch(modeltype){
		case 1 :
			// lambda phi^4
		    output.printlog("Running Quintessence model with lambda phi^4 potential.");
		    output.printvalue("Mass", mass);
		    output.printvalue("lambda", lambda);
			break;
		case 2 :
			// Exponential
			output.printlog("Running Quintessence model with exponential potential.");
            output.printvalue("alpha", alpha);
            output.printvalue("beta", beta);
			break;
		case 3 :
			// User defined potential
			output.printlog("Running Quintessence model with user-defined potential: linear for sequestering.");
            output.printvalue("Mass", mass);
			break;
		case 0 :
		default : // Catch all
			// mass term only
            output.printlog("Running Quintessence model with mass potential.");
            output.printvalue("Mass", mass);
	}
    output.printvalue("phi0", data[1]);
    output.printvalue("phidot0", data[2]);

	// Success!
	return 0;

}

/*
 * The derivatives routine is given the state of the system (a, phi, \dot{\phi} and H)
 * as well as the parameters of the system, and returns the derivatives
 * \dot{a}, \dot{\phi}, \ddot{\phi} and \dot{H}
 */
int Quintessence::derivatives(const double data[], double derivs[], Parameters &params) {

	// The first section here is model independent

	// Extract data for easier reading of the code
	double a = data[0];
	double a2 = a * a;   // a^2
	double phi = data[1];
	double phidot = data[2];
	double hubble = data[3];

	// Computing \dot{a}, knowing the hubble rate
	derivs[0] = hubble * a;
	// Computing \dot{\phi}
	derivs[1] = phidot;

	// Model dependent code below

	// Computing \ddot{\phi}. This requires solving the scalar equation of motion for \ddot{\phi}.
	derivs[2] = - 2.0 * hubble * phidot - a2 * potentialprime(phi);

	// Computing \dot{H}. This requires solving the acceleration equation for \dot{H}.
	// Note that pressure does not depend on \dot{H} in this model,
	// so we pass in 0 for \dot{H} when calculating pressure
	double press = pressure(data, 0.0);
    derivs[3] = 0.5 * (- params.rhoR() / a2 - 3.0 * a2 * press - hubble * hubble + params.rhoK());

	return GSL_SUCCESS;
}

// Returns the ratio rho_Q/rho_0
// where rho_0 = 3 H_0^2 m_P^2 (the factor of H_0^2 m_P^2 is already factored out in the action)
double Quintessence::energydensity(const double data[]){
	// Extract data for easier reading of the code
	double a = data[0];
	double phi = data[1];
	double phidot = data[2];

	return (0.5 * phidot * phidot / a / a + potential(phi)) / 3;
	// The factor of 1/3 is correct. Note that a cosmological constant will contribute
	// 8 pi G Lambda / 3 H_0^2 = Lambda / \rho_c = Omega_Lambda.
}
// Returns the ratio P_Q/rho_0
double Quintessence::pressure(const double data[], const double hdot){
	// Extract data for easier reading of the code
	double a = data[0];
	double phi = data[1];
	double phidot = data[2];

	return (0.5 * phidot * phidot / a / a - potential(phi)) / 3;
}
