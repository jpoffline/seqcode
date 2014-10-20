/*
 * model.h
 *
 * This defines an abstract "Model" class. Any dark energy model can be built off this class. At the very least, it must
 * implement the equations of motion in the derivatives function. Fancier implementations can also compute the speed of sound
 * and check for ghosts.
 *
 * Also available is a state function, which can compute a variety of (useful?) information from the present state.
 *
 */

#ifndef MODEL_H_
#define MODEL_H_

#include "params.h"
#include "inireader.h"
#include "output.h"
#include <cmath>
// stringstreams
#include <iostream>
#include <sstream>
#include <string>
// GSL error numbers
#include <gsl/gsl_errno.h>

// This defines the Model abstract class.
// Individual models will inherit this class and implement the appropriate functions.
class Model {
	public:
		// The derivatives routine is given the state of the system (a, phi and \dot{\phi} as well as H in some cases)
		//  as well as the parameters of the system, and returns the derivatives (\dot{a}, \dot{\phi}, \ddot{\phi} (and \dot{H}, if necessary))
		virtual int derivatives(const double data[], double derivs[], Parameters &params) = 0;
		// This function returns information on the state of the model
		// A mostly-model independent implementation is provided
		virtual void getstate(const double data[], const double time, double info[], Parameters &params);

		// A function to allow the model to initialize itself
		// In particular, this function will have to calculate H from the Friedmann equation as an initial condition
		// if the model is evolving H.
		// It can also read in information from the ini file if desired.
		// The return value is a string to be output to the log
		virtual int init(double data[], const double time, Parameters &params, IniReader &init, Output &output) = 0;

		// The speedofsound2 returns the speed of sound squared, given the state of the system
		virtual double speedofsound2(const double data[]) {return 0;}
		// The implementsSOS function returns whether or not a class actually implements the speedofsound2 function
		virtual bool implementsSOS() {return false;}

		// The isghost function is given the state of the system
		// and returns whether or not the theory has become ghostlike
		virtual bool isghost(const double data[]) {return false;}
		// The implementsghost function returns whether or not a class actually implements the isghost function
		virtual bool implementsghost() {return false;}

		// The speedoftensor2 returns the speed of tensor perturbations squared, given the state of the system
		virtual double speedoftensor2(const double data[]) {return 0;}
		// The implementsSOT function returns whether or not a class actually implements the speedoftensor2 function
		virtual bool implementsSOT() {return false;}

		// The istensorghost function is given the state of the system
		// and returns whether or not the tensor perturbations have become ghostlike
		virtual bool istensorghost(const double data[]) {return false;}
		// The implementstensorghost function returns whether or not a class actually implements the istensorghost function
		virtual bool implementstensorghost() {return false;}

		// This function is used to report whether or not the system is in a valid configuration
		// The most notable example of a system in an invalid configuration is F(R) where F'' = 0
		// If this function returns false, the evolution is aborted
		virtual bool isvalidconfig(const double data[]) = 0;

		// Virtual destructor
		virtual ~Model() {return;}

		// Call to return the class name
		std::string classname() {return section;}

	protected:
		// Functions to return the energy density (given the data), and the pressure (given data and \dot{H})
		virtual double energydensity(const double data[]) = 0;
		virtual double pressure(const double data[], const double hdot) = 0;
		// Name of the class. Should be set in the init routine. Used as the section name in the ini file.
		std::string section;
};

#endif /* MODEL_H_ */
