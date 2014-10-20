#include "simplecheck.h"

// Checks the state of the system after every timestep
void SimpleCheck::checkstate (const double data[], const double time, IntParams &params, Output &output, const double status[]){
	// We want to check the following conditions:
	// Error in Friedmann equation sufficiently small
	// Equation of state is not phantom
	// Speed of sound is not superluminal
	// Speed of sound is not imaginary
	// Perturbations are not ghost-like
	// Tensor speed is not superluminal
	// Tensor speed is not imaginary
	// Tensors are not ghost-like
	// Energy density is not negative

	// Check the error in the Friedmann equation
	if (abs(status[16]) > 1e-9) {
		std::stringstream message;
		message << "Warning: Relative error in Friedmann equation is larger than 1e-9, time t = " << time;
		output.printlog(message.str());
		friederror = true;
	}

	// Check for phantom equation of state for DE
	if (status[15] < -1) {
		std::stringstream message;
		message << "Warning: Dark energy EOS is phantom, time t = " << time;
		output.printlog(message.str());
		phantom = true;
	}

	// Check for speed of sound
	if (params.getmodel().implementsSOS()) {
		double speed2 = params.getmodel().speedofsound2(data);
		if (speed2 > 1) {
			std::stringstream message;
			message << "Warning: Speed of sound is superluminal, time t = " << time;
			output.printlog(message.str());
			superluminal = true;
		} else if (speed2 < 0) {
			std::stringstream message;
			message << "Warning: Speed of sound is imaginary, time t = " << time;
			output.printlog(message.str());
			laplacian = true;
		}
	}

	// Check for ghosts
	if (params.getmodel().implementsghost()) {
		if (params.getmodel().isghost(data)) {
			std::stringstream message;
			message << "Warning: Perturbations are ghostlike, time t = " << time;
			output.printlog(message.str());
			ghost = true;
		}
	}

	// Check for speed of sound (tensors)
	if (params.getmodel().implementsSOT()) {
		double speed2 = params.getmodel().speedoftensor2(data);
		if (speed2 > 1) {
			std::stringstream message;
			message << "Warning: Tensor speed is superluminal, time t = " << time;
			output.printlog(message.str());
			tsuperluminal = true;
		} else if (speed2 < 0) {
			std::stringstream message;
			message << "Warning: Tensor speed is imaginary, time t = " << time;
			output.printlog(message.str());
			tlaplacian = true;
		}
	}

	// Check for ghosts (tensors)
	if (params.getmodel().implementstensorghost()) {
		if (params.getmodel().istensorghost(data)) {
			std::stringstream message;
			message << "Warning: Tensor perturbations are ghostlike, time t = " << time;
			output.printlog(message.str());
			tghost = true;
		}
	}

	// Check for negative energy density in DE
	// Usually we'll crash from a div/0 error before getting here though!
	if (status[13] < 0) {
		std::stringstream message;
		message << "Warning: Dark energy has negative energy density, time t = " << time;
		output.printlog(message.str());
		negenergy = true;
	}

}

void SimpleCheck::checkfinal (const double data[], const double time, IntParams &params, Output &output, const double status[]){
	// We want to check the following conditions:
	// Equation of state of dark energy is very close to -1

	bool finaleos = false;

	// Check for EOS of DE (warn if the EOS is greater than -0.9)
	if (status[15] > -0.9) {
		output.printlog("Warning: final equation of state of dark energy is > -0.9");
		finaleos = true;
	}

	// Report on any errors received throughout evolution
	if (finaleos || friederror || phantom || superluminal || laplacian || ghost || tsuperluminal || tlaplacian || tghost || negenergy) {
		output.printlog("");
		output.printlog("Warning summary:");
		if (finaleos) output.printvalue("FinalEOS", 1);
		if (friederror) output.printvalue("FriedError", 1);
		if (phantom) output.printvalue("Phantom", 1);
		if (superluminal) output.printvalue("Superluminal", 1);
		if (laplacian) output.printvalue("Laplacian", 1);
		if (ghost) output.printvalue("Ghost", 1);
		if (tsuperluminal) output.printvalue("TensorSuperluminal", 1);
		if (tlaplacian) output.printvalue("TensorLaplacian", 1);
		if (tghost) output.printvalue("TensorGhost", 1);
		if (negenergy) output.printvalue("NegativeEnergy", 1);
		output.printlog("");
	}

}
