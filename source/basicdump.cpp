#include "basicdump.h"

using std::cout;
using std::endl;
using std::scientific;
using std::setprecision;

// Function to print data before the run starts
void BasicDump::printinfo(const double data[], Parameters &params) {

	// Require scientific notation, 8 digits
	myLog->precision(8);

	// Print the header
	*myLog << "# Cosmological parameters" << endl;
	*myLog << "OmegaMh2 = " << params.rhoM() << endl;
	*myLog << "OmegaBh2 = " << params.rhoB() << endl;
	*myLog << "OmegaKh2 = " << params.rhoK() << endl;
    *myLog << "OmegaRh2 = " << params.rhoR() << endl;
	*myLog << "Tgamma = " << params.Tgamma() << endl; // K
	*myLog << "zinit = " << params.z0() << endl;
	*myLog << "Neff = " << params.Neff() << endl << endl;

}

// Function to print a header before output starts (e.g., column headings)
void BasicDump::printheading() {

	// Dumping everything. Make nice headers.
	//*myData << "time, a, redshift, H, Hdot, phi, phidot, phiddot, Omega_m, Omega_r, Omega_k, Omega_Q, w_total, rho_Q/rho_c, P_Q/rho_c, w_Q, Error" << endl;
	*myLog << "# Beginning run" << endl;
	*myLog << "# Columns in data file are as follows:"
		   << endl
		   << "# time, a, redshift, H, Hdot, phi, phidot, phiddot, Omega_m, Omega_r, Omega_k, Omega_Q, w_total, rho_Q/rho_c, P_Q/rho_c, w_Q, Error, R"
		   << endl << endl;

	// Set up the output form
	*myData << scientific;
	myData->precision(15);

}

// Function to print information after each timestep
void BasicDump::printstep(const double data[], const double time, const double status[]) {

	// We are just dumping everything here, except for consistency flags
	for (int i = 0; i < 16; i++)
		*myData << status[i] << ", ";
//	*myData << status[16] << endl;
	// Make it also print out Ricci scalar
	*myData << status[16] << ", " << 6.0 * (status[4] + status[3] * status[3]) / status[1] / status[1] << endl;

}

// Function to print a header before postprocessed output starts (e.g., column headings)
void BasicDump::postprintheading() {

    if (!doingpp) return; // We shouldn't be here if we were told that postprocessing wasn't happening

	// Dumping everything. Make nice headers.
	*myLog << "# Beginning postprocessing" << endl;
	*myLog << "# Columns in post-processed data file are as follows:"
		   << endl
		   << "# redshift, DC, DM, DA, DL, mu"
		   << endl
		   << "# DC, DM, DA and DL are all in Mpc."
		   << endl << endl;

	// Set up the output form
	*myPostData << scientific;
	myPostData->precision(15);

	// Print headings
    *myPostData << "# redshift, DC, DM, DA, DL, mu" << endl;
}

// Function to print information after each timestep
void BasicDump::postprintstep(const double z, const double H, const double DC, const double DM, const double DA, const double DL, const double mu) {

    if (!doingpp) return; // We shouldn't be here if we were told that postprocessing wasn't happening

    // We are just dumping everything here
	*myPostData << z << ", "
			    << DC << ", "
			    << DM << ", "
			    << DA << ", "
			    << DL << ", "
			    << mu << endl;

}

// Constructor
BasicDump::BasicDump(bool postprocess, const std::string &filename, const std::string &postname) {
    // Are we postprocessing? If not, don't open the postprocess file.
    doingpp = postprocess;

	// Construct filenames
	string logfile = filename + ".log";
	string datfile = filename + ".dat";
	string postdatfile = filename + postname + ".dat";

	// Open the log files
	myLog = new std::ofstream(logfile.c_str());
	myData = new std::ofstream(datfile.c_str());
	if (doingpp) myPostData = new std::ofstream(postdatfile.c_str());

}

// Destructor
BasicDump::~BasicDump() {

	// Close the log files...
	myLog->close();
	myData->close();
	if (doingpp) myPostData->close();

	// ... and release memory
	delete myLog;
	delete myData;
	if (doingpp) delete myPostData;

}

// Did the files open correctly?
bool BasicDump::filesready() {

	if (myLog->is_open() && myData->is_open()) {
	    if (doingpp) {
	        if (myPostData->is_open())
	            return true;
	        else
	            return false;
	    }
	    else
	        return true;
	}
	else
		return false;

}

// Print a "We're done!" message
void BasicDump::printfinish(const double time) {
	*myLog << setprecision(4) << endl << "# Evolution complete in " << time << " milliseconds." << endl;
}

// Print a line to the log file
void BasicDump::printlog(const std::string &output) {
	if (output == "") {
		*myLog << endl;
	} else {
		*myLog << "# " << output << endl;
	}
}

// Print a key and string to the log file
void BasicDump::printvalue(const std::string &key, const std::string &value) {
	*myLog << key << " = " << value << endl;
}

// Print a key and double to the log file
void BasicDump::printvalue(const std::string &key, const double value) {
    *myLog << key << " = " << setprecision(8) << value << endl;
}

// Print a key and integer to the log file
void BasicDump::printvalue(const std::string &key, const int value) {
    *myLog << key << " = " << value << endl;
}
