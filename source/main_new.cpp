/*
 * main.cpp
 *
 * The purpose of this software is to evolve scalar field dark energy models through cosmological time.
 *
 * This wrapper performs a single evolution. It outputs a detailed log of its data and various integrated distance measures.
 *
 * This software requires the GSL libraries and the C++ BOOST libraries (see www.boost.org)
 * These are most easily installed using a package manager (libboost-all-dev on ubuntu)
 *
 * Jolyon K. Bloomfield and Jonathan A. Pearson, March 2014
 *
 */

// Include all supporting machinery
#include "main.h"

using std::cout;
using std::endl;
using std::setprecision;
using namespace std;
// Our program entry point
// This entry point is just a wrapper around the evolution routines.
// It sets up the input parameters as well as the output file, and otherwise just calls the routines to run the evolution.
int main(int argc, char* argv[]) {

    int result; // Just a number for return values

    //******************//
    // Input parameters //
    //******************//

	// Read the input parameters file, named "params.ini" by default.
	// If there is a command line argument, assume that it is the filename for the parameters file.
	IniReader inifile;
	if (argc > 1)
		inifile.read(argv[1]);
	else
		inifile.read("params.ini");

	// Values in the ini file can be set using the following command:
    // inifile.setparam("desiredh", "Cosmology", 1);
	// First parameter is the key name, the second is the section name, and the third is the value, either an integer, string or double

    // Set up the cosmological parameters
    Parameters myParams(inifile);

    // Figure out if we're going to postprocess or not
    bool dopostprocess = inifile.getiniBool("postprocess", false, "Single");

	// Turn off post-processing for now
	if(dopostprocess)
		dopostprocess = false;

    //**************//
    // Output class //
    //**************//

	/*

    // Set up the filenames to output
    string outputdir = inifile.getiniString("logdir", "logs", "Function");
    string basename = inifile.getiniString("runname", "run", "Function");
    string postname = inifile.getiniString("postname", "d", "Function"); // Only used if dopostprocess is true

    // Go and find our appropriate file name (using 4 digit numbers as the default)
    string outputname = getfilename(outputdir, basename, postname, inifile.getiniInt("numberpad", 4, "Function"), dopostprocess);

    // Set up the output class
    std::string parsestring = inifile.getiniString("outputclass", "BasicDump", "Single");
    Output *myOutput;
    if (parsestring == "BasicDump")
        myOutput = new BasicDump(dopostprocess, outputname, postname);
    else
        myOutput = new BasicDump(dopostprocess, outputname, postname);    // BasicDump is the default

    // Check that output is a go
    if (!myOutput->filesready()) {
        // End gracefully if not
        cout << "Unable to open files for output." << endl;
        delete myOutput;
        return -1;
    }
*/
	Print2Memory myOutput;
    // Load SN1a data
    // Not so much here in a single run, but when doing multiple runs, it makes sense to store the SN1a data in memory
    // rather than reading the file each time it's needed.
    // Here, we load up the SN1a data (if postprocessing is true)
    vector<vector<double> > SN1adata;
    if (dopostprocess) {
        string sn1afile = inifile.getiniString("union21", "SCPUnion2.1_mu_vs_z.txt", "Function");
        if (loadSN1adata(sn1afile, SN1adata) != 0) {
            // Could not find data file
            myOutput.printlog("Warning: cannot find Union2.1 SN1a data file.");
            cout << "Warning: cannot find Union2.1 SN1a data file." << endl;
        }
    }


    //*******************//
    // Do the evolution! //
    //*******************//

    // Print some stuff to the screen
    cout << "Beginning evolution." << endl;
    //cout << "Outputting to " << outputname << endl;

    // Start timing!
    boost::timer::cpu_timer myTimer;

	// parameter name to be dialed
	string paramname = "phi0"; 
	// Section name of dialed parameter
	string paramsection = "Cosmology"; 
	// Value of the parameter which will be returned as the sequestering-compatible parameter
	double param_seq;
	// Start value of the parameter to be dialed
	double paramstart = inifile.getiniDouble("paramstart", 0.0, "Seq");
	double paramend = inifile.getiniDouble("paramend", 2.5, "Seq");
	// Required accuracy of the sequestering constraint
	int maxSteps = inifile.getiniInt("maxSteps", 2000, "Seq");
	double DeltaParam = ( paramend - paramstart) / maxSteps;
	
	int Steps = 0;
	int NRflag = 0;
	double param = paramstart;
	double histIntRicciScalar;

	vector<double> results;
	results.push_back(0.0); // "result"
	results.push_back(0.0); // "<R>"
	results.push_back(0.0); // "R_end"
	results.push_back(0.0); // "a_end"	
	ofstream dumpconv;
	dumpconv.open( inifile.getiniString("dumpdir", "testconv", "Seq") + "/" + inifile.getiniString("dumpfilename", "test.dat", "Seq") );
	/*
	while(true){
		
		inifile.setparam(paramname,paramsection,param);
	    results = doEvolution(inifile, myParams, myOutput, SN1adata, dopostprocess);
		result = int(results[0]);
		histIntRicciScalar = results[1];
		
		if( Steps > maxSteps ){
			NRflag = 1;
			break;
		}
		
		dumpconv << param << " " << results[1] << " " << results[2] << std::endl;
		param = param + DeltaParam;
		Steps++;
		
	}
	*/
	double p1_min = inifile.getiniDouble("p1_min", 0.0, "Seq");
	double p1_max = inifile.getiniDouble("p1_max", 0.0, "Seq");
	double delta_p1 = (p1_max - p1_min)/inifile.getiniDouble("p1_n", 10.0, "Seq");
	double p2_min = inifile.getiniDouble("p2_min", 0.0, "Seq");
	double p2_max = inifile.getiniDouble("p2_max", 0.0, "Seq");
	double delta_p2 = (p2_max - p2_min)/inifile.getiniDouble("p2_n", 10.0, "Seq");		 
		
	std::cout << "p1 " << p1_min << " " << p1_max << " " << delta_p1 << std::endl;
	std::cout << "p2 " << p2_min << " " << p2_max << " " << delta_p2 << std::endl;
	for(double p1 = p1_min; p1 < p1_max; p1 += delta_p1){
		inifile.setparam("smallesta","Seq",p1);
		for(double p2 = p2_min; p2 < p2_max; p2 += delta_p2){
			inifile.setparam("phi0","Cosmology",p2);
			results = doEvolution(inifile, myParams, myOutput, SN1adata, dopostprocess);
			dumpconv << p1 << " " << p2 << " " << results[1] << " " << results[2] << " " << results[3] << std::endl;
		}
		dumpconv << std::endl;
	}
	
	dumpconv.close();
	
	if(NRflag == 0)
		std::cout << "histIntRicciScalar = " << histIntRicciScalar << endl;
	
	
    // Stop timing
    myTimer.stop();

    // Interpret the result of the evolution
    if (result == 0) {
        // Success! Print a nice message
        myOutput.printfinish(myTimer.elapsed().wall / 1e6);
        cout << setprecision(4) << "Evolution complete in " << myTimer.elapsed().wall / 1e6 << " milliseconds." << endl;
    }
    else if (result == -1) {
        // Initialization error
        cout << "Error initializing model; terminating." << endl;
    }
    else if (result == 1) {
        // Integration error
        cout << "Integration error; terminating." << endl;
    }
    else if (result == 2) {
        // NAN error
        cout << "NAN error; terminating." << endl;
    }
    else if (result == 3) {
        // Did not get t a = 1 error
        cout << "Did not evolve to a = 1 within appropriate conformal time; terminating." << endl;
    }
    else if (result == 4) {
        // Invalid state error
        cout << "Invalid state reached; terminating." << endl;
    }
    else if (result == 5) {
        // Invalid state error
        cout << "Error integrating distance measures; terminating." << endl;
    }


    //**********//
    // Clean up //
    //**********//

    // No memory leaks!
//    delete *myOutput;

	// Exit gracefully
	return 0;
}
