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
	
    // Load SN1a data
    // Not so much here in a single run, but when doing multiple runs, it makes sense to store the SN1a data in memory
    // rather than reading the file each time it's needed.
    // Here, we load up the SN1a data (if postprocessing is true)
    vector<vector<double> > SN1adata;
    if (dopostprocess) {
        string sn1afile = inifile.getiniString("union21", "SCPUnion2.1_mu_vs_z.txt", "Function");
        if (loadSN1adata(sn1afile, SN1adata) != 0) {
            // Could not find data file
            myOutput->printlog("Warning: cannot find Union2.1 SN1a data file.");
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

	
	double HIRl_p, HIRl_c = 1000, Omk_p, tmax, w0, amax, a_now, tarmfrac;
	double HIRl_t = inifile.getiniDouble("HIRt", 0.1, "Seq");
	double p1_min = inifile.getiniDouble("p1_min", 0.0, "Seq");
	double p1_max = inifile.getiniDouble("p1_max", 0.0, "Seq");
	double p2_min = inifile.getiniDouble("p2_min", 0.0, "Seq");
	double p2_max = inifile.getiniDouble("p2_max", 0.0, "Seq");
	double p2_sigma = (p2_max - p2_min) / 100.0;
	double Omk = p2_min;
	vector<double> results2;
	Omk = p2_min + (p2_max - p2_min) * 0.5 + NormalRand() * p2_sigma;
	HIRl_c = 1000;
	cout << "Finding a sequestering solution..." << endl;
	while(true){
		
		// Pick a new Omk that is inside the prior range
		while(true){
			Omk_p = Omk + NormalRand() * p2_sigma;
			if(Omk_p <= p2_max && Omk_p >= p1_min)
				break;
		}

		// Now that we have a sensible choice of Omk,
		// run the evolver...
		inifile.setparam("evtoend","Seq",true);
		inifile.setparam("Omegakh2","Cosmology",Omk_p);
		Parameters myParams_1(inifile);
		results2 = doEvolution(inifile, myParams_1, *myOutput, SN1adata, false);
		// ... and pull out the value of <R>
		HIRl_p = log10(abs(results2[1]));
		
		// If <R> for this choice of Omk is smaller than the previous one,
		// then keep it.
		if(HIRl_p < HIRl_c){
			Omk = Omk_p;
			HIRl_c = HIRl_p;
		}
		
		// If <R> is smaller than te desired threshold "error",
		//	then we have found a sequestering solution.	
		if(HIRl_c < HIRl_t){
			amax = results2[4];
			w0 = results2[5];
			tmax = results2[6];
			tarmfrac = results2[7];
			a_now = results2[8];
			break;
		}
	}
	cout << "A sequestering solution has been found with Omkh2 = " << Omk << endl;
	inifile.setparam("evtoend","Seq",false);
	Parameters myParams_2(inifile);
	
	
	vector<double> results = doEvolution(inifile, myParams_2, *myOutput, SN1adata, dopostprocess);
	//vector<double> results = doEvolution(inifile, myParams, *myOutput, SN1adata, dopostprocess);
	
	cout << endl;
	cout << "outputname = " << outputname << endl;	
	cout << "phi0 = " << inifile.getiniDouble("phi0", 0.0, "Cosmology") << endl;
	cout << "phidot0 = " << inifile.getiniDouble("phidot0", 0.0, "Cosmology") << endl;	
	cout << "mass3 = " << inifile.getiniDouble("mass", 0.0, "Quintessence") << endl;	
	cout << "Omega_k h^2 = " << inifile.getiniDouble("Omegakh2", 0.0, "Cosmology") << endl;
	cout << "<R> = " << results2[1] << endl;
	cout << "R_end = " << results2[2] << endl;
	cout << "a_end = " << results2[3] << endl;
	cout << "a_max = " << results2[4] << endl;
	cout << "1 + w0 = " << 1.0 + results[5] << endl;
	cout << "t_final = " << results2[6] << endl;
	cout << "t_final / t0 = " << results2[7] << endl;
	cout << endl;
	
    // Stop timing
    myTimer.stop();
	result = results[0];
    // Interpret the result of the evolution
    if (result == 0) {
        // Success! Print a nice message
        myOutput->printfinish(myTimer.elapsed().wall / 1e6);
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
    delete myOutput;

	// Exit gracefully
	return 0;
}
