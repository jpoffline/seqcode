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
	/*
	ofstream dumpconv;
	dumpconv.open( inifile.getiniString("dumpdir", "testconv", "Seq") + "/" + inifile.getiniString("dumpfilename", "test.dat", "Seq") );
	string p1_name = "mass";
	string p1_section = "Quintessence";
	string p2_name = "Omegakh2";
	string p2_section = "Cosmology";
	double p1_min = inifile.getiniDouble("p1_min", 0.0, "Seq");
	double p1_max = inifile.getiniDouble("p1_max", 0.0, "Seq");
	double delta_p1 = (p1_max - p1_min)/inifile.getiniDouble("p1_n", 10.0, "Seq");
	double p2_min = inifile.getiniDouble("p2_min", 0.0, "Seq");
	double p2_max = inifile.getiniDouble("p2_max", 0.0, "Seq");
	double delta_p2 = (p2_max - p2_min)/inifile.getiniDouble("p2_n", 10.0, "Seq");		 
	std::cout << p1_name << " = " << p1_min << " " << p1_max << " " << delta_p1 << "(logspace)" << std::endl;
	std::cout << p2_name << " = " << p2_min << " " << p2_max << " " << delta_p2 << std::endl;
	
	for(double p1 = p1_min; p1 < p1_max; p1 += delta_p1){
		inifile.setparam(p1_name,p1_section,pow(10.0,p1));
		for(double p2 = p2_min; p2 < p2_max; p2 += delta_p2){
			inifile.setparam("evtoend","Seq",true);
			inifile.setparam(p2_name,p2_section,p2);
			Parameters myParams_1(inifile);
			results = doEvolution(inifile, myParams_1, myOutput, SN1adata, dopostprocess);
			// Now do a run to get the likelihood
			inifile.setparam("evtoend","Seq",false);
			Parameters myParams_2(inifile);
			doEvolution(inifile, myParams_2, myOutput, SN1adata, true);
			dumpconv << p1 << " " << p2 << " ";
			for(int n = 0; n < results.size(); n++)
				dumpconv << results[n] << " ";
			dumpconv << myOutput.getvalue("combinationchi", 0.0);
			dumpconv << std::endl;
		}
		dumpconv << std::endl;
	}
	*/
	
	
	// HOPEFULLY, A ROUTINE TO FIND A SEQUESTERING SOLUTION VIA SOME MCMC-ISH METHOD.
	
	
	
	ofstream dumpconv;
	dumpconv.open( inifile.getiniString("dumpdir", "testconv", "Seq") + "/" + inifile.getiniString("dumpfilename", "test.dat", "Seq") );
	
	double HIRl_p, HIRl_c = 1000, Omk_p, tmax, w0, amax;
	double HIRl_t = inifile.getiniDouble("HIRt", 0.1, "Seq");

	double p1_min = inifile.getiniDouble("p1_min", 0.0, "Seq");
	double p1_max = inifile.getiniDouble("p1_max", 0.0, "Seq");
	string p1_name = "mass";
	string p1_section = "Quintessence";
	double delta_p1 = (p1_max - p1_min)/inifile.getiniDouble("p1_n", 10.0, "Seq");
	double p2_min = inifile.getiniDouble("p2_min", 0.0, "Seq");
	double p2_max = inifile.getiniDouble("p2_max", 0.0, "Seq");
	double p2_sigma = (p2_max - p2_min) / 100.0;
	double Omk = p2_min;
	double tarmfrac, a_now;
	for(double p1 = p1_min; p1 < p1_max; p1 += delta_p1){
		
		inifile.setparam(p1_name,p1_section,pow(10.0,p1));
		Omk = p2_min + (p2_max - p2_min) * 0.5 + NormalRand() * p2_sigma;
		HIRl_c = 1000;
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
			results = doEvolution(inifile, myParams_1, myOutput, SN1adata, false);
			// ... and pull out the value of <R>
			HIRl_p = log10(abs(results[1]));
			
			// If <R> for this choice of Omk is smaller than the previous one,
			// then keep it.
			if(HIRl_p < HIRl_c){
				Omk = Omk_p;
				HIRl_c = HIRl_p;
			}
	
			// If <R> is smaller than te desired threshold "error",
			//	then we have found a sequestering solution.	
			if(HIRl_c < HIRl_t){
				amax = results[4];
				w0 = results[5];
				tmax = results[6];
				tarmfrac = results[7];
				a_now = results[8];
				break;
			}
		}
		inifile.setparam("evtoend","Seq",false);
		Parameters myParams_2(inifile);
		doEvolution(inifile, myParams_2, myOutput, SN1adata, true);
		dumpconv << p1 << " " <<  Omk << " " << HIRl_c << " ";
		dumpconv << exp(-0.5 * myOutput.getvalue("combinationchi", 0.0) ) << " ";
		dumpconv << amax << " " << tmax << " " << w0 << " " << tarmfrac << " " << a_now << endl;
		/*
		1: mass
		2: Omk
		3: <R>
		4: L
		5: amax
		6: tmax
		7: w0
		
		*/
	}
	
	
	
	dumpconv.close();
	if(NRflag == 0)
		std::cout << "histIntRicciScalar = " << histIntRicciScalar << endl;	
	
    // Stop timing
    myTimer.stop();
	result = (int)results[0];
	
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
