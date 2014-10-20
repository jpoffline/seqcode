/*
 * main.cpp
 *
 * This is the program entry point for the program.
 *
 * The purpose of this software is to evolve scalar field dark energy models through cosmological time.
 *
 * This wrapper performs a sweep of evolutions across one or two parameters, as specified in the ini file.
 * It outputs likelihood information on each parameter combination.
 *
 * This software requires the GSL libraries and the C++ BOOST libraries (see www.boost.org)
 * These are most easily installed using a package manager (libboost-all-dev on ubuntu)
 *
 * Jolyon K. Bloomfield and Jonathan A. Pearson, March 2014
 *
 */

#include "main.h"

// Structure for storing results in a vector
typedef struct expresults {
    double data[7];
} expresults;
const int numchisqds = 7;

using std::cout;
using std::endl;
using std::setprecision;
using std::scientific;

// Our program entry point
// This entry point is just a wrapper around the evolution routines.
// It sets up the input parameters as well as the output file, and otherwise just calls the routines.
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

    // Progressbar stuff
    bool showprogress = inifile.getiniBool("progress", true, "Sweep");
    float progress = 0.0;
    int barcount = 0;


    //**************//
    // Output class //
    //**************//

    // Set up the filenames to output
    string outputdir = inifile.getiniString("logdir", "logs", "Function");
    string basename = inifile.getiniString("runname", "run", "Function");
    string postname = inifile.getiniString("postname", "d", "Function");

    // Go and find our appropriate file name (using 4 digit numbers as the default)
    string outputname = getfilename(outputdir, basename, "", inifile.getiniInt("numberpad", 4, "Function"), true) + ".dat";

    // Set up the output class. Print2Memory is used for sweeps
    Print2Memory myOutput;

    // Check that output is a go
    std::ofstream outputstream(outputname.c_str());
    if (!outputstream.is_open()) {
        // End gracefully if not
        cout << "Unable to open file for output. Terminating." << endl;
        return -1;
    }


    //******************//
    // Set up the sweep //
    //******************//

    int numparams = inifile.getiniInt("numparams", 1, "Sweep");
    if (numparams != 2) numparams = 1;

    // First parameter
    string section1 = inifile.getiniString("section", "Cosmology", "Sweep");
    string param1 = inifile.getiniString("param", "phi0", "Sweep");
    double lower1 = inifile.getiniDouble("lower", -1.0, "Sweep");
    double upper1 = inifile.getiniDouble("upper", 1.0, "Sweep");
    int numsteps1 = inifile.getiniDouble("steps", 20, "Sweep");
    if (numsteps1 < 1) numsteps1 = 1;
    double stepper1;
    double stepsize1 = (upper1 - lower1) / static_cast<double> (numsteps1);
    bool logsteps1 = inifile.getiniBool("logsteps", false, "Sweep");
    double ratio1 = 0;
    if (logsteps1 && (lower1 < 0 || upper1 < 0)) {
        // Return an error; we can't take logs of negative parameter values
        cout << "Error: when using log steps, cannot have negative parameter bounds. Terminating." << endl;
        return -1;
    }
    if (logsteps1) ratio1 = log10(upper1 / lower1) / numsteps1;

    // Second parameter
    string section2 = inifile.getiniString("sectionb", "Cosmology", "Sweep");
    string param2 = inifile.getiniString("paramb", "phi0dot", "Sweep");
    double lower2 = inifile.getiniDouble("lowerb", -1.0, "Sweep");
    double upper2 = inifile.getiniDouble("upperb", 1.0, "Sweep");
    int numsteps2 = inifile.getiniDouble("stepsb", 20, "Sweep");
    if (numsteps2 < 1) numsteps2 = 1;
    double stepsize2 = (upper2 - lower2) / static_cast<double> (numsteps2);
    double stepper2;
    bool logsteps2 = inifile.getiniBool("logstepsb", false, "Sweep");
    double ratio2 = 0;
    if (logsteps2 && (lower2 < 0 || upper2 < 0)) {
        // Return an error; we can't take logs of negative parameter values
        cout << "Error: when using log steps, cannot have negative parameter bounds. Terminating." << endl;
        return -1;
    }
    if (logsteps2) ratio2 = log10(upper2 / lower2) / numsteps2;

    // Report what we're doing to screen
    if (logsteps1)
        cout << "Sweeping over " << param1 << " from " << lower1 << " to " << upper1 << " in " << numsteps1 << " log steps." << endl;
    else
        cout << "Sweeping over " << param1 << " from " << lower1 << " to " << upper1 << " in " << numsteps1 << " steps." << endl;
    if (numparams == 2) {
        if (logsteps2)
            cout << "Sweeping over " << param2 << " from " << lower2 << " to " << upper2 << " in " << numsteps2 << " log steps." << endl;
        else
            cout << "Sweeping over " << param2 << " from " << lower2 << " to " << upper2 << " in " << numsteps2 << " steps." << endl;
    }

    // Note that we increase the number of steps by one, as this is the number of samples that we'll actually be doing
    numsteps1++;
    numsteps2++;
    if (numparams == 1) numsteps2 = 1;
    int totnumsteps = numsteps1 * numsteps2;
    cout << "Total number of samples: " << totnumsteps << endl;
    cout << "Outputting to " << outputname << endl;

    // Storage for various values
    vector<double> parameter1;
    vector<double> parameter2;
    vector<expresults> chisquareds;
    expresults filling;
    // Allocate storage space
    parameter1.reserve(totnumsteps);
    parameter2.reserve(totnumsteps);
    chisquareds.reserve(totnumsteps);

    // Load SN1a data
    vector<vector<double> > SN1adata;
    string sn1afile = inifile.getiniString("union21", "SCPUnion2.1_mu_vs_z.txt", "Function");
    if (loadSN1adata(sn1afile, SN1adata) != 0) {
        // Could not find data file
        cout << "Warning: cannot find Union2.1 SN1a data file." << endl;
    }


    //*******************//
    // Perform the sweep //
    //*******************//

    // Start timing!
    boost::timer::cpu_timer myTimer;

    // Loop through the parameter space
    // Second parameter first (only runs through once if not being used)
    for (int i2 = 0; i2 < numsteps2; i2++) {
        // If we are doing a two parameter scan, update the second parameter here
        if (numparams == 2) {
            if (!logsteps2) {
                stepper2 = lower2 + stepsize2 * i2;
                inifile.setparam(param2, section2, stepper2);
            } else {
                stepper2 = lower2 * pow(10, ratio2 * i2);
                inifile.setparam(param2, section2, log10(stepper2));
            }
        }

        // First parameter next
        for (int i1 = 0; i1 < numsteps1; i1++) {
            // Update the first parameter here
            if (!logsteps1) {
                stepper1 = lower1 + stepsize1 * i1;
                inifile.setparam(param1, section1, stepper1);
            } else {
                stepper1 = lower1 * pow(10, ratio1 * i1);
                inifile.setparam(param1, section1, log10(stepper1));
            }

            // Set up the cosmological parameters (done here in case something significant changed in the sweeping)
            Parameters myParams(inifile);

            // Do the evolution
            result = doEvolution(inifile, myParams, myOutput, SN1adata, true);

            // Interpret the result of the evolution
            if (result == 0) {
                // Everything was successful. Now we can save the results!

                // Store the parameter values
                parameter1.push_back(stepper1);
                if (numparams == 2) parameter2.push_back(stepper2);

                // Extract the chi^2 values
                filling.data[0] = myOutput.getvalue("WMAPchi", 0.0);
                filling.data[1] = myOutput.getvalue("PLANCKchi", 0.0);
                filling.data[2] = myOutput.getvalue("SNchi", 0.0);
                filling.data[3] = myOutput.getvalue("Hubblechi", 0.0);
                filling.data[4] = myOutput.getvalue("BAOtotalchi", 0.0);
                filling.data[5] = myOutput.getvalue("BAOtotalrchi", 0.0);
                filling.data[6] = myOutput.getvalue("combinationchi", 0.0);
                // Plop that on the stack too!
                chisquareds.push_back(filling);

            }

            // Update the progress bar every 20 cycles
            if (showprogress && ++barcount >= 20) {
                barcount = 0;
                progress = (float) (i2 * numsteps1 + i1 + 1) / (float) totnumsteps;
                updateprogress(progress);
            }

        }

        // Update the progress bar (to make sure it shows 100% at the end)
        barcount = 0;
        progress = (float) ((i2 + 1) * numsteps1) / (float) totnumsteps;
        if (showprogress) updateprogress(progress);
    }
    // Clear the progressbar
    if (showprogress) std::cout << std::endl;

    //****************//
    // Postprocessing //
    //****************//

    // We have the chi^2 values for all experiments for each value of the parameters that we've scanned over
    // What we want to do now is to convert these chi^2 values into likelihood values

    // Storage
    vector<expresults> likelihoods;
    totnumsteps = chisquareds.size(); // Just in case some evolution failed and the actual number is less than we expected
    likelihoods.reserve(totnumsteps);
    // Iterate over everything, calculating the likelihood L = e^{-chi2/2}
    for (int i = 0; i < totnumsteps; i++) {
        filling = chisquareds[i];
        for (int j = 0; j < numchisqds; j++) {
            filling.data[j] = exp(- 0.5 * filling.data[j]);
        }
        likelihoods.push_back(filling);
    }

    // Find the highest likelihood values, so that we can normalize using them
    for (int j = 0; j < numchisqds; j++) {
        filling.data[j] = 0;
    }
    for (int i = 0; i < totnumsteps; i++) {
        for (int j = 0; j < numchisqds; j++) {
            if (likelihoods[i].data[j] > filling.data[j])
                filling.data[j] = likelihoods[i].data[j];
        }
    }

    // Normalize the likelihoods
    for (int i = 0; i < totnumsteps; i++)
        for (int j = 0; j < numchisqds; j++)
            likelihoods[i].data[j] /= filling.data[j];

    // Output the likelihood data to file!
    // Print the heading
    outputstream << scientific << setprecision(8) << "# " << param1;
    if (numparams == 2) outputstream << "\t" << param2;
    outputstream << "\tWMAP\tPLANCK\tSN\tHubble\tBAO (SDSS)\tBAO (SDSSR)\tCombined" << endl;
    // Loop over all results, printing them too
    for (int i = 0; i < totnumsteps; i++) {
		// Print a new line every time parameter2 increments -- for gnuplot purposes
		if( i != 0 )
			if( numparams == 2 && parameter2[i] != parameter2[i-1] )
				outputstream << endl;
		
        outputstream << parameter1[i];
        if (numparams == 2) outputstream << "\t" << parameter2[i];
        for (int j = 0; j < numchisqds; j++)
            outputstream << "\t" << likelihoods[i].data[j];
        outputstream << endl;
    }
    outputstream.close();


    //**********//
    // Clean up //
    //**********//

    // Stop timing
    myTimer.stop();
    double ms = myTimer.elapsed().wall / 1e6;
    if (ms < 1e3)
        cout << setprecision(4) << "Sweep complete in " << ms << " milliseconds.";
    else
        cout << setprecision(4) << "Sweep complete in " << ms / 1000.0 << " seconds.";
    cout << endl;

	// Exit gracefully
	return 0;
}
