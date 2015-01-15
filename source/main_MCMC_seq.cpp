/*
 * main.cpp
 *
 * This is the program entry point for the program.
 *
 * It acts as a wrapper around the evolution routines. All that this wrapper does is to set
 * up the appropriate input and output objects.
 *
 * The purpose of this software is to evolve scalar field models through cosmological time.
 *
 * This software requires the GSL libraries.
 *
 * This software requires the C++ BOOST libraries (see www.boost.org)
 * These are most easily installed using a package manager (libboost-all-dev on ubuntu)
 *
 * Jolyon K. Bloomfield and Jonathan A. Pearson, March 2014
 *
 */

#include "main.h"

using namespace std;

double ComputeLikelihood(IniReader& inifile, Print2Memory& output, string *names, string *sections, double *parameters, int numparams, bool usingSN1a, vector<vector<double> > &SN1adata);
void GetProposedParameters(double *priors, double *current, double *proposed, bool *logs, int numparams);

// Our program entry point
int main(int argc, char* argv[]) {

    // Start timing!
    boost::timer::cpu_timer myTimer;

    //******************//
    // Input parameters //
    //******************//

    // Read the input parameters file, named "params.ini" by default.
    // If there is a command line argument, assume that it is the filename for the parameters file.
    IniReader inifile;
    string paramsfile;
    if (argc > 1)
        paramsfile = argv[1];
    else
        paramsfile = "params.ini";
    inifile.read(paramsfile.c_str());

    //**************//
    // Declarations //
    //**************//

    // Print2Memory is used for MCMC outputting
    Print2Memory myOutput;

    // Seed random number generator
    RNGtool.seed(time(NULL));

    // Output directories and filenames
    string chaindir = inifile.getiniString("chaindir", "chains", "MCMC");
    string chainsubdir = inifile.getiniString("chainsubdir", "run1", "MCMC");
    string outputdir = chaindir + "/" + chainsubdir;
    string basename = inifile.getiniString("chainprefix", "chain", "MCMC");
    string priorsfile = inifile.getiniString("priorsfile", "priors.txt", "MCMC");

    // If the priors file doesn't exist, abort
    // Make sure that the directories exist
    if (!boost::filesystem::exists(priorsfile)) {
        cout << "Priors file (" << priorsfile << ") does not exist. Aborting." << endl;
        return -1;
    }

    // Make sure that the directories exist
    if (!boost::filesystem::exists(chaindir + "/")) {
        // Directory doesn't exist. Make it.
        boost::filesystem::create_directory(chaindir);
        std::cout << "Creating directory " << chaindir << "/" << std::endl;
    }
    if (!boost::filesystem::exists(outputdir + "/")) {
        // Directory doesn't exist. Make it.
        boost::filesystem::create_directory(outputdir);
        std::cout << "Creating directory " << outputdir << "/" << std::endl;
    }


    // Copy the params.ini and priors file being used into the subdir
    if (boost::filesystem::exists(paramsfile.c_str())) {
        boost::filesystem::copy_file(paramsfile, outputdir + "/params.ini", boost::filesystem::copy_option::overwrite_if_exists);
        boost::filesystem::copy_file(priorsfile, outputdir + "/priors.txt", boost::filesystem::copy_option::overwrite_if_exists);
    }

    // MCMC parameters
    // Number of MCMC steps to burn
    int MCMCburninsteps = inifile.getiniInt("MCMCburninsteps", 1000, "MCMC");
    // Number of MCMC steps to take
    int MCMCnumsteps = inifile.getiniInt("MCMCnumsteps", 40000, "MCMC");
    // Number of MCMC chains
    int MCMCnumchains = inifile.getiniInt("numchains", 5, "MCMC");

    // Various variables
    int MCMCstep, MCMCaccept_counter;
    double lower, upper;
    bool usingSN1a;

    // Variables for holding the current and proposed likelihoods, and their ratio
    double L_current, L_proposed, LikelihoodRatio;

    // Progressbar stuff
    bool showprogress = inifile.getiniBool("progress", true, "MCMC");
    float progress = 0.0;
    int barcount = 0;


    //****************//
    // Initialization //
    //****************//

    // Get the priors
	ifstream priorsin;
    string line;
	priorsin.open(priorsfile.c_str());
	vector<PARAMPRIORS> spriors;
	string dummy;
	if(priorsin){
		while(!priorsin.eof()){
		    // Extract the line
		    getline(priorsin, line);
		    // Check whether the line is a comment
		    if (line[0] != '#' && line.length() > 0) {
		        if (line[0] == 'L' && line[1] == ' ') {
		            // Log parameter
                    // Convert the string into a stringstream for extraction
                    stringstream stream(line);
                    PARAMPRIORS temp;
                    stream >> dummy >> temp.section >> temp.name >> temp.lower >> temp.upper >> temp.sigma;
                    temp.logparam = true;
                    // Make sure that the upper and lower bounds are positive
                    if (temp.upper < 0 or temp.lower < 0) {
                        cout << "Error: bounds for " << temp.name << " must be positive in order to investigate parameter in log space. Terminating." << endl;
                        return -1;
                    }
                    spriors.push_back(temp);
		        } else {
                    // Convert the string into a stringstream for extraction
                    stringstream stream(line);
                    PARAMPRIORS temp;
                    stream >> temp.section >> temp.name >> temp.lower >> temp.upper >> temp.sigma;
                    temp.logparam = false;
                    spriors.push_back(temp);
		        }
		    }
		}
		priorsin.close();
	}
	
	// Number of parameters
	int numparams = spriors.size();
	
	// Report priors info to screen
	cout << "The priors are:" << endl;
	for(int n = 0; n < numparams; n++) {
	    if (spriors[n].logparam) cout << "Log10 ";
		cout << spriors[n].section << " " << spriors[n].name
			 << " :: " << spriors[n].lower << "\t" << spriors[n].upper << "\t" << spriors[n].sigma << endl;
	}
	
	// Initialize containers to store various parameters
	// Current values of the parameters
	double *current = new double[numparams];
	// Proposed values of the parameters
	double *proposed = new double[numparams];
	// Array to hold the prior info
	double *priors = new double[3 * numparams];
	// Names of the parameters
	string *names = new string[numparams];
	// Sections that the parameters are in
	string *sections = new string[numparams];	
    // Whether or not parameters are logarithmic
    bool *logs = new bool[numparams];
	
	// Populate prior array from input prior struct
	for(int n = 0; n < numparams; n++){
		names[n] = spriors[n].name;
		sections[n] = spriors[n].section;
		logs[n] = spriors[n].logparam;
		// priors contains the lower, upper, and sigma values
		priors[n] = spriors[n].lower;
		priors[n + numparams] = spriors[n].upper;
		priors[n + 2 * numparams] = spriors[n].sigma;
	}

    // Load SN1a data
    vector<vector<double> > SN1adata;
    string sn1afile = inifile.getiniString("union21", "SCPUnion2.1_mu_vs_z.txt", "Function");
    if (loadSN1adata(sn1afile, SN1adata) != 0) {
        // Could not find data file
        cout << "Warning: cannot find Union2.1 SN1a data file." << endl;
    }

	double Omk, Omk_p;
	double Omk_min = inifile.getiniDouble("p2_min", 0.0, "Seq");
	double Omk_max = inifile.getiniDouble("p2_max", 0.0, "Seq");
	double Omk_sigma = (Omk_max - Omk_min) / 100.0;
	double HIRl_c, HIR_p, amax, w0, tmax, tarmfrac;
	double HIRl_t = inifile.getiniDouble("HIRt", 0.1, "Seq");
	int seq_find_step, seq_find_steps = 200;
	bool seqfound;
	bool findsequester = inifile.getiniBool("findsequester", false, "MCMC");
	
	if(findsequester)
		std::cout << "MCMC with the sequestering solution" << std::endl;
	else
		std::cout << "MCMC without sequestering constraint" << std::endl;
	
	vector<double> results;
	//******************//
    // Begin the chains //
    //******************//

	for(int chain = 1; chain < MCMCnumchains + 1; chain++){
		
        // Go and find an appropriate file name (using 4 digit numbers as the default)
        string outputname = getfilename(outputdir, basename, "", inifile.getiniInt("numberpad", 4, "Function"), false) + ".dat";
        // Open up file to dump chain info
        ofstream MCMCchainfile(outputname.c_str());
        // Report chain number to screen
        cout << "Chain #" << chain << " outputting to " << outputname << endl;
        // Make sure the file actually opened
        if (!MCMCchainfile.is_open()) {
            // Big problem!
            cout << "Error: could not open output file! Terminating." << endl;
            break;
        }

        // Zero the MCMC step number
        MCMCstep = 0;
        // Zero the MCMC acceptance counter
        MCMCaccept_counter = 0;

        // Print some information on the chain to file
        MCMCchainfile << "# Chain " << chain << " for model " << inifile.getiniString("Model", "LambdaCDM", "Cosmology")
                      << " using following datasets:" << endl;
        // Extract the bitmask and report the experiments being used
        int bitmaskcheck = inifile.getiniInt("chicombo", 29, "Function");
        if (bitmaskcheck < 1) bitmaskcheck = 29;
        if (bitmaskcheck > 127) bitmaskcheck = 29;
        unsigned int bitmask = bitmaskcheck;
        if ((bitmask & 1) > 0) {
            MCMCchainfile << "# Type 1a Supernovae" << endl;
            usingSN1a = true;
        } else {usingSN1a = false;}
        if ((bitmask & 2) > 0) MCMCchainfile << "# BAO data (SDSS)" << endl;
        if ((bitmask & 4) > 0) MCMCchainfile << "# BAO data (SDSSR)" << endl;
        if ((bitmask & 8) > 0) MCMCchainfile << "# Hubble prior" << endl;
        if ((bitmask & 16) > 0) MCMCchainfile << "# WMAP distance posteriors" << endl;
        if ((bitmask & 32) > 0) MCMCchainfile << "# Planck distance posteriors" << endl;
        // Print the priors
        MCMCchainfile << "# Priors:" << endl;
        for(int n = 0; n < numparams; n++) {
            MCMCchainfile << "#\t" << spriors[n].section << "\t";
            if (spriors[n].logparam) MCMCchainfile << "Log10 ";
            MCMCchainfile << spriors[n].name << "\t" << spriors[n].lower << "\t" << spriors[n].upper << "\t" << spriors[n].sigma << endl;
        }
        MCMCchainfile << "# ";
        for(int n = 0; n < numparams; n++)
            MCMCchainfile << spriors[n].name << "\t";
        MCMCchainfile << "Likelihood (unnormalised)" << endl;
        // Set up precision outputting
        MCMCchainfile << scientific << setprecision(15);

        // Start off at a random position in parameter space
        // But, make sure that it has a nonzero likelihood (otherwise it can meander almost indefinitely!)
        L_current = 0;
        while (L_current < 1e-200) {
            for(int param = 0; param < numparams; param++){
                lower = priors[param];
                upper = priors[param + numparams];
                current[param] = lower + UnitRand() * (upper - lower);
            }
            // Calculate the likelihood of the initial guess
            L_current = ComputeLikelihood(inifile, myOutput, names, sections, current, numparams, usingSN1a, SN1adata);
        }

		 
        // Reset progress counter
        barcount = 0;


		// Start off with a randomly chosen Omk
		// This will be refined to find the sequestering solution
		while(true && findsequester){
			Omk = Omk_min + UnitRand() * (Omk_max - Omk_min);
			if(Omk >= Omk_min || Omk <= Omk_max)
					break;
		}

		// Start the sampling
		while(true){
		
            // Increment our step counter
            MCMCstep++;

			/////////////////////////////////////////////////////////
			//// BEGIN ALGORITHM TO FIND SEQUESTERING SOLUTION
			//
			// Get some proposed parameters: need to dial Omk until sequestering solution is found
			// for rest of the values of the parameters.
			/////////////////////////////////////////////////////////
			seqfound = false;
			while(true){
				seq_find_step = 0;
				// (1) Get a set of parameters 
				GetProposedParameters(priors, current, proposed, logs, numparams);
				// (1.1) Set values of the parameters
				for(int n = 0; n < numparams; n++)	
					inifile.setparam(names[n], sections[n], proposed[n]);
			
				// (2) Dial Omega_k h^2 until a sequestering solution is found
				// (2.0.0) The starting value of Omk is inherited... this just makes finding
				//         a sequestering solution "easier".
				// (2.0.1) Start off the current value of <R> as something "quite large" 
				//		   ... the code will find something better than this
				HIRl_c = 1000;
				// (2.0.2) Set the evolver to run to Armageddon
				inifile.setparam("evtoend","Seq",true);
				while(true){
					// We only dial Omk to find sequestering solutions
					if(findsequester){
						// (2.1.1) Pick a new Omk that is inside the prior range
						while(true){
							Omk_p = Omk + NormalRand() * Omk_sigma;
							if(Omk_p <= Omk_max && Omk_p >= Omk_min)
								break;
						}
		
						// (2.1.2) Now that we have a sensible choice of Omk, run the evolver...
						// (2.1.2.1) Set the Omegakh2 value in inifile to be the choice found above
						inifile.setparam("Omegakh2","Cosmology",Omk_p);
					}
					// (2.1.2.2) Create a myParams_1 to be fed into the doEvolution routine
					Parameters myParams_1(inifile);
					// (2.1.2.3) Do the evolution
					results = doEvolution(inifile, myParams_1, myOutput, SN1adata, false);
					// ... and pull out the value of log10(|<R>|)
					HIR_p = log10(abs(results[1]));
		
					// (2.1.3) If <R> for this choice of Omk is smaller than the previous one, then keep it.
					if(HIR_p < HIRl_c){
						Omk = Omk_p;
						HIRl_c = HIR_p;
					}

					// (2.1.4) If <R> is smaller than the desired threshold "error", then we have found a sequestering solution.	
					if(HIRl_c < HIRl_t || !findsequester){
						// Remember some useful information about this run
						amax = results[4];
						w0 = results[5];
						tmax = results[6];
						tarmfrac = results[7];
						seqfound = true;
						break;
					}
					// If its taking too long to find a sequestering solution,
					// with this given set of parameters, then stop, and get a new set.
					if(seq_find_step > seq_find_steps && !seqfound)
						break;
					else
						seq_find_step++;
				}
				// If we've found a sequestering solution, then break out of this "find sequestering"
				// loop, and continue to compute the likelihood.
				if(seqfound)
					break;
				
			}
			
			// Make sure we only run the evolver to a = 1 (the current time)
			inifile.setparam("evtoend","Seq",false);
			/////////////////////////////////////////////////////////
			//// END ALGORITHM TO FIND SEQUESTERING SOLUTION
			/////////////////////////////////////////////////////////
						
			// Get the value of the likelihood with these proposed parameters		
			L_proposed = ComputeLikelihood(inifile, myOutput, names, sections, proposed, numparams, usingSN1a, SN1adata);
			// Compute likelihood ratio	
			LikelihoodRatio = L_proposed / L_current;
			// Decide whether to accept the proposed parameters.
			if( L_proposed >= L_current || UnitRand() < LikelihoodRatio ){
				// Increment acceptance counter
				MCMCaccept_counter++;
				// Store proposed parameters
				memcpy(current, proposed, numparams*sizeof(double));
				L_current = L_proposed;
			}
		
			// Dump to file after burn-in
			if(MCMCstep > MCMCburninsteps){
                for(int n = 0; n < numparams; n++) {
                    if (!spriors[n].logparam) {
                        MCMCchainfile << current[n] << "\t";
                    } else {
                        MCMCchainfile << log10(current[n]) << "\t";
                    }
                }
                MCMCchainfile << L_current;
				MCMCchainfile << "\t" << Omk << "\t" << amax << "\t" <<  tmax << "\t" << tarmfrac << "\t" << HIRl_c << endl;
			}
			else{
			    // Only start the acceptance counter after the burn-in period has ended
				MCMCaccept_counter = 1; // Set to one because we accept the first point we count
			}
		
            // Update the progress bar every 20 cycles
            if (showprogress && ++barcount >= 20) {
                barcount = 0;
                progress = (float) MCMCstep / (float) MCMCnumsteps;
                updateprogress(progress);
            }

			// Halt if we've taken enough steps
			if(MCMCstep >= MCMCnumsteps) break;
			
		} // END sampling while(){}
	
	    // Clear the progressbar
	    if (showprogress) {
	        updateprogress(1.0);
	        std::cout << std::endl;
	    }

	    // Write final counts for this chain
		MCMCchainfile << "# Number of samples (after burn-in):\t" << MCMCstep - MCMCburninsteps << endl
		               << "# Number of acceptances:\t" << MCMCaccept_counter << endl;
		MCMCchainfile.close();
		cout << setprecision(4) << "Chain complete. Acceptance rate = "
			 << 100 * MCMCaccept_counter / (float) (MCMCstep - MCMCburninsteps) << "%" << endl;

	} // END chain-loop


    //**********//
    // Clean up //
    //**********//

	delete current;
	delete proposed;
	delete priors;
	delete logs;

    // Stop timing
    myTimer.stop();
    double ms = myTimer.elapsed().wall / 1e6;
    if (ms < 1e3)
        cout << setprecision(4) << "Run complete in " << ms << " milliseconds.";
    else
        cout << setprecision(4) << "Run complete in " << ms / 1000.0 << " seconds.";
    cout << endl;

    // Exit gracefully
    return 0;

}

// Computes the likelihood of given parameters
double ComputeLikelihood(IniReader& inifile, Print2Memory& output, string *names, string *sections, double *parameters, int numparams, bool usingSN1a, vector<vector<double> > &SN1adata){
	
    // Make a clean slate for collecting data
    output.printfinish(0.0);

	// Set values of the parameters (note: need to get name & section for inifile)
	for(int n = 0; n < numparams; n++)	
		inifile.setparam(names[n], sections[n], parameters[n]);
	
	inifile.setparam("evtoend","Seq",false);
	
    // Set up the cosmological parameters 
    Parameters myParams(inifile);
	// Do the evolution, and return the likelihood for the data combination
	// defined by "combinationchi"
    double result;
	
	vector<double> results = doEvolution(inifile, myParams, output, SN1adata, true);
	
	// Now check if original run (which was to a = 1, today) is sensible
    if (results[0] == 0)
    if (usingSN1a) {
        // If the supernovae data is being used, then that has a base chi^2 of around 560
        // This translates to very small likelihoods, as e^{-500} is a very small number
        // In order to prevent underflow errors, we add an offset to the chi^2 (which doesn't change anything)
        result = exp( - 0.5 * (output.getvalue("combinationchi", 0.0) - 550.0));
    } else {
        result = exp( - 0.5 * output.getvalue("combinationchi", 0.0));
    }
	else
		result = 0;

    return result;
} // END ComputeLikelihood()
 

// Proposes a new point in parameterspace to jump to
void GetProposedParameters(double *priors, double *current, double *proposed, bool *logs, int numparams){
	
	// Get the upper and lower bounds, and sigma on the prior for this parameter
	double upper, lower, sigma;
	// Holding variable of the proposed parameter value
	double prop;
	// Holding variable of the current parameter value
	double thisval;
	// Is this a log parameter
	bool logparam;
	// Loop over all parameters
	for(int param = 0; param <  numparams; param++){
		
		thisval = current[param];
		lower = priors[param];
		upper = priors[param + numparams];
		sigma = priors[param + 2 * numparams];
		logparam = logs[param];

		// This process may choose parameter values which live outside the prior range,
		// and so must repeat until a parameter is found which is inside prior range.
		while( true ){
			
			// Jump around to a nearby parameter with normal distribution with the given sigma
			if (!logparam) {
			    prop = thisval + NormalRand() * sigma;
			} else {
			    prop = thisval * pow(10, NormalRand() * sigma);
			}
			
			// Make sure the proposal is inside the prior range before getting out
			if(prop <= upper && prop >= lower) break;
			
		}
		// Place the new proposal into the array
		proposed[param] = prop;
	}
	
} // END GetProposedParameters() 
