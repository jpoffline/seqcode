/*
 * evolve.cpp
 *
 * This source file contains the important functions for this program.
 *
 * If the routines of this program are to be used in another program, these routines are the ones to call.
 *
 */
#include "evolve.h"

// Function that actually does the evolution
vector<double> BeginEvolution(Integrator&, IntParams&, double*, const double, const double, Output&, Consistency&, vector<double>&, vector<double>&, double, double,Parameters&);

// Function that the integrator calls to obtain derivatives
static int intfunc(double, const double*, double*, void*);


// This is the workhorse of the code.
// It takes in the input parameters and an output object, and runs everything from there.
// This modularization is set up so that anybody can call the evolution routines.
//
// Input parameters:
// inifile: the input parameters
// params: class containing cosmological parameters (somewhat degenerate with inifile, but these values won't be read from inifile)
// output: the outputting class
// postprocess: whether or not to perform postprocessing
//
// Return values:
// 0: success
// -1: error in initialization
// 1: Integration error
// 2: NAN error
// 3: Did not get to a = 1 in allotted time
// 4: Model reports invalid state
// 5: Error in postprocessing integration of distance measures
vector<double> doEvolution(IniReader& inifile, Parameters& params, Output& output, vector<vector<double> > &SN1adata, bool postprocess) {

    //****************//
    // Initialization //
    //****************//

    // Set up the integrator
    Integrator myIntegrator;

    // Set up the model class
    Model *myModel;
    std::string parsestring = inifile.getiniString("model", "LambdaCDM", "Cosmology");
    if (parsestring == "Quintessence")
        myModel = new Quintessence();
    else if (parsestring == "LinearW")
        myModel = new LinearW();
    else if (parsestring == "Kessence")
        myModel = new Kessence();
    else if (parsestring == "KGB")
        myModel = new KGB();
    else if (parsestring == "Fr")
        myModel = new Fr();
    else
        myModel = new LambdaCDM();    // LambdaCDM is the default

    // Load the model and parameters into a class to pass into the integration routine
    IntParams myIntParams(params, *myModel);

    // Set up the consistency check class
    parsestring = inifile.getiniString("consistencyclass", "None", "Function");
    Consistency *myChecker;
    if (parsestring == "SimpleCheck")
        myChecker = new SimpleCheck();
    else
        myChecker = new Consistency();  // Default option, which has no checking

    // Create vectors for redshift and hubble
    vector<double> redshift;
    vector<double> hubble;

    int result = 0; // For information coming back from functions
	
	
    //********************//
    // Initial Conditions //
    //********************//

    // The initial value of a is extracted from the starting redshift
    // The initial values of \phi and \dot{\phi} are read from the input parameters
    double phi0 = inifile.getiniDouble("phi0", 0.0, "Cosmology");
    double phidot0 = inifile.getiniDouble("phidot0", 0.0, "Cosmology");
    double data[4] = { 1.0 / (1.0 + params.z0()), phi0, phidot0, 0.0 };
    // The data array stores a, \phi, \dot{phi} and H through the evolution
    // H is calculated in the initialization of the model

    // Start and end times
    double starttime = inifile.getiniDouble("starttime", 0.0, "Function");
    double endtime = starttime + inifile.getiniDouble("maxtime", 10.0, "Function");

    // Get the output class to write out information on the run
    output.printinfo(data, params);

    // Allow the model to initialize itself
    result = myModel->init(data, starttime, params, inifile, output);
    if (result != 0){
        delete myChecker;
        delete myModel;
        result = -1;
    }

    // Check that the model has an internally consistent state with the initial data (it should have aborted already if so, but be safe)
    if (myModel->isvalidconfig(data) == false) {
        delete myChecker;
        delete myModel;
        result = 4;
    }

    // Write the model name to the output log
    output.printvalue("Model", myModel->classname());
    output.printlog(""); // Whitespace for prettiness


    //***********//
    // Evolution //
    //***********//
	
	
	double RT = inifile.getiniDouble("RicciThresh", 0.0, "Seq");
	double sma = inifile.getiniDouble("smallesta", 0.0, "Seq");

    // Do the evolution!
	// Vector to hold returned values
	// entry 1 = "result", entry 2 = "histintRicci"
	vector<double> returns = BeginEvolution(myIntegrator, myIntParams, data, starttime, endtime, output, *myChecker, hubble, redshift, RT, sma,params);
	
	result = int(returns[0]);

    //*********************************//
    // Clean up of Hubble and Redshift //
    //*********************************//

    if (result == 0) { // when the evolution was a success

        // We have H and z starting with high z going to z = 0. We want these reversed.
        reverse(hubble.begin(),hubble.end());
        reverse(redshift.begin(),redshift.end());

        // Update the parameters class with the new hubble value
        params.seth(hubble[0]);

        // Now go and report the actual density fractions and hubble parameter to the logs
        output.printlog("Values from evolution of model are as follows:");
        output.printvalue("modelh", params.geth());
        output.printvalue("modelOmegaR", params.getOmegaR());
        output.printvalue("modelOmegaM", params.getOmegaM());
        output.printvalue("modelOmegaB", params.getOmegaB());
        output.printvalue("modelOmegaK", params.getOmegaK());
        output.printvalue("modelOmegaLambda", 1 - params.getOmegaK() - params.getOmegaR() - params.getOmegaM());
        output.printlog("");

    }


    //*********************************//
    // Post-processing: initialization //
    //*********************************//

    // Construct vectors for distance measures
    int numrows = hubble.size();
    vector<double> DC;
    vector<double> DM;
    vector<double> DA;
    vector<double> DL;
    vector<double> mu;
    // Other distances that are computed
    double rs; // sound horizon at z_CMB
    double rd; // sound horizon at z_drag

    //*******************//
    // Distance measures //
    //*******************//

    if (result == 0 && postprocess) { // when the evolution was a success AND we're doing the postprocessing

        // Get the output class to print any headings for postprocessing
        output.postprintheading();

        // Allocate space for the vectors
        DC.reserve(numrows);
        DM.reserve(numrows);
        DA.reserve(numrows);
        DL.reserve(numrows);
        mu.reserve(numrows);

        // Calculate distance measurements
        result = PostProcessingDist(hubble, redshift, DC, DM, DA, DL, mu, rs, rd, params, output);

        // Flag any errors
        if (result == 1) result = 5;

    }


    //**************//
    // chi^2 values //
    //**************//

    if (result == 0 && postprocess) { // when the distance calculations were a success AND we're doing the postprocessing

        double BAOchi, BAOrchi, WMAPchi, Planckchi, Hubblechi, SNchi;

        // Put in some whitespace!
        output.printlog("");

        // First, do chi^2 of WMAP and Planck distance posteriors
        chi2CMB(redshift, DA, rs, output, params, WMAPchi, Planckchi);

        // Next, do chi^2 of SN1a
        SNchi = chi2SN1a(redshift, mu, output, inifile, SN1adata);

        // Do chi^2 for hubble value
        Hubblechi = chi2hubble(params, inifile.getiniDouble("desiredh", 0.7, "Cosmology"), inifile.getiniDouble("sigmah", 0.03, "Cosmology"), output);

        // Finally, do chi^2 of BAO measurements
        chi2BAO(rd, redshift, hubble, DA, params, output, BAOchi, BAOrchi);

        // Calculate the combination chi^2 from the bitmask
        double combination = 0;
        unsigned int bitmask = params.getbitmask();
        // 1 = SN1a
        // 2 = BAO (using SDSS)
        // 4 = BAO (using SDSSR)
        // 8 = Hubble
        // 16 = WMAP
        // 32 = Planck
        combination += (bitmask & 1) * SNchi;
        combination += ((bitmask & 2) >> 1) * BAOchi;
        combination += ((bitmask & 4) >> 2) * BAOrchi;
        combination += ((bitmask & 8) >> 3) * Hubblechi;
        combination += ((bitmask & 16) >> 4) * WMAPchi;
        combination += ((bitmask & 32) >> 5) * Planckchi;
        output.printvalue("combinationchi", combination);
        output.printvalue("chicombo", (int) bitmask);

    }
	


    //*********//
    // Cleanup //
    //*********//

    // All these classes have now done their job, and can go away
    delete myChecker;
    delete myModel;

    // Return the result of the evolution (a vector called "returns")
    return returns;

}


//
// Helper routines
//

vector<double> BeginEvolution(Integrator &integrator, IntParams &params, double data[],
        const double starttime, const double endtime, Output &output, Consistency &check,
        vector<double>& hubble, vector<double>& redshift, double RicciThreshold, double smallesta, Parameters& modp) {
    // This routine takes in a number of parameters, and performs the cosmological background evolution
    // integrator is the class that handles integration steps
    // params is the class that stores the cosmological parameters
    // data is an array of four elements that stores a, \phi, \dot{\phi} and H. This is the data that gets evolved.
    // starttime stores the starting value of conformal time for the evolution
    // endtime stores a value of conformal time after which we should terminate evolution
    // output is a class that writes things to the screen and output files as appropriate
    // check is a class that checks the data for internal self-consistency (e.g., ghosts, superluminality, etc)
    // hubble and redshift are vectors used for storing this data for the purpose of postprocessing

	vector<double> returns;
	int nreturns = 10;
	for(int nr = 0; nr < nreturns; nr++)
		returns.push_back(0.0);
		
    // We need our own time variable to step forwards
    double time = starttime;
    // The result from the integrator. It returns GSL_SUCCESS (0) when everything works
    int result;
    // An array to hold the status information
    double status[18];
    // And a double to hold the stepsize
    double stepsize;

    // Variables for storing the current state, in case we overshoot.
    double oldtime;
    double olddata[4];

	// BEGIN :: JAP: sequestering variables
	double a, RicciScalar;
	double Seq_vol = 0.0; // zero the measure
	double Seq_num = 0.0; // zero the numerator integrand
	double Seq_den = 0.0; // zero the denominator integrand
	double amax = 0.0;
	double wnow = 0.0;
	// END :: JAP: sequestering variables

	
    // Write our any header information to the logs
    output.printheading();
    // Extract the initial state from the model
    params.getmodel().getstate(data, time, status, params.getparams());
    // Add the starting hubble and redshift values to the vectors
    hubble.push_back(status[3]);
    redshift.push_back(status[2]);
    // Write the initial conditions to the log
	output.printstep(data, time, status);
    // Do a consistency check on the initial conditions
    check.checkstate(data, time, params, output, status);

	bool inFuture = false;
	bool akill = false;
	if( smallesta != 1 )
		akill = true;
    // This is the loop that takes successive integration steps. Halt if we go past the maximum evolution time.
	int aflag = 0;	
	
    while (time < endtime) {

        // Before taking a timestep, check for a valid model configuration
        if (params.getmodel().isvalidconfig(data) == false) {
            // The model reports it is internally inconsistent. Abort.
            output.printlog("Model reports inconsistent state. Terminating.");
            output.printvalue("FatalError", 1);
            returns[0] = 4; // invalid state error

        }

        // Store old data before taking a new step
        oldtime = time;
        olddata[0] = data[0];
        olddata[1] = data[1];
        olddata[2] = data[2];
        olddata[3] = data[3];

        // Take a step
        result = integrator.dointstep(intfunc, params, data, time, endtime);

        // If the step failed, return an integration error
        if (result != GSL_SUCCESS) {
           // output.printlog("Integration routine failed.");
            //output.printvalue("FatalError", 1);
            returns[0] = 1; // integration error
        }
		
        // Extract the state from the model
        params.getmodel().getstate(data, time, status, params.getparams());

        // Add the hubble and redshift values to the vectors
        hubble.push_back(status[3]);
        redshift.push_back(status[2]);
		
		// Get the scale factor
		a = status[1];
		
		// Compute Ricci Scalar in conformal time,
		RicciScalar = 6.0 * ( status[4] + status[3] * status[3] - modp.rhoK() * modp.getconH0() * modp.getconH0() ) / a / a;

		// Compute measure in conformal time,
		Seq_vol = a * a * a * a;
		
		// Dump the Ricci scalar into the last slot of "status" -- mainly to be printed out
		status[17] = RicciScalar*Seq_vol;
		
        // Get the output class to write out the state of the system
		output.printstep(data, time, status);

        // Take a look at the consistency of the data
        check.checkstate(data, time, params, output, status);

        // Make sure that nothing has become NaN
        if (check.checknan(data, time, status)) {
            // Something has become not-a-number
			output.printlog("A quantity has become NaN. Terminating.");
			output.printvalue("FatalError", 1);
            returns[0] = 2; // NAN error
        }
		
       
		/*
		// KEEP THIS COMMENTED OUT!!!
		// JAP
		This stuff will help for steeply varying a near the crunch...
	    stepsize = integrator.getstepsize();
        if (data[0] + 4.0 * data[3] * stepsize > 20) {
            // Reduce the stepsize
            // Calculate the exact amount that the stepsize will need to be in order to get to a = 0
            // in the linear approximation
            double temp = (20 - data[0]) / 4.0 / data[3];
            integrator.setstepsize(0.001 * temp);
			std::cout << "refining step-size on the upwards stretch before crunch" << std::endl;
        }
		*/
		
		if( a > amax)
			amax = a;
		
		// Refine step-size near crunch singularity
        stepsize = integrator.getstepsize();
        if (data[0] + 20.0 * data[3] * stepsize < 0 ) {
            // Reduce the stepsize
            // Calculate the exact amount that the stepsize will need to be in order to get to a = 0
            // in the linear approximation
            double temp =  - data[0] / 20.0 / data[3];
            integrator.setstepsize( 0.00000000001 * temp );
        }
		
		// Check to find out if we are "in the future"
		if( a > 1 )
			inFuture = true;
		if(!inFuture)
			wnow = status[15];
		
		// When nearing the crunch, start to refine the stepsize
		if(data[0] < 0.08 && inFuture)
			integrator.setstepsize( 0.0000000001 * integrator.getstepsize() );
		
		if(data[0] < 0.01 && inFuture)
			integrator.setstepsize( 0.0000000000001 * integrator.getstepsize() );
				
		
		// Compute the numerator-term of the constraint for sequestering
		Seq_num += RicciScalar * Seq_vol;
		
		// Compute the denominator-term of the constraint for sequestering
		Seq_den += Seq_vol;
			
		// Reasons to stop...

		if( a < smallesta && inFuture && akill ){
			aflag = 1;
			break;
		}
				
		// If we are about to drop below a < 0, kill it (else serious issue ensue...)
		//if( data[0] + 0.01 * data[3] * stepsize < 0 )
		//	break;
		
		
		// \JAP
				
		if(  (RicciScalar > 1E10 || RicciScalar < -1E10) && inFuture ){
			break;
		}
		
    }
		
	// Construct value of the historic integral of the Ricci Scalar, 
	// to be sent back to calling routine.	
	double histintRicci = Seq_num / Seq_den;

    // Take a look at the consistency of the final state
    check.checkfinal(data, time, params, output, status);

    // Return success!
	if (returns[0] == 0)
		returns[0] = 0;
	
	returns[1] = histintRicci; // <R>
	returns[2] = RicciScalar; // R_end
	returns[3] = data[0]; // a_end
	returns[4] = amax;	// a_max
	returns[5] = wnow;
	returns[6] = time;
	return returns;
	
}

static int intfunc(double t, const double data[], double derivs[], void *params) {
    // This routine calculates the derivatives for a, phi, \dot{\phi} and \dot{H}
    // It is called by the integration routine.

    // Extract parameters
    IntParams myParams = *(IntParams *) params;

    // Call the derivatives routine in the model to calculate the derivatives appropriately
    // Note that they don't depend on time
    return myParams.getmodel().derivatives(data, derivs, myParams.getparams());

}
