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
static int BeginEvolution(Integrator&, IntParams&, double*, const double, const double, Output&, Consistency&, vector<double>&, vector<double>&);

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
int doEvolution(IniReader& inifile, Parameters& params, Output& output, vector<vector<double> > &SN1adata, bool postprocess) {

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
        return -1;
    }

    // Check that the model has an internally consistent state with the initial data (it should have aborted already if so, but be safe)
    if (myModel->isvalidconfig(data) == false) {
        delete myChecker;
        delete myModel;
        return 4;
    }

    // Write the model name to the output log
    output.printvalue("Model", myModel->classname());
    output.printlog(""); // Whitespace for prettiness


    //***********//
    // Evolution //
    //***********//

    // Do the evolution!
    result = BeginEvolution(myIntegrator, myIntParams, data, starttime, endtime, output, *myChecker, hubble, redshift);


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

    // Return the result of the evolution
    return result;

}


//
// Helper routines
//

static int BeginEvolution(Integrator &integrator, IntParams &params, double data[],
        const double starttime, const double endtime, Output &output, Consistency &check,
        vector<double>& hubble, vector<double>& redshift) {
    // This routine takes in a number of parameters, and performs the cosmological background evolution
    // integrator is the class that handles integration steps
    // params is the class that stores the cosmological parameters
    // data is an array of four elements that stores a, \phi, \dot{\phi} and H. This is the data that gets evolved.
    // starttime stores the starting value of conformal time for the evolution
    // endtime stores a value of conformal time after which we should terminate evolution
    // output is a class that writes things to the screen and output files as appropriate
    // check is a class that checks the data for internal self-consistency (e.g., ghosts, superluminality, etc)
    // hubble and redshift are vectors used for storing this data for the purpose of postprocessing

    // We need our own time variable to step forwards
    double time = starttime;
    // The result from the integrator. It returns GSL_SUCCESS (0) when everything works
    int result;
    // An array to hold the status information
    double status[17];
    // And a double to hold the stepsize
    double stepsize;

    // Variables for storing the current state, in case we overshoot.
    double oldtime;
    double olddata[4];

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

    // This is the loop that takes successive integration steps. Halt if we go past the maximum evolution time.
    while (time < endtime) {

        // Before taking a timestep, check for a valid model configuration
        if (params.getmodel().isvalidconfig(data) == false) {
            // The model reports it is internally inconsistent. Abort.
            output.printlog("Model reports inconsistent state. Terminating.");
            output.printvalue("FatalError", 1);
            return 4; // invalid state error

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
            output.printlog("Integration routine failed.");
            output.printvalue("FatalError", 1);
            return 1; // integration error
        }
		
        // If we've shot past a = 1, then interpolate back to a = 1.
        if (data[0] > 1.0) {
            // Calculate the time at which a = 0
            double olda = olddata[0];
            double newa = data[0];
            double time1 = oldtime + (time - oldtime) * (1 - olda) / (newa - olda);

            // Calculate the value of H
            double oldH = olddata[3];
            double newH = data[3];
            double H1 = oldH + (time1 - oldtime) / (time - oldtime) * (newH - oldH);

            // Calculate the value of phi
            double oldphi = olddata[1];
            double newphi = data[1];
            double phi1 = oldphi + (time1 - oldtime) / (time - oldtime) * (newphi - oldphi);

            // Calculate the value of phidot
            double oldphid = olddata[2];
            double newphid = data[2];
            double phid1 = oldphid + (time1 - oldtime) / (time - oldtime) * (newphid - oldphid);

            // Put it all back into the data
            time = time1;
            data[0] = 1.0;
            data[1] = phi1;
            data[2] = phid1;
            data[3] = H1;

        }
		
        // Extract the state from the model
        params.getmodel().getstate(data, time, status, params.getparams());

        // Add the hubble and redshift values to the vectors
        hubble.push_back(status[3]);
        redshift.push_back(status[2]);

        // Get the output class to write out the state of the system
        output.printstep(data, time, status);

        // Take a look at the consistency of the data
        check.checkstate(data, time, params, output, status);

        // Make sure that nothing has become NaN
        if (check.checknan(data, time, status)) {
            // Something has become not-a-number
            output.printlog("A quantity has become NaN. Terminating.");
            output.printvalue("FatalError", 1);
            return 2; // NAN error
        }
		
        // Get out if we're done
        if (data[0] >= 1.0)
            break;
		
		 
		
        // If we're nearing a = 1, be careful about overshooting
        // Estimated step size in a is H * stepsize
        stepsize = integrator.getstepsize();
        if (data[0] + 2.0 * data[3] * stepsize > 1.0) {
            // Reduce the stepsize
            // Calculate the exact amount that the stepsize will need to be in order to get to a = 1
            // in the linear approximation
            double temp = (1.0 - data[0]) / 2.0 / data[3];
            integrator.setstepsize(0.9 * temp);
            // Eventually we'll run into the minimum step size and we'll cross the finishline
            // Also note that Hdot is usually negative, so it will typically take a little bit more than temp
            // to cross the finish line.
            // This tends to take about 20 steps to hit a = 1, which should be good enough to have
            // derivatives near a = 1 under control.
        }
		
				
    }

    // We reached here either because we got to the max time, or a = 1

    // Return the error if we didn't get to a = 1
    if (data[0] < 1.0) {
        output.printlog("Did not reach a=1 during expected evolution time.");
        output.printvalue("Incomplete", 1);
        return 3; // "Did not complete" error
    }

    // Take a look at the consistency of the final state
    check.checkfinal(data, time, params, output, status);

    // Return success!
    return 0;
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
