/*
 * process.h
 *
 * Describes the routines used for performing post-processing
 */

#ifndef PROCESS_H_
#define PROCESS_H_

#include "params.h"
#include "output.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_integration.h>
#include <boost/filesystem.hpp> // Used for making sure sn1a file exists
#include <iomanip>

using std::vector;

// We need a data type to pass a spline as well as the spline accelerator into the integration function
struct splinetools {
    gsl_interp_accel *acc;
    gsl_spline *spline;
    double param; // Also allow an extra parameter to go through to the integrator
};

// Routine that performs all postprocessing
int PostProcessingDist(vector<double>&, vector<double>&, vector<double>&, vector<double>&,
		vector<double>&, vector<double>&, vector<double>&, double &rs, double &rd, Parameters&, Output&);

// Routines that return the integrand for various integrals
int PPintfunc(double, const double*, double*, void*);
double rsintfunc(double z, void *params);
double rsintfuncinf(double z, void *params);

// Routine to compute chi^2 values for supernovae data
double chi2SN1a(vector<double>& redshift, vector<double>& mu, Output &output, IniReader &init, vector<vector<double> > &SN1adata);

// Routine to compute chi^2 values for CMB data
void chi2CMB(vector<double>& redshift, vector<double>& mu, double rs, Output &output, Parameters &params, double& WMAP, double& Planck);

// Routine to compute chi^2 values for BAO data
void chi2BAO(double rdrag, vector<double>& redshift, vector<double>& hubble, vector<double>& DA, Parameters &params, Output &output, double& BAO, double& BAOr);

// Routine to compute chi^2 value for Hubble
double chi2hubble(Parameters &params, double desired, double sigma, Output &output);

// Routines to compute chi^2 for WMAP and Planck distance posteriors
double chi2WMAP (double lA, double R, double z);
double chi2Planck (double lA, double R, double omegaB);

#endif /* PROCESS_H_ */
