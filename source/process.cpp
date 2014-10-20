/*
 * process.cpp
 *
 * Contains the routines to perform post-processing of data after the fact.
 *
 * At the moment, just computes the various distance measures used in cosmology.
 *
 */

#include "process.h"

using namespace std;

static inline double getzCMB(Parameters &params) {
	// Calculate redshift of recombination (CMB formation)
	// Computed using Eq E.1 from astro-ph/9510117v2 (note that Omega_0 = Omega_m) for percent level accuracy
	double OmegaBh2 = params.rhoB(); // Note that this is OmegaB * h^2
	double g1 = 0.0783 * pow( OmegaBh2, -0.238) / (1.0 + 39.5 * pow(OmegaBh2, 0.763));
	double g2 = 0.560 / (1.0 + 21.1 * pow(OmegaBh2, 1.81));
	return 1048 * (1.0 + 0.00124 * pow(OmegaBh2, -0.738)) * (1.0 + g1 * pow(params.rhoM(), g2));
}

static inline double getzdrag(Parameters &params) {
	// Calculate redshift of drag epoch
	// Computed using fitting forms from astro-ph/9709112 (note that Omega_0 = Omega_m) for percent level accuracy
    double OmegaBh2 = params.rhoB(); // Note that this is OmegaB * h^2
	double OmegaMh2 = params.rhoM(); // Note that this is OmegaM * h^2
	double b1 = 0.313 * pow(OmegaMh2, -0.419) * (1 + 0.607 * pow(OmegaMh2, 0.674));
	double b2 = 0.238 * pow(OmegaMh2, 0.223);
	return 1291 * pow(OmegaMh2, 0.251) * (1 + b1 * pow(OmegaBh2, b2)) / (1 + 0.659 * pow(OmegaMh2, 0.828));
}

int PostProcessingDist(vector<double>& hubble,
		vector<double>& redshift,
		vector<double>& DC,
		vector<double>& DM,
		vector<double>& DA,
		vector<double>& DL,
		vector<double>& mu,
		double &rs,
		double &rd,
		Parameters &params, Output &output) {
	// Computes distance measures

	// Takes in a vector of hubble values and z values (in ascending z order)
	// The rest of the vectors are assumed to be empty
	// rs will be given the sound horizon scale at z_CMB
	// rd will be given the sound horizon scale at z_drag
	// Also takes in some basic parameters about the cosmology
	// as well as the output class
	// Note that all parameters have now been scaled to the model, rather than to the input values

	// Store the number of rows
	int numrows = hubble.size();

	// Step 1: Construct an interpolation of H and z
	// Trick: Get a pointer to the array for H and z
	double* pH = &hubble[0];
	double* pz = &redshift[0];

	// Initialize the spline
	splinetools myspline;
	myspline.acc = gsl_interp_accel_alloc();
	myspline.spline = gsl_spline_alloc(gsl_interp_cspline, numrows);

	// Construct the spline
	gsl_spline_init (myspline.spline, pz, pH, numrows);

	// To get the value from the spline, use the following command
	// H(z) = gsl_spline_eval (myspline.spline, z, myspline.acc);

	// Step 2: Integrate
	// Now that we have our spline, it's time to integrate to obtain the appropriate quantities
	// We compute the quantity D_C(z) = c \int_0^z dz'/H(z')/(1 + z') where this is the full H (with units)
    // We begin by computing \int_0^z dz'/\tilde{H}(z')/(1 + z') where \tilde{H} = H / (100 km/s/Mpc)

	// Set up the integration stuff
	gsl_odeiv2_step *step;
	gsl_odeiv2_control *control;
	gsl_odeiv2_evolve *evolve;

	// Integration method as supplied by GSL
	const gsl_odeiv2_step_type *Type = gsl_odeiv2_step_rk8pd;

	// This is the number of elements that are being integrated
	size_t numelements = 1;

	// Initialize GSL
	step = gsl_odeiv2_step_alloc(Type, numelements);
	control = gsl_odeiv2_control_yp_new(0.0, 1e-12); // Only care about relative error (this is quite a stringent tolerance)
	evolve = gsl_odeiv2_evolve_alloc(numelements);

	// This is the stepsize that has been recommended by the integrator
	double stepsize;

	// Set the initial stepsize to be 10% of the first quantity
	stepsize = redshift[1] / 10.0;

	// Sets up the system to integrate, including the function that specifies the derivatives, NULL for the Jacobian, the number of elements
	// being integrated, and the parameters to pass through to the function
	gsl_odeiv2_system sys = { PPintfunc, NULL, numelements, &myspline };

	// We need to store the current value of the function we're integrating as an array
	double data[1] = {0.0};
	double currentz = 0.0;
	double zfinish;

	int status = GSL_SUCCESS;
	// We want to integrate from z = 0 to each z value in the data
	for (int i = 0; i < numrows; ++i) {
		// Read in the next z value to get to
		zfinish = redshift[i];

		// Integrate forwards in z until we get to the point we want
		while (currentz < zfinish) {

		    // Limit the step size in order to improve accuracy
		    if ((redshift[i] - redshift[i-1]) / 10.0 > stepsize)
		        stepsize = (redshift[i] - redshift[i-1]) / 10.0;

			// Do the integration
			status = gsl_odeiv2_evolve_apply(evolve, control, step, &sys, &currentz, zfinish, &stepsize, data);

			// Check for problems
			if (status != GSL_SUCCESS)
				break;

		}

		// Check for problems
		if (status != GSL_SUCCESS)
			break;

		// We're up to the next value of z; add it to the vector
		DC.push_back(data[0]);

	}

	// Release the integrator memory
	gsl_odeiv2_evolve_free (evolve);
	gsl_odeiv2_control_free (control);
	gsl_odeiv2_step_free (step);

	// Get out with an error response if necessary
	if (status != GSL_SUCCESS) {
		// Error message
		std::stringstream erroroutput;
		erroroutput << "Error in integrating distance measures: " << status << std::endl;
		output.printlog(erroroutput.str());
		output.printvalue("DistanceError", 1);

		// Release the spline memory
		gsl_spline_free (myspline.spline);
		gsl_interp_accel_free (myspline.acc);
		return 1;
    }

	// We have the quantity \int_0^z dz'/\tilde{H}(z')/(1 + z') where \tilde{H} = H / (100 km/s/Mpc)
	// To obtain D_C, we need to multiply this by c/H0
	// Store this quantity
	double conH0 = params.getconH0(); // in Mpc
	// Multiply all values by this in order to get results in Mpc
    for (int i = 0; i < numrows; i++) {
        DC[i] *= conH0;
    }


	// Step 3: Compute the other distance measures
	// Begin by calculating DH
    double DH = params.getDH();

	// Extract OmegaK from the parameters
	double OmegaK = params.getOmegaK();

	// The computation of DM depends on whether OmegaK is 0, positive or negative
	if (OmegaK > 0) {
		// k < 0
		double rootk = sqrt(OmegaK);
		for (int i = 0; i < numrows; i++) {
			DM.push_back(DH * sinh(rootk * DC[i] / DH) / rootk);
		}
	}
	else if (OmegaK < 0) {
		// k > 0
		double rootk = sqrt(-OmegaK);
		for (int i = 0; i < numrows; i++) {
			DM.push_back(DH * sin(rootk * DC[i] / DH) / rootk);
		}
	}
	else {
		// OmegaK == 0
		DM = DC;
	}

	// Next, compute DA and DL
	DA = DM;
	DL = DM;
	for (int i = 0; i < numrows; i++) {
		currentz = 1 + redshift[i];
		DA[i] /= currentz;
		DL[i] *= currentz;
	}

	// Finally, compute mu. As DL is already in Mpc, this is straightforward.
	for (int i = 0; i < numrows; i++) {
		mu.push_back(25 + 5 * log10 (DL[i]));
	}


    // Step 4: Output all data.
	// At z = 0, all distance measures are zero, H = 1, and mu = -infinity.
	// This is a problem for mu, which becomes infinite.
	// For this reason, we remove the first entry of each distance measurement when reporting data.
	for (int i = 1; i < numrows; i++) {
		output.postprintstep(redshift[i], hubble[i], DC[i], DM[i], DA[i], DL[i], mu[i]);
	}


	// Step 5: Calculate the sound horizon distance at recombination and drag epochs
	// This one is a little more complicated.
	// The expression we need is the following.
    // r_s = c \frac{1}{\sqrt{3}} \int^{\infty}_{z_{CMB}} \frac{dz'}{(1 + z') H(z')} \frac{1}{\sqrt{1 + 3 \Omega_B / 4 \Omega_\gamma (1 + z')}}
	//     = c / H_0 \frac{1}{\sqrt{3}} \int^{\infty}_{z_{CMB}} \frac{dz'}{(1 + z') \tilde{\ch}(z')} \frac{1}{\sqrt{1 + 3 \Omega_B / 4 \Omega_\gamma (1 + z')}}
	// Start by computing the quantity that depends on OmegaB and OmegaR, and storing it in the integration helper
	double alpha = 3 * params.rhoB() / 4 / params.rhoGamma();
	myspline.param = alpha;
	// Note that this integral goes from zCMB to infinity. Our evolution doesn't go to infinite redshift however.
	// We need to evaluate the integral over our span of redshift, as well as over the infinite portion.
	// As we don't have information outside our span, we use a LambdaCDM prediction for that period, which should be accurate for most purposes.

	// As we are only computing one value here, rather than a value at various redshifts, we can use a specialized integration routine.
	// Set up the integration workspace
	gsl_integration_workspace * workspace = gsl_integration_workspace_alloc (1000);

	// First, we integrate over the data that we have
	double intresult, error;
	gsl_function intFunc;
	intFunc.function = &rsintfunc;
	intFunc.params = &myspline;

	// Perform the integration: z_CMB
	gsl_integration_qag (&intFunc, getzCMB(params), params.z0(), 1e-8, 0, 1000, 6, workspace, &intresult, &error);
	// Store the result
	rs = intresult;

	// Perform the integration: z_drag
	gsl_integration_qag (&intFunc, getzdrag(params), params.z0(), 1e-8, 0, 1000, 6, workspace, &intresult, &error);
	// Store the result
	rd = intresult;

	// Next, we integrate over the infinite part of the integral. Note that we need to pass in all the parameters to calculate this integrand.
	intFunc.function = &rsintfuncinf;
	intFunc.params = &params;

	// Perform the integration in two steps. Once, to redshift 1e5, using the more accurate integrator. Then, to infinite redshift.
	// First, finite part
	gsl_integration_qag (&intFunc, params.z0(), 1e5, 1e-8, 0, 1000, 6, workspace, &intresult, &error);
	// Store the result
	rs += intresult;
	rd += intresult;

	// Infinite part.
	gsl_integration_qagiu (&intFunc, 1e5, 1e-8, 0, 1000, workspace, &intresult, &error);
	// Store the result
	rs += intresult;
	rd += intresult;

	// Include the multiplicative factor of 1/sqrt(3) and the factor of c / H0
	const double rootm3 = 1 / sqrt(3.0);
	rs *= rootm3 * conH0;
	rd *= rootm3 * conH0;

	// Release the memory from the integrator
	gsl_integration_workspace_free (workspace);

	// Release the spline memory
	gsl_spline_free (myspline.spline);
	gsl_interp_accel_free (myspline.acc);

	// Report the result
    output.printvalue("SoundScaleCMB", rs); // Sound horizon scale in Mpc at z_CMB
    output.printvalue("SoundScaleDrag", rd); // Sound horizon scale in Mpc at z_drag

	// Success!
	return 0;
}

int PPintfunc(double z, const double data[], double derivs[], void *params) {
	// This routine returns the derivative of the function we wish to integrate.
	// In this case, we're integrating \int_0^z dz'/H(z')/(1 + z')
	// Thus, the derivative is 1/H(z')/(1 + z')

	// Extract parameters
	splinetools myParams = *(splinetools *) params;

	// Calculate H(z)
	double H = gsl_spline_eval (myParams.spline, z, myParams.acc);

	// Calculate the derivative
	derivs[0] = 1.0 / (1.0 + z) / H;

	// Return success!
	return GSL_SUCCESS;

}

double rsintfunc(double z, void *params) {
	// This routine returns the derivative of the function we wish to integrate.
	// In this case, we're integrating \int dz / (1 + z') H(z') sqrt{1 + alpha / (1 + z')}
	// Thus, the derivative is 1 / (1 + z') H(z') sqrt{1 + alpha / (1 + z')}

	// Extract parameters
	splinetools myParams = *(splinetools *) params;

	// Calculate H(z)
	double H = gsl_spline_eval (myParams.spline, z, myParams.acc);

	// Calculate the derivative
	return 1 / (1.0 + z) / H / sqrt(1 + myParams.param / (1.0 + z));

}

double rsintfuncinf(double z, void *params) {
	// This routine returns the derivative of the function we wish to integrate.
	// In this case, we're integrating \int dz / (1 + z') H(z') sqrt{1 + alpha / (1 + z')}
	// Thus, the derivative is 1 / (1 + z') H(z') sqrt{1 + alpha / (1 + z')}
	// However, here we do not have the benefit of a spline

	// Extract parameters
    Parameters myParams = *(Parameters *) params;

	// Calculate H(z) based on LambdaCDM FRW evolution with radiation, matter and curvature
	double a = 1.0 + z;
	double H = sqrt(a * a * myParams.rhoR() + a * myParams.rhoM() + myParams.rhoK());
	double alpha = 3 * myParams.rhoB() / 4 / myParams.rhoR();

	// Calculate the derivative
	return 1 / a / H / sqrt(1.0 + alpha / a);

}

double chi2SN1a(vector<double>& redshift, vector<double>& mu, Output &output, IniReader &init, vector<vector<double> > &SN1adata) {
	// This routine computes the chi^2 value for supernovae measurements
    // Returns the chi^2 value, as well as sending it to the output class

	// Step 1: We will need to know mu at specific redshifts. Construct the required interpolater.

	// Trick: Get pointers to the various arrays
	// At z = 0, mu = -infinity.
	// For this reason, we remove the first entry when constructing the interpolation.
	double* pz = &redshift[1];
	double* pmu = &mu[1];
	int numrows = redshift.size() - 1;
    double chi2 = 0;

	// Initialize the splines
	// Distance modulus is used for the SN1a data
	splinetools muspline;
	muspline.acc = gsl_interp_accel_alloc();
	muspline.spline = gsl_spline_alloc(gsl_interp_cspline, numrows);

	// Construct the spline
	gsl_spline_init (muspline.spline, pz, pmu, numrows);

	// Step 2: Read in all of the data from the SN1a data

    // SN1adata contains a vector of doubles. Each row is the data from a supernovae.
    // Each row contains four pieces of information: redshift, distance modulus, error on distance modulus, and a piece of
    // information that is not of use to us.

    // Iterate through each row, and construct the chi^2
    double val = 0;
    int numsn = SN1adata.size();
    // If the file failed to load, then numsn = 0.
    for (int i = 0; i < numsn; i++) {
        // Add (mu - mu(z))^2 / sigma^2 to the chi^2
        val = (SN1adata[i][1] - gsl_spline_eval (muspline.spline, SN1adata[i][0], muspline.acc)) / SN1adata[i][2];
        chi2 += val * val;
    }

    // Having gotten here, present the results
    output.printvalue("SNchi", chi2);
    // The minimum chi^2 for LambdaCDM for this data set is ~562

	// Release the spline memory
	gsl_spline_free (muspline.spline);
	gsl_interp_accel_free (muspline.acc);

	return chi2;

}

void chi2CMB(vector<double>& redshift, vector<double>& DA, double rs, Output &output, Parameters &params, double &WMAP, double&Planck) {
	// This routine computes the chi^2 value for CMB distance posteriors
    // Returns the chi^2 values in WMAP and Planck

	// We use the WMAP distance posteriors on the acoustic scale and the shift parameter
	// See http://lambda.gsfc.nasa.gov/product/map/dr3/pub_papers/fiveyear/cosmology/wmap_5yr_cosmo_reprint.pdf
	// This needs the angular diameter distance at recombination, so we'll need a new interpolater
	// Extract the appropriate data
	double* pz = &redshift[0];
	double* pDA = &DA[0];
	int numrows = redshift.size();

	// Construct the spline
	splinetools DAspline;
	DAspline.acc = gsl_interp_accel_alloc();
	DAspline.spline = gsl_spline_alloc(gsl_interp_cspline, numrows);
	gsl_spline_init (DAspline.spline, pz, pDA, numrows);

	// Get the redshift of the CMB and drag epochs
	double zCMB = getzCMB(params);
	double zdrag = getzdrag(params);

	// Calculate l_A
	double la = (1.0 + zCMB) * gsl_spline_eval (DAspline.spline, zCMB, DAspline.acc) * M_PI / rs;

	// Calculate R
	double R = sqrt(params.getOmegaM()) * (1.0 + zCMB) * gsl_spline_eval (DAspline.spline, zCMB, DAspline.acc) / params.getDH();

	// Release the spline memory
	gsl_spline_free (DAspline.spline);
	gsl_interp_accel_free (DAspline.acc);

	// Finally, compute the chi^2 values for WMAP
	double chi2wmap = chi2WMAP(la, R, zCMB);

	// For Planck, we use the results from 1304.4514v2. Note that this uses Omega_B h^2 instead of zCMB
	double chi2p = chi2Planck(la, R, params.rhoB());

	// Output everything
    output.printvalue("AcousticScale", la);
    output.printvalue("ShiftR", R);
    output.printvalue("zCMB", zCMB);
    output.printvalue("zdrag", zdrag);
    output.printlog(""); // Whitespace!
    output.printvalue("WMAPchi", chi2wmap);
    output.printvalue("PLANCKchi", chi2p);

    WMAP = chi2wmap;
    Planck = chi2p;

}

// Computes the chi^2 value for the CMB distance posteriors from WMAP 9 data
// See Table 11 in arxiv.org/pdf/1212.5226v3.pdf and numbers below
double chi2WMAP (double lA, double R, double z) {
	// Construct deltas
	double deltal = lA - 302.4;
	double deltaR = R - 1.7246;
	double deltaz = z - 1090.88;

	double chi2 = 0;
	// Add contributions from diagonals first
	chi2 += deltal * deltal * 3.182;
	chi2 += deltaR * deltaR * 11887.879;
	chi2 += deltaz * deltaz * 4.556;
	// Add contributions from off-diagonal terms next
	chi2 += 2 * deltal * deltaR * 18.253;
	chi2 += - 2 * deltal * deltaz * 1.429;
	chi2 += - 2 * deltaR * deltaz * 193.808;

	// Return the result
	return chi2;
}

// Computes the chi^2 value for the CMB distance posteriors from Planck 1 data using 1304.4514v2 results
double chi2Planck (double lA, double R, double omegaB) {
	// Construct deltas
	double deltal = lA - 301.57;
	double deltaR = R - 1.7407;
	double deltaB = omegaB - 0.02228;
	// Divide the deltas by the standard deviations (we're given a correlation matrix instead of a covariance matrix)
	deltal /= 0.18;
	deltaR /= 0.0094;
	deltaB /= 0.00030;

	double chi2 = 0;
	// We drop the n_s row/col from the paper's matrix in Eq 13, and invert the matrix to obtain the inverse covariance matrix.
	// Add contributions from diagonals first
	chi2 += deltal * deltal * 1.39378;
	chi2 += deltaR * deltaR * 2.19775;
	chi2 += deltaB * deltaB * 1.93992;
	// Add contributions from off-diagonal terms next
	chi2 += - 2 * deltal * deltaR * 0.620578;
	chi2 +=   2 * deltal * deltaB * 0.160517;
	chi2 +=   2 * deltaR * deltaB * 1.25913;

	// Return the result
	return chi2;
}

// Routine to compute a chi^2 value for hubble, based on some desired value
double chi2hubble(Parameters &params, double desired, double sigma, Output &output){

	if (sigma == 0.0) return 0;

	double x = (params.geth() - desired) / sigma;
	if (fabs(x) < 1e-7) x = 0;
	// If the chi^2 is less than 10^-14, then just report it as zero

	output.printvalue("Hubblechi", x * x);

	return x * x;

}

// Helper routine that computes DV given various values
static inline double getDV(double z, double DA, double DH) {
	// DV = ((1 + z)^2 DA^2 z DH)^(1/3)
	return pow((1 + z) * (1 + z) * DA * DA * z * DH, 1.0 / 3.0);
}

// Helper routine that calculates the chi^2 for SDSS
static inline double getSDSSresult(double d2, double d35) {
	double x1 = d2 - 0.1905;
	double x2 = d35 - 0.1097;

	return x1 * x1 * 30124 + x2 * x2 * 86977 - 2 * x1 * x2 * 17227;
}

// Helper routine that calculates the chi^2 for WiggleZ
static inline double getWiggleZresult(double A44, double A6, double A73) {
	double x1 = A44 - 0.474;
	double x2 = A6 - 0.442;
	double x3 = A73 - 0.424;

	double chi2 = 0;
	// Add contributions from diagonals first
	chi2 += x1 * x1 * 1040.3;
	chi2 += x2 * x2 * 3720.3;
	chi2 += x3 * x3 * 2914.9;
	// Add contributions from off-diagonal terms next
	chi2 += - 2 * x1 * x2 * 807.5;
	chi2 += 2 * x1 * x3 * 336.8;
	chi2 += - 2 * x2 * x3 * 1551.9;

	return chi2;
}

// Routine to compute chi^2 values for BAO observations
void chi2BAO(double rdrag, vector<double>& redshift, vector<double>& hubble, vector<double>& DA, Parameters &params, Output &output, double &chitotal, double &chirtotal) {
	// Takes in the sound horizon at the drag epoch (in Mpc), the redshift, hubble and angular diameter distance measures,
	// the evolution parameters (we're going to need Omega_M), and the output class for logging purposes.
    // Returns the sum of BAO chi^2 values in chitotal (SDSS) and chirtotal (SDSSr)

	// Each separate experiment seems to have their own way of presenting results, so we'll have to do them one by one, carefully
	// We'll need to interpolate to obtain H(z) and DA(z) values
	double Hval;
	double DAval;
	double Hval2;
	double DAval2;
	double Hval3;
	double DAval3;
	// Note that what we actually want is c / H(z)
	// What we have is the dimensionless H. So, we need c / H_0 / \tilde{H}
	double conH0 = params.getconH0();
	// Also recall that we have conformal H, not cosmic time H. To obtain cosmic time H, divide by a (or multiply by 1+z)
	double DV;
	double DV2;
	double DV3;
	double result;
	// Zero the running totals
	chitotal = 0; // Using SDSS
	chirtotal = 0; // Using SDSSR

	// Extract the appropriate data
	double* pz = &redshift[0];
	double* pDA = &DA[0];
	double* pH = &hubble[0];
	int numrows = redshift.size();

	// Construct the splines
	// DA
	splinetools DAspline;
	DAspline.acc = gsl_interp_accel_alloc();
	DAspline.spline = gsl_spline_alloc(gsl_interp_cspline, numrows);
	gsl_spline_init (DAspline.spline, pz, pDA, numrows);
	// H
	splinetools Hspline;
	Hspline.acc = gsl_interp_accel_alloc();
	Hspline.spline = gsl_spline_alloc(gsl_interp_cspline, numrows);
	gsl_spline_init (Hspline.spline, pz, pH, numrows);


	// 6dFGS: z = 0.106
	DAval = gsl_spline_eval (DAspline.spline, 0.106, DAspline.acc);
	Hval = gsl_spline_eval (Hspline.spline, 0.106, Hspline.acc) * 1.106;
	DV = getDV(0.106, DAval, conH0 / Hval);
	result = (rdrag / DV - 0.336) / 0.015;

	output.printvalue("6dFGSchi", result * result);
    chitotal += result * result;
    chirtotal += result * result;


	// SDSS: z = 0.2, z = 0.35
	DAval = gsl_spline_eval (DAspline.spline, 0.2, DAspline.acc);
	Hval = gsl_spline_eval (Hspline.spline, 0.2, Hspline.acc) * 1.2;
	DAval2 = gsl_spline_eval (DAspline.spline, 0.35, DAspline.acc);
	Hval2 = gsl_spline_eval (Hspline.spline, 0.35, Hspline.acc) * 1.35;
	DV = getDV(0.2, DAval, conH0 / Hval);
	DV2 = getDV(0.35, DAval2, conH0 / Hval2);

	result = getSDSSresult(rdrag / DV, rdrag / DV2);
	output.printvalue("SDSSchi", result);
    chitotal += result;

	// SDSSR: z = 0.35
	// Use DV2 from above
	result = (DV2 / rdrag - 8.88) / 0.17;
	output.printvalue("SDSSRchi", result * result);
    chirtotal += result * result;


	// WiggleZ: z = 0.44, z = 0.6, z = 0.73
	DAval = gsl_spline_eval (DAspline.spline, 0.44, DAspline.acc);
	Hval = gsl_spline_eval (Hspline.spline, 0.44, Hspline.acc) * 1.44;
	DV = getDV(0.44, DAval, conH0 / Hval);

	DAval2 = gsl_spline_eval (DAspline.spline, 0.6, DAspline.acc);
	Hval2 = gsl_spline_eval (Hspline.spline, 0.6, Hspline.acc) * 1.6;
	DV2 = getDV(0.6, DAval2, conH0 / Hval2);

	DAval3 = gsl_spline_eval (DAspline.spline, 0.73, DAspline.acc);
	Hval3 = gsl_spline_eval (Hspline.spline, 0.73, Hspline.acc) * 1.73;
	DV3 = getDV(0.73, DAval3, conH0 / Hval3);

	// Convert these DV values into A values
	double fact = sqrt(params.getOmegaM()) / params.getDH();
	double A1 = DV * fact / 0.44;
	double A2 = DV2 * fact / 0.6;
	double A3 = DV3 * fact / 0.73;

	result = getWiggleZresult(A1, A2, A3);
	output.printvalue("WiggleZchi", result);
    chitotal += result;
    chirtotal += result;


	// BOSS DR9: z = 0.57
	DAval = gsl_spline_eval (DAspline.spline, 0.57, DAspline.acc);
	Hval = gsl_spline_eval (Hspline.spline, 0.57, Hspline.acc) * 1.57;
	DV = getDV(0.57, DAval, conH0 / Hval);
	result = (DV / rdrag - 13.67) / 0.22;

	output.printvalue("BOSSDR9chi", result * result);
    chitotal += result * result;
    chirtotal += result * result;


	// BOSS DR11: z = 2.34
	DAval = gsl_spline_eval (DAspline.spline, 2.34, DAspline.acc);
	Hval = gsl_spline_eval (Hspline.spline, 2.34, Hspline.acc) * 3.34;
	double alphapar = conH0 / Hval / rdrag / 8.708;
	double alphaperp = DAval / rdrag / 11.59;
	result = pow(alphapar, 0.7) * pow(alphaperp, 0.3);
	double chi = (result - 1.025) / 0.021;

	output.printvalue("BOSSDR11chi", chi * chi);
    chitotal += chi * chi;
    chirtotal += chi * chi;


    // Totals
    output.printvalue("BAOtotalchi", chitotal);
    output.printvalue("BAOtotalrchi", chirtotal);


	// Release the spline memory
	gsl_spline_free (DAspline.spline);
	gsl_interp_accel_free (DAspline.acc);
	gsl_spline_free (Hspline.spline);
	gsl_interp_accel_free (Hspline.acc);

}
