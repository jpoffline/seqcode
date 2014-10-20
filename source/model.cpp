#include "model.h"

/* The state routine is given the state of the system as well as the parameters of the model,
   and returns information in the info array. The return values are as follows:

   * 0 time
   * 1 a
   * 2 Redshift
   * 3 H = \dot{a}/a
   * 4 \dot{H}
   * 5 phi
   * 6 \dot{phi}
   * 7 \ddot{phi}
   * 8 Omega_matter (present value)
   * 9 Omega_radiation (present value)
   * 10 Omega_k (present value)
   * 11 Omega_Q (present value)
   * 12 w_total
   * 13 rho_Q / rho_c
   * 14 P_Q / rho_c
   * 15 w_Q
   * 16 Error

   This routine should be largely model independent, and as such is included in the abstract Model class.

 */
void Model::getstate(const double data[], const double time, double info[], Parameters &params) {

	// Extract data for easier reading of the code
	double a = data[0];
	double phi = data[1];
	double phidot = data[2];
	// Calculate a^2, a^4
	double a2 = a * a;
	double a4 = a2 * a2;
	// Go and compute the derivatives
	double derivs[4];
	derivatives(data, derivs, params);
	// Hubble parameter H = \dot{a}/a
	double hubble = derivs[0] / a;
	double hubble2 = hubble * hubble;
	// Energy density and pressure of dark energy
	double energy = energydensity(data);
	double press = pressure(data, derivs[3]);

	// time
	info[0] = time;
	// a
	info[1] = a;
	// Redshift
	info[2] = 1.0 / a - 1.0;
	// Hubble
	info[3] = hubble;
	// \dot{H}
	info[4] = derivs[3];
	// phi
	info[5] = phi;
	// \dot{phi}
	info[6] = phidot;
	// \ddot{phi}
	info[7] = derivs[2];
	// Omega_matter (present value)
	info[8] = params.rhoM() / a / hubble2;
	// Omega_radiation (present value)
	info[9] = params.rhoR() / a2 / hubble2;
	// Omega_k (present value)
	info[10] = params.rhoK() / hubble2;
	// Omega_Q (present value)
	info[11] = a2 * energy / hubble2;
	// w_total
	info[12] = (params.rhoR() / 3 + a4 * press) / (a * params.rhoM() + params.rhoR() + a4 * energy);
	// rho_Q / rho_0
	info[13] = energy;
	// P_Q / rho_0
	info[14] = press;
	// w_Q
	info[15] = press / energy;

	/* The Friedmann equation reads
 	   H^2 = \Omega_m / a + \Omega_r / a^2 + \Omega_k + a^2 \rho_Q / \rho_0
 	   The components on the right may require H for their calculation, but that's ok; they're all been calculated
 	   The LHS is then the exact Hubble parameter, given the constituents on the RHS
 	   Dividing this equation by ~H, the Hubble parameter that we have, we obtain
 	   H^2 / ~H^2 = \Omega_m / a ~H^2 + \Omega_r / a^2 ~H^2 + \Omega_k / ~H^2 + a^2 \rho_Q / \rho_c ~H^2
	   The components on the RHS have already been computed (in info[8-11]), so we can obtain that simply.
	   Taking the square root of that quantity, we get H / ~H.
	   In computing the relative error in H, we want to calculate
	   (~H - H) / ~H = 1 - H / ~H, which is easily obtained from what we have.
	   I believe this is model independent.
	 */

	info[16] = 1 - sqrt(info[8] + info[9] + info[10] + info[11]);

}
