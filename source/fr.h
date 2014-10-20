/*
 * fr.h
 *
 * Outlines the F(R) class
 *
 * Uses an action S = \int d^4x \sqrt{-g} m_P^2/2 (R + F(R)) + matter
 *
 * Translated into a scalar field theory, this becomes the following:
 *
 * S = \int d^4x \sqrt{-g} m_P^2/2 (R + F[phi] + (R - phi) F'[phi]) + matter
 */

#ifndef FR_H_
#define FR_H_

#include "model.h"

class Fr : public Model {

    public:
        // Here are the functions that are overridden by the F(R) class
        int derivatives(const double data[], double derivs[], Parameters &params);
        int init(double data[], const double time, Parameters &params, IniReader &init, Output &output);

        // We actually need to test for a valid configuration here
        bool isvalidconfig(const double data[]);

        // The speed of sound and speed of tensor propagation is unity in F(R) models.
        double speedofsound2(const double data[]) {return 1;}
        double speedoftensor2(const double data[]) {return 1;}
        bool implementsSOS() {return true;}
        bool implementsSOT() {return true;}

        // The no-ghost condition for both tensor and scalar modes is the same; we only need to implement one function
        bool implementsghost() {return true;}
        bool implementstensorghost() {return true;}
        bool isghost(const double data[]);
        bool istensorghost(const double data[]) {return isghost(data);}


    private:
        // Here are some overridden internal functions. They're pretty self-explanatory.
        double pressure(const double data[], const double hdot);
        double energydensity(const double data[]);

        // Just some parameters
        double alpha;
        double n;

        // Function to calculate the Lagrangian and all its appropriate derivatives:
        // F, F', F'', F'''
        int computelagrangian(const double data[]);
        // Each time the stuff is calculated, store both the data and it, so as not to waste computation time
        double storeddata[4];
        double valF;
        double valFp;
        double valFpp;
        double valFppp;

};

#endif /* FR_H_ */
