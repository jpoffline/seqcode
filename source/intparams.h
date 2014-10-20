/*
 * intparams.h
 *
 * This defines the IntParams class, which is essentially a holding class for a Parameters class
 * and a Model class. It is intended that this class is passed into the integration routine.
 *
 */

#ifndef INTPARAMS_H_
#define INTPARAMS_H_

#include <iostream> // for the NULL constant
#include "params.h"
#include "model.h"

class IntParams {

	public:
		// Constructor. Point the class members to the two objects that were passed in.
		IntParams(Parameters &params, Model &model){
			myParams = &params;
			myModel = &model;
		}

		// Destructor. Remove references to objects that now no longer need to be pointed to.
		~IntParams() {
			myParams = NULL;
			myModel = NULL;
		}

		// Getters for the stored classes
		// Return the model class that was stored
		Model &getmodel() {
			return *myModel;
		}
		// Return the parameters class that was stored
		Parameters &getparams() {
			return *myParams;
		}

	private:
		// Internal storage for the classes
		Model *myModel;
		Parameters *myParams;

};

#endif /* INTPARAMS_H_ */
