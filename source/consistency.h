/*
 * consistency.h
 *
 * This header defines a class for a module that checks consistency conditions on the data as its being evolved.
 * Note that this is not an abstract class; it can be used as is if no consistency checks are desired.
 *
 */

#ifndef CONSISTENCY_H_
#define CONSISTENCY_H_

#include "intparams.h"
#include "output.h"

class Consistency {
	public:
		// These functions take as arguments: data, time, parameters, and output class
		// Check the state of the data in the middle of a run
		virtual void checkstate (const double*, const double, IntParams&, Output&, const double*) {}
		// Check the state of the data at the end of a run
		virtual void checkfinal (const double*, const double, IntParams&, Output&, const double*) {}

		// Check for NaN results. Return true if any are found.
		bool checknan (const double data[], const double time, const double status[]) {
			if (isnan(time))
				return true;

			for (int i = 0; i < 4; i++) {
				if (isnan(data[i]))
					return true;
			}

			for (int i = 0; i < 17; i++) {
				if (isnan(status[i]))
					return true;
			}

			return false;
		}

		virtual ~Consistency() {} // Destructor
};


#endif /* CONSISTENCY_H_ */
