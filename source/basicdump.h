/*
 * basicdump.h
 *
 * This is an output module. It dumps basically everything to file, creating a log file that can be machine read as a parameter file.
 *
 */

#ifndef BASICDUMP_H_
#define BASICDUMP_H_

#include "output.h"
#include <iostream>
#include <fstream>
#include <iomanip>

class BasicDump : public Output {
	public:
		// Functions that are overridden from the Output class
		void printinfo(const double data[], Parameters&);
		void printheading();
		void printstep(const double data[], const double time, const double status[]);
		void postprintheading();
		void postprintstep(const double z, const double H, const double DC, const double DM, const double DA, const double DL, const double mu);
		void printfinish(const double time);
		bool filesready();
		void printlog(const std::string&);
		void printvalue(const std::string&, const std::string&);
        void printvalue(const std::string&, const double);
        void printvalue(const std::string&, const int);

		// Constructor
		BasicDump(bool postprocess, const std::string &filename = "run", const std::string &postname = "d");
		// Destructor
		~BasicDump(); // Overrides Output class

	private:
		// File output information
		std::ofstream *myLog;
		std::ofstream *myData;
		std::ofstream *myPostData;

		// Postprocessing?
		bool doingpp;

};

#endif /* BASICDUMP_H_ */
