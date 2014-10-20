/*
 * main.h
 *
 * This is a general header that provides supporting routines to programs implementing the code. As such, this should be considered
 * the one-stop file to include to use all of the available tools in the code.
 *
 * Mostly, it contains includes for the rest of the code, as well as a few routines.
 *
 */

#ifndef MAIN_H_
#define MAIN_H_

#include "inireader.h"
#include "params.h"
#include "output.h"
#include "basicdump.h"
#include "print2memory.h"
#include "evolve.h"

#include <ostream>
#include <string>
#include <cmath>
#include <iomanip>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <boost/filesystem.hpp>
#include <boost/timer/timer.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/normal_distribution.hpp>


// Structure used to store priors
struct PARAMPRIORS{
    string section;
    string name;
    double lower;
    double upper;
    double sigma;
    bool logparam;
};


// Random number generation tools
boost::random::mt19937 RNGtool;
boost::random::uniform_real_distribution<> RNGunit(0.0, 1.0);
boost::random::normal_distribution<> RNGnormal(0.0, 1.0);
// Usage:
// double x = RNGunit(RNGtool) // random real from [0,1)
// double x = RNGnormal(RNGtool) // random real from real distribution with mean 0 and standard deviation 1

// Return random number from unit interval
inline double UnitRand(){return RNGunit(RNGtool);}

// Function to return number from N(0,1): unit normal distribution
inline double NormalRand(){return RNGnormal(RNGtool);}


// Function that finds an appropriate filename (padding is number of characters in the number)
std::string getfilename(const std::string &dir, const std::string &filebase, const std::string &postbase, const int padding = 4, bool dopostprocess = false) {
    // This routine takes in a directory and an output name
    // It goes and finds the first available filename of the form dir / filebase 00001 etc
    // eg., dir/run00001.log and dir/run00001.dat
    // It checks that both files are free
    // If dopostprocess is true, then also checks for dir/run00001d.dat, where the d is whatever is in postbase
    // Note that even if the output is going to screen, this routine won't make anything bad happen

    using namespace boost::filesystem;

    // Firstly, make sure that the directory exists
    if (!exists(dir + "/")) {
        // Directory doesn't exist. Make it.
        create_directory(dir);
        std::cout << "Creating directory " << dir << "/" << std::endl;
    }

    // Secondly, find a unique filename
    for (int counter = 1; ; counter++) {
        // Construct the file number
        string filenum;
        std::ostringstream convert;
        convert << counter;
        filenum = convert.str();
        // Pad the file number with the appropriate number of zeroes
        int len = filenum.length();
        for (int i = 0; i < padding - len; i++)
            filenum = "0" + filenum;

        // Check for the files
        if (exists(dir + "/" + filebase + filenum + ".log"))
            continue;
        if (exists(dir + "/" + filebase + filenum + ".dat"))
            continue;
        if (dopostprocess && exists(dir + "/" + filebase + filenum + postbase + ".dat"))
            continue;

        // If we got to here, we have a unique filename; return it
        return dir + "/" + filebase + filenum;
    }

    // We really shouldn't get here, but parsers like making sure there's a return
    return "error";

}


// Loads the SN1a data from the file into the vector SN1adata
// Returns -1 if file could not be found, 0 for success
int loadSN1adata(const std::string &file, vector<vector<double> > &SN1adata) {

    // Check that the file exists before further processing
    if (boost::filesystem::exists(file)) {
        // Begin by reading in the data file
        std::ifstream f(file.c_str());
        string l;

        // Read in the file one line at a time into the string l
        while(getline(f, l)) {
            // If the first character is a # (comment), move onto the next line
            if (l[0] == '#')
                continue;
            // Convert the string l into a stringstream s
            std::stringstream s(l);
            string extract;
            vector<double> row;
            // Extract entries one at a time, using a tab as a delimiter
            // Ignore the first entry, which is the supernovae name
            bool first = true;
            while(getline(s, extract, '\t')) {
                if (first) {
                    first = false;
                } else {
                    row.push_back(atof(extract.c_str()));
                }
            }
            SN1adata.push_back(row);
        }

        // Return success
        return 0;
    }
    else {
        // Could not find file
        return -1;
    }

}


// A routine that prints a progressbar to screen
// Make sure to call std::cout << endl; to clear the bar after finishing
void updateprogress(float progress) {
    int barWidth = 70;
    std::cout << "[";
    int pos = barWidth * progress;
    for (int i = 0; i < barWidth; ++i) {
        if (i < pos) std::cout << "=";
        else if (i == pos) std::cout << ">";
        else std::cout << " ";
    }
    std::cout << "] " << int(progress * 100.0) << " %\r";
    std::cout.flush();
}

#endif /* MAIN_H_ */
