/*
 * print2memory.h
 *
 * This defines an output class whereby the log data is saved in memory
 * Note that this does not include the evolution data, just the input parameters and final results
 *
 * This class works by storing information in a property tree. Whenever a value is sent to be printed, it stores it.
 * It also stores some information from the parameters of the evolution.
 *
 * The "printfinal" routine instructs the class clears all stored information.
 *
 * There is an extra routine beyond the base class that allows for properties to be read.
 *
 */

#ifndef PRINT2MEMORY_H_
#define PRINT2MEMORY_H_

#include "output.h"
#include <boost/property_tree/ptree.hpp>

using std::string;

class Print2Memory : public Output {
    public:
        // Function to extract data from the evolution
        void printinfo(const double data[], Parameters&);

        // Routine to store information as it comes in
        void printvalue(const std::string&, const std::string&);
        void printvalue(const std::string&, const double);
        void printvalue(const std::string&, const int);

        // Function to clear the information from the store
        virtual void printfinish(const double time);

        // Overloaded function to extract saved information
        string getvalue (const string &key, const string &def = "");
        double getvalue (const string &key, const double &def = 0.0);
        double getvalue (const string &key, const int &def = 0);

    private:
        // The property tree, where everything is stored
        boost::property_tree::ptree data;

};

#endif /* PRINT2MEMORY_H_ */
