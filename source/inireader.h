/*
 * inireader.h
 *
 * This file reads in the initialization file, and provides methods to extract entries from it.
 * Heavily relies upon the BOOST C++ libraries
 *
 */

#ifndef INIREADER_H_
#define INIREADER_H_

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <boost/filesystem.hpp>
#include <iostream>
#include <sstream>

using std::string;

class IniReader {
    public:
        void read (const std::string &filename);

        // These two routines allow the entire parameter structure to be imported or exported
        // Returns the full property tree, typically so that setinidata can be used
        boost::property_tree::ptree getdata() {return inifile;}
        // Instead of reading an ini file, this allows the data to be submitted directly,
        // allowing a programmatic setting of the ini data
        void setdata(boost::property_tree::ptree data) {inifile = data;}

        // Routines that allow data to be written directly to the stored parameters
        void setparam(const std::string &key, const std::string &section, const std::string &value);
        void setparam(const std::string &key, const std::string &section, const double value);
        void setparam(const std::string &key, const std::string &section, const int value);

        // The following overloaded routines are all identical, but written for different data types.
        // They all return information from the ini file, of the specified data type.
        // The parameters are identical for each:
        // key is the name of the key whose value is desired
        // def is the default to be returned if no value is present
        // section is an optional argument for the section of the ini file (defaults to no section)
        string getiniString (const string &key, const string &def = "", const string &section = "");
        double getiniDouble (const string &key, const double &def = 0.0, const string &section = "");
        int getiniInt (const string &key, const int &def = 0, const string &section = "");
        bool getiniBool (const string &key, const bool &def = false, const string &section = "");

    private:
        // Stores the content of the ini file
        boost::property_tree::ptree inifile;

        // Returns the appropriate section from the ini file
        bool getsection(const std::string &section, boost::property_tree::ptree &result);

};

#endif /* INIREADER_H_ */
