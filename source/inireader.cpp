#include "inireader.h"

using namespace boost::property_tree;

void IniReader::read (const string &filename) {
    // Reads the ini file

    // If no params.ini file exists, use default parameters
    if (!boost::filesystem::exists(filename)) {
        std::cout << "Warning: " << filename << " does not exist. Using default parameters." << std::endl;
    } else {
        ini_parser::read_ini(filename, inifile);
    }

}

// Routines that allow data to be written directly to the stored parameters
void IniReader::setparam(const string &key, const string &section, const string &value) {
    std::stringstream keyname;
    if (section == "") {
        keyname << key;
    } else {
        keyname << section << "." << key;
    }
    inifile.put(keyname.str(), value);
}
void IniReader::setparam(const string &key, const string &section, const double value) {
    std::stringstream keyname;
    if (section == "") {
        keyname << key;
    } else {
        keyname << section << "." << key;
    }
    inifile.put(keyname.str(), value);
}
void IniReader::setparam(const string &key, const string &section, const int value) {
    std::stringstream keyname;
    if (section == "") {
        keyname << key;
    } else {
        keyname << section << "." << key;
    }
    inifile.put(keyname.str(), value);
}


// The following overloaded routines are all identical, but written for different data types.
// They all return information from the ini file, of the specified data type.
// The parameters are identical for each:
// key is the name of the key whose value is desired
// def is the default to be returned if no value is present
// section is an optional argument for the section of the ini file (defaults to no section)
string IniReader::getiniString (const string &key, const string &def, const string &section) {

    ptree usetree;
    if (getsection(section, usetree) == true) {
        // The section exists (or there was no section), go and pull out the data
        return usetree.get<string>(key, def);
    } else {
        // The section doesn't exist, return the default
        return def;
    }

}

double IniReader::getiniDouble (const string &key, const double &def, const string &section) {

    ptree usetree;
    if (getsection(section, usetree) == true) {
        // The section exists (or there was no section), go and pull out the data
        return usetree.get<double>(key, def);
    } else {
        // The section doesn't exist, return the default
        return def;
    }

}

int IniReader::getiniInt (const string &key, const int &def, const string &section) {

    ptree usetree;
    if (getsection(section, usetree) == true) {
        // The section exists (or there was no section), go and pull out the data
        return usetree.get<int>(key, def);
    } else {
        // The section doesn't exist, return the default
        return def;
    }

}

bool IniReader::getiniBool (const string &key, const bool &def, const string &section) {
    // Note that 1 is true, anything else is false (I think?)

    ptree usetree;
    if (getsection(section, usetree) == true) {
        // The section exists (or there was no section), go and pull out the data
        return usetree.get<bool>(key, def);
    } else {
        // The section doesn't exist, return the default
        return def;
    }

}

// Returns the appropriate section from the ini file
bool IniReader::getsection(const string &section, ptree &result) {

    // If we are not in a section, return the tree
    if (section == "") {
        result = inifile;
        return true;
    } else {
        // If we are in a section, try to construct the tree for that section
        try {
            result = inifile.get_child(section);
        }
        catch (...){
            // Likely got here because that section doesn't exist
            return false;
        }
    }
    return true;

}
