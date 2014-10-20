/*
 * print2memory.cpp
 *
 * Contains code supporting the Print2Memory output class
 */

#include "print2memory.h"

// Function to print data before the run starts
void Print2Memory::printinfo(const double data[], Parameters& params){

    // Save the cosmological parameters
    printvalue("OmegaMh2", params.rhoM());
    printvalue("OmegaBh2", params.rhoB());
    printvalue("OmegaKh2", params.rhoK());
    printvalue("OmegaRh2", params.rhoR());
    printvalue("Tgamma", params.Tgamma());
    printvalue("zinit", params.z0());
    printvalue("Neff", params.Neff());

}

// Routine to print a key and value to the log file
// Three overloaded versions
void Print2Memory::printvalue(const std::string& key, const std::string& value){
    data.put(key, value);
}
void Print2Memory::printvalue(const std::string& key, const double value){
    data.put(key, value);
}
void Print2Memory::printvalue(const std::string& key, const int value){
    data.put(key, value);
}

// Function to clear all current data
void Print2Memory::printfinish(const double time) {
    // Clears the data tree
    data.clear();
}

// Overloaded function to extract saved information
string Print2Memory::getvalue (const string &key, const string &def) {
    return data.get<string>(key, def);
}
double Print2Memory::getvalue (const string &key, const double &def) {
    return data.get<double>(key, def);
}
double Print2Memory::getvalue (const string &key, const int &def) {
    return data.get<int>(key, def);
}
