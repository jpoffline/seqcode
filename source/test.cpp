
#include "main.h"

using namespace std;
double SEQ_evolveit(IniReader& inifile, Print2Memory& output, bool usingSN1a, vector<vector<double> > &SN1adata, string paramname, string paramsection, double param);
double NewtRaph_find(IniReader& inifile, Print2Memory& output, bool usingSN1a, vector<vector<double> > &SN1adata, string paramname, string paramsection, double paramstart, double DeltaParam, double seq_acc, int NRstepsMAX, double param_seq);

int main(int argc, char* argv[]) {

    // Start timing!
    boost::timer::cpu_timer myTimer;

    // Read the input parameters file, named "params.ini" by default.
    // If there is a command line argument, assume that it is the filename for the parameters file.
    IniReader inifile;
    string paramsfile;
    if (argc > 1)
        paramsfile = argv[1];
    else
        paramsfile = "params.ini";
    inifile.read(paramsfile.c_str());

    // Print2Memory is used for outputting
    Print2Memory myOutput;
	
	bool usingSN1a;
	
    // Load SN1a data
    vector<vector<double> > SN1adata;
    string sn1afile = inifile.getiniString("union21", "SCPUnion2.1_mu_vs_z.txt", "Function");
    if (loadSN1adata(sn1afile, SN1adata) != 0) {
        // Could not find data file
        cout << "Warning: cannot find Union2.1 SN1a data file." << endl;
    }
	
	// parameter name to be dialed
	string paramname = "phi0"; 
	// Section name of dialed parameter
	string paramsection = "Cosmology"; 
	// Value of the parameter which will be returned as the sequestering-compatible parameter
	double param_seq;
	// Start value of the parameter to be dialed
	double paramstart = 1.0;
	// Size of jumps to be made for finite differences
	double DeltaParam = 0.01;
	// Required accuracy of Newton-Raphson 
	double seq_acc = 1E-10; 
	// Maximum number of NR iterations before stopping
	int NRstepsMAX = 1E5;
	// Use Newton Raphson to obtain the value of "paramname" which is consistent with a sequestering solution 
	bool isfound = NewtRaph_find(inifile, myOutput, usingSN1a, SN1adata, paramname, paramsection, paramstart, DeltaParam, seq_acc, NRstepsMAX, param_seq);
	
	if(isfound){
		cout << "sequestering-compatible parameter found :: " + paramname + " = " << param_seq << endl;
	}
	else
		cout << "NR maxsteps reached" << endl;
	
    // Stop timing
    myTimer.stop();
    double ms = myTimer.elapsed().wall / 1e6;
    if (ms < 1e3)
        cout << setprecision(4) << "Run complete in " << ms << " milliseconds.";
    else
        cout << setprecision(4) << "Run complete in " << ms / 1000.0 << " seconds.";
    cout << endl;

    // Exit gracefully
    return 0;
	
} // END main()


double NewtRaph_find(IniReader& inifile, Print2Memory& output, bool usingSN1a, vector<vector<double> > &SN1adata, string paramname, string paramsection, double paramstart, double DeltaParam, double seq_acc, int NRstepsMAX, double param_seq){

	/*
		This routine uses a Newton-Raphson method to
		obtain sequestering solutions.
	*/

	// This is the parameters we will be dialing.
	// Start it off at some fiducial value -- hopefully, close to the
	// sequestering solution. That needs some trial & error to find.
	double param = paramstart;
	double HIR0, HIRf, HIRb, HIRn;
	double retval;
	
	// Counter for the number of Newton-Raphson steps
	int NRsteps = 0;
	
	while(true){
		NRsteps++;
		// Main NR algorithm: 
		// Compute HIR = <R> at different parameter values
		// n.b., HIR = historic integral of Ricci Scalar
		// (1) Get value of HIR at current parameter value
		HIR0 = SEQ_evolveit(inifile, output, usingSN1a, SN1adata, paramname, paramsection, param);
		
		// (2) Get HIR either side of current parameter value (for use in finite difference)
		// (2a) forwards
		HIRf = SEQ_evolveit(inifile, output, usingSN1a, SN1adata, paramname, paramsection, param + DeltaParam);
		// (2b) backwards
		HIRb = SEQ_evolveit(inifile, output, usingSN1a, SN1adata, paramname, paramsection, param - DeltaParam);
		
		// Get "next value" of the parameter using Newton-Raphson (NR) method,
		// and finite difference derivative (FD) approximation:
		// NR :: p_{n+1} = p_n - HIR(p_n) / HIR'(p),
		// FD :: HIR'(p) = ( HIR(p+DeltaParam) - HIR(p-DeltaParam) ) / 2 * DeltaParam
		param = param - 2.0 * DeltaParam * HIR0 / ( HIRf - HIRb );
	
		// Compute value of HIR at this new parameter value
		HIRn = SEQ_evolveit(inifile, output, usingSN1a, SN1adata, paramname, paramsection, param);

		cout << HIR0 << endl;
		cout << HIRf << endl;
		cout << HIRb << endl;
		cout << HIRn << endl;		
		// If |HIR| < accuracy required, then we have found a parameter thats good for sequestering.
		if( abs( HIRn ) < seq_acc ){
			retval = param_seq;
			cout << NRsteps << endl;
			return true;
			break;
		}
		
		// If we are just going "indefinitely", exit
		if( NRsteps > NRstepsMAX ){
			return false;
			break;
		}

	} // END NR while-loop
	


} // END NewtRaph_find()

double SEQ_evolveit(IniReader& inifile, Print2Memory& output, bool usingSN1a, vector<vector<double> > &SN1adata, string paramname, string paramsection, double param){
	
	/*
		This routine sets up an instance of deevolve with a given
		value of the parameter which is used to dial to obtain
		sequestering solutions
	*/
	cout << param << endl;
	// Set the value of the parameter which we are dialing
	inifile.setparam(paramname,paramsection,param);
	
	// Setup myParams
	Parameters myParams(inifile);
	
	// Do the evolution, with the given value of the parameter
	// Use "false" here, since we dont want to do any post-processing (data stuff)
	doEvolution(inifile, myParams, output, SN1adata, false);
	
	// Return the value of the historic integral of the Ricci scalar
	return output.getvalue("histintRicci",0.0);
	
} // END evolveit()




// EOF