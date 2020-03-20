#!/usr/bin/env python
"""This is a simple example for tuning C code with the kernel tuner"""

from kernel_tuner import tune_kernel
from collections import OrderedDict

kernel_string = """ 
#include <omp.h>
#include <iostream>
#include <cstdlib>
#include <cstring>
#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;
using namespace std;

extern "C" float solver_tuner(void) {
	// get benchmark suite
	string directory="/home/muhos/tuning_probs/";
  //string directory="/home/muhos/solver_tuner/";
	vector<string> cnfs;
	for (auto& p : fs::directory_iterator(directory)) {
		string file = p.path().string();
		string fileName = p.path().filename().string();
		string ext = ".cnf";
		if (fileName.find(ext) != -1)  {
			cnfs.push_back(file);	
		}
	}
	// start tuning
	// parameters playground
	string argIncSmall       = " --inc-small="       + to_string(incSmall);
	string argIncBig         = " --inc-big="         + to_string(incBig);
	string argLBDFrozen      = " --lbd-frozen="      + to_string(LBDFrozen);
	string argLBDMinSize     = " --lbd-min-size="    + to_string(LBDMinSize);
	string argLBDMin         = " --lbd-min="         + to_string(LBDMin);
	string argRestRate       = " --rf-rate="         + to_string(restRate);
	string argRestDefuse     = " --rb-rate="         + to_string(restDefuse);
	string argVarDecay       = " --var-decay="       + to_string(varDecay);
	string argVarDecayRate   = " --var-decay-r="     + to_string(varDecayRate);
	string argCBTDist        = " --cbt-dist="        + to_string(CBTDist);
	string argCBTConf        = " --cbt-conf="        + to_string(CBTConf);
	//============================
	// solver executable binary
	float sumTime = 0.0;
  string exec =  "./parafrost ";
  string fixedArgs = " -q -no-perf --timeout=1000 --pdm=3 --pre-delay=50";
  fixedArgs += argIncSmall + argIncBig;
	fixedArgs += argLBDFrozen + argLBDMinSize + argLBDMin;
	fixedArgs += argRestRate + argRestDefuse;
	fixedArgs += argVarDecay + argVarDecayRate;
	fixedArgs += argCBTDist + argCBTConf;
	string allArgs = fixedArgs;
 cout << "== Running the tuner on " << cnfs.size() << " files." << endl;
 cout << "== Arguments: " << allArgs << endl;  
	for (string& f: cnfs) {
		// exec solver with arguments
		string solver = exec + f + allArgs;
		double start = omp_get_wtime();
		int status = system(solver.c_str());
		sumTime += (float)(omp_get_wtime() - start);
	}
	cout << "== Runtime = " << sumTime << " sec." << endl;
  return sumTime;
}
"""

args = []

tune_params = OrderedDict()
tune_params["incSmall"] = [300, 500, 600]
tune_params["incBig"] = [1000, 2000, 3000]
tune_params["LBDFrozen"] = [20, 30, 40]
tune_params["LBDMinSize"] = [20, 30, 40]
tune_params["LBDMin"] = [5, 8, 10, 20]
tune_params["restRate"] = [0.8, 0.85, 0.9]
tune_params["restDefuse"] = [1.5, 2, 3]
tune_params["varDecay"] = [0.5, 0.6, 0.7, 0.8]
tune_params["varDecayRate"] = [0.001, 0.005, 0.01]
tune_params["CBTDist"] = [100, 500]
tune_params["CBTConf"] = [5000, 10000]

tune_kernel("solver_tuner", kernel_string, 1, args, tune_params, iterations=1, strategy="basinhopping", strategy_options={"method":"BFGS","maxiter":100,"T":1.0}, answer=[], compiler_options=['-O3','--std=c++11','-lstdc++fs'])
