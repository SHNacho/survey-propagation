
#include <iostream>
#include "SPSolver.h"

using namespace std;

int main(){
	string file = "data/cnf_1.txt";
	FactorGraph fg(file);
    SPSolver solver(&fg, 0.04);
    if(solver.surveyInspiredDecimation()){
        cout << "OK" << endl;
    } else {
        cout << ":-(" << endl; 
    }
    return 0;
}
