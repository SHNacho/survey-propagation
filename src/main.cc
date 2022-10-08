
#include <iostream>
#include "SPSolver.h"

using namespace std;

int main(){
	string file = "data/cnf_1.txt";
    int solved = 0;
    for (int i = 1; i <= 50; ++i){
        file = "data/S-4.21/cnf_" + to_string(i) + ".txt";
	    sp::FactorGraph fg(file);
        sp::SPSolver solver(&fg, 0.02);
        if(solver.surveyInspiredDecimation()){
            solved++;
            cout << "OK" << endl;
        } else {
            cout << ":-(" << endl; 
        }
    }

    cout << "Porcentaje de instancias resueltas: " << solved*1.0 / 50.0 * 100.0 << endl;
    return 0;
}
