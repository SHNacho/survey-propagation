#ifndef FACTOR_GRAPH_H_
#define FACTOR_GRAPH_H_

#include <iostream>
#include <vector>


using namespace std;

class Literal;

class Clause{
public:
	vector<Literal*> literals;
	int unassigned_literals;
	bool satisfied;
	
	// Methods:
	Clause();
};


class Variable{
public:
	vector<Literal*> literals; //

	int value; //-1, 0, 1 -- 0=unfixed
	
	double p; //product of (1-eta) of literals of same sign
	double m; //product of (1-eta) of literals of oposite sign
	
	int pzero; //equal to 1 if p equal to 0
	int mzero; //equal to 1 if m equal to 0

	// Biases -- equations (28, 29, 30)
	double wp; //plus
	double wz; //zero
	double wm; //minus
	
	// Methods:
};

/**
* Un literal representa el estado de una variable (positivo o negativo)
* en una cláusula.
* Se puede ver como una arista del grafo ya que une una variable con una
* cláusula
*/
class Literal{
public:
	Variable* var;
	Clause* cl;
	int type; // negative: -1; positive: 1
	bool enabled;
	double survey;
	
	// Methods:
	Literal(Variable* var, Clause* cl, int type);
	void Disable();
};

class FactorGraph{
public:
	vector<Literal*> literals;
	vector<Variable*> variables;
	vector<Clause*> clauses;

	int unassigned_vars;

	// Methods:
	bool simplify(Variable* var);
	bool fix(Variable* var, int val);
	bool fixUnitClause(Clause* c);
};

#endif // FACTOR_GRAPH_H_
