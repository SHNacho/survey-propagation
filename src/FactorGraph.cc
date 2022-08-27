#include "FactorGraph.h"
#include <iostream>
#include <string>
#include <unordered_set>
#include <fstream>
#include <sstream>
#include <vector>

using namespace std;

Variable::Variable(int id){
	this->value = 0;
	this->id = id;
	this->p = 1;
	this->m = 1;
	this->mzero = 0;
	this->pzero = 0;
	wp = wm = wz = 0;
}

Literal::Literal(Variable* var, Clause* cl, int type){
	this->cl = cl;
	this->var = var;
	this->type = type;
	cl->literals.push_back(this);
	var->literals.push_back(this);
	this->enabled = true;
}

Clause::Clause(){
	this->satisfied = false;
	this->unassigned_literals = 0;
}

FactorGraph::FactorGraph(string file){
	this->unassigned_vars = 0;

	fstream fs;
	fs.open(file);
	
	string line;
	getline(fs, line);
	while(line[0] == 'c'){
		getline(fs, line);
	};

	if(line[0] == 'p'){
		vector<string> tokens = splitString(line);
		int numVariables = stoi(tokens[2]);
		int numClauses = stoi(tokens[3]);

		for (int i = 0; i < numVariables; ++i){
			variables.push_back(new Variable(i+1));
			unassigned_vars++;
		}

		for (int i = 0; i < numClauses; ++i){
			clauses.push_back(new Clause);
		}
	}

	int cl = 0;
	while (getline(fs, line)){
		vector<string> tokens = splitString(line);
		for (string t : tokens){
			int var = stoi(t);
			if (var != 0){
				int type = var > 0 ? 1 : -1;
				int id = abs(var);
				literals.push_back(new Literal(variables[id - 1], clauses[cl], type));
				clauses[cl]->unassigned_literals++;
			}
		}
		cl++;
	}

	fs.close();
}

vector<string> FactorGraph::splitString(string s) {
	const char delim = ' ';
	stringstream stream(s);

	vector<string> tokens;

	string token;
	while (getline(stream, token, delim)) {
		tokens.push_back(token);
	}

	return tokens;
}

bool FactorGraph::fix(Variable* var, int val){
	if(var->value) // si es 1 o -1, ya está asignado
		return true;

	var->value = val;
	unassigned_vars--;
	
	for(Literal* l : var->literals){
		l->cl->unassigned_literals--;
	}
	
	return simplify(var);
}

/**
* Simpify after fixing var
*/
bool FactorGraph::simplify(Variable* var){
	// Para cada cláusula de la variable
	for(Literal* l : var->literals){
		// Se deshabilita la arista
		l->enabled = false;

		Clause* c = l->cl;

		if(!c->satisfied){
			// Se comprueba si la asignación satisface la cláusula
			if(l->type == var->value)
				c->satisfied = true;
			// Si no la satisface se comprueba si es unitaria
			else if(c->unassigned_literals == 1){
				// Se asigna la variable de la cláusula
				if(!fixUnitClause(c)) return false;
			}
			// Si no está satisfecha y no contiene literales, contradicción
			else if(c->unassigned_literals == 0){
				cout << "CONTRADICTION" << endl;
				return false;
			}
		}
	}
	return true;
}

/**
* Asigna un valor al literal de la cláusula unitaria
*/
bool FactorGraph::fixUnitClause(Clause* c){
	// Se busca el literal que queda sin asignar
	for(Literal* l : c->literals){
		if(l->enabled)
			// Se asigna
			return fix(l->var, l->type);
	}
	return false;
}















