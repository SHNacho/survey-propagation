#include "FactorGraph.h"
#include <iostream>
#include <string>
#include <unordered_set>
#include <fstream>
#include <sstream>
#include <vector>
#include <set>

using namespace std;

namespace sp{
	//// VARIABLE ////
	Variable::Variable(int id){
		this->value = 0;
		this->id = id;
		this->p = 1;
		this->m = 1;
		this->mzero = 0;
		this->pzero = 0;
		this->ws_break = 0;
		this->ws = false;
		wp = wm = wz = 0;
	}

	Variable::Variable(const Variable &copy){
			this->value = copy.value;
			this->id = copy.id;
			this->p = copy.p;
			this->m = copy.m;
			this->mzero = copy.mzero;
			this->pzero = copy.pzero;
			this->wp = copy.wp;
			this->wm = copy.wm;
			this->wz = copy.wz;
			this->literals = copy.literals;
			this->ws_break = copy.ws_break;
			this->ws = copy.ws;
	}

	vector<Clause*> Variable::TLC(){
		vector<Clause*> vTLC;
		for(Literal* l : literals) if(l->enabled){
			Clause* c = l->cl;
			if(!c->satisfied){ // Quitamos las que ya se hayan satisfecho por SP
				if(this->value == l->type){
					vTLC.push_back(c);
				}
			}
		}
		return vTLC;
	}

	//// LITERAL ////

	Literal::Literal(Variable* var, Clause* cl, int type){
		this->cl = cl;
		this->var = var;
		this->type = type;
		this->cl->literals.push_back(this);
		this->var->literals.push_back(this);
		this->enabled = true;
	}

	Literal::Literal(const Literal &copy){
		this->var = copy.var;
		this->cl = copy.cl;
		this->type = copy.type;
		this->enabled = copy.enabled;
		this->survey = copy.survey;
	}
	//// CLAUSE ////

	Clause::Clause(){
		this->satisfied = false;
		this->unassigned_literals = 0;
		this->ws_unvisited = true;
	}

	Clause::Clause(const Clause &copy){
		this->literals = copy.literals;
		this->unassigned_literals = copy.unassigned_literals;
		this->satisfied = copy.satisfied;
		this->ws_unvisited = true;
	}

	bool Clause::isSAT(){
		for(Literal* l : literals){
			if(l->type == l->var->value)
				return true;
		}

		return false;
	}

	int Clause::NT(){
		int count = 0;
		for(Literal* l : literals) if(l->enabled){
			Variable* v = l->var;
			if(l->type == v->value) count++;		
		}
		return count;
	}

	//// FACTOR GRAPH //

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

	FactorGraph::FactorGraph(vector<vector<int>> &cls){
		set<Variable*, comp> set_vars;

		for (int i = 0; i < cls.size(); ++i){
			clauses.push_back(new Clause);

			for(int j = 0; j < cls[i].size(); ++j){
				int var = cls[i][j];
				int type = var > 0 ? 1 : -1;
				int id = abs(var);

				auto ret = set_vars.insert(new Variable(abs(var)));
				literals.push_back(new Literal(*(ret.first), clauses[i], type));
				clauses[i]->unassigned_literals++;
			}
		}

		this->variables = vector<Variable*>(set_vars.begin(), set_vars.end());
		this->unassigned_vars = variables.size();
	}

	FactorGraph::FactorGraph(const FactorGraph &copy){
		this->literals = copy.literals;
		this->variables = copy.variables;
		this->clauses = copy.clauses;
	}

	FactorGraph::~FactorGraph(){
		for(Literal* l : literals) if(l != nullptr)
			delete l;
		for(Variable* v : variables) if(v != nullptr)
			delete v;
		for(Clause* c : clauses) if(c != nullptr)
			delete c;
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

	bool FactorGraph::isSAT(){
		for(Clause* c : clauses){
			if(!c->isSAT())
				return false;
		}
		return true;
	}

	bool FactorGraph::fix(Variable* var, int val){
		if(var->value) // si es 1 o -1, ya está asignado
			return true;

		var->value = val;
		unassigned_vars--;


		return simplify(var);
	}

	/**
	* Simpify after fixing var
	*/
	bool FactorGraph::simplify(Variable* var){
		// Para cada cláusula de la variable
		for(Literal* l : var->literals) if(l->enabled){
			// Se deshabilita la arista
			l->enabled = false;

			Clause* c = l->cl;
			c->unassigned_literals--;

			if(!c->satisfied){
				// Se comprueba si la asignación satisface la cláusula
				if(l->type == var->value){
					c->satisfied = true;
					for (Literal* l : c->literals)
						l->enabled = false;
				}
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
}
