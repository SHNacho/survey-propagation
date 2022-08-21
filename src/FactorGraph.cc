#include "FactorGraph.h"
#include <iostream>
#include <unordered_set>

using namespace std;

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
	for(Literal* l : var->literals){
		// Se deshabilita la arista
		l->enabled = false;

		Clause* c = l->cl;
		c->unassigned_literals--;

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















