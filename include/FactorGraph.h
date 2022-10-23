#ifndef FACTOR_GRAPH_H_
#define FACTOR_GRAPH_H_

#include <iostream>
#include <vector>


using namespace std;

namespace sp{
	class Literal;
	
	class Clause{
	public:
		vector<Literal*> literals;
		int unassigned_literals;
		bool satisfied;
	
		// WalkSAT variables
		bool ws_unvisited;
		
		// Methods:
		Clause();
		Clause(const Clause &copy);
		bool isSAT();
	
		//WalkSAT methods
		/**
		 * @brief Number of true literals in the clause
		 * 
		 * @return int Number of true literals in the clause
		 */
		int NT();
	};
	
	
	class Variable{
	public:
		int id;
		int value; //-1, 0, 1 -- 0=unfixed
		
		double p; //product of (1-eta) of positive literals
		double m; //product of (1-eta) of negative literals
		
		int pzero; //equal to 1 if p equal to 0
		int mzero; //equal to 1 if m equal to 0
	
		// Biases -- equations (28, 29, 30)
		double wp; //plus
		double wz; //zero
		double wm; //minus
		
		vector<Literal*> literals; //
	
		// WalkSAT variables:
		bool ws; // Indica si la variable se usa en WS
		int ws_break;
		
		// Methods:
		Variable(int id);
		Variable(const Variable &copy);
	
		// WalkSAT Methods:
		/**
		 * @brief Obtiene las cláusualas que contienen un literal de la variable positivo
		 * 
		 * @return vector<Clause*> 
		 */
		vector<Clause*> TLC();
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
		Literal(const Literal &copy);
		void Disable();
	};
	
	class FactorGraph{
	public:
		vector<Literal*> literals;
		vector<Variable*> variables;
		vector<Clause*> clauses;
	
		int unassigned_vars;
	
		// Methods:
		FactorGraph(string file);
		FactorGraph(vector<vector<int>> &clauses);
		FactorGraph(const FactorGraph &copy);
		~FactorGraph();
		vector<string> splitString(string str);
		bool simplify(Variable* var);
		bool fix(Variable* var, int val, bool sp);
		bool fixUnitClause(Clause* c);
		bool isSAT();
	
	};
	
	struct comp
	{
	    bool operator()(Variable* lhs, Variable* rhs) const
	    {
			return lhs->id < rhs->id;
	    }
	};
}

#endif // FACTOR_GRAPH_H_
