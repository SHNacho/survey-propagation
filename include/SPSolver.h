#ifndef SP_SOLVER_H_
#define SP_SOLVER_H_

#include <iostream>
#include <random>
#include <vector>

#include "FactorGraph.h"

using namespace std;

namespace sp{
	class SPSolver{
	public:
	
		// SP constants
		const float PARAMAGNET = 0.01;
		const float EPSILON = 0.01;
		const double EPS = 1.0e-16;
		const int ITERATIONS = 1000;
	
		// WalkSAT constants
		const double WS_NOISE = 0.567;
		const int WS_MAX_TRIES = 10;
		const int WS_MAX_STEPS = 1000;
	
		default_random_engine rng;
	
		FactorGraph* fg;
		vector< pair<int, int> > sol;
	
		float alpha;	// Fracción de variables a asignar
		int SPIter;		// Número de iteraciones de SP
	
		// Methods:
		SPSolver(FactorGraph* fg, float alpha);
		/**
		 * @brief Itera sobre todas las cláusulas, actualizando
		 * las surveys de cada una de sus variables
		 *
		 * @return Valor de la survey más alta
		*/
		double iterate();
		/**
		 * @brief Actualiza todas las surveys y productos de (1 - survey)
		 * de una cláusula
		 *
		 * @return Valor de la mayor diferencia entre la survey anterior y la
		 * nueva calculada. max(|new_survey - old_survey|)
		*/
		double updateSurvey(Clause* c); //update_eta 
		/**
		 * @brief Actualiza de forma iteratica las surveys hasta que estas convergen a EPSILON
		 * o hasta que se alcanza el número máximo de iteraciones
		 * 
		 * @return true Si convergen 
		 * @return false Si no convergen
		 */
		bool surveyPropagation(); //Converge
		/**
		 * @brief Calcula los productos (1 - survey) de todas las variables
		 */
		void computeSubProducts(); //compute_pi
		/**
		 * @brief Calcula los sesgos de una variable
		 * 
		 * @param var Variable sobre la que calcular los sesgos
		 */
		void computeBias(Variable* var); //compute_fields
		/**
		 * @brief Comprueba y satisface las cláusulas unitarias
		 */
		void unitPropagation();
		/**
		 * @brief 
		 * 
		 * @return true 
		 * @return false 
		 */
		bool surveyInspiredDecimation();
		/**
		 * @brief Asigna un conjunto de variables
		 * 
		 * @param quant Cantidad de variables a asignar
		 * @return true Si se han asignado correctamente
		 * @return false Si ha habido alguna contradicción
		 */
		bool fixChunk(int quant, vector<Variable*> vars);
		/**
		 * @brief Ejecuta el algoritmo WalkSat sobre la fórmula
		 * 
		 */
		bool WalkSat();
	
		/**
		 * @brief Selecciona una variable de una cláusula a la que 
		 * modificar su valor asignado
		 * 
		 * @param c Cláusula sobre la que se escoge la variable
		 * @return Variable* Puntero a la variable escogida
		 */
		Variable* pickVar(Clause* c);
	};
	
	bool biasComparator(Variable* v1, Variable* v2);
}


#endif // SP_SOLVER_H_

