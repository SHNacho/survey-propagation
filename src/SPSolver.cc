#include <algorithm>
#include <iostream>
#include "SPSolver.h"

using namespace std;

SPSolver::SPSolver(float alpha){
	rng = default_random_engine {};
	this->alpha = alpha;
}

double SPSolver::iterate(){
// Update surveys of clauses in a random permutation
	double eps=0.0, maxeps=0.0;
	
	vector<Clause*> clauses = this->fg->clauses;
	shuffle(clauses.begin(), clauses.end(), this->rng);

	for(Clause* c : clauses){
		eps = updateSurvey(c);	
		if(eps > maxeps)
			maxeps = eps;
	}

	return maxeps;
}

double SPSolver::updateSurvey(Clause* c){
	double m, p;
	double wt, wn;
	vector<double> prods;
	double allprod;
	int zeroes = 0;
	double eps;

	for(Literal* l : c->literals){
		Variable* var = l->var;
		if(var->value==0){
			if(l->type < 0){
				m = var->mzero?0:var->m;
				if(var->pzero==0)
					p = var->p / (1 - l->survey);
				else if(var->pzero == 1 && (1 - l->survey) < EPS)
					p = var->p;
				else
					p = 0;
				wn = p * (1 - m);
				wt = m;
			} else {
				p = var->pzero ? 0 : var->p;
				if(var->mzero == 0)
					m = var->m / (1 - l->survey);
				else if (var->mzero == 1 && (1 - l->survey) < EPS)
					m = var->m;
				else
					m = 0;
				wn = m * (1 - p);
				wt = p;
			}	
			
			double prod = wn / (wn + wt); 
			prods.push_back(prod);

			if(prod < EPS){
				if(++zeroes == 2)
					break;
			} else {
				allprod *= prod;
			}
		}
	}

	int i = 0;
	eps = 0.0;
	for(Literal* l : c->literals){
		Variable* var = l->var;
		if(var->value==0){
			double newsurvey = 0;
			if(!zeroes){
				newsurvey = allprod/prods[i];
			} else if(zeroes == 1 && prods[i] < EPS) {
				newsurvey = allprod;
			} else {
				newsurvey = 0;
			}

			if(l->type < 0){
				if(1-l->survey>EPS) {
					if(1-newsurvey>EPS) {
						var->p*=(1-newsurvey)/(1-l->survey);
					} else {
						var->p/=(1-l->survey);
						var->pzero++;
					} 
				} else {
					if(1-newsurvey>EPS) {
						var->p*=(1-newsurvey);
						var->pzero--;
					} 
				}
			} else {
				if(1 - l->survey > EPS){
					if(1 - newsurvey > EPS){
						var->m *= (1 - newsurvey) / (1 - l->survey);
					} else {
						var->m /= (1 - l->survey);
						var->mzero++;
					}
				} else {
					if(1 - newsurvey > EPS){
						var->m *= (1 - newsurvey);
						var->mzero--;
					}
				}
			}
			double diff = abs(l->survey - newsurvey);
      		if (eps < diff)
        		eps = abs(diff);

      		l->survey = newsurvey;
      		i++;
		}
	}
	return eps;
}

bool SPSolver::surveyPropagation(){
	double eps = 0;
	int iter = 0;
	computeSubProducts();
	do{
		eps = iterate();
	} while((eps > EPSILON) && (iter++ < ITERATIONS));
	if(eps <= EPSILON){
		return true;
	} else {
		return false;
	}
}


void SPSolver::computeSubProducts(){
	for (Variable* v : fg->variables) if (v->value == 0){
		v->p = 1;
		v->m = 1;
		v->pzero = 0;
		v->mzero = 0;

		for (Literal* l : v->literals) if(!l->cl->satisfied){
			if(l->type < 0){
				if(1 - l->survey > EPS){
					v->p *= 1 - l->survey;
				} else {
					v->pzero++;
				}
			} else {
				if(1 - l->survey > EPS){
					v->m *= 1 - l->survey;
				} else {
					v->mzero++;
				}
			}

		}
	}
}


void SPSolver::computeBias(Variable* var){
	double p, m;
	p = var->pzero ? 0 : var->p;
	m = var->mzero ? 0 : var->m;

	var->wz = p*m;
	var->wp = m - var->wz;
	var->wm = p - var->wz;

	// Normalización de los sesgos
	double norm = var->wm + var->wz + var->wp;
	var->wm /= norm;
	var->wp /= norm;
	var->wz /= norm;
}


void SPSolver::unitPropagation(){
	for (Clause* c : fg->clauses){
		if(c->unassigned_literals == 1 && !c->satisfied){
			fg->fixUnitClause(c);
		}
	}
}


bool SPSolver::surveyInspiredDecimation(){
	// Asignación aleatoria de las surveys
	std::uniform_real_distribution<float> distribution(0.0, 1.0);
	for(Literal* l : fg->literals){
		if(l->enabled){
			l->survey = distribution(rng);
		}
	}
	// Se cualculan cuantas variables se asignarán en cada paso
	int fixPerStep = fg->unassigned_vars * alpha;
	// Se realiza unitPropagation por si hay cláusulas unitarias
	unitPropagation();
	// Mientras que survey propagation llegue a un estado de convergencia
	while(surveyPropagation() && fg->unassigned_vars){
		double summag = 0;
		double maxmag;
		int unasignedVars = 0;
		
		// Se calcula el sesgo de cada la variable
		for(Variable* v : fg->variables) if(v->value == 0){
			computeBias(v);
			maxmag=v->wp > v->wm ? v->wp : v->wm;
			summag += maxmag;
			unasignedVars++;
		}
		
		// Se ordenan las variables en función del sesgo
		vector<Variable*> sortedVars = fg->variables;
		sort(sortedVars.begin(), sortedVars.end(), biasComparator);

		// Si se llega a un estado paramagnético, se resuelve por walksat
		if(summag/unasignedVars < PARAMAGNET){
			//WALKSAT
			cout << "WALKSAT" << endl;
			return true;
		}

		// Asignación de las variables ordenadas por bias
		while(fg->unassigned_vars && fixPerStep--){
			vector<Variable*>::iterator it = sortedVars.begin();
			while((*it)->value)
				it++;
			computeBias(*it);
			int val = (*it)->wp > (*it)->wm ? -1 : 1;
			if (!fg->fix(*it, val))
				return false;
		}
	}
	
	return true;
}


bool biasComparator(Variable* v1, Variable* v2){
	double b1 = fabs(v1->wm - v1->wp);
	double b2 = fabs(v2->wm - v2->wp);
	return (b1 > b2);
}












