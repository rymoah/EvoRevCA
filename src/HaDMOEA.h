#include "ecf/Algorithm.h"
#include "ecf/ECF.h"
#include "ecf/ECF_base.h"
#include <iostream>
#include <fstream>

#pragma once

class HaDMOEA : public Algorithm
{
public:
	std::vector <IndividualP> *parentPop;
	boost::shared_ptr<std::vector <std::vector <IndividualP>>> fronts;
	SelectionOperatorP selRandomOp, selWorstOp;
	
public:
	HaDMOEA();

	int numberOfNeighbours;
	bool initialize(StateP state);

	bool advanceGeneration(StateP state, DemeP deme);

	// it is always sorted from the lowest value to the highest,
    // regardless of whether in this particular case it is about maximization or minimization
	void quickSort(std::vector <IndividualP> *group, int left, int right, std::string prop, int objective);

	// sorts the population according to the fitness of the individual with respect to one selected target function
	void sortBasedOnProperty(std::vector <IndividualP>* deme, double* fMin, double* fMax, std::string prop, int objective);

	// checks the dominance between two multi-objective fitness
    // returns -1 if fitness1 dominates fitness2
    // returns 0 if neither dominates the other
    // returns 1 if fitness2 dominates fitness1
    // doesn't take into account any other things like 'nc' or 'crowding-distance', only the values of the goal functions are looked at 
	int checkDominance(MOFitnessP fitness1, MOFitnessP fitness2);

	// updates the rank field in the MOFitness of the individual 
	void nonDomSorting(boost::shared_ptr<std::vector <IndividualP>> pool, int N, boost::shared_ptr<std::vector <std::vector <IndividualP>>> fronts);

	// updates the crowding_distance field in the MOFitness of the individual
	void crowdedDistanceEst(StateP state, std::vector <IndividualP> *deme);
	void crowdedDistanceEst1(StateP state, std::vector <IndividualP> *deme);
	void crowdedDistanceEst3(StateP state, std::vector <IndividualP> *deme);


	void makeNewPop(StateP state, DemeP deme);

};
typedef boost::shared_ptr<HaDMOEA> HaDMOEAP;
