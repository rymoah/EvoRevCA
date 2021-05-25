#include "ecf/ECF.h"
#include "HaDMOEA.h"
#include "ecf/SelRandomOp.h"
#include "ecf/SelWorstOp.h"


HaDMOEA::HaDMOEA() {
	name_ = "HaDMOEA";
	// create selection operators needed
	// in this case, SelRandomOp and SelWorstOp
	selRandomOp = static_cast<SelectionOperatorP> (new SelRandomOp);

	selWorstOp = static_cast<SelectionOperatorP> (new SelWorstOp);
	this->parentPop = new std::vector <IndividualP>();
	this->fronts = boost::shared_ptr<std::vector <std::vector <IndividualP>>> (new std::vector <std::vector <IndividualP>>());

}


bool HaDMOEA::initialize(StateP state)
{
	// initialize all operators
	selRandomOp->initialize(state);
	selWorstOp->initialize(state);
	numberOfNeighbours = 5;
	return true;
}

void HaDMOEA::quickSort(std::vector <IndividualP> *group, int left, int right, std::string prop, int objective = -1) {
	int i = left, j = right;
	
	IndividualP tmp;
	MOFitnessP pivotMO =  boost::static_pointer_cast<MOFitness> (group->at((left + right) / 2)->fitness);
	double pivot =  pivotMO->getProperty(prop, objective);
 
	/* partition */
	while (i <= j) {
		MOFitnessP moFitnessI =  boost::static_pointer_cast<MOFitness> (group->at(i)->fitness);
		while (moFitnessI->getProperty(prop, objective) <  pivot) {
			i++;
			moFitnessI =  boost::static_pointer_cast<MOFitness> (group->at(i)->fitness);
		}
		MOFitnessP moFitnessJ =  boost::static_pointer_cast<MOFitness> (group->at(j)->fitness);
		while (pivot < moFitnessJ->getProperty(prop, objective)) {				
			j--;
			moFitnessJ =  boost::static_pointer_cast<MOFitness> (group->at(j)->fitness);
		}

		if (i <= j) {
			tmp = group->at(i);
			group->at(i) = group->at(j);
			group->at(j) = tmp;
			i++;
			j--;
		}
	};
 
	/* recursion */
	if (left < j)
		quickSort(group, left, j, prop, objective);
	if (i < right)
		quickSort(group, i, right, prop, objective);
}



void HaDMOEA::sortBasedOnProperty(std::vector <IndividualP>* deme, double* fMin, double* fMax, std::string prop, int objective = -1) {
	int left = 0;
	int right = deme->size()-1;
	quickSort(deme, left, right, prop, objective);

	MOFitnessP fitness = boost::static_pointer_cast<MOFitness> (deme->at(left)->fitness);
	*fMin = fitness->getProperty(prop, objective);
	fitness = boost::static_pointer_cast<MOFitness> (deme->at(right)->fitness);
	*fMax = fitness->getProperty(prop, objective);
}

int HaDMOEA::checkDominance(MOFitnessP fitness1, MOFitnessP fitness2) {
	double eps = 1E-9;
	uint size = fitness1->size();
	if (size != fitness2->size()) {
		// it must not happen 
	}
	int dominance = 0;
	for (uint i = 0; i<size; i++) {
		FitnessP f1 = fitness1->at(i);
		FitnessP f2 = fitness2->at(i);
		if (abs(f1->getValue() - f2->getValue()) > DBL_EPSILON) {
			if (dominance == 0) {
				if (f1->isBetterThan(f2)) {
					dominance = -1;
				} else {
					dominance = 1;
				}
			} else if (dominance == -1) {
				if (f1->isBetterThan(f2)) {
					// dominance stays -1 
				} else {
					return 0;
				}
			} else {
				if (f1->isBetterThan(f2)) {
					return 0;
				} else {
					// dominance stays -1 
				}
			}
		}
	}
	return dominance;
}


// O(M * N^2),
// M => number of goal functions 
// N population size for sorting 
void HaDMOEA::nonDomSorting(boost::shared_ptr<std::vector <IndividualP>> pool, int N, boost::shared_ptr<std::vector <std::vector <IndividualP>>> fronts)
{
	for(int i = 0; i<fronts->size(); i++){
		fronts->at(i).clear();
	}
	fronts->clear();

	std::vector <IndividualP> Q /*= *(new std::vector <IndividualP>())*/;
	int collectedSoFar = 0;
	int p = 1; // lowest (best) rank 

	// calculate nc and Sp for each individual in the population 
	// and must go all the way to the end (not to the pool-> size () - 1) because 
	// in case the whole population is in the first front, the last one must be included in the first front 
	for (uint i = 0; i < pool->size(); i++) {
		IndividualP ind = pool->at(i);
		MOFitnessP fitnessI =  boost::static_pointer_cast<MOFitness> (ind->fitness);

		if (i == 0) {
			fitnessI->nc = 0;
			fitnessI->Sp = new std::vector<IndividualP>();
			fitnessI->rank = 0;
		}

		for (uint j = i+1; j < pool->size(); j++) {
			IndividualP other = pool->at(j);
			MOFitnessP fitnessJ = boost::static_pointer_cast<MOFitness> (other->fitness);

			if (i == 0) {
				fitnessJ->nc = 0;
				fitnessJ->Sp = new std::vector<IndividualP>();
				fitnessJ->rank = 0;
			}

			int dominance = checkDominance(fitnessI, fitnessJ);
			if (dominance == -1) {
				fitnessI->Sp->push_back(other);
				fitnessJ->nc++;
			} else if (dominance == 1) {
				fitnessI->nc++;
				fitnessJ->Sp->push_back(ind);
			}
		}
		
		// initially, pareto optimal solutions enter the set Q
		// all solutions that have a 'domination count' (nc) equal to 0 
		if (fitnessI->nc == 0) {
			fitnessI->rank = p;
			Q.push_back(ind);
			collectedSoFar++;
		}
	}




	// odrediti rankove rjesenjima
	//while (Q.size() != 0 && collectedSoFar < N) {
	while (Q.size() != 0) {
		fronts->push_back(Q);

		
		p++;
		// if the set Q is empty, it means that we have assigned a rank to all solutions 
		// in the set newQ the solutions of the next front of non-domination are collected 
		// means only those solutions that are dominated exclusively by solutions from the set Q 
		std::vector <IndividualP> newQ; 
		for (uint i=0; i<Q.size(); i++) {
			//pool->push_back(Q.at(i));
			MOFitnessP fitnessI = boost::static_pointer_cast<MOFitness> (Q.at(i)->fitness);

			for (uint j=0; j<fitnessI->Sp->size(); j++) {
				// for each solution we go through a set of solutions over which it dominates and update nc 
				IndividualP dominated = fitnessI->Sp->at(j);
				MOFitnessP fitnessJ = boost::static_pointer_cast<MOFitness> (dominated->fitness);
				fitnessJ->nc--;
				if (fitnessJ->nc == 0) {
					fitnessJ->rank = p;
					newQ.push_back(dominated);
					collectedSoFar++;
				}
			}
		}
		Q.clear();
		Q = newQ;
	}

// O(M * N * logN)
void HaDMOEA::crowdedDistanceEst(StateP state, std::vector <IndividualP> *deme) {
	MOFitnessP fitness =  boost::static_pointer_cast<MOFitness> (deme->at(0)->fitness);
	int objCount = fitness->size();
	for (uint i = 0; i<objCount; i++) {

		double fMin;
		double fMax;

		// we sort the entire population of solutions according to one goal 'i'
		// population is sorted so that the results of the objective function 'and' are ranked from smallest to largest,
		// whether the objective of eandf is about maximization or minimization
		sortBasedOnProperty(deme, &fMin, &fMax, "objective", i);

		// now the solutions are sorted, the zero index is the best 
		// the last one is the worst 
		for (uint j = 0; j<deme->size(); j++) {
			fitness =  boost::static_pointer_cast<MOFitness> (deme->at(j)->fitness);


			int l = j-1,r=j+1;
			double fl, fr;
			MOFitnessP fitnessLeft, FitnessRight;

			std::vector<double> distances;

			for(int k = 0; k<numberOfNeighbours; k++){
				if(l<0){
					fl = std::numeric_limits<double>::infinity();
				} else{
					fitnessLeft =  boost::static_pointer_cast<MOFitness> (deme->at(l)->fitness);
					fl = fitnessLeft->getValueOfObjective(i);
					//l--;
				}

				if(r>=deme->size()){
					fr = std::numeric_limits<double>::infinity();
				} else {
					FitnessRight =  boost::static_pointer_cast<MOFitness> (deme->at(r)->fitness);
					fr = FitnessRight->getValueOfObjective(i);
					//r++;
				}

				double dl = abs(fitness->getValueOfObjective(i)-fl);
				double dr = abs(fitness->getValueOfObjective(i)-fr);
				//std::cout<<"dl: " <<dl;
				//std::cout<<"dr: " <<dr;
				double distance;

				if(dl<dr){
					l--;
					distance = dl;
				} else {
					r++;
					distance = dr;
				}

				distance /= fMax - fMin;
				distances.push_back(distance);
			}


			double harmonicSum = 0;
			for(int k=0; k<distances.size(); k++){
				//if(distances.at(k)==0){
				//	harmonicSum += 10000;
				//} else 
				harmonicSum += 1/distances.at(k);
			}
			
		//	std::cout<<"Harmonic sum: " <<numberOfNeighbours/harmonicSum<< std::endl;

			fitness->crowding_distance += numberOfNeighbours/harmonicSum;
		}

	}

}


// O(M * N * logN)
void HaDMOEA::crowdedDistanceEst3(StateP state, std::vector <IndividualP> *deme) {
	MOFitnessP fitness =  boost::static_pointer_cast<MOFitness> (deme->at(0)->fitness);
	MOFitnessP fitnessj =  boost::static_pointer_cast<MOFitness> (deme->at(0)->fitness);	
	int objCount = fitness->size();

	std::vector<double> mins;
	std::vector<double> maxs;

	for(int k=0; k<objCount; k++){
		MOFitnessP fit =  boost::static_pointer_cast<MOFitness> (deme->at(0)->fitness);
		mins.push_back(fit->getValueOfObjective(k));
		maxs.push_back(fit->getValueOfObjective(k));
	}
	

	for(int i =0; i< deme->size(); i++){
		MOFitnessP fit =  boost::static_pointer_cast<MOFitness> (deme->at(i)->fitness);
		for(int k = 0; k<objCount; k++){
			if(fit->getValueOfObjective(k)<mins.at(k)){
				mins[k] = fit->getValueOfObjective(k);
			}
			if(fit->getValueOfObjective(k)>maxs.at(k)){
				maxs[k] = fit->getValueOfObjective(k);
			}
		}
	}
	

	
	std::vector<std::vector<double>> distanceMatrix;

	//calculate distance matrix
	
	for(int i = 0; i < deme->size(); i++){
		std::vector<double> row;
		fitness =  boost::static_pointer_cast<MOFitness> (deme->at(i)->fitness);

		for(int j = 0; j < deme->size(); j++){
			fitnessj =  boost::static_pointer_cast<MOFitness> (deme->at(j)->fitness);
			
			double distance = 0;

			if(i==j){
				row.push_back(distance);
				continue;
			}

			for(int k=0; k<objCount; k++){
				distance+= pow((fitness->getValueOfObjective(k)-fitnessj->getValueOfObjective(k))/(maxs[k]-mins[k]),2);
			}
			distance = sqrt(distance);
			row.push_back(distance);

		}
		distanceMatrix.push_back(row);

	}

	for(int i = 0; i < deme->size(); i++){
		fitness =  boost::static_pointer_cast<MOFitness> (deme->at(i)->fitness);

		std::vector<double> row = distanceMatrix.at(i);
		std::sort(row.begin(), row.end());
		
		double harmonicSum = 0;

		for(int k = 0; k<numberOfNeighbours; k++){
			harmonicSum += 1.0/(row.at(k));
		}	

		fitness->crowding_distance = (float)numberOfNeighbours/harmonicSum;
	}
	
	

}

// performs selection, crossover, mutation and creates a new population of the same size
// it is currently a steadyStateTournament with a size 3 tournament
void HaDMOEA::makeNewPop(StateP state, DemeP deme) {
	for(uint iIter = 0; iIter < deme->size(); iIter++) {
		
		ECF_LOG(state, 5, "Individuals in tournament: ");

		std::vector<IndividualP> tournament = *(new std::vector<IndividualP>());
		for (uint i = 0; i < 2; ++i) {
			// select a random individual for the tournament
			tournament.push_back(selRandomOp->select(*deme));
			ECF_LOG(state, 5, uint2str(i) + ": " + tournament[i]->toString());
		}

		// select the worst from the tournament
		IndividualP worst = selWorstOp->select(tournament);
		ECF_LOG(state, 5, "The worst from the tournament: " + worst->toString());

		// remove pointer to 'worst' individual from vector 'tournament'
		removeFrom(worst, tournament);

		IndividualP myMate = selRandomOp->select(*deme);

		// crossover the first two (remaining) individuals in the tournament
		mate(tournament[0], myMate, worst);

		// perform mutation on new individual
		mutate(worst);

		// create new fitness
		evaluate(worst);
		ECF_LOG(state, 5, "New individual: " + worst->toString());
	}

}


bool HaDMOEA::advanceGeneration(StateP state, DemeP deme)
{
	static bool firstGen = true;

	// ne prvi put
	if(!firstGen)
		this->makeNewPop(state, deme);
	firstGen = false;

	int N = deme->size();
	bool initialGeneration = parentPop->size() == 0;
	for (uint i = 0; i<parentPop->size(); i++) {
		deme->push_back((IndividualP) parentPop->at(i)->copy());
	}

	fronts->clear();
	//int lastFront;
	//nonDomSorting(deme, N, fronts, &lastFront);
	nonDomSorting(deme, N, fronts);
	//crowdedDistanceEst(state, deme, lastFront);
	deme->clear();

	uint i = 0;
	uint size = 0;
	while (size + fronts->at(i).size() < N) {
		//crowdedDistanceEst(state, &(fronts->at(i)));
		for (uint j = 0; j<fronts->at(i).size(); j++) {
			deme->push_back(fronts->at(i).at(j));
		}
		size += fronts->at(i).size();
		i++;
		std::cout<<"i: "<<i<<std::endl;
	}

	if (!initialGeneration) {
		double fMin;
		double fMax;
		std::vector<IndividualP> combined;

		for(int front = 0; front<= i; front++){
			for(int ind =0; ind<fronts->at(front).size(); ind++){
				combined.push_back(fronts->at(front).at(ind));
			}
		}

		crowdedDistanceEst3(state, &combined);
		sortBasedOnProperty(&(fronts->at(i)), &fMin, &fMax, "crowding_distance");

		int howMany = N - deme->size();
		int lastFrontSize = fronts->at(i).size();
		for (int j = lastFrontSize-1; j >= (lastFrontSize - howMany); j--) {
			deme->push_back(fronts->at(i).at(j));
		}
	}

	// we save the current population to the parent population 
	parentPop->clear();
	for (uint i = 0; i<deme->size(); i++) {
		this->parentPop->push_back((IndividualP) deme->at(i)->copy());
	}

	// moved to the beginning of the generation 
	//this->makeNewPop(state, deme);


	return true;
}
