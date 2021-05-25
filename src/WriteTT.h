#ifndef Write_TT
#define Write_TT


//
// auxiliary operator to print the best individual
//

bool evaluateVerbose;


class WriteTT : public Operator
{
private:
	StateP state_;
	uint freq, maxGen, eval;
	double fitnessVal;
	bool termEval, termFitness, termStagn, termGen;

public:

	bool initialize(StateP state)
	{
		termEval = termFitness = termStagn = termGen = false;

		if(state->getRegistry()->isModified("term.eval"))
			termEval = true;
		if(state->getRegistry()->isModified("term.stagnation"))
			termStagn = true;
		if(state->getRegistry()->isModified("term.fitnessval"))
			termFitness = true;
		if(state->getRegistry()->isModified("term.maxgen"))
			termGen = true;

		state_ = state;
		voidP sptr = state->getRegistry()->getEntry("term.stagnation");
		freq = *((uint*)sptr.get());

		sptr = state->getRegistry()->getEntry("term.maxgen");
		maxGen = *((uint*)sptr.get());

		sptr = state->getRegistry()->getEntry("term.eval");
		eval = *((uint*)sptr.get());

		sptr = state->getRegistry()->getEntry("term.fitnessval");
		fitnessVal = *((double*)sptr.get());

		return true;
	}


	bool operate(StateP state)
	{
		bool output = false;
		uint currentGen = state->getGenerationNo();
		double minFitness = state->getStats()->getFitnessMin();
		double maxFitness = state->getStats()->getFitnessMax();

		if(termGen && currentGen >= maxGen)
			output = true;

		if(termStagn && (currentGen - state->getPopulation()->getHof()->getLastChange()) > freq)
			output = true;

		if(termEval && state->getEvaluations() >= eval)
			output = true;

		if(termFitness && (fitnessVal >= minFitness && fitnessVal <= maxFitness))
			output = true;


		if (output) {
			//std::cout << state->getGenerationNo() << endl;
			//std::cout << state->getTerminateCond() << endl;
			//std::cout << state->getEvaluations() << endl;

			std::vector<IndividualP> hof = state->getHoF()->getBest();
			IndividualP ind = hof[0];

			evaluateVerbose = true;
			state->getAlgorithm()->evaluate(ind);
			evaluateVerbose = false;

			ECF_LOG(state, 1, ind->toString());
		}

		return true;
	}
};

#endif
