#ifndef EvalOp_h
#define EvalOp_h

// multiobjective define directive
//#define MOEA


class EvalOp : public EvaluateOp
{
protected:
	StateP state_;

public:
	void registerParameters(StateP);
	FitnessP evaluate(IndividualP individual);
	std::vector<int> getTT(IndividualP individual, int& vars, uint genId = 0);
	bool initialize(StateP);

	uint nVariables;
	int variant;
	uint omega;
	int bestScore;
	std::vector< std::vector<bool> > inputMap;
	std::vector< std::string > inputNames;
	std::vector<bool> results;
	std::vector<int> tt;
	bool showTruth;
	enum Genotype { bitstring, tree };
	Genotype genotype;

};
typedef boost::shared_ptr<EvalOp> EvalOpP;

#endif
