#include<iostream>
#include "ecf/ECF.h"
#include "EvalOp.h"
#include "Primitives.cpp"
#include "WriteTT.h"



#ifndef MOEA

int main(int argc, char** argv)
{
	StateP state (new State);

	// crate evaluation operator
	EvalOp* e = new EvalOp;
	
	state->setEvalOp(e);

	// create tree genotype
	TreeP tree (new Tree::Tree);

	// create new functions and add them to function set
	Tree::PrimitiveP ifl (new If);
	tree->addFunction(ifl);
	Tree::PrimitiveP orp (new Or);
	tree->addFunction(orp);
	Tree::PrimitiveP andp (new And);
	tree->addFunction(andp);
	Tree::PrimitiveP notp (new Not);
	tree->addFunction(notp);
	Tree::PrimitiveP xorp (new Xor);
	tree->addFunction(xorp);
	Tree::PrimitiveP and2p (new And2);
	tree->addFunction(and2p);
	Tree::PrimitiveP xnorp (new XNor);
	tree->addFunction(xnorp);


	// custom type terminals
	for(uint i = 0; i < 20; i++) {
		Tree::PrimitiveP myTerm = (Tree::PrimitiveP) new BoolV;
		std::string name = "v" + uint2str(i);
		myTerm->setName(name);
		tree->addTerminal(myTerm);
	}

	// register genotype with our primitives
	state->addGenotype(tree);

	WriteTT *b = new WriteTT;
	state->addOperator((OperatorP) b);

	// initialize and start evaluation
	if(!state->initialize(argc, argv))
		return 1;


	// read individual from file in 3rd agrument
	// (to check constructions on different number of variables)
	if(argc == 3) {
		XMLNode xDeme = XMLNode::parseFile(argv[2], "Deme");
		uint iNode = 0;
		while(1) {
			XMLNode xInd = xDeme.getChildNode("Individual", iNode);
			if(xInd.isEmpty())
				break;
			iNode++;

			IndividualP ind = (IndividualP) new Individual(state);
			ind->read(xInd);

			evaluateVerbose = true;
			e->showTruth = true;
			ind->fitness = state->getAlgorithm()->evalOp_->evaluate(ind);
			//cout << ind->toString();
		}

		//for(uint eval = 0; eval < 0; eval++) {
		//	EvalOp* ev = (EvalOp*) (e->evaluators[eval]->getAlgorithm()->evalOp_.get());
		//	ev->showTruth = true;
		//	ev->evaluate(ind);
		//}
		return 0;
	}

	state->run();


	// after the evolution: show best
	//{
	//	std::vector<IndividualP> hof = state->getHoF()->getBest();
	//	IndividualP ind = hof[0];
	//	std::cout << ind->toString();
	//	e->showTruth = true;
	//	state->getEvalOp()->evaluate(ind);
	//}


	// load default logfile, read tree, add TT to same file 
	if(argc < 3)
		return 0;
	XMLNode xHof = XMLNode::parseFile(argv[2], "HallOfFame");
	XMLNode xInd = xHof.getChildNode();
	IndividualP ind = (IndividualP) new Individual(state);
	ind->read(xInd);
	state->getAlgorithm()->evalOp_->evaluate(ind);
	std::stringstream ss;
	//for(uint i = 0; i < e->tt.size(); i++)
	//	ss << e->tt[i];// << " ";
	//ss << std::endl;
	std::ofstream out(argv[2], std::ios_base::app);
	out << ss.str();

	return 0;
}

#endif



// multiobjective version
#ifdef MOEA

#include "HaDMOEA.h"

int main(int argc, char **argv)
{
	StateP state (new State);

	EvalOp* e = new EvalOp;
	state->setEvalOp(e);

	// create tree genotype
	TreeP tree (new Tree::Tree);

	// create new functions and add them to function set
	Tree::PrimitiveP ifl (new If);
	tree->addFunction(ifl);
	Tree::PrimitiveP orp (new Or);
	tree->addFunction(orp);
	Tree::PrimitiveP andp (new And);
	tree->addFunction(andp);
	Tree::PrimitiveP notp (new Not);
	tree->addFunction(notp);
	Tree::PrimitiveP xorp (new Xor);
	tree->addFunction(xorp);
	Tree::PrimitiveP and2p (new And2);
	tree->addFunction(and2p);
	Tree::PrimitiveP xnorp (new XNor);
	tree->addFunction(xnorp);


	// custom type terminals
	for(uint i = 0; i < 20; i++) {
		Tree::PrimitiveP myTerm = (Tree::PrimitiveP) new BoolV;
		std::string name = "v" + uint2str(i);
		myTerm->setName(name);
		tree->addTerminal(myTerm);
	}

	// custom type functions
	for(uint i = 0; i < 20; i++) {
		Tree::PrimitiveP myTerm = (Tree::PrimitiveP) new BoolV;
		std::string name = "f" + uint2str(i);
		myTerm->setName(name);
		tree->addTerminal(myTerm);
	}

	// register genotype with our primitives
	state->addGenotype(tree);

	WriteTT *b = new WriteTT;
	state->addOperator((OperatorP) b);

	// initialize and start evaluation
	if(!state->initialize(argc, argv))
		return 1;

	state->run();

	std::stringstream dirName;
	dirName << "mkdir pareto";
	system(dirName.str().c_str());
		
	//
	// print the first Pareto front in a separate file:
	//
	if(true){
		ECF_LOG(state, 1, "Best in " + uint2str(state->getGenerationNo()));
		//state->getAlgorithm()
		HaDMOEAP alg(new HaDMOEA());
		
		alg->nonDomSorting(state->getPopulation()->getLocalDeme(), state->getPopulation()->getLocalDeme()->size(),alg->fronts);
		double mi, ma;
		alg->sortBasedOnProperty(&alg->fronts->at(0),&mi, &ma,"objective",0);
			
		std::ofstream frontValues;
		frontValues.open (".\\pareto\\paretoValues.txt");
		std::ofstream front;
		front.open(".\\pareto\\paretoFront.txt");
		for(int i=0; i< alg->fronts->at(0).size(); i++){
				front << alg->fronts->at(0).at(i)->toString();
				MOFitnessP fitness = boost::static_pointer_cast<MOFitness> (alg->fronts->at(0).at(i)->fitness);
				for (uint f = 0; f < fitness->size(); f++)
						frontValues << fitness->at(f)->getValue() << "\t";
				frontValues << endl;
		}
		front.close();
		frontValues.close();
	}

	return 0;
}

#endif
