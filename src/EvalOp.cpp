#include <cmath>
#include <ecf/ECF.h>
#include "EvalOp.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <iostream>
#include <string>
using namespace std;

// from WriteTT.h
extern bool evaluateVerbose;


/**
    * Check the compatibility of two landscapes L1 and L2, having the same
    * width d and the same center omega. L1 and L2 are incompatible if there
    * exists a position i between 0 and d-1 such that L1 has a 0 and L2 as a 1,
    * or vice versa.
    * @param l1    First landscape, represented as an int array
    * @param l2    Second landscape, represented as an int array
    * @return      true if the l1 and l2 are compatible, false otherwise
    */
bool checkCompatibility(std::vector<int>& l1, std::vector<int>& l2) {
        
    int d = l1.size();
    bool compatible=true;
    int i=0;
        
    while((i<d) && (compatible)) {
            
        //The condition "if(l1[i]=0 AND l2[i]=1) OR (l1[i]=1 AND l2[i]=1)
        //can be compactly described as "if(l1[i]+l2[1]=1)"
        if(l1[i]+l2[i] == 1) {
                
            //A mismatch has been found, hence the landscapes are incompatible
            compatible = false;
                
        } else {
                
            i++;
                
        }
            
    }
        
    return compatible;
}


/**
    * Check the compatibility of a landscape against a set of other landscapes
    * with the same width d and center omega.
    * 
    * @param landscape     an int array representing a landscape
    * @param setland       an int matrix representing a set of landscapes
    * @return              The number of elements of setland that landscape is
    *                      compatible with
    */
int checkCompatibilitySet(std::vector<int>& landscape, std::vector< std::vector<int> >& setland) {
        
    int ncomp = 0;
        
    for(int i=0; i<(int)setland.size(); i++) {
            
        if(checkCompatibility(landscape, setland[i])) {
            ncomp++;
        }
            
    }
        
    return ncomp;
}


/**
    * Build the neighborhood landscapes of a given landscape L.
    * 
    * @param L     an array of int representing a complementing landscape
    * @param omega the center of the landscape L
    * @return      an int matrix of size d*d where each row represents a
    *              neighborhood landscape of L. In particular, the landscape
    *              in position omega corresponds to L.
    * 
    */
void tabulateLandscape(std::vector<int>& L, int omega, std::vector< std::vector<int> >& M) {
        
    int d = L.size();           //diameter of the local rule/width of the landscape
    // original in Java: int[][] M = new int[d][d];  //matrix holding all neighborhood landscapes
	M.resize(d);
	for(int i = 0; i < d; i++)
		M[i].resize(d);
        
    //Initialization: Copy L in the row omega of M, and fill the upper and
    //lower triangular parts of the matrix with don't cares
    // original in Java: System.arraycopy(L, 0, M[omega], 0, d);
	for(int i = 0; i < d; i++)
		M[omega][i] = L[i];
        
    for(int i=0; i<omega; i++) {
            
        for(int j=0; j<omega-i; j++) {
                
            M[i][j] = 2;
                
        }
            
    }
        
    for(int i=d-1; i>omega; i--) {
            
        for(int j=d-1; j>(d-1+omega-i); j--) {
                
            M[i][j] = 2;
                
        }
    }
        
    //1st phase: build the matrix above row omega by shifting L to the right
    for(int i=1; i<=omega; i++) {
            
        for(int j=i; j<d; j++) {
                
            //There are three cases: if we are in column omega, then we have
            //to fill the position with the origin symbol *. Else, if we
            //are in column omega+i, we have to fill the position with a
            //don't care -, since that is the position aligned with origin
            //of L. Otherwise, we fill the position with the symbol of L
            //shifted i places on the left
            if(j==omega) {
                    
                M[omega-i][j] = 3;
                    
            } else {
                    
                if(j==omega+i) {
                        
                    M[omega-i][j] = 2;
                        
                } else {
                        
                    M[omega-i][j] = L[j-i];
                        
                }
                        
            }
                
        }
            
    }
        
    //2nd phase: build the matrix below row omega by shifting L to the left
    for(int i=1; i<d-omega; i++) {
            
        for(int j=0; j<d-i; j++) {
                
            //Again, there are three cases: for j=omega we put the origin
            //symbol, for j=omega-i we put the don't care, and in all other
            //cases we shift the value of L by i places on the right
            if(j==omega) {
                    
                M[omega+i][j] = 3;
                    
            } else {
                    
                if(j==omega-i) {
                        
                    M[omega+i][j] = 2;
                        
                } else {
                        
                    M[omega+i][j] = L[j+i];
                        
                }
                    
            }
                
        }
            
    }
        
    //return M;
}


/**
    * Compute the compatibility fitness of a set S of landscapes with the same
    * width d and center omega. For each element S[i] of S, build the neighborhood
    * landscapes with the tabulateLandscapes() method. Then, check the compatibility
    * of each neigborhood landscape (except the one in position omega, which
    * corresponds to S[i]) against set S with checkCompatibilitySet(), and add
    * the result of checkCompatibilitySet() to the fitness value.
    * @param setland   an int matrix representing a set of landscapes
    * @param omega     the center of the landscapes
    * @return          The compatibility fitness value of the landscapes in
    *                  setland. A fitness value of 0 means that the landscapes
    *                  are mutually incompatible (hence, they give rise to a
    *                  locally invertible marker CA).
    */
int compFitnessLandscapes(std::vector< std::vector<int> >& setland, int omega) {
        
    int fitness = 0;
        
    for(int i=0; i<(int)setland.size(); i++) {
            
        //Build the neighborhood landscapes matrix of the landscape setland[i]
        static std::vector< std::vector<int> > nbrlandscapes;
		tabulateLandscape(setland[i], omega, nbrlandscapes);
            
        //Check compatibility of each neighborhood landscape (except row
        //omega) against setland, and add the result to the fitness value
        for(int j=0; j<(int)nbrlandscapes.size(); j++) {
                
            if(j != omega) {
                    
                fitness += checkCompatibilitySet(nbrlandscapes[j], setland);
                    
            }
                
        }
            
    }
        
    return fitness;
        
}


void bin2Int(std::vector<bool>& binstring, std::vector<int>& toRet) {
       
    toRet.resize(binstring.size());
       
    for(int i=0; i<(int)binstring.size(); i++) {
           
        if(binstring[i]) {
            toRet[i] = 1;
        } else {
            toRet[i] = 0;
        }
           
    }
       
    //return toRet;
}


/**
    * Converts a decimal number in a binary string (int version).
    * 
    * @param   dNum    a decimal number
    * @param   length  the length of the binary string necessary to hold dNum
    * @return  bNum    the conversion of dNum as a binary string
    */
void dec2BinMod(int dNum, int length, std::vector<bool>& bNum) {
        
    bNum.resize(length);
    int temp = dNum;
    int i = 0;

    while(temp!=0) {

        int mod = temp%2;
        temp = temp/2;

        if(mod==1) {
            bNum[i] = true;
        }

        i++;
    }

    //return bNum;
}


/**
    * Given the truth table of a generating function, build the corresponding
    * set of complementing landscapes defining the marker CA rule. This is
    * accomplished by first computing the support of the generating function.
    * Then, each vector in the support is converted to the corresponding
    * landscape by inserting the origin symbol * in position omega.
    * 
    * @param ttable    A boolean array specifying the truth table of a
    *                  generating function
    * @param d         diameter of the CA local rule (the generating function
    *                  has d-1 variables)
    * @param omega     center of the landscapes
    * @return          An int matrix representing the set of complementing
    *                  landscapes that specifies the CA local rule
    */
void convertGenFuncToLandscapes(std::vector<bool>& ttable, int d, int omega, std::vector< std::vector<int> >& setland) {
        
    //Step 1: get the support of the generating function
    std::vector<int> nonzeropos;
    for(int i=0; i<(int)ttable.size(); i++) {
        if(ttable[i]) {
            nonzeropos.push_back(i);
        }
    }
    //nonzeropos.trimToSize();
    // original in Java: boolean[][] support = new boolean[nonzeropos.capacity()][];
	std::vector< std::vector<bool> > support(nonzeropos.size());
    for(int i=0; i<(int)support.size(); i++) {
        dec2BinMod(nonzeropos[i], d-1, support[i]);
    }
        
    //Step 2: convert each vector of the support in the corresponding
    //landscape with the origin symbol * in position omega
    // original in Java: int[][] setland = new int[support.length][d];
	setland.resize(support.size());
	for(uint i = 0; i < setland.size(); i++)
		setland[i].resize(d);

    for(int i=0; i<(int)setland.size(); i++) {
            
        //Step 2a: convert the boolean vector in a int vector of 0s and 1s
        std::vector<int> intvect;
		bin2Int(support[i], intvect);
            
        //Step 2b: copy the values of the int vector in the i-th landscape,
        //except for position omega which is set to the origin symbol
        for(int j=0; j<d; j++) {
            if(j<omega) {
                setland[i][j] = intvect[j];
            } else {
                if(j==omega) {
                    setland[i][j] = 3;
                } else {
                    setland[i][j] = intvect[j-1];
                }
            }
        }
            
    }
        
    //return setland;
}



void EvalOp::registerParameters(StateP state)
{
	state->getRegistry()->registerEntry("variables", (voidP) (new uint(1)), ECF::UINT);
	state->getRegistry()->registerEntry("variant", (voidP) (new int(0)), ECF::INT, "0 - SOGP, 1 - lexGP");
	state->getRegistry()->registerEntry("omega", (voidP) (new uint(1)), ECF::UINT);
	state->getRegistry()->registerEntry("genotype", (voidP) new string("bitstring"), ECF::STRING, "genotype variant (bitstring, tree)");
}


bool EvalOp::initialize(StateP state)
{
	state_ = state;
	showTruth = false;
	inputNames.clear();
	inputMap.clear();

	bestScore = 1;

	voidP sptr = state->getRegistry()->getEntry("variables");
	nVariables = *((uint*)sptr.get());

	sptr = state->getRegistry()->getEntry("variant");
	variant = *((int*)sptr.get());

	sptr = state->getRegistry()->getEntry("omega");
	omega = *((uint*)sptr.get());

	std::string genString = *((std::string*) state->getRegistry()->getEntry("genotype").get());
	genotype = bitstring;
	if (genString == "bitstring")
		genotype = bitstring;
	else if (genString == "tree")
		genotype = tree;

	if (genotype == tree) {
		// build variable names 
		for (uint var = 0; var < nVariables; var++)
			inputNames.push_back("v" + uint2str(var));

		// build a complete table of combinations of input variables 
		inputMap.resize(nVariables);
		bool val; //, t = true, f = false;
		uint ttSize = (uint)pow(2., (int)nVariables);

		results.resize(ttSize);
		tt.resize(ttSize);

		for (uint count = 0; count < ttSize; count++) {
			uint mask = count;
			for (uint var = 0; var < nVariables; var++) {
				if (mask % 2)
					val = true;
				else
					val = false;
				mask /= 2;
				inputMap[var].push_back(val);
			}
		}

		// set values to all variables (do not change during evolution!) 
		Tree::Tree* tree = (Tree::Tree*) state->getGenotypes().at(0).get();
		for (uint var = 0; var < nVariables; var++) {
			tree->setTerminalValue(inputNames[var], &inputMap[var]);
		}
	}

	return true;
}


#ifndef MOEA

FitnessP EvalOp::evaluate(IndividualP individual)
{
	FitnessP fitness (new FitnessMin);

	// populate bitVector depending on genotype
	if (genotype == bitstring) {
		BitString::BitString* bitstr = (BitString::BitString*) individual->getGenotype().get();
		results = bitstr->bits;
	}
	else if (genotype == tree) {
		// get tree from the individual
		Tree::Tree* tree = (Tree::Tree*) individual->getGenotype().get();
		// build TT
		tree->execute(&results);
	}

	// Hamming weight factor
	int hwf = results.size();
	for(uint i = 0; i < results.size(); i++)
		if(results[i] == true)
			hwf--;

	// REMARK: hwf actually contains the number of 0s in results, not of 1s.
	// Hence it is the complement of the Hamming weight
	if(hwf == results.size()) {
		fitness->setValue(1. * hwf * hwf * hwf);
		return fitness;
	}

	static std::vector< std::vector<int> > setland;
	convertGenFuncToLandscapes(results, nVariables + 1, omega, setland);

	int score = compFitnessLandscapes(setland, omega);

	switch (variant) {
	// SOGP, just minimize score
	case 0: {
		break; }

	// lexi GP: reach min score 0, then keep 0 and optimize hwf
	case 1: {
		// try to diversify solutions with min fitness
		if (score == 0)
			score -= (results.size() - hwf);

		if (score < bestScore) {
			bestScore = score;
			ofstream output("solutions.txt", ios_base::app);
			for (uint i = 0; i < results.size(); i++)
				output << results[i];
			output << endl;
			output.close();
			output.open("hw_values.txt", ios_base::app);
			output << results.size() - hwf << endl;
			output.close();
		}
		break; }
	}

	fitness->setValue(score);

	if(evaluateVerbose) {
		ECF_LOG(state_, 1, "Truth table:");
		stringstream ss;
		for(uint i = 0; i < results.size(); i++)
			ss << results[i];// << " ";
		ss << endl;
		ss << "score: " << score << endl;
		ss << "hwf: " << hwf << endl;
		ECF_LOG(state_, 1, ss.str());

		if(showTruth) {
			ofstream out("tt.txt", ios_base::app);
			for(uint i = 0; i < results.size(); i++)
				out << results[i];// << " ";
			out << "\n";
			out.close();
		}
	}

	return fitness;
}

#endif



#ifdef MOEA

// MOEA verzija
FitnessP EvalOp::evaluate(IndividualP individual)
{
	MOFitnessP fitness = static_cast<MOFitnessP> (new MOFitness);
	fitness->push_back((FitnessP) new FitnessMin);
	fitness->push_back((FitnessP) new FitnessMax);

	// get tree from the individual
	Tree::Tree* tree = (Tree::Tree*) individual->getGenotype().get();

	// build TT
	tree->execute(&results);

	// Hamming weight factor
	int hwf = results.size();
	for(uint i = 0; i < results.size(); i++)
		if(results[i] == true)
			hwf--;

	// REMARK: hwf actually contains the number of 0s in results, not of 1s.
	// Hence it is the complement of the Hamming weight
	if(hwf == results.size()) {
		fitness->at(0)->setValue(1. * hwf * hwf * hwf);
		fitness->at(1)->setValue(0);
		return fitness;
	}

	static std::vector< std::vector<int> > setland;
	convertGenFuncToLandscapes(results, nVariables + 1, omega, setland);

	int score = compFitnessLandscapes(setland, omega);

	fitness->at(0)->setValue(score);
	fitness->at(1)->setValue(results.size() - hwf);

	if(evaluateVerbose) {
		ECF_LOG(state_, 1, "Truth table:");
		stringstream ss;
		for(uint i = 0; i < results.size(); i++)
			ss << results[i];// << " ";
		ss << endl;
		ss << "score: " << score << endl;
		ss << "hwf: " << hwf << endl;
		ECF_LOG(state_, 1, ss.str());
	}

	return fitness;
}

#endif


std::vector<int> EvalOp::getTT(IndividualP individual, int& vars, uint genId)
{
	// get tree from the individual (genotype with genId)
	Tree::Tree* tree = (Tree::Tree*) individual->getGenotype(genId).get();

	// set values to all variables (do not change during evolution!) 
	static bool zadano = false;
	if(!zadano) {
		zadano = true;
		for(uint var = 0; var < nVariables; var++) {
			tree->setTerminalValue(inputNames[var], &inputMap[var]);
		}
	}


	//int score = 0;

	tree->execute(&results);

	uint tries = (uint) tt.size();
	for(uint i = 0; i < tries; i++)
		tt[i] = results[i];

	vars = nVariables;
	return tt;
}
