#include "ecf/ECF.h"



// custom terminal class - boolean vector
class BoolV : public Tree::Primitives::Primitive
{
public:
	std::vector<bool> value_;

	BoolV()
	{
		nArguments_ = 0;
	}
	void execute(void* result, Tree::Tree& tree)
	{
		std::vector<bool>& res = *(std::vector<bool>*)result;
		res = value_;
	}
	void setValue(void* value)
	{
		value_ = *(std::vector<bool>*)value;
	}
	~BoolV()
	{	}
};



//
// Boolean function primitives
//
class Or: public Tree::Primitives::Primitive
{
public:
	Or()
	{
		nArguments_ = 2;
		name_ = "OR";
	}

	void execute(void* result, Tree::Tree& tree)
	{
		std::vector<bool>& vOr = *(std::vector<bool>*)result;
		uint size = (uint) vOr.size();

		std::vector<bool> arg1(size), arg2(size);

		getNextArgument(&arg1, tree);
		getNextArgument(&arg2, tree);

		for(uint i = 0; i < size; i++)
			vOr[i] = arg1[i] || arg2[i];
	}

	~Or()
	{	}
};

class Nor : public Tree::Primitives::Primitive
{
public:
	Nor()
	{
		nArguments_ = 2;
		name_ = "NOR";
	}

	void execute(void* result, Tree::Tree& tree)
	{
		std::vector<bool>& vOr = *(std::vector<bool>*)result;
		uint size = (uint)vOr.size();

		std::vector<bool> arg1(size), arg2(size);

		getNextArgument(&arg1, tree);
		getNextArgument(&arg2, tree);

		for (uint i = 0; i < size; i++)
			//vOr[i] = (!arg1[i]) || (!arg2[i]);
			vOr[i] = (!arg1[i]) && (!arg2[i]);
	}

	~Nor()
	{	}
};


class Xor: public Tree::Primitives::Primitive
{
public:
	Xor()
	{
		nArguments_ = 2;
		name_ = "XOR";
	}

	void execute(void* result, Tree::Tree& tree)
	{
		std::vector<bool>& vXor = *(std::vector<bool>*)result;
		uint size = (uint) vXor.size();

		std::vector<bool> arg1(size), arg2(size);

		getNextArgument(&arg1, tree);
		getNextArgument(&arg2, tree);

		for(uint i = 0; i < size; i++)
			vXor[i] = (arg1[i] && !arg2[i]) || (!arg1[i] && arg2[i]);
	}

	~Xor()
	{	}
};


class XNor: public Tree::Primitives::Primitive
{
public:
	XNor()
	{
		nArguments_ = 2;
		name_ = "XNOR";
	}

	void execute(void* result, Tree::Tree& tree)
	{
		std::vector<bool>& vXNor = *(std::vector<bool>*)result;
		uint size = (uint) vXNor.size();

		std::vector<bool> arg1(size), arg2(size);

		getNextArgument(&arg1, tree);
		getNextArgument(&arg2, tree);

		for(uint i = 0; i < size; i++)
			vXNor[i] = (!(arg1[i] && !arg2[i]) || (!arg1[i] && arg2[i]));
	}

	~XNor()
	{	}
};


class And: public Tree::Primitives::Primitive
{
public:
	And()
	{
		nArguments_ = 2;
		name_ = "AND";
	}

	void execute(void* result, Tree::Tree& tree)
	{
		std::vector<bool>& vAnd = *(std::vector<bool>*)result;
		uint size = (uint) vAnd.size();

		std::vector<bool> arg1(size), arg2(size);

		getNextArgument(&arg1, tree);
		getNextArgument(&arg2, tree);

		for(uint i = 0; i < size; i++)
			vAnd[i] = arg1[i] && arg2[i];
	}

	~And()
	{	}
};


class Nand : public Tree::Primitives::Primitive
{
public:
	Nand()
	{
		nArguments_ = 2;
		name_ = "NAND";
	}

	void execute(void* result, Tree::Tree& tree)
	{
		std::vector<bool>& vAnd = *(std::vector<bool>*)result;
		uint size = (uint)vAnd.size();

		std::vector<bool> arg1(size), arg2(size);

		getNextArgument(&arg1, tree);
		getNextArgument(&arg2, tree);

		for (uint i = 0; i < size; i++)
			vAnd[i] = (!arg1[i]) || (!arg2[i]);
	}

	~Nand()
	{	}
};

class And2: public Tree::Primitives::Primitive
{
public:
	And2()
	{
		nArguments_ = 2;
		name_ = "AND2";
	}

	void execute(void* result, Tree::Tree& tree)
	{
		std::vector<bool>& vAnd2 = *(std::vector<bool>*)result;
		uint size = (uint) vAnd2.size();

		std::vector<bool> arg1(size), arg2(size);

		getNextArgument(&arg1, tree);
		getNextArgument(&arg2, tree);

		for(uint i = 0; i < size; i++)
			vAnd2[i] = arg1[i] && (!arg2[i]);
	}

	~And2()
	{	}
};


class Not: public Tree::Primitives::Primitive
{
public:
	Not()
	{
		nArguments_ = 1;
		name_ = "NOT";
	}

	void execute(void* result, Tree::Tree& tree)
	{
		std::vector<bool>& vNot = *(std::vector<bool>*)result;
		uint size = (uint) vNot.size();

		std::vector<bool> arg1(size);
		getNextArgument(&arg1, tree);

		for(uint i = 0; i < size; i++)
			vNot[i] = !arg1[i];
	}

	~Not()
	{	}
};


class If: public Tree::Primitives::Primitive
{
public:
	If()
	{
		nArguments_ = 3;
		name_ = "IF";
	}

	void execute(void* result, Tree::Tree& tree)
	{
		std::vector<bool>& res = *(std::vector<bool>*)result;
		uint size = (uint) res.size();

		std::vector<bool> arg(size), res1(size), res2(size);
		getNextArgument(&arg, tree);
		getNextArgument(&res1, tree);
		getNextArgument(&res2, tree);

		for(uint i = 0; i < size; i++)
			if(arg[i]) {
				res[i] = res1[i];
			} else {
				res[i] = res2[i];
			}
	}

	~If()
	{	}
};
