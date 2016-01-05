#pragma once

#include "CSR.h"
#include "GridInfomation.h"

#include "LinearEquationSolver.h"

class NSEquationSolver
{
public:
	GridInfo* GridNode;
	GridInfo* GridCell;

	NSEquationSolver();
	~NSEquationSolver();

private:

};

NSEquationSolver::NSEquationSolver()
{
}

NSEquationSolver::~NSEquationSolver()
{
}

