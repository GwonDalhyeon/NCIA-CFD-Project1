#include "CSR.h"
#include "NSEquation.h"
#include "PoissonEquation.h"
#include "LevelSet.h"


void main()
{
	GridInfo testGrid1 = GridInfo(0.0, 1.0, 101);

	PoissonEquationSolver testPoisson = PoissonEquationSolver(testGrid1);
	testPoisson.solvePoissonEquationJumpCondi(1);
}