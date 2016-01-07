#include "CSR.h"
#include "NSEquation.h"
#include "PoissonEquation.h"
#include "LevelSet.h"


void main()
{
	//GridInfo testGrid1d = GridInfo(0.0, 1.0, 101);
	//PoissonEquationSolver testPoisson1d = PoissonEquationSolver(testGrid1d);
	//testPoisson1d.solvePoissonEquationJumpCondi(1);

	GridInfo testGrid2d = GridInfo(0.0, 1.0, 101, 0.0, 1.0, 101);
	PoissonEquationSolver testPoisson2d = PoissonEquationSolver(testGrid2d);
	testPoisson2d.solvePoissonEquationJumpCondi(2);
}