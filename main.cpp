#include "PoissonSolver.h"

#include "AdvectionMethod2D.h"
#include "SurfaceReconstruction.h"
#include "LevelSetReinitializationProblem.h"

void main()
{

	//SurfaceReconst<double>::alpha = min(testGrid2d.dx, testGrid2d.dy);

	//SurfaceReconst<double> surface;
	//surface.surfaceReconstructionSolver(2);
	
	Reinitialzation reinitial;
	reinitial. reinitializationSolver(4);
	
	//PoissonSolver testPoisson2d;
	//testPoisson2d.solvePoissonJumpCondi(2, testGrid2d);
	
	//GridInfo testGrid1d(0.0, 1.0, 101);
	//LaplaceEquationSolver testLaplace = LaplaceEquationSolver(testGrid1d);
	//testLaplace.solveLaplaceEquationJumpCondi(1);

	//GridInfo testGrid1d(0.0, 1.0, 101);
	//PoissonEquationSolver testPoisson1d = PoissonEquationSolver(testGrid1d);
	//testPoisson1d.solvePoissonEquationJumpCondi(1);

	//GridInfo testGrid2d(0.0, 1.0, 101, 0.0, 1.0, 101);
	//PoissonEquationSolver testPoisson2d = PoissonEquationSolver(testGrid2d);
	//testPoisson2d.solvePoissonEquationJumpCondi(2);
}