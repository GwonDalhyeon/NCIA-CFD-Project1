#include "Field2D.h"
#include "VectorND.h"
#include "CSR.h"

#include "PoissonSolver.h"

#include "AdvectionMethod2D.h"
#include "SurfaceReconstruction.h"

double SurfaceReconst<double>::alpha;
void main()
{
	Grid2D testGrid2d(0.0, 1.0, 101, 0.0, 1.0, 101);

	SurfaceReconst<double>::alpha = min(testGrid2d.dx, testGrid2d.dy);

	double x;

	x = SurfaceReconst<double>::heaviside(0.1);
	cout << x << endl;
	//Grid2D grid(0,1,10, 0,1,20);

	
	//PoissonSolver testPoisson2d;
	//testPoisson2d.solvePoissonJumpCondi(2, testGrid2d);
	

	//CSR<double> testCSR();
	//cout << testCSR << endl;
	//cout <<1 << endl;
	
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