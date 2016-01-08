#include "LaplaceEquation.h"
#include "PoissonEquation.h"

void main()
{

	Array2D<int> testArray(2, 4);
	testArray = 2;
	cout << testArray(0,0) << "\n";
	Array2D<int> testArray2 = testArray;
	testArray2 = testArray;
	
	Vector2D<double> testVector(1,3);
	Vector2D<double> testVector1 = testVector;


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