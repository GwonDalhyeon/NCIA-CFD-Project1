#include "LaplaceEquation.h"
#include "PoissonEquation.h"
#include "Field2D.h"
#include "VectorND.h"
#include "CSR.h"

void main()
{
	Grid2D grid(1, 2, 100, 3, 5, 400);
	//Field2D<double> field(grid);
	//cout << field(0, 0) << endl;

	Array2D<double> testArray(2, 4);
	testArray(0, 0) = 1;
	testArray(0, 1) = 4;
	testArray(1, 3) = 7;
	CSR<double> testCSR(testArray);
	cout << testCSR << endl;
	cout <<1 << endl;

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