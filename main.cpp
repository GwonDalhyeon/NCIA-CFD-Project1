#include "LaplaceEquation.h"
#include "PoissonEquation.h"
#include "Vector2D.h"
#include "Array2D.h"



void main()
{
	Vector2D<double> testVector1(-3, 4);
	//cout << testVector1 << "\n";
	//cout << testVector1.magnitude() << "\n";
	Vector2D<double> testVector2 = normalize(testVector1);
	Vector2D<double> testVector3 = testVector1;
	testVector3.normalize();
	cout << testVector1 << "\n";
	cout << testVector2 << "\n";
	cout << testVector3 << "\n";

	//testVector2 =100.0 + testVector1;
	//Vector2D<double> testVector3 = 1.0 /testVector1;
	//testVector1 += 2;
	//cout << testVector1.x << " " << testVector1.y << "\n";
	//cout << testVector2.x << " " << testVector2.y << "\n";
	//cout << testVector3.x << " " << testVector3.y << "\n";


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