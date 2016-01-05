#include "CSR.h"
#include "NSEquation.h"
//#include "PoissonEquation.h"
#include "LevelSet.h"


void main()
{
	double* a = new double[10];
	for (size_t i = 0; i < 10; i++)
	{
		a[i] = 1;
	}
	double* b = NULL;

	delete a;
	delete b;

	GridInfo testGrid1 = GridInfo(0.0, 1.0, 101);
	
	LevelSet testLevelSet = LevelSet(testGrid1);
	cout << testGrid1.deltaX << "\n";
	cout << testLevelSet.test << "\n";
	if (testLevelSet.grid == nullptr)
	{
		cout << "null" << "\n";
		cout << testLevelSet.grid << "\n";
	}
	else
	{
		cout << testLevelSet.grid->deltaX << "\n";
	}
	
}