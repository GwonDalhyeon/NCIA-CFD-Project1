#include "CSR.h"
#include "NSEquation.h"

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
	//cout << testGrid1.x[0] << "\n";
	//cout << testGrid1.x[10] << "\n";
	GridInfo testGrid2 = GridInfo(testGrid1);
	//testGrid2 = testGrid1;

	testGrid2.X0 = 2;
	cout << testGrid1.X0 << "\n";
	cout << testGrid2.X0 << "\n";
	// change.
	
}