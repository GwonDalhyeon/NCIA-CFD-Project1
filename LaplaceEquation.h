#pragma once
#include "CSR.h"
#include "GridInfomation.h"
#include "LevelSet.h"

#include "LinearEquationSolver.h"

class LaplaceEquationSolver
{
public:
	double* laplaceMatrix;
	double* laplaceVector;
	double* f;
	double* jCondition1;
	double* jCondition2;
	double* solution;
	double* tempSol;
	double leftBdry, rightBdry;

	GridInfo grid;
	LevelSet levelSet;

	LaplaceEquationSolver();
	~LaplaceEquationSolver();

	int index(int i, int j);
	int index(int i, int j, int k);
	int indexVec(int i, int j);
	int indexVec(int i, int j, int k);
	//int indexMat(int i, int j);
	//int indexMat(int i, int j, int k);

	// A Boundary Condition Capturing Method
	void generateJumpCondi(int example, GridInfo inputGrid);
	void generateLaplaceMatrixJumpCondi();
	void generateLaplaceVectorJumpCondi();
	void solveLaplaceEquationJumpCondi(int example, GridInfo inputGrid);
	void outputResult();
private:

};

LaplaceEquationSolver::LaplaceEquationSolver()
{
	if (grid.dimension == 1)
	{
		laplaceMatrix = new double[grid.numMatX*grid.numMatX];
		laplaceVector = new double[grid.numMatX];
		f = new double[grid.numMatX];
		jCondition1 = new double[grid.numX];
		jCondition2 = new double[grid.numX];
		solution = new double[grid.numX];
		for (int i = 0; i < grid.numMatX*grid.numMatX; i++)
		{
			laplaceMatrix[i] = 0;
		}
		for (int i = 0; i < grid.numMatX; i++)
		{
			laplaceVector[i] = 0;
			f[i] = 0;
		}
		for (int i = 0; i < grid.numX; i++)
		{
			jCondition1[i] = 0;
			jCondition2[i] = 0;
			solution[i] = 0;
		}
	}
	if (grid.dimension == 2)
	{
		laplaceMatrix = new double[grid.numMatX*grid.numMatX*grid.numMatY*grid.numMatY];
		laplaceVector = new double[grid.numMatX*grid.numMatY];
		f = new double[grid.numMatX*grid.numMatY];
		jCondition1 = new double[grid.numX*grid.numY];
		jCondition2 = new double[grid.numX*grid.numY];
		solution = new double[grid.numX*grid.numY];
		for (int i = 0; i < grid.numMatX*grid.numMatX*grid.numMatY*grid.numMatY; i++)
		{
			laplaceMatrix[i] = 0;
		}
		for (int i = 0; i < grid.numMatX*grid.numMatY; i++)
		{
			laplaceVector[i] = 0;
			f[i] = 0;
		}
		for (int i = 0; i < grid.numX*grid.numY; i++)
		{
			jCondition1[i] = 0;
			jCondition2[i] = 0;
			solution[i] = 0;
		}
	}
}

LaplaceEquationSolver::~LaplaceEquationSolver()
{
	delete laplaceMatrix, laplaceVector, f, jCondition1, jCondition2;
}

inline int LaplaceEquationSolver::index(int i, int j)
{
	return i + j*grid.numX;
}

inline int LaplaceEquationSolver::index(int i, int j, int k)
{
	return i + j*grid.numX + k*grid.numX*grid.numY;
}

inline int LaplaceEquationSolver::indexVec(int i, int j)
{
	return i + j*grid.numMatX;
}

inline int LaplaceEquationSolver::indexVec(int i, int j, int k)
{
	return i + j*grid.numMatX + k*grid.numMatX*grid.numMatY;
}

inline void LaplaceEquationSolver::generateJumpCondi(int example, GridInfo inputGrid)
{
	/////////////////////////////////////////////////////////////////////////////
	//
	//     Laplace equation. Example 1
	//
	/////////////////////////////////////////////////////////////////////////////


	//double X0, X1, Y0, Y1, Z0, Z1;
	//double deltaX, deltaY, deltaZ;
	//int numX, numY, numZ;
	//X0 = 0; X1 = 1; numX = 101;
	//Y0 = 0; Y1 = 1; numY = 101;
	//Z0 = 0; Z1 = 1; numZ = 101;

	grid = inputGrid;
	levelSet = LevelSet(grid);
	
	for (int i = 0; i < grid.numX; i++)
	{
		levelSet.phi[i] = grid.x[i] - 0.5 - grid.deltaX / 2;
	}

	leftBdry = 0, rightBdry = 1;
	f[0] = leftBdry / (grid.deltaX*grid.deltaX);
	f[grid.numMatX - 1] = rightBdry / (grid.deltaX*grid.deltaX);
	for (int i = 1; i < grid.numMatX - 1; i++)
	{
		f[i] = 0;
	}

	int centIndex, rightIndex;
	for (int i = 0; i < grid.numX - 1; i++)
	{
		centIndex = i;
		rightIndex = i + 1;
		if ((levelSet.phi[centIndex] <= 0 && levelSet.phi[rightIndex]>0) || (levelSet.phi[centIndex]>0 && levelSet.phi[rightIndex] <= 0))
		{
			jCondition1[centIndex] = 1;
			jCondition1[rightIndex] = 1;
			jCondition2[centIndex] = 0;
			jCondition2[rightIndex] = 0;
		}
	}

}

inline void LaplaceEquationSolver::generateLaplaceMatrixJumpCondi()
{
	if (grid.dimension == 1)
	{
		for (int i = 0; i < grid.numMatX; i++)
		{
			if (i>0 && i<grid.numMatX - 1)
			{
				laplaceMatrix[i*grid.numMatX + i - 1]	= -1 / (grid.deltaX*grid.deltaX);
				laplaceMatrix[i*grid.numMatX + i]		= 2 / (grid.deltaX*grid.deltaX);
				laplaceMatrix[i*grid.numMatX + i + 1]	= -1 / (grid.deltaX*grid.deltaX);
			}
			else if (i == 0)
			{
				laplaceMatrix[i*grid.numMatX + i]		= 2 / (grid.deltaX*grid.deltaX);
				laplaceMatrix[i*grid.numMatX + i + 1]	= -1 / (grid.deltaX*grid.deltaX);
			}
			else
			{
				laplaceMatrix[i*grid.numMatX + i - 1]	= -1 / (grid.deltaX*grid.deltaX);
				laplaceMatrix[i*grid.numMatX + i]		= 2 / (grid.deltaX*grid.deltaX);
			}
		}
	}
	else if (grid.dimension == 2)
	{

	}
}

inline void LaplaceEquationSolver::generateLaplaceVectorJumpCondi()
{
	double aGamma = 0;
	double bGamma = 0;
	double theta = 0;

	if (grid.dimension == 1)
	{
		for (int i = 1; i < grid.numMatX - 1; i++)
		{
			if ((levelSet.phi[i + 1] <= 0 && levelSet.phi[i + 1 + 1]>0) || (levelSet.phi[i + 1]>0 && levelSet.phi[i + 1 + 1] <= 0))
			{
				theta				= abs(levelSet.phi[i + 1]) / (abs(levelSet.phi[i + 1]) + abs(levelSet.phi[i + 1 + 1]));
				aGamma				= (jCondition1[i + 1] * abs(levelSet.phi[i + 1 + 1]) + jCondition1[i + 1 + 1] * abs(levelSet.phi[i + 1])) / (abs(levelSet.phi[i]) + abs(levelSet.phi[i + 1 + 1]));
				bGamma				= (jCondition2[i + 1] * abs(levelSet.phi[i + 1 + 1]) + jCondition2[i + 1 + 1] * abs(levelSet.phi[i + 1])) / (abs(levelSet.phi[i]) + abs(levelSet.phi[i + 1 + 1]));
				laplaceVector[i]	= laplaceVector[i] - (aGamma / (grid.deltaX*grid.deltaX) + bGamma*(1 - theta) / grid.deltaX);
				
				i++;
				laplaceVector[i]	= laplaceVector[i] - (-aGamma / (grid.deltaX*grid.deltaX) + bGamma*theta / grid.deltaX);
			}
		}
	}
	else if (grid.dimension == 2)
	{

	}
}

inline void LaplaceEquationSolver::solveLaplaceEquationJumpCondi(int example, GridInfo inputGrid)
{

}

inline void LaplaceEquationSolver::outputResult()
{
	ofstream solutionFile;
	solutionFile.open("E:\Data/laplace.txt");

	for (int i = 0; i < grid.numX; i++)
	{
		solutionFile << i << " " << grid.x[i] << " " << solution[i] << endl;
		//solutionFile<<grid.x[i]<< endl;
	}
	solutionFile.close();
}

