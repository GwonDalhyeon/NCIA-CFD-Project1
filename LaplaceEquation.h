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
	double* beta;
	double* f;
	double* jCondition1;
	double* jCondition2;

	GridInfo grid;
	LevelSet levelSet;

	LaplaceEquationSolver();
	~LaplaceEquationSolver();
	void generateLaplaceMatrix();
	void generateLaplaceVector();

private:

};

LaplaceEquationSolver::LaplaceEquationSolver()
{
	if (grid.dimension == 1)
	{
		poissonMatrix = new double[grid.numMatX*grid.numMatX];
		for (int i = 0; i < grid.numMatX*grid.numMatX; i++)
		{
			poissonMatrix[i] = 0;
		}
	}
	if (grid.dimension == 2)
	{
		poissonMatrix = new double[grid.numMatX*grid.numMatX*grid.numMatY*grid.numMatY];

		for (int i = 0; i < grid.numMatX*grid.numMatX*grid.numMatY*grid.numMatY; i++)
		{
			poissonMatrix[i] = 0;
		}
	}
}

LaplaceEquationSolver::~LaplaceEquationSolver()
{
	delete laplaceMatrix, laplaceVector, beta, f, jCondition1, jCondition2;
}

inline void LaplaceEquationSolver::generateLaplaceMatrix()
{
}

inline void LaplaceEquationSolver::generateLaplaceVector()
{
}

double* laplaceMatrix(GridInfo& domainInfo, int dimension)
{
	if (dimension == 1)
	{
		double* A = new double[domainInfo.numMatX*domainInfo.numMatX];
		for (int i = 0; i < domainInfo.numMatX*domainInfo.numMatX; i++)
		{
			A[i] = 0;
		}
		for (int i = 0; i < domainInfo.numMatX; i++)
		{
			if (i>0 && i<domainInfo.numMatX - 1)
			{
				A[i*domainInfo.numMatX + i - 1] = -1 / (domainInfo.deltaX*domainInfo.deltaX);
				A[i*domainInfo.numMatX + i] = 2 / (domainInfo.deltaX*domainInfo.deltaX);
				A[i*domainInfo.numMatX + i + 1] = -1 / (domainInfo.deltaX*domainInfo.deltaX);
			}
			else if (i == 0)
			{
				A[i*domainInfo.numMatX + i] = 2 / (domainInfo.deltaX*domainInfo.deltaX);
				A[i*domainInfo.numMatX + i + 1] = -1 / (domainInfo.deltaX*domainInfo.deltaX);
			}
			else
			{
				A[i*domainInfo.numMatX + i - 1] = -1 / (domainInfo.deltaX*domainInfo.deltaX);
				A[i*domainInfo.numMatX + i] = 2 / (domainInfo.deltaX*domainInfo.deltaX);
			}
		}
		return A;
	}
	else if (dimension == 2)
	{

	}
}

double* laplaceVector(double* levelSet, double* f, double* jCondition1, double* jCondition2, GridInfo& domainInfo, int dimension)
{
	double aGamma = 0;
	double bGamma = 0;
	double theta = 0;

	if (dimension == 1)
	{
		double* b = new double[domainInfo.numMatX];
		for (int i = 0; i < domainInfo.numMatX; i++)
		{
			b[i] = f[i];
		}

		for (int i = 1; i < domainInfo.numMatX - 1; i++)
		{
			if ((levelSet[i + 1] <= 0 && levelSet[i + 1 + 1]>0) || (levelSet[i + 1]>0 && levelSet[i + 1 + 1] <= 0))
			{
				theta = abs(levelSet[i + 1]) / (abs(levelSet[i + 1]) + abs(levelSet[i + 1 + 1]));
				aGamma = (jCondition1[i + 1] * abs(levelSet[i + 1 + 1]) + jCondition1[i + 1 + 1] * abs(levelSet[i + 1])) / (abs(levelSet[i]) + abs(levelSet[i + 1 + 1]));
				bGamma = (jCondition2[i + 1] * abs(levelSet[i + 1 + 1]) + jCondition2[i + 1 + 1] * abs(levelSet[i + 1])) / (abs(levelSet[i]) + abs(levelSet[i + 1 + 1]));
				b[i] = b[i] - (aGamma / (domainInfo.deltaX*domainInfo.deltaX) + bGamma*(1 - theta) / domainInfo.deltaX);
				//cout<<b[i]<<endl;
				//cout<<i+1<< endl;
				i++;
				b[i] = b[i] - (-aGamma / (domainInfo.deltaX*domainInfo.deltaX) + bGamma*theta / domainInfo.deltaX);
				//cout<<b[i]<<endl;
			}
		}
		//cout<<endl;
		return b;
	}
	else if (dimension == 2)
	{
		double* b = new double[domainInfo.numX*domainInfo.numY];


		return b;
	}
}

