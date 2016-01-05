#pragma once
#include "CSR.h"
#include "GridInfomation.h"
#include "LevelSet.h"

#include "LinearEquationSolver.h"

class PoissonEquationSolver
{
public:
	double* poissonMatrix;
	double* poissonVector;

	GridInfo* grid;

	PoissonEquationSolver();
	~PoissonEquationSolver();

private:

};

PoissonEquationSolver::PoissonEquationSolver()
{
}

PoissonEquationSolver::~PoissonEquationSolver()
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

double* poissonMatrix(double* levelSet, double* beta, GridInfo& domainInfo, int dimension)
{
	double tempBeta = 0;
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
				if ((levelSet[i + 1]>0 && levelSet[i + 1 + 1] <= 0) || (levelSet[i + 1] <= 0 && levelSet[i + 1 + 1]>0))
				{
					tempBeta = beta[i + 1] * beta[i + 1 + 1] * (abs(levelSet[i + 1]) + abs(levelSet[i + 1 + 1])) / (beta[i + 1 + 1] * abs(levelSet[i + 1]) + beta[i + 1] * abs(levelSet[i + 1 + 1]));
					A[i*domainInfo.numMatX + i - 1] = -1 / (domainInfo.deltaX*domainInfo.deltaX)*(beta[i + 1] + beta[i]) / 2;
					A[i*domainInfo.numMatX + i] = 1 / (domainInfo.deltaX*domainInfo.deltaX)*(beta[i + 1] + beta[i]) / 2 + 1 / (domainInfo.deltaX*domainInfo.deltaX)*(tempBeta);
					A[i*domainInfo.numMatX + i + 1] = -1 / (domainInfo.deltaX*domainInfo.deltaX)*tempBeta;
					//cout<<i<<" "<< tempBeta<<endl;
				}
				else if ((levelSet[i]>0 && levelSet[i + 1] <= 0) || (levelSet[i] <= 0 && levelSet[i + 1]>0))
				{
					tempBeta = beta[i] * beta[i + 1] * (abs(levelSet[i]) + abs(levelSet[i + 1])) / (beta[i + 1] * abs(levelSet[i]) + beta[i] * abs(levelSet[i + 1]));
					A[i*domainInfo.numMatX + i - 1] = -1 / (domainInfo.deltaX*domainInfo.deltaX)*tempBeta;
					A[i*domainInfo.numMatX + i] = 1 / (domainInfo.deltaX*domainInfo.deltaX)*(beta[i + 1] + beta[i + 1 + 1]) / 2 + 1 / (domainInfo.deltaX*domainInfo.deltaX)*(tempBeta);
					A[i*domainInfo.numMatX + i + 1] = -1 / (domainInfo.deltaX*domainInfo.deltaX)*(beta[i + 1] + beta[i + 1 + 1]) / 2;
					//cout<<i<<" "<< tempBeta<<endl;
				}
				else
				{
					A[i*domainInfo.numMatX + i - 1] = -1 / (domainInfo.deltaX*domainInfo.deltaX)*(beta[i + 1] + beta[i]) / 2;
					A[i*domainInfo.numMatX + i] = 1 / (domainInfo.deltaX*domainInfo.deltaX)*(beta[i + 1] + beta[i]) / 2 + 1 / (domainInfo.deltaX*domainInfo.deltaX)*(beta[i + 1] + beta[i + 1 + 1]) / 2;
					A[i*domainInfo.numMatX + i + 1] = -1 / (domainInfo.deltaX*domainInfo.deltaX)*(beta[i + 1] + beta[i + 1 + 1]) / 2;
				}

			}
			else if (i == 0)
			{
				A[i*domainInfo.numMatX + i] = 1 / (domainInfo.deltaX*domainInfo.deltaX)*beta[i] + 1 / (domainInfo.deltaX*domainInfo.deltaX)*(beta[i] + beta[i + 1]) / 2;
				A[i*domainInfo.numMatX + i + 1] = -1 / (domainInfo.deltaX*domainInfo.deltaX)*(beta[i + 1] + beta[i + 1 + 1]) / 2;
			}
			else
			{
				A[i*domainInfo.numMatX + i - 1] = -1 / (domainInfo.deltaX*domainInfo.deltaX)*(beta[i + 1] + beta[i]) / 2;
				A[i*domainInfo.numMatX + i] = 1 / (domainInfo.deltaX*domainInfo.deltaX)*(beta[i + 1] + beta[i]) / 2 + 1 / (domainInfo.deltaX*domainInfo.deltaX)*beta[i + 1];
			}
		}
		return A;
	}
	else if (dimension == 2)
	{
		double* A = new double[domainInfo.numMatX*domainInfo.numMatX*domainInfo.numMatY*domainInfo.numMatY];
		int matIndex, matLeftIndex, matRightIndex, matTopIndex, matBottomIndex;
		int centIndex, leftIndex, rightIndex, topIndex, bottomIndex;
		for (int i = 0; i < domainInfo.numMatX*domainInfo.numMatX*domainInfo.numMatY*domainInfo.numMatY; i++)
		{
			A[i] = 0;
		}

		for (int j = 0; j < domainInfo.numMatY; j++)
		{
			//A[i*domainInfo.numMatX*domainInfo.numMatY + i + j*domainInfo.numMatX*domainInfo.numMatX*domainInfo.numMatY + j*domainInfo.numMatX ] = 0;
			//A[i*domainInfo.numMatX*domainInfo.numMatY + i + j*domainInfo.numMatX*(domainInfo.numMatX*domainInfo.numMatY + 1) ] = 0;
			for (int i = 0; i < domainInfo.numMatX; i++)
			{
				matIndex = i*domainInfo.numMatX*domainInfo.numMatY + i + j*domainInfo.numMatX*(domainInfo.numMatX*domainInfo.numMatY + 1);
				matLeftIndex = i*domainInfo.numMatX*domainInfo.numMatY + i - 1 + j*domainInfo.numMatX*(domainInfo.numMatX*domainInfo.numMatY + 1);
				matRightIndex = i*domainInfo.numMatX*domainInfo.numMatY + i + 1 + j*domainInfo.numMatX*(domainInfo.numMatX*domainInfo.numMatY + 1);
				matBottomIndex = i*domainInfo.numMatX*domainInfo.numMatY + i + j*domainInfo.numMatX*domainInfo.numMatX*domainInfo.numMatY + (j - 1)*domainInfo.numMatX;
				matTopIndex = i*domainInfo.numMatX*domainInfo.numMatY + i + j*domainInfo.numMatX*domainInfo.numMatX*domainInfo.numMatY + (j + 1)*domainInfo.numMatX;

				centIndex = i + 1 + (j + 1)*domainInfo.numX;
				leftIndex = i + (j + 1)*domainInfo.numX;
				rightIndex = i + 1 + 1 + (j + 1)*domainInfo.numX;
				bottomIndex = i + 1 + j*domainInfo.numX;
				topIndex = i + 1 + (j + 1 + 1)*domainInfo.numX;

				if (j>0 && j<domainInfo.numMatY - 1)
				{
					if ((levelSet[centIndex]>0 && levelSet[topIndex] <= 0) || (levelSet[centIndex] <= 0 && levelSet[topIndex]>0))
					{
						tempBeta = beta[centIndex] * beta[topIndex] * (abs(levelSet[centIndex]) + abs(levelSet[topIndex])) / (beta[topIndex] * abs(levelSet[centIndex]) + beta[centIndex] * abs(levelSet[topIndex]));
						A[matBottomIndex] = -1 / (domainInfo.deltaY*domainInfo.deltaY)*(beta[centIndex] + beta[bottomIndex]) / 2;
						A[matIndex] = A[matIndex] + 1 / (domainInfo.deltaY*domainInfo.deltaY)*(beta[centIndex] + beta[bottomIndex]) / 2 + 1 / (domainInfo.deltaY*domainInfo.deltaY)*(tempBeta);
						A[matTopIndex] = -1 / (domainInfo.deltaY*domainInfo.deltaY)*tempBeta;
						//cout<<i<<" "<< tempBeta<<endl;
						//cout<<i<<" " << j<<" "<<A[matBottomIndex]<<" "<< A[matIndex]<<" " <<A[matTopIndex] <<endl;
						//cout<<"";

					}
					else if ((levelSet[bottomIndex]>0 && levelSet[centIndex] <= 0) || (levelSet[bottomIndex] <= 0 && levelSet[centIndex]>0))
					{
						tempBeta = beta[bottomIndex] * beta[centIndex] * (abs(levelSet[bottomIndex]) + abs(levelSet[centIndex])) / (beta[centIndex] * abs(levelSet[bottomIndex]) + beta[bottomIndex] * abs(levelSet[centIndex]));
						A[matBottomIndex] = -1 / (domainInfo.deltaY*domainInfo.deltaY)*tempBeta;
						A[matIndex] = A[matIndex] + 1 / (domainInfo.deltaY*domainInfo.deltaY)*(beta[centIndex] + beta[topIndex]) / 2 + 1 / (domainInfo.deltaY*domainInfo.deltaY)*(tempBeta);
						A[matTopIndex] = -1 / (domainInfo.deltaY*domainInfo.deltaY)*(beta[centIndex] + beta[topIndex]) / 2;
						//cout<<i<<" "<< tempBeta<<endl;
						//cout<<i<<" " << j<<" "<<A[matBottomIndex]<<" "<< A[matIndex]<<" " <<A[matTopIndex] <<endl;
						//cout<<"";
					}
					else
					{
						A[matBottomIndex] = -1 / (domainInfo.deltaY*domainInfo.deltaY)*(beta[centIndex] + beta[bottomIndex]) / 2;
						A[matIndex] = A[matIndex] + 1 / (domainInfo.deltaY*domainInfo.deltaY)*(beta[centIndex] + beta[topIndex]) / 2 + 1 / (domainInfo.deltaY*domainInfo.deltaY)*(beta[centIndex] + beta[bottomIndex]) / 2;
						A[matTopIndex] = -1 / (domainInfo.deltaY*domainInfo.deltaY)*(beta[centIndex] + beta[topIndex]) / 2;

						//cout<<i<<" " << j<<" "<<A[matBottomIndex]<<" "<< A[matIndex]<<" " <<A[matTopIndex] <<endl;
						//cout<<"";
					}

					if (i>0 && i<domainInfo.numMatX - 1)
					{
						if ((levelSet[centIndex]>0 && levelSet[rightIndex] <= 0) || (levelSet[centIndex] <= 0 && levelSet[rightIndex]>0))
						{
							tempBeta = beta[centIndex] * beta[rightIndex] * (abs(levelSet[centIndex]) + abs(levelSet[rightIndex])) / (beta[rightIndex] * abs(levelSet[centIndex]) + beta[centIndex] * abs(levelSet[rightIndex]));
							A[matLeftIndex] = -1 / (domainInfo.deltaX*domainInfo.deltaX)*(beta[centIndex] + beta[leftIndex]) / 2;
							A[matIndex] = A[matIndex] + 1 / (domainInfo.deltaX*domainInfo.deltaX)*(beta[centIndex] + beta[leftIndex]) / 2 + 1 / (domainInfo.deltaX*domainInfo.deltaX)*(tempBeta);
							A[matRightIndex] = -1 / (domainInfo.deltaX*domainInfo.deltaX)*tempBeta;
							//cout<<i<<" " << j<<" "<<A[matLeftIndex]<<" "<< A[matIndex]<<" " <<A[matRightIndex] <<endl;
							//cout<<"";
						}
						else if ((levelSet[leftIndex]>0 && levelSet[centIndex] <= 0) || (levelSet[leftIndex] <= 0 && levelSet[centIndex]>0))
						{
							tempBeta = beta[leftIndex] * beta[centIndex] * (abs(levelSet[leftIndex]) + abs(levelSet[centIndex])) / (beta[centIndex] * abs(levelSet[leftIndex]) + beta[leftIndex] * abs(levelSet[centIndex]));
							A[matLeftIndex] = -1 / (domainInfo.deltaX*domainInfo.deltaX)*tempBeta;
							A[matIndex] = A[matIndex] + 1 / (domainInfo.deltaX*domainInfo.deltaX)*(beta[centIndex] + beta[rightIndex]) / 2 + 1 / (domainInfo.deltaX*domainInfo.deltaX)*(tempBeta);
							A[matRightIndex] = -1 / (domainInfo.deltaX*domainInfo.deltaX)*(beta[centIndex] + beta[rightIndex]) / 2;
							//cout<<i<<" "<< tempBeta<<endl;
							//cout<<i<<" " << j<<" "<<A[matLeftIndex]<<" "<< A[matIndex]<<" " <<A[matRightIndex] <<endl;
							//cout<<"";
						}
						else
						{
							A[matLeftIndex] = -1 / (domainInfo.deltaX*domainInfo.deltaX)*(beta[centIndex] + beta[leftIndex]) / 2;
							A[matIndex] = A[matIndex] + 1 / (domainInfo.deltaX*domainInfo.deltaX)*(beta[centIndex] + beta[rightIndex]) / 2 + 1 / (domainInfo.deltaX*domainInfo.deltaX)*(beta[centIndex] + beta[leftIndex]) / 2;
							A[matRightIndex] = -1 / (domainInfo.deltaX*domainInfo.deltaX)*(beta[centIndex] + beta[rightIndex]) / 2;
							//cout<<i<<" " << j<<" "<<A[matLeftIndex]<<" "<< A[matIndex]<<" " <<A[matRightIndex] <<endl;
							//cout<<"";
						}
					}
					else if (i == 0)
					{
						A[matIndex] = A[matIndex] + 1 / (domainInfo.deltaX*domainInfo.deltaX)*(beta[centIndex] + beta[rightIndex]) / 2 + 1 / (domainInfo.deltaX*domainInfo.deltaX)*(beta[centIndex] + beta[leftIndex]) / 2;
						A[matRightIndex] = -1 / (domainInfo.deltaX*domainInfo.deltaX)*(beta[centIndex] + beta[rightIndex]) / 2;
						//cout<<i<<" " << j<<" "<<0<<" "<< A[matIndex]<<" " <<A[matRightIndex] <<endl;
						//cout<<"";
					}
					else
					{
						A[matLeftIndex] = -1 / (domainInfo.deltaX*domainInfo.deltaX)*(beta[centIndex] + beta[leftIndex]) / 2;
						A[matIndex] = A[matIndex] + 1 / (domainInfo.deltaX*domainInfo.deltaX)*(beta[centIndex] + beta[rightIndex]) / 2 + 1 / (domainInfo.deltaX*domainInfo.deltaX)*(beta[centIndex] + beta[leftIndex]) / 2;
						//cout<<i<<" " << j<<" "<<A[matLeftIndex]<<" "<< A[matIndex]<<" " <<0<<endl;
						//cout<<"";
					}
				}
				else if (j == 0)
				{
					A[matIndex] = A[matIndex] + 1 / (domainInfo.deltaY*domainInfo.deltaY)*(beta[centIndex] + beta[topIndex]) / 2 + 1 / (domainInfo.deltaY*domainInfo.deltaY)*(beta[centIndex] + beta[bottomIndex]) / 2;
					A[matTopIndex] = -1 / (domainInfo.deltaY*domainInfo.deltaY)*(beta[centIndex] + beta[topIndex]) / 2;
					//cout<<i<<" " << j<<" "<<0<<" "<< A[matIndex]<<" " <<A[matTopIndex] <<endl;
					//cout<<"";
					if (i>0 && i<domainInfo.numMatX - 1)
					{
						if ((levelSet[centIndex]>0 && levelSet[rightIndex] <= 0) || (levelSet[centIndex] <= 0 && levelSet[rightIndex]>0))
						{
							tempBeta = beta[centIndex] * beta[rightIndex] * (abs(levelSet[centIndex]) + abs(levelSet[rightIndex])) / (beta[rightIndex] * abs(levelSet[centIndex]) + beta[centIndex] * abs(levelSet[rightIndex]));
							A[matLeftIndex] = -1 / (domainInfo.deltaX*domainInfo.deltaX)*(beta[centIndex] + beta[leftIndex]) / 2;
							A[matIndex] = A[matIndex] + 1 / (domainInfo.deltaX*domainInfo.deltaX)*(beta[centIndex] + beta[leftIndex]) / 2 + 1 / (domainInfo.deltaX*domainInfo.deltaX)*(tempBeta);
							A[matRightIndex] = -1 / (domainInfo.deltaX*domainInfo.deltaX)*tempBeta;
							//cout<<i<<" "<< tempBeta<<endl;
							//cout<<i<<" " << j<<" "<<A[matLeftIndex]<<" "<< A[matIndex]<<" " <<A[matRightIndex] <<endl;
							//cout<<"";
						}
						else if ((levelSet[leftIndex]>0 && levelSet[centIndex] <= 0) || (levelSet[leftIndex] <= 0 && levelSet[centIndex]>0))
						{
							tempBeta = beta[leftIndex] * beta[centIndex] * (abs(levelSet[leftIndex]) + abs(levelSet[centIndex])) / (beta[centIndex] * abs(levelSet[leftIndex]) + beta[leftIndex] * abs(levelSet[centIndex]));
							A[matLeftIndex] = -1 / (domainInfo.deltaX*domainInfo.deltaX)*tempBeta;
							A[matIndex] = A[matIndex] + 1 / (domainInfo.deltaX*domainInfo.deltaX)*(beta[centIndex] + beta[rightIndex]) / 2 + 1 / (domainInfo.deltaX*domainInfo.deltaX)*(tempBeta);
							A[matRightIndex] = -1 / (domainInfo.deltaX*domainInfo.deltaX)*(beta[centIndex] + beta[rightIndex]) / 2;
							//cout<<i<<" "<< tempBeta<<endl;
							//cout<<i<<" " << j<<" "<<A[matLeftIndex]<<" "<< A[matIndex]<<" " <<A[matRightIndex] <<endl;
							//cout<<"";
						}
						else
						{
							A[matLeftIndex] = -1 / (domainInfo.deltaX*domainInfo.deltaX)*(beta[centIndex] + beta[leftIndex]) / 2;
							A[matIndex] = A[matIndex] + 1 / (domainInfo.deltaX*domainInfo.deltaX)*(beta[centIndex] + beta[rightIndex]) / 2 + 1 / (domainInfo.deltaX*domainInfo.deltaX)*(beta[centIndex] + beta[leftIndex]) / 2;
							A[matRightIndex] = -1 / (domainInfo.deltaX*domainInfo.deltaX)*(beta[centIndex] + beta[rightIndex]) / 2;
							//cout<<i<<" " << j<<" "<<A[matLeftIndex]<<" "<< A[matIndex]<<" " <<A[matRightIndex] <<endl;
							//cout<<"";
						}
					}
					else if (i == 0)
					{
						A[matIndex] = A[matIndex] + 1 / (domainInfo.deltaX*domainInfo.deltaX)*(beta[centIndex] + beta[rightIndex]) / 2 + 1 / (domainInfo.deltaX*domainInfo.deltaX)*(beta[centIndex] + beta[leftIndex]) / 2;
						A[matRightIndex] = -1 / (domainInfo.deltaX*domainInfo.deltaX)*(beta[centIndex] + beta[rightIndex]) / 2;
						//cout<<i<<" " << j<<" "<<0<<" "<< A[matIndex]<<" " <<A[matRightIndex] <<endl;
						//cout<<"";
					}
					else
					{
						A[matLeftIndex] = -1 / (domainInfo.deltaX*domainInfo.deltaX)*(beta[centIndex] + beta[leftIndex]) / 2;
						A[matIndex] = A[matIndex] + 1 / (domainInfo.deltaX*domainInfo.deltaX)*(beta[centIndex] + beta[rightIndex]) / 2 + 1 / (domainInfo.deltaX*domainInfo.deltaX)*(beta[centIndex] + beta[leftIndex]) / 2;
						//cout<<i<<" " << j<<" "<<A[matLeftIndex]<<" "<< A[matIndex]<<" " <<0<<endl;
						//cout<<"";
					}
				}
				else
				{
					A[matBottomIndex] = -1 / (domainInfo.deltaY*domainInfo.deltaY)*(beta[centIndex] + beta[bottomIndex]) / 2;
					A[matIndex] = A[matIndex] + 1 / (domainInfo.deltaY*domainInfo.deltaY)*(beta[centIndex] + beta[topIndex]) / 2 + 1 / (domainInfo.deltaY*domainInfo.deltaY)*(beta[centIndex] + beta[bottomIndex]) / 2;
					//cout<<i<<" " << j<<" "<<A[matBottomIndex]<<" "<< A[matIndex]<<" " <<0 <<endl;
					//cout<<"";
					if (i>0 && i<domainInfo.numMatX - 1)
					{
						if ((levelSet[centIndex]>0 && levelSet[rightIndex] <= 0) || (levelSet[centIndex] <= 0 && levelSet[rightIndex]>0))
						{
							tempBeta = beta[centIndex] * beta[rightIndex] * (abs(levelSet[centIndex]) + abs(levelSet[rightIndex])) / (beta[rightIndex] * abs(levelSet[centIndex]) + beta[centIndex] * abs(levelSet[rightIndex]));
							A[matLeftIndex] = -1 / (domainInfo.deltaX*domainInfo.deltaX)*(beta[centIndex] + beta[leftIndex]) / 2;
							A[matIndex] = A[matIndex] + 1 / (domainInfo.deltaX*domainInfo.deltaX)*(beta[centIndex] + beta[leftIndex]) / 2 + 1 / (domainInfo.deltaX*domainInfo.deltaX)*(tempBeta);
							A[matRightIndex] = -1 / (domainInfo.deltaX*domainInfo.deltaX)*tempBeta;
							//cout<<i<<" "<< tempBeta<<endl;
							//cout<<i<<" " << j<<" "<<A[matLeftIndex]<<" "<< A[matIndex]<<" " <<A[matRightIndex] <<endl;
							//cout<<"";
						}
						else if ((levelSet[leftIndex]>0 && levelSet[centIndex] <= 0) || (levelSet[leftIndex] <= 0 && levelSet[centIndex]>0))
						{
							tempBeta = beta[leftIndex] * beta[centIndex] * (abs(levelSet[leftIndex]) + abs(levelSet[centIndex])) / (beta[centIndex] * abs(levelSet[leftIndex]) + beta[leftIndex] * abs(levelSet[centIndex]));
							A[matLeftIndex] = -1 / (domainInfo.deltaX*domainInfo.deltaX)*tempBeta;
							A[matIndex] = A[matIndex] + 1 / (domainInfo.deltaX*domainInfo.deltaX)*(beta[centIndex] + beta[rightIndex]) / 2 + 1 / (domainInfo.deltaX*domainInfo.deltaX)*(tempBeta);
							A[matRightIndex] = -1 / (domainInfo.deltaX*domainInfo.deltaX)*(beta[centIndex] + beta[rightIndex]) / 2;
							//cout<<i<<" "<< tempBeta<<endl;
							//cout<<i<<" " << j<<" "<<A[matLeftIndex]<<" "<< A[matIndex]<<" " <<A[matRightIndex] <<endl;
							//cout<<"";
						}
						else
						{
							A[matLeftIndex] = -1 / (domainInfo.deltaX*domainInfo.deltaX)*(beta[centIndex] + beta[leftIndex]) / 2;
							A[matIndex] = A[matIndex] + 1 / (domainInfo.deltaX*domainInfo.deltaX)*(beta[centIndex] + beta[rightIndex]) / 2 + 1 / (domainInfo.deltaX*domainInfo.deltaX)*(beta[centIndex] + beta[leftIndex]) / 2;
							A[matRightIndex] = -1 / (domainInfo.deltaX*domainInfo.deltaX)*(beta[centIndex] + beta[rightIndex]) / 2;
							//cout<<i<<" " << j<<" "<<A[matLeftIndex]<<" "<< A[matIndex]<<" " <<A[matRightIndex] <<endl;
							//cout<<"";
						}
					}
					else if (i == 0)
					{
						A[matIndex] = A[matIndex] + 1 / (domainInfo.deltaX*domainInfo.deltaX)*(beta[centIndex] + beta[rightIndex]) / 2 + 1 / (domainInfo.deltaX*domainInfo.deltaX)*(beta[centIndex] + beta[leftIndex]) / 2;
						A[matRightIndex] = -1 / (domainInfo.deltaX*domainInfo.deltaX)*(beta[centIndex] + beta[rightIndex]) / 2;
						//cout<<i<<" " << j<<" "<<0<<" "<< A[matIndex]<<" " <<A[matRightIndex] <<endl;
						//cout<<"";
					}
					else
					{
						A[matLeftIndex] = -1 / (domainInfo.deltaX*domainInfo.deltaX)*(beta[centIndex] + beta[leftIndex]) / 2;
						A[matIndex] = A[matIndex] + 1 / (domainInfo.deltaX*domainInfo.deltaX)*(beta[centIndex] + beta[rightIndex]) / 2 + 1 / (domainInfo.deltaX*domainInfo.deltaX)*(beta[centIndex] + beta[leftIndex]) / 2;
						//cout<<i<<" " << j<<" "<<A[matLeftIndex]<<" "<< A[matIndex]<<" " <<0 <<endl;
						//cout<<"";
					}
				}
			}
		}

		return A;
	}
}

double* poissonVector(double* levelSet, double* f, double* beta, double* jCondition1, double* jCondition2, GridInfo& domainInfo, int dimension)
{
	double aGamma = 0;
	double bGamma = 0;
	double theta = 0;
	double tempBeta = 0;

	if (dimension == 1)
	{
		double* fL = new double[domainInfo.numMatX];
		double* fR = new double[domainInfo.numMatX];
		double* b = new double[domainInfo.numMatX];

		double normalLeft, normalCenter, normalRight;

		for (int i = 0; i < domainInfo.numMatX; i++)
		{
			normalLeft = unitNormal1D(levelSet, i, domainInfo);
			normalCenter = unitNormal1D(levelSet, i + 1, domainInfo);
			normalRight = unitNormal1D(levelSet, i + 1 + 1, domainInfo);

			if (levelSet[i]>0 && levelSet[i + 1] <= 0)
			{
				theta = abs(levelSet[i]) / (abs(levelSet[i]) + abs(levelSet[i + 1]));
				aGamma = (jCondition1[i] * abs(levelSet[i + 1]) + jCondition1[i + 1] * abs(levelSet[i])) / (abs(levelSet[i]) + abs(levelSet[i + 1]));
				bGamma = (jCondition2[i] * normalLeft*abs(levelSet[i + 1]) + jCondition2[i + 1] * normalCenter*abs(levelSet[i])) / (abs(levelSet[i]) + abs(levelSet[i + 1]));
				tempBeta = beta[i] * beta[i + 1] * (abs(levelSet[i]) + abs(levelSet[i + 1])) / (beta[i + 1] * abs(levelSet[i]) + beta[i] * abs(levelSet[i + 1]));
				fL[i] = tempBeta*aGamma / (domainInfo.deltaX*domainInfo.deltaX) - tempBeta*bGamma*theta / (beta[i] * domainInfo.deltaX);
				//cout<<levelSet[i]<<" "<<levelSet[i+1]<<endl;
				//cout<<i<<" " <<i+1 <<endl;
				//cout << i<< " "<< fL[i]<<endl;
				//cout<<"theta"<<theta<<endl;
				//cout<<"aGamma : "<<aGamma<<endl;
				//cout<<"bGamma : "<<bGamma<<endl;
				//cout<<"tempBeta : "<<tempBeta<<endl;
				//cout<<endl;
			}
			else if (levelSet[i] <= 0 && levelSet[i + 1]>0)
			{
				theta = abs(levelSet[i]) / (abs(levelSet[i]) + abs(levelSet[i + 1]));
				aGamma = (jCondition1[i] * abs(levelSet[i + 1]) + jCondition1[i + 1] * abs(levelSet[i])) / (abs(levelSet[i]) + abs(levelSet[i + 1]));
				bGamma = (jCondition2[i] * normalLeft*abs(levelSet[i + 1]) + jCondition2[i + 1] * normalCenter*abs(levelSet[i])) / (abs(levelSet[i]) + abs(levelSet[i + 1]));
				tempBeta = beta[i] * beta[i + 1] * (abs(levelSet[i]) + abs(levelSet[i + 1])) / (beta[i + 1] * abs(levelSet[i]) + beta[i] * abs(levelSet[i + 1]));
				fL[i] = -tempBeta*aGamma / (domainInfo.deltaX*domainInfo.deltaX) + tempBeta*bGamma*theta / (beta[i] * domainInfo.deltaX);
				//cout<<levelSet[i]<<" "<<levelSet[i+1]<<endl;
				//cout<<i<<" " <<i+1 <<endl;
				//cout << i<< " "<< fL[i]<<endl;
				//cout<<"theta"<<theta<<endl;
				//cout<<"aGamma : "<<aGamma<<endl;
				//cout<<"bGamma : "<<bGamma<<endl;
				//cout<<"tempBeta : "<<tempBeta<<endl;
				//cout<<endl;
			}
			else
			{
				fL[i] = 0;
			}

			if (levelSet[i + 1] <= 0 && levelSet[i + 1 + 1]>0)
			{
				theta = abs(levelSet[i + 1 + 1]) / (abs(levelSet[i + 1]) + abs(levelSet[i + 1 + 1]));
				aGamma = (jCondition1[i + 1] * abs(levelSet[i + 1 + 1]) + jCondition1[i + 1 + 1] * abs(levelSet[i + 1])) / (abs(levelSet[i + 1]) + abs(levelSet[i + 1 + 1]));
				bGamma = (jCondition2[i + 1] * normalCenter*abs(levelSet[i + 1 + 1]) + jCondition2[i + 1 + 1] * normalRight*abs(levelSet[i + 1])) / (abs(levelSet[i + 1]) + abs(levelSet[i + 1 + 1]));
				tempBeta = beta[i + 1] * beta[i + 1 + 1] * (abs(levelSet[i + 1]) + abs(levelSet[i + 1 + 1])) / (beta[i + 1 + 1] * abs(levelSet[i + 1]) + beta[i + 1] * abs(levelSet[i + 1 + 1]));
				fR[i] = tempBeta*aGamma / (domainInfo.deltaX*domainInfo.deltaX) + tempBeta*bGamma*theta / (beta[i + 1 + 1] * domainInfo.deltaX);
				//cout<<levelSet[i+1]<<" "<<levelSet[i+1+1]<<endl;
				//cout<<i+1<<" " <<i+1+1 <<endl;
				//cout << i<< " "<< fR[i]<<endl;
				//cout<<"theta"<<theta<<endl;
				//cout<<"aGamma : "<<aGamma<<endl;
				//cout<<"bGamma : "<<bGamma<<endl;
				//cout<<"tempBeta : "<<tempBeta<<endl;
				//cout<<endl;

			}
			else if (levelSet[i + 1]>0 && levelSet[i + 1 + 1] <= 0)
			{
				theta = abs(levelSet[i + 1 + 1]) / (abs(levelSet[i + 1]) + abs(levelSet[i + 1 + 1]));
				aGamma = (jCondition1[i + 1] * abs(levelSet[i + 1 + 1]) + jCondition1[i + 1 + 1] * abs(levelSet[i + 1])) / (abs(levelSet[i + 1]) + abs(levelSet[i + 1 + 1]));
				bGamma = (jCondition2[i + 1] * normalCenter*abs(levelSet[i + 1 + 1]) + jCondition2[i + 1 + 1] * normalRight*abs(levelSet[i + 1])) / (abs(levelSet[i + 1]) + abs(levelSet[i + 1 + 1]));
				tempBeta = beta[i + 1] * beta[i + 1 + 1] * (abs(levelSet[i + 1]) + abs(levelSet[i + 1 + 1])) / (beta[i + 1 + 1] * abs(levelSet[i + 1]) + beta[i + 1] * abs(levelSet[i + 1 + 1]));
				fR[i] = -tempBeta*aGamma / (domainInfo.deltaX*domainInfo.deltaX) - tempBeta*bGamma*theta / (beta[i + 1 + 1] * domainInfo.deltaX);
				//cout<<levelSet[i+1]<<" "<<levelSet[i+1+1]<<endl;
				//cout<<i+1<<" " <<i+1+1 <<endl;
				//cout << i<< " "<< fR[i]<<endl;
				//cout<<"theta"<<theta<<endl;
				//cout<<"aGamma : "<<aGamma<<endl;
				//cout<<"bGamma : "<<bGamma<<endl;
				//cout<<"tempBeta : "<<tempBeta<<endl;
				//cout<<endl;
			}
			else
			{
				fR[i] = 0;
			}

			b[i] = -f[i] - fR[i] - fL[i];

		}

		delete[] fR, fL;
		return b;
	}
	else if (dimension == 2)
	{
		double* fL = new double[domainInfo.numMatX*domainInfo.numMatY];
		double* fR = new double[domainInfo.numMatX*domainInfo.numMatY];
		double* fB = new double[domainInfo.numMatX*domainInfo.numMatY];
		double* fT = new double[domainInfo.numMatX*domainInfo.numMatY];
		double* b = new double[domainInfo.numMatX*domainInfo.numMatY];

		double* normalLeft = new double[2];
		double* normalRight = new double[2];
		double* normalCenter = new double[2];
		double* normalBottom = new double[2];
		double* normalTop = new double[2];

		int matIndex;
		int centIndex, leftIndex, rightIndex, topIndex, bottomIndex;

		for (int j = 0; j < domainInfo.numMatY; j++)
		{
			for (int i = 0; i < domainInfo.numMatX; i++)
			{
				matIndex = i + j*domainInfo.numMatX;

				centIndex = i + 1 + (j + 1)*domainInfo.numX;
				leftIndex = i + (j + 1)*domainInfo.numX;
				rightIndex = i + 1 + 1 + (j + 1)*domainInfo.numX;
				bottomIndex = i + 1 + j*domainInfo.numX;
				topIndex = i + 1 + (j + 1 + 1)*domainInfo.numX;

				unitNormal2D(levelSet, i - 1, j, domainInfo, normalLeft);
				unitNormal2D(levelSet, i + 1, j, domainInfo, normalRight);
				unitNormal2D(levelSet, i, j, domainInfo, normalCenter);
				unitNormal2D(levelSet, i, j - 1, domainInfo, normalBottom);
				unitNormal2D(levelSet, i, j + 1, domainInfo, normalTop);

				if (levelSet[leftIndex]>0 && levelSet[centIndex] <= 0)
				{
					theta = abs(levelSet[leftIndex]) / (abs(levelSet[leftIndex]) + abs(levelSet[centIndex]));
					aGamma = (jCondition1[leftIndex] * abs(levelSet[centIndex]) + jCondition1[centIndex] * abs(levelSet[leftIndex])) / (abs(levelSet[leftIndex]) + abs(levelSet[centIndex]));
					bGamma = (jCondition2[leftIndex] * normalLeft[0] * abs(levelSet[centIndex]) + jCondition2[centIndex] * normalCenter[0] * abs(levelSet[leftIndex])) / (abs(levelSet[leftIndex]) + abs(levelSet[centIndex]));
					tempBeta = beta[leftIndex] * beta[centIndex] * (abs(levelSet[leftIndex]) + abs(levelSet[centIndex])) / (beta[centIndex] * abs(levelSet[leftIndex]) + beta[leftIndex] * abs(levelSet[centIndex]));
					fL[matIndex] = tempBeta*aGamma / (domainInfo.deltaX*domainInfo.deltaX) - tempBeta*bGamma*theta / (beta[leftIndex] * domainInfo.deltaX);
					//cout<<levelSet[i]<<" "<<levelSet[i+1]<<endl;
					//cout<<i<<" " <<i+1 <<endl;
					//cout << i<< " "<< fL[i]<<endl;
					//cout<<"theta"<<theta<<endl;
					//cout<<"aGamma : "<<aGamma<<endl;
					//cout<<"bGamma : "<<bGamma<<endl;
					//cout<<"tempBeta : "<<tempBeta<<endl;
					//cout<<endl;
				}
				else if (levelSet[leftIndex] <= 0 && levelSet[centIndex]>0)
				{
					theta = abs(levelSet[leftIndex]) / (abs(levelSet[leftIndex]) + abs(levelSet[centIndex]));
					aGamma = (jCondition1[leftIndex] * abs(levelSet[centIndex]) + jCondition1[centIndex] * abs(levelSet[leftIndex])) / (abs(levelSet[leftIndex]) + abs(levelSet[centIndex]));
					bGamma = (jCondition2[leftIndex] * normalLeft[0] * abs(levelSet[centIndex]) + jCondition2[centIndex] * normalCenter[0] * abs(levelSet[leftIndex])) / (abs(levelSet[leftIndex]) + abs(levelSet[centIndex]));
					tempBeta = beta[leftIndex] * beta[centIndex] * (abs(levelSet[leftIndex]) + abs(levelSet[centIndex])) / (beta[centIndex] * abs(levelSet[leftIndex]) + beta[leftIndex] * abs(levelSet[centIndex]));
					fL[matIndex] = -tempBeta*aGamma / (domainInfo.deltaX*domainInfo.deltaX) + tempBeta*bGamma*theta / (beta[leftIndex] * domainInfo.deltaX);
					//cout<<levelSet[i]<<" "<<levelSet[i+1]<<endl;
					//cout<<i<<" " <<i+1 <<endl;
					//cout << i<< " "<< fL[i]<<endl;
					//cout<<"theta"<<theta<<endl;
					//cout<<"aGamma : "<<aGamma<<endl;
					//cout<<"bGamma : "<<bGamma<<endl;
					//cout<<"tempBeta : "<<tempBeta<<endl;
					//cout<<endl;
				}
				else
				{
					fL[matIndex] = 0;
				}

				if (levelSet[centIndex] <= 0 && levelSet[rightIndex]>0)
				{
					theta = abs(levelSet[rightIndex]) / (abs(levelSet[centIndex]) + abs(levelSet[rightIndex]));
					aGamma = (jCondition1[centIndex] * abs(levelSet[rightIndex]) + jCondition1[rightIndex] * abs(levelSet[centIndex])) / (abs(levelSet[centIndex]) + abs(levelSet[rightIndex]));
					bGamma = (jCondition2[centIndex] * normalCenter[0] * abs(levelSet[rightIndex]) + jCondition2[rightIndex] * normalRight[0] * abs(levelSet[centIndex])) / (abs(levelSet[centIndex]) + abs(levelSet[rightIndex]));
					tempBeta = beta[centIndex] * beta[rightIndex] * (abs(levelSet[centIndex]) + abs(levelSet[rightIndex])) / (beta[rightIndex] * abs(levelSet[centIndex]) + beta[centIndex] * abs(levelSet[rightIndex]));
					fR[matIndex] = tempBeta*aGamma / (domainInfo.deltaX*domainInfo.deltaX) + tempBeta*bGamma*theta / (beta[rightIndex] * domainInfo.deltaX);
					//cout<<levelSet[i+1]<<" "<<levelSet[i+1+1]<<endl;
					//cout<<i+1<<" " <<i+1+1 <<endl;
					//cout << i<< " "<< fR[i]<<endl;
					//cout<<"theta"<<theta<<endl;
					//cout<<"aGamma : "<<aGamma<<endl;
					//cout<<"bGamma : "<<bGamma<<endl;
					//cout<<"tempBeta : "<<tempBeta<<endl;
					//cout<<endl;

				}
				else if (levelSet[centIndex]>0 && levelSet[rightIndex] <= 0)
				{
					theta = abs(levelSet[rightIndex]) / (abs(levelSet[centIndex]) + abs(levelSet[rightIndex]));
					aGamma = (jCondition1[centIndex] * abs(levelSet[rightIndex]) + jCondition1[rightIndex] * abs(levelSet[centIndex])) / (abs(levelSet[centIndex]) + abs(levelSet[rightIndex]));
					bGamma = (jCondition2[centIndex] * normalCenter[0] * abs(levelSet[rightIndex]) + jCondition2[rightIndex] * normalRight[0] * abs(levelSet[centIndex])) / (abs(levelSet[centIndex]) + abs(levelSet[rightIndex]));
					tempBeta = beta[centIndex] * beta[rightIndex] * (abs(levelSet[centIndex]) + abs(levelSet[rightIndex])) / (beta[rightIndex] * abs(levelSet[centIndex]) + beta[centIndex] * abs(levelSet[rightIndex]));
					fR[matIndex] = -tempBeta*aGamma / (domainInfo.deltaX*domainInfo.deltaX) - tempBeta*bGamma*theta / (beta[rightIndex] * domainInfo.deltaX);
					//cout<<levelSet[i+1]<<" "<<levelSet[i+1+1]<<endl;
					//cout<<i+1<<" " <<i+1+1 <<endl;
					//cout << i<< " "<< fR[i]<<endl;
					//cout<<"theta"<<theta<<endl;
					//cout<<"aGamma : "<<aGamma<<endl;
					//cout<<"bGamma : "<<bGamma<<endl;
					//cout<<"tempBeta : "<<tempBeta<<endl;
					//cout<<endl;
				}
				else
				{
					fR[matIndex] = 0;
				}

				if (levelSet[bottomIndex]>0 && levelSet[centIndex] <= 0)
				{
					theta = abs(levelSet[bottomIndex]) / (abs(levelSet[bottomIndex]) + abs(levelSet[centIndex]));
					aGamma = (jCondition1[bottomIndex] * abs(levelSet[centIndex]) + jCondition1[centIndex] * abs(levelSet[bottomIndex])) / (abs(levelSet[bottomIndex]) + abs(levelSet[centIndex]));
					bGamma = (jCondition2[bottomIndex] * normalBottom[1] * abs(levelSet[centIndex]) + jCondition2[centIndex] * normalCenter[1] * abs(levelSet[bottomIndex])) / (abs(levelSet[bottomIndex]) + abs(levelSet[centIndex]));
					tempBeta = beta[bottomIndex] * beta[centIndex] * (abs(levelSet[bottomIndex]) + abs(levelSet[centIndex])) / (beta[centIndex] * abs(levelSet[bottomIndex]) + beta[bottomIndex] * abs(levelSet[centIndex]));
					fB[matIndex] = tempBeta*aGamma / (domainInfo.deltaX*domainInfo.deltaX) - tempBeta*bGamma*theta / (beta[bottomIndex] * domainInfo.deltaX);
					//cout<<levelSet[i]<<" "<<levelSet[i+1]<<endl;
					//cout<<i<<" " <<i+1 <<endl;
					//cout << i<< " "<< fL[i]<<endl;
					//cout<<"theta"<<theta<<endl;
					//cout<<"aGamma : "<<aGamma<<endl;
					//cout<<"bGamma : "<<bGamma<<endl;
					//cout<<"tempBeta : "<<tempBeta<<endl;
					//cout<<endl;
				}
				else if (levelSet[bottomIndex] <= 0 && levelSet[centIndex]>0)
				{
					theta = abs(levelSet[bottomIndex]) / (abs(levelSet[bottomIndex]) + abs(levelSet[centIndex]));
					aGamma = (jCondition1[bottomIndex] * abs(levelSet[centIndex]) + jCondition1[centIndex] * abs(levelSet[bottomIndex])) / (abs(levelSet[bottomIndex]) + abs(levelSet[centIndex]));
					bGamma = (jCondition2[bottomIndex] * normalBottom[1] * abs(levelSet[centIndex]) + jCondition2[centIndex] * normalCenter[1] * abs(levelSet[bottomIndex])) / (abs(levelSet[bottomIndex]) + abs(levelSet[centIndex]));
					tempBeta = beta[bottomIndex] * beta[centIndex] * (abs(levelSet[bottomIndex]) + abs(levelSet[centIndex])) / (beta[centIndex] * abs(levelSet[bottomIndex]) + beta[bottomIndex] * abs(levelSet[centIndex]));
					fB[matIndex] = -tempBeta*aGamma / (domainInfo.deltaX*domainInfo.deltaX) + tempBeta*bGamma*theta / (beta[bottomIndex] * domainInfo.deltaX);
					//cout<<levelSet[i]<<" "<<levelSet[i+1]<<endl;
					//cout<<i+1<<" " <<j+1 <<endl;
					//cout <<fB[matIndex]<<endl;
					//cout<<"theta"<<theta<<endl;
					//cout<<"aGamma : "<<aGamma<<endl;
					//cout<<"bGamma : "<<bGamma<<endl;
					//cout<<"tempBeta : "<<tempBeta<<endl;
					//cout<<endl;
				}
				else
				{
					fB[matIndex] = 0;
				}

				if (levelSet[centIndex] <= 0 && levelSet[topIndex]>0)
				{
					theta = abs(levelSet[topIndex]) / (abs(levelSet[centIndex]) + abs(levelSet[topIndex]));
					aGamma = (jCondition1[centIndex] * abs(levelSet[topIndex]) + jCondition1[topIndex] * abs(levelSet[centIndex])) / (abs(levelSet[centIndex]) + abs(levelSet[topIndex]));
					bGamma = (jCondition2[centIndex] * normalCenter[1] * abs(levelSet[topIndex]) + jCondition2[topIndex] * normalTop[1] * abs(levelSet[centIndex])) / (abs(levelSet[centIndex]) + abs(levelSet[topIndex]));
					tempBeta = beta[centIndex] * beta[topIndex] * (abs(levelSet[centIndex]) + abs(levelSet[topIndex])) / (beta[topIndex] * abs(levelSet[centIndex]) + beta[centIndex] * abs(levelSet[topIndex]));
					fT[matIndex] = tempBeta*aGamma / (domainInfo.deltaX*domainInfo.deltaX) + tempBeta*bGamma*theta / (beta[topIndex] * domainInfo.deltaX);
					//cout<<levelSet[i+1]<<" "<<levelSet[i+1+1]<<endl;
					//cout<<i+1<<" " <<j+1 <<endl;
					//cout << fT[matIndex]<<endl;
					//cout<<"theta"<<theta<<endl;
					//cout<<"aGamma : "<<aGamma<<endl;
					//cout<<"bGamma : "<<bGamma<<endl;
					//cout<<"tempBeta : "<<tempBeta<<endl;
					//cout<<endl;

				}
				else if (levelSet[centIndex]>0 && levelSet[topIndex] <= 0)
				{
					theta = abs(levelSet[topIndex]) / (abs(levelSet[centIndex]) + abs(levelSet[topIndex]));
					aGamma = (jCondition1[centIndex] * abs(levelSet[topIndex]) + jCondition1[topIndex] * abs(levelSet[centIndex])) / (abs(levelSet[centIndex]) + abs(levelSet[topIndex]));
					bGamma = (jCondition2[centIndex] * normalCenter[1] * abs(levelSet[topIndex]) + jCondition2[topIndex] * normalTop[1] * abs(levelSet[centIndex])) / (abs(levelSet[centIndex]) + abs(levelSet[topIndex]));
					tempBeta = beta[centIndex] * beta[topIndex] * (abs(levelSet[centIndex]) + abs(levelSet[topIndex])) / (beta[topIndex] * abs(levelSet[centIndex]) + beta[centIndex] * abs(levelSet[topIndex]));
					fT[matIndex] = -tempBeta*aGamma / (domainInfo.deltaX*domainInfo.deltaX) - tempBeta*bGamma*theta / (beta[topIndex] * domainInfo.deltaX);
					//cout<<levelSet[i+1]<<" "<<levelSet[i+1+1]<<endl;
					//cout<<i+1<<" " <<j+1 <<endl;
					//cout <<fT[matIndex]<<endl;
					//cout<<"theta"<<theta<<endl;
					//cout<<"aGamma : "<<aGamma<<endl;
					//cout<<"bGamma : "<<bGamma<<endl;
					//cout<<"tempBeta : "<<tempBeta<<endl;
					//cout<<endl;
				}
				else
				{
					fT[matIndex] = 0;
				}

				b[matIndex] = -f[matIndex] - fR[matIndex] - fL[matIndex] - fB[matIndex] - fT[matIndex];
			}

		}

		//ofstream ffl, ffr, ffb, fft, bbb, fff;
		//ffl.open("E:\Data/fl.txt");
		//ffr.open("E:\Data/fr.txt");
		//ffb.open("E:\Data/fb.txt");
		//fft.open("E:\Data/ft.txt");
		//bbb.open("E:\Data/b.txt");
		//fff.open("E:\Data/f.txt");

		//for (int j = 0; j < domainInfo.numMatY; j++)
		//{
		//	for (int i = 0; i < domainInfo.numMatX; i++)
		//	{
		//		ffl<<fL[i+j*domainInfo.numMatX]<< " ";
		//		ffr<<fB[i+j*domainInfo.numMatX]<< " ";
		//		ffb<<fB[i+j*domainInfo.numMatX]<< " ";
		//		fft<<fT[i+j*domainInfo.numMatX]<< " ";
		//		bbb<<-b[i+j*domainInfo.numMatX]<< " ";
		//		fff<<f[i+j*domainInfo.numMatX]<< " ";
		//	}
		//		ffl<< endl;
		//		ffr<< endl;
		//		ffb<< endl;
		//		fft<< endl;
		//		bbb<< endl;
		//		fff<< endl;
		//}
		//ffl.close();
		//ffr.close();
		//ffb.close();
		//fft.close();
		//bbb.close();
		//fff.close();

		delete[] fR, fL, fT, fB, normalBottom, normalLeft, normalRight, normalTop;
		return b;
	}
}