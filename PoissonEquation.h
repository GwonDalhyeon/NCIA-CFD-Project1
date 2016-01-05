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
	double* beta;
	double* f;
	double* jCondition1;
	double* jCondition2;

	GridInfo grid;
	LevelSet levelSet;

	PoissonEquationSolver();
	~PoissonEquationSolver();
	void generatePoissonMatrix();
	void generatePoissonVector();

private:

};

PoissonEquationSolver::PoissonEquationSolver()
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

PoissonEquationSolver::~PoissonEquationSolver()
{
	delete poissonMatrix, poissonVector, beta, f, jCondition1, jCondition2;
}

inline void PoissonEquationSolver::generatePoissonMatrix()
{
	double tempBeta = 0;
	if (grid.dimension == 1)
	{
		for (int i = 0; i < grid.numMatX; i++)
		{
			if (i>0 && i<grid.numMatX - 1)
			{
				if ((levelSet.phi[i + 1]>0 && levelSet.phi[i + 1 + 1] <= 0) || (levelSet.phi[i + 1] <= 0 && levelSet.phi[i + 1 + 1]>0))
				{
					tempBeta = beta[i + 1] * beta[i + 1 + 1] * (abs(levelSet.phi[i + 1]) + abs(levelSet.phi[i + 1 + 1])) / (beta[i + 1 + 1] * abs(levelSet.phi[i + 1]) + beta[i + 1] * abs(levelSet.phi[i + 1 + 1]));
					poissonMatrix[i*grid.numMatX + i - 1] = -1 / (grid.deltaX*grid.deltaX)*(beta[i + 1] + beta[i]) / 2;
					poissonMatrix[i*grid.numMatX + i] = 1 / (grid.deltaX*grid.deltaX)*(beta[i + 1] + beta[i]) / 2 + 1 / (grid.deltaX*grid.deltaX)*(tempBeta);
					poissonMatrix[i*grid.numMatX + i + 1] = -1 / (grid.deltaX*grid.deltaX)*tempBeta;
					//cout<<i<<" "<< tempBeta<<endl;
				}
				else if ((levelSet.phi[i]>0 && levelSet.phi[i + 1] <= 0) || (levelSet.phi[i] <= 0 && levelSet.phi[i + 1]>0))
				{
					tempBeta = beta[i] * beta[i + 1] * (abs(levelSet.phi[i]) + abs(levelSet.phi[i + 1])) / (beta[i + 1] * abs(levelSet.phi[i]) + beta[i] * abs(levelSet.phi[i + 1]));
					poissonMatrix[i*grid.numMatX + i - 1] = -1 / (grid.deltaX*grid.deltaX)*tempBeta;
					poissonMatrix[i*grid.numMatX + i] = 1 / (grid.deltaX*grid.deltaX)*(beta[i + 1] + beta[i + 1 + 1]) / 2 + 1 / (grid.deltaX*grid.deltaX)*(tempBeta);
					poissonMatrix[i*grid.numMatX + i + 1] = -1 / (grid.deltaX*grid.deltaX)*(beta[i + 1] + beta[i + 1 + 1]) / 2;
					//cout<<i<<" "<< tempBeta<<endl;
				}
				else
				{
					poissonMatrix[i*grid.numMatX + i - 1] = -1 / (grid.deltaX*grid.deltaX)*(beta[i + 1] + beta[i]) / 2;
					poissonMatrix[i*grid.numMatX + i] = 1 / (grid.deltaX*grid.deltaX)*(beta[i + 1] + beta[i]) / 2 + 1 / (grid.deltaX*grid.deltaX)*(beta[i + 1] + beta[i + 1 + 1]) / 2;
					poissonMatrix[i*grid.numMatX + i + 1] = -1 / (grid.deltaX*grid.deltaX)*(beta[i + 1] + beta[i + 1 + 1]) / 2;
				}

			}
			else if (i == 0)
			{
				poissonMatrix[i*grid.numMatX + i] = 1 / (grid.deltaX*grid.deltaX)*beta[i] + 1 / (grid.deltaX*grid.deltaX)*(beta[i] + beta[i + 1]) / 2;
				poissonMatrix[i*grid.numMatX + i + 1] = -1 / (grid.deltaX*grid.deltaX)*(beta[i + 1] + beta[i + 1 + 1]) / 2;
			}
			else
			{
				poissonMatrix[i*grid.numMatX + i - 1] = -1 / (grid.deltaX*grid.deltaX)*(beta[i + 1] + beta[i]) / 2;
				poissonMatrix[i*grid.numMatX + i] = 1 / (grid.deltaX*grid.deltaX)*(beta[i + 1] + beta[i]) / 2 + 1 / (grid.deltaX*grid.deltaX)*beta[i + 1];
			}
		}
	}
	else if (grid.dimension == 2)
	{
		int matIndex, matLeftIndex, matRightIndex, matTopIndex, matBottomIndex;
		int centIndex, leftIndex, rightIndex, topIndex, bottomIndex;

		for (int j = 0; j < grid.numMatY; j++)
		{
			//poissonMatrix[i*grid.numMatX*grid.numMatY + i + j*grid.numMatX*grid.numMatX*grid.numMatY + j*grid.numMatX ] = 0;
			//poissonMatrix[i*grid.numMatX*grid.numMatY + i + j*grid.numMatX*(grid.numMatX*grid.numMatY + 1) ] = 0;
			for (int i = 0; i < grid.numMatX; i++)
			{
				matIndex = i*grid.numMatX*grid.numMatY + i + j*grid.numMatX*(grid.numMatX*grid.numMatY + 1);
				matLeftIndex = i*grid.numMatX*grid.numMatY + i - 1 + j*grid.numMatX*(grid.numMatX*grid.numMatY + 1);
				matRightIndex = i*grid.numMatX*grid.numMatY + i + 1 + j*grid.numMatX*(grid.numMatX*grid.numMatY + 1);
				matBottomIndex = i*grid.numMatX*grid.numMatY + i + j*grid.numMatX*grid.numMatX*grid.numMatY + (j - 1)*grid.numMatX;
				matTopIndex = i*grid.numMatX*grid.numMatY + i + j*grid.numMatX*grid.numMatX*grid.numMatY + (j + 1)*grid.numMatX;

				centIndex = i + 1 + (j + 1)*grid.numX;
				leftIndex = i + (j + 1)*grid.numX;
				rightIndex = i + 1 + 1 + (j + 1)*grid.numX;
				bottomIndex = i + 1 + j*grid.numX;
				topIndex = i + 1 + (j + 1 + 1)*grid.numX;

				if (j>0 && j<grid.numMatY - 1)
				{
					if ((levelSet.phi[centIndex]>0 && levelSet.phi[topIndex] <= 0) || (levelSet.phi[centIndex] <= 0 && levelSet.phi[topIndex]>0))
					{
						tempBeta = beta[centIndex] * beta[topIndex] * (abs(levelSet.phi[centIndex]) + abs(levelSet.phi[topIndex])) / (beta[topIndex] * abs(levelSet.phi[centIndex]) + beta[centIndex] * abs(levelSet.phi[topIndex]));
						poissonMatrix[matBottomIndex] = -1 / (grid.deltaY*grid.deltaY)*(beta[centIndex] + beta[bottomIndex]) / 2;
						poissonMatrix[matIndex] = poissonMatrix[matIndex] + 1 / (grid.deltaY*grid.deltaY)*(beta[centIndex] + beta[bottomIndex]) / 2 + 1 / (grid.deltaY*grid.deltaY)*(tempBeta);
						poissonMatrix[matTopIndex] = -1 / (grid.deltaY*grid.deltaY)*tempBeta;
						//cout<<i<<" "<< tempBeta<<endl;
						//cout<<i<<" " << j<<" "<<poissonMatrix[matBottomIndex]<<" "<< poissonMatrix[matIndex]<<" " <<poissonMatrix[matTopIndex] <<endl;
						//cout<<"";

					}
					else if ((levelSet.phi[bottomIndex]>0 && levelSet.phi[centIndex] <= 0) || (levelSet.phi[bottomIndex] <= 0 && levelSet.phi[centIndex]>0))
					{
						tempBeta = beta[bottomIndex] * beta[centIndex] * (abs(levelSet.phi[bottomIndex]) + abs(levelSet.phi[centIndex])) / (beta[centIndex] * abs(levelSet.phi[bottomIndex]) + beta[bottomIndex] * abs(levelSet.phi[centIndex]));
						poissonMatrix[matBottomIndex] = -1 / (grid.deltaY*grid.deltaY)*tempBeta;
						poissonMatrix[matIndex] = poissonMatrix[matIndex] + 1 / (grid.deltaY*grid.deltaY)*(beta[centIndex] + beta[topIndex]) / 2 + 1 / (grid.deltaY*grid.deltaY)*(tempBeta);
						poissonMatrix[matTopIndex] = -1 / (grid.deltaY*grid.deltaY)*(beta[centIndex] + beta[topIndex]) / 2;
						//cout<<i<<" "<< tempBeta<<endl;
						//cout<<i<<" " << j<<" "<<poissonMatrix[matBottomIndex]<<" "<< poissonMatrix[matIndex]<<" " <<poissonMatrix[matTopIndex] <<endl;
						//cout<<"";
					}
					else
					{
						poissonMatrix[matBottomIndex] = -1 / (grid.deltaY*grid.deltaY)*(beta[centIndex] + beta[bottomIndex]) / 2;
						poissonMatrix[matIndex] = poissonMatrix[matIndex] + 1 / (grid.deltaY*grid.deltaY)*(beta[centIndex] + beta[topIndex]) / 2 + 1 / (grid.deltaY*grid.deltaY)*(beta[centIndex] + beta[bottomIndex]) / 2;
						poissonMatrix[matTopIndex] = -1 / (grid.deltaY*grid.deltaY)*(beta[centIndex] + beta[topIndex]) / 2;

						//cout<<i<<" " << j<<" "<<poissonMatrix[matBottomIndex]<<" "<< poissonMatrix[matIndex]<<" " <<poissonMatrix[matTopIndex] <<endl;
						//cout<<"";
					}

					if (i>0 && i<grid.numMatX - 1)
					{
						if ((levelSet.phi[centIndex]>0 && levelSet.phi[rightIndex] <= 0) || (levelSet.phi[centIndex] <= 0 && levelSet.phi[rightIndex]>0))
						{
							tempBeta = beta[centIndex] * beta[rightIndex] * (abs(levelSet.phi[centIndex]) + abs(levelSet.phi[rightIndex])) / (beta[rightIndex] * abs(levelSet.phi[centIndex]) + beta[centIndex] * abs(levelSet.phi[rightIndex]));
							poissonMatrix[matLeftIndex] = -1 / (grid.deltaX*grid.deltaX)*(beta[centIndex] + beta[leftIndex]) / 2;
							poissonMatrix[matIndex] = poissonMatrix[matIndex] + 1 / (grid.deltaX*grid.deltaX)*(beta[centIndex] + beta[leftIndex]) / 2 + 1 / (grid.deltaX*grid.deltaX)*(tempBeta);
							poissonMatrix[matRightIndex] = -1 / (grid.deltaX*grid.deltaX)*tempBeta;
							//cout<<i<<" " << j<<" "<<poissonMatrix[matLeftIndex]<<" "<< poissonMatrix[matIndex]<<" " <<poissonMatrix[matRightIndex] <<endl;
							//cout<<"";
						}
						else if ((levelSet.phi[leftIndex]>0 && levelSet.phi[centIndex] <= 0) || (levelSet.phi[leftIndex] <= 0 && levelSet.phi[centIndex]>0))
						{
							tempBeta = beta[leftIndex] * beta[centIndex] * (abs(levelSet.phi[leftIndex]) + abs(levelSet.phi[centIndex])) / (beta[centIndex] * abs(levelSet.phi[leftIndex]) + beta[leftIndex] * abs(levelSet.phi[centIndex]));
							poissonMatrix[matLeftIndex] = -1 / (grid.deltaX*grid.deltaX)*tempBeta;
							poissonMatrix[matIndex] = poissonMatrix[matIndex] + 1 / (grid.deltaX*grid.deltaX)*(beta[centIndex] + beta[rightIndex]) / 2 + 1 / (grid.deltaX*grid.deltaX)*(tempBeta);
							poissonMatrix[matRightIndex] = -1 / (grid.deltaX*grid.deltaX)*(beta[centIndex] + beta[rightIndex]) / 2;
							//cout<<i<<" "<< tempBeta<<endl;
							//cout<<i<<" " << j<<" "<<poissonMatrix[matLeftIndex]<<" "<< poissonMatrix[matIndex]<<" " <<poissonMatrix[matRightIndex] <<endl;
							//cout<<"";
						}
						else
						{
							poissonMatrix[matLeftIndex] = -1 / (grid.deltaX*grid.deltaX)*(beta[centIndex] + beta[leftIndex]) / 2;
							poissonMatrix[matIndex] = poissonMatrix[matIndex] + 1 / (grid.deltaX*grid.deltaX)*(beta[centIndex] + beta[rightIndex]) / 2 + 1 / (grid.deltaX*grid.deltaX)*(beta[centIndex] + beta[leftIndex]) / 2;
							poissonMatrix[matRightIndex] = -1 / (grid.deltaX*grid.deltaX)*(beta[centIndex] + beta[rightIndex]) / 2;
							//cout<<i<<" " << j<<" "<<poissonMatrix[matLeftIndex]<<" "<< poissonMatrix[matIndex]<<" " <<poissonMatrix[matRightIndex] <<endl;
							//cout<<"";
						}
					}
					else if (i == 0)
					{
						poissonMatrix[matIndex] = poissonMatrix[matIndex] + 1 / (grid.deltaX*grid.deltaX)*(beta[centIndex] + beta[rightIndex]) / 2 + 1 / (grid.deltaX*grid.deltaX)*(beta[centIndex] + beta[leftIndex]) / 2;
						poissonMatrix[matRightIndex] = -1 / (grid.deltaX*grid.deltaX)*(beta[centIndex] + beta[rightIndex]) / 2;
						//cout<<i<<" " << j<<" "<<0<<" "<< poissonMatrix[matIndex]<<" " <<poissonMatrix[matRightIndex] <<endl;
						//cout<<"";
					}
					else
					{
						poissonMatrix[matLeftIndex] = -1 / (grid.deltaX*grid.deltaX)*(beta[centIndex] + beta[leftIndex]) / 2;
						poissonMatrix[matIndex] = poissonMatrix[matIndex] + 1 / (grid.deltaX*grid.deltaX)*(beta[centIndex] + beta[rightIndex]) / 2 + 1 / (grid.deltaX*grid.deltaX)*(beta[centIndex] + beta[leftIndex]) / 2;
						//cout<<i<<" " << j<<" "<<poissonMatrix[matLeftIndex]<<" "<< poissonMatrix[matIndex]<<" " <<0<<endl;
						//cout<<"";
					}
				}
				else if (j == 0)
				{
					poissonMatrix[matIndex] = poissonMatrix[matIndex] + 1 / (grid.deltaY*grid.deltaY)*(beta[centIndex] + beta[topIndex]) / 2 + 1 / (grid.deltaY*grid.deltaY)*(beta[centIndex] + beta[bottomIndex]) / 2;
					poissonMatrix[matTopIndex] = -1 / (grid.deltaY*grid.deltaY)*(beta[centIndex] + beta[topIndex]) / 2;
					//cout<<i<<" " << j<<" "<<0<<" "<< poissonMatrix[matIndex]<<" " <<poissonMatrix[matTopIndex] <<endl;
					//cout<<"";
					if (i>0 && i<grid.numMatX - 1)
					{
						if ((levelSet.phi[centIndex]>0 && levelSet.phi[rightIndex] <= 0) || (levelSet.phi[centIndex] <= 0 && levelSet.phi[rightIndex]>0))
						{
							tempBeta = beta[centIndex] * beta[rightIndex] * (abs(levelSet.phi[centIndex]) + abs(levelSet.phi[rightIndex])) / (beta[rightIndex] * abs(levelSet.phi[centIndex]) + beta[centIndex] * abs(levelSet.phi[rightIndex]));
							poissonMatrix[matLeftIndex] = -1 / (grid.deltaX*grid.deltaX)*(beta[centIndex] + beta[leftIndex]) / 2;
							poissonMatrix[matIndex] = poissonMatrix[matIndex] + 1 / (grid.deltaX*grid.deltaX)*(beta[centIndex] + beta[leftIndex]) / 2 + 1 / (grid.deltaX*grid.deltaX)*(tempBeta);
							poissonMatrix[matRightIndex] = -1 / (grid.deltaX*grid.deltaX)*tempBeta;
							//cout<<i<<" "<< tempBeta<<endl;
							//cout<<i<<" " << j<<" "<<poissonMatrix[matLeftIndex]<<" "<< poissonMatrix[matIndex]<<" " <<poissonMatrix[matRightIndex] <<endl;
							//cout<<"";
						}
						else if ((levelSet.phi[leftIndex]>0 && levelSet.phi[centIndex] <= 0) || (levelSet.phi[leftIndex] <= 0 && levelSet.phi[centIndex]>0))
						{
							tempBeta = beta[leftIndex] * beta[centIndex] * (abs(levelSet.phi[leftIndex]) + abs(levelSet.phi[centIndex])) / (beta[centIndex] * abs(levelSet.phi[leftIndex]) + beta[leftIndex] * abs(levelSet.phi[centIndex]));
							poissonMatrix[matLeftIndex] = -1 / (grid.deltaX*grid.deltaX)*tempBeta;
							poissonMatrix[matIndex] = poissonMatrix[matIndex] + 1 / (grid.deltaX*grid.deltaX)*(beta[centIndex] + beta[rightIndex]) / 2 + 1 / (grid.deltaX*grid.deltaX)*(tempBeta);
							poissonMatrix[matRightIndex] = -1 / (grid.deltaX*grid.deltaX)*(beta[centIndex] + beta[rightIndex]) / 2;
							//cout<<i<<" "<< tempBeta<<endl;
							//cout<<i<<" " << j<<" "<<poissonMatrix[matLeftIndex]<<" "<< poissonMatrix[matIndex]<<" " <<poissonMatrix[matRightIndex] <<endl;
							//cout<<"";
						}
						else
						{
							poissonMatrix[matLeftIndex] = -1 / (grid.deltaX*grid.deltaX)*(beta[centIndex] + beta[leftIndex]) / 2;
							poissonMatrix[matIndex] = poissonMatrix[matIndex] + 1 / (grid.deltaX*grid.deltaX)*(beta[centIndex] + beta[rightIndex]) / 2 + 1 / (grid.deltaX*grid.deltaX)*(beta[centIndex] + beta[leftIndex]) / 2;
							poissonMatrix[matRightIndex] = -1 / (grid.deltaX*grid.deltaX)*(beta[centIndex] + beta[rightIndex]) / 2;
							//cout<<i<<" " << j<<" "<<poissonMatrix[matLeftIndex]<<" "<< poissonMatrix[matIndex]<<" " <<poissonMatrix[matRightIndex] <<endl;
							//cout<<"";
						}
					}
					else if (i == 0)
					{
						poissonMatrix[matIndex] = poissonMatrix[matIndex] + 1 / (grid.deltaX*grid.deltaX)*(beta[centIndex] + beta[rightIndex]) / 2 + 1 / (grid.deltaX*grid.deltaX)*(beta[centIndex] + beta[leftIndex]) / 2;
						poissonMatrix[matRightIndex] = -1 / (grid.deltaX*grid.deltaX)*(beta[centIndex] + beta[rightIndex]) / 2;
						//cout<<i<<" " << j<<" "<<0<<" "<< poissonMatrix[matIndex]<<" " <<poissonMatrix[matRightIndex] <<endl;
						//cout<<"";
					}
					else
					{
						poissonMatrix[matLeftIndex] = -1 / (grid.deltaX*grid.deltaX)*(beta[centIndex] + beta[leftIndex]) / 2;
						poissonMatrix[matIndex] = poissonMatrix[matIndex] + 1 / (grid.deltaX*grid.deltaX)*(beta[centIndex] + beta[rightIndex]) / 2 + 1 / (grid.deltaX*grid.deltaX)*(beta[centIndex] + beta[leftIndex]) / 2;
						//cout<<i<<" " << j<<" "<<poissonMatrix[matLeftIndex]<<" "<< poissonMatrix[matIndex]<<" " <<0<<endl;
						//cout<<"";
					}
				}
				else
				{
					poissonMatrix[matBottomIndex] = -1 / (grid.deltaY*grid.deltaY)*(beta[centIndex] + beta[bottomIndex]) / 2;
					poissonMatrix[matIndex] = poissonMatrix[matIndex] + 1 / (grid.deltaY*grid.deltaY)*(beta[centIndex] + beta[topIndex]) / 2 + 1 / (grid.deltaY*grid.deltaY)*(beta[centIndex] + beta[bottomIndex]) / 2;
					//cout<<i<<" " << j<<" "<<poissonMatrix[matBottomIndex]<<" "<< poissonMatrix[matIndex]<<" " <<0 <<endl;
					//cout<<"";
					if (i>0 && i<grid.numMatX - 1)
					{
						if ((levelSet.phi[centIndex]>0 && levelSet.phi[rightIndex] <= 0) || (levelSet.phi[centIndex] <= 0 && levelSet.phi[rightIndex]>0))
						{
							tempBeta = beta[centIndex] * beta[rightIndex] * (abs(levelSet.phi[centIndex]) + abs(levelSet.phi[rightIndex])) / (beta[rightIndex] * abs(levelSet.phi[centIndex]) + beta[centIndex] * abs(levelSet.phi[rightIndex]));
							poissonMatrix[matLeftIndex] = -1 / (grid.deltaX*grid.deltaX)*(beta[centIndex] + beta[leftIndex]) / 2;
							poissonMatrix[matIndex] = poissonMatrix[matIndex] + 1 / (grid.deltaX*grid.deltaX)*(beta[centIndex] + beta[leftIndex]) / 2 + 1 / (grid.deltaX*grid.deltaX)*(tempBeta);
							poissonMatrix[matRightIndex] = -1 / (grid.deltaX*grid.deltaX)*tempBeta;
							//cout<<i<<" "<< tempBeta<<endl;
							//cout<<i<<" " << j<<" "<<poissonMatrix[matLeftIndex]<<" "<< poissonMatrix[matIndex]<<" " <<poissonMatrix[matRightIndex] <<endl;
							//cout<<"";
						}
						else if ((levelSet.phi[leftIndex]>0 && levelSet.phi[centIndex] <= 0) || (levelSet.phi[leftIndex] <= 0 && levelSet.phi[centIndex]>0))
						{
							tempBeta = beta[leftIndex] * beta[centIndex] * (abs(levelSet.phi[leftIndex]) + abs(levelSet.phi[centIndex])) / (beta[centIndex] * abs(levelSet.phi[leftIndex]) + beta[leftIndex] * abs(levelSet.phi[centIndex]));
							poissonMatrix[matLeftIndex] = -1 / (grid.deltaX*grid.deltaX)*tempBeta;
							poissonMatrix[matIndex] = poissonMatrix[matIndex] + 1 / (grid.deltaX*grid.deltaX)*(beta[centIndex] + beta[rightIndex]) / 2 + 1 / (grid.deltaX*grid.deltaX)*(tempBeta);
							poissonMatrix[matRightIndex] = -1 / (grid.deltaX*grid.deltaX)*(beta[centIndex] + beta[rightIndex]) / 2;
							//cout<<i<<" "<< tempBeta<<endl;
							//cout<<i<<" " << j<<" "<<poissonMatrix[matLeftIndex]<<" "<< poissonMatrix[matIndex]<<" " <<poissonMatrix[matRightIndex] <<endl;
							//cout<<"";
						}
						else
						{
							poissonMatrix[matLeftIndex] = -1 / (grid.deltaX*grid.deltaX)*(beta[centIndex] + beta[leftIndex]) / 2;
							poissonMatrix[matIndex] = poissonMatrix[matIndex] + 1 / (grid.deltaX*grid.deltaX)*(beta[centIndex] + beta[rightIndex]) / 2 + 1 / (grid.deltaX*grid.deltaX)*(beta[centIndex] + beta[leftIndex]) / 2;
							poissonMatrix[matRightIndex] = -1 / (grid.deltaX*grid.deltaX)*(beta[centIndex] + beta[rightIndex]) / 2;
							//cout<<i<<" " << j<<" "<<poissonMatrix[matLeftIndex]<<" "<< poissonMatrix[matIndex]<<" " <<poissonMatrix[matRightIndex] <<endl;
							//cout<<"";
						}
					}
					else if (i == 0)
					{
						poissonMatrix[matIndex] = poissonMatrix[matIndex] + 1 / (grid.deltaX*grid.deltaX)*(beta[centIndex] + beta[rightIndex]) / 2 + 1 / (grid.deltaX*grid.deltaX)*(beta[centIndex] + beta[leftIndex]) / 2;
						poissonMatrix[matRightIndex] = -1 / (grid.deltaX*grid.deltaX)*(beta[centIndex] + beta[rightIndex]) / 2;
						//cout<<i<<" " << j<<" "<<0<<" "<< poissonMatrix[matIndex]<<" " <<poissonMatrix[matRightIndex] <<endl;
						//cout<<"";
					}
					else
					{
						poissonMatrix[matLeftIndex] = -1 / (grid.deltaX*grid.deltaX)*(beta[centIndex] + beta[leftIndex]) / 2;
						poissonMatrix[matIndex] = poissonMatrix[matIndex] + 1 / (grid.deltaX*grid.deltaX)*(beta[centIndex] + beta[rightIndex]) / 2 + 1 / (grid.deltaX*grid.deltaX)*(beta[centIndex] + beta[leftIndex]) / 2;
						//cout<<i<<" " << j<<" "<<poissonMatrix[matLeftIndex]<<" "<< poissonMatrix[matIndex]<<" " <<0 <<endl;
						//cout<<"";
					}
				}
			}
		}
	}
}

inline void PoissonEquationSolver::generatePoissonVector()
{
	double aGamma = 0;
	double bGamma = 0;
	double theta = 0;
	double tempBeta = 0;

	if (grid.dimension == 1)
	{
		double* fL = new double[grid.numMatX];
		double* fR = new double[grid.numMatX];
		double* b = new double[grid.numMatX];

		double normalLeft, normalCenter, normalRight;

		for (int i = 0; i < grid.numMatX; i++)
		{
			normalLeft = levelSet.unitNormal(i);
			normalCenter = levelSet.unitNormal(i + 1);
			normalRight = levelSet.unitNormal(i + 1 + 1);

			if (levelSet.phi[i]>0 && levelSet.phi[i + 1] <= 0)
			{
				theta = abs(levelSet.phi[i]) / (abs(levelSet.phi[i]) + abs(levelSet.phi[i + 1]));
				aGamma = (jCondition1[i] * abs(levelSet.phi[i + 1]) + jCondition1[i + 1] * abs(levelSet.phi[i])) / (abs(levelSet.phi[i]) + abs(levelSet.phi[i + 1]));
				bGamma = (jCondition2[i] * normalLeft*abs(levelSet.phi[i + 1]) + jCondition2[i + 1] * normalCenter*abs(levelSet.phi[i])) / (abs(levelSet.phi[i]) + abs(levelSet.phi[i + 1]));
				tempBeta = beta[i] * beta[i + 1] * (abs(levelSet.phi[i]) + abs(levelSet.phi[i + 1])) / (beta[i + 1] * abs(levelSet.phi[i]) + beta[i] * abs(levelSet.phi[i + 1]));
				fL[i] = tempBeta*aGamma / (grid.deltaX*grid.deltaX) - tempBeta*bGamma*theta / (beta[i] * grid.deltaX);
				//cout<<levelSet.phi[i]<<" "<<levelSet.phi[i+1]<<endl;
				//cout<<i<<" " <<i+1 <<endl;
				//cout << i<< " "<< fL[i]<<endl;
				//cout<<"theta"<<theta<<endl;
				//cout<<"aGamma : "<<aGamma<<endl;
				//cout<<"bGamma : "<<bGamma<<endl;
				//cout<<"tempBeta : "<<tempBeta<<endl;
				//cout<<endl;
			}
			else if (levelSet.phi[i] <= 0 && levelSet.phi[i + 1]>0)
			{
				theta = abs(levelSet.phi[i]) / (abs(levelSet.phi[i]) + abs(levelSet.phi[i + 1]));
				aGamma = (jCondition1[i] * abs(levelSet.phi[i + 1]) + jCondition1[i + 1] * abs(levelSet.phi[i])) / (abs(levelSet.phi[i]) + abs(levelSet.phi[i + 1]));
				bGamma = (jCondition2[i] * normalLeft*abs(levelSet.phi[i + 1]) + jCondition2[i + 1] * normalCenter*abs(levelSet.phi[i])) / (abs(levelSet.phi[i]) + abs(levelSet.phi[i + 1]));
				tempBeta = beta[i] * beta[i + 1] * (abs(levelSet.phi[i]) + abs(levelSet.phi[i + 1])) / (beta[i + 1] * abs(levelSet.phi[i]) + beta[i] * abs(levelSet.phi[i + 1]));
				fL[i] = -tempBeta*aGamma / (grid.deltaX*grid.deltaX) + tempBeta*bGamma*theta / (beta[i] * grid.deltaX);
				//cout<<levelSet.phi[i]<<" "<<levelSet.phi[i+1]<<endl;
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

			if (levelSet.phi[i + 1] <= 0 && levelSet.phi[i + 1 + 1]>0)
			{
				theta = abs(levelSet.phi[i + 1 + 1]) / (abs(levelSet.phi[i + 1]) + abs(levelSet.phi[i + 1 + 1]));
				aGamma = (jCondition1[i + 1] * abs(levelSet.phi[i + 1 + 1]) + jCondition1[i + 1 + 1] * abs(levelSet.phi[i + 1])) / (abs(levelSet.phi[i + 1]) + abs(levelSet.phi[i + 1 + 1]));
				bGamma = (jCondition2[i + 1] * normalCenter*abs(levelSet.phi[i + 1 + 1]) + jCondition2[i + 1 + 1] * normalRight*abs(levelSet.phi[i + 1])) / (abs(levelSet.phi[i + 1]) + abs(levelSet.phi[i + 1 + 1]));
				tempBeta = beta[i + 1] * beta[i + 1 + 1] * (abs(levelSet.phi[i + 1]) + abs(levelSet.phi[i + 1 + 1])) / (beta[i + 1 + 1] * abs(levelSet.phi[i + 1]) + beta[i + 1] * abs(levelSet.phi[i + 1 + 1]));
				fR[i] = tempBeta*aGamma / (grid.deltaX*grid.deltaX) + tempBeta*bGamma*theta / (beta[i + 1 + 1] * grid.deltaX);
				//cout<<levelSet.phi[i+1]<<" "<<levelSet.phi[i+1+1]<<endl;
				//cout<<i+1<<" " <<i+1+1 <<endl;
				//cout << i<< " "<< fR[i]<<endl;
				//cout<<"theta"<<theta<<endl;
				//cout<<"aGamma : "<<aGamma<<endl;
				//cout<<"bGamma : "<<bGamma<<endl;
				//cout<<"tempBeta : "<<tempBeta<<endl;
				//cout<<endl;

			}
			else if (levelSet.phi[i + 1]>0 && levelSet.phi[i + 1 + 1] <= 0)
			{
				theta = abs(levelSet.phi[i + 1 + 1]) / (abs(levelSet.phi[i + 1]) + abs(levelSet.phi[i + 1 + 1]));
				aGamma = (jCondition1[i + 1] * abs(levelSet.phi[i + 1 + 1]) + jCondition1[i + 1 + 1] * abs(levelSet.phi[i + 1])) / (abs(levelSet.phi[i + 1]) + abs(levelSet.phi[i + 1 + 1]));
				bGamma = (jCondition2[i + 1] * normalCenter*abs(levelSet.phi[i + 1 + 1]) + jCondition2[i + 1 + 1] * normalRight*abs(levelSet.phi[i + 1])) / (abs(levelSet.phi[i + 1]) + abs(levelSet.phi[i + 1 + 1]));
				tempBeta = beta[i + 1] * beta[i + 1 + 1] * (abs(levelSet.phi[i + 1]) + abs(levelSet.phi[i + 1 + 1])) / (beta[i + 1 + 1] * abs(levelSet.phi[i + 1]) + beta[i + 1] * abs(levelSet.phi[i + 1 + 1]));
				fR[i] = -tempBeta*aGamma / (grid.deltaX*grid.deltaX) - tempBeta*bGamma*theta / (beta[i + 1 + 1] * grid.deltaX);
				//cout<<levelSet.phi[i+1]<<" "<<levelSet.phi[i+1+1]<<endl;
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
	}
	else if (grid.dimension == 2)
	{
		double* fL = new double[grid.numMatX*grid.numMatY];
		double* fR = new double[grid.numMatX*grid.numMatY];
		double* fB = new double[grid.numMatX*grid.numMatY];
		double* fT = new double[grid.numMatX*grid.numMatY];
		double* b = new double[grid.numMatX*grid.numMatY];

		double* normalLeft = new double[2];
		double* normalRight = new double[2];
		double* normalCenter = new double[2];
		double* normalBottom = new double[2];
		double* normalTop = new double[2];

		int matIndex;
		int centIndex, leftIndex, rightIndex, topIndex, bottomIndex;

		for (int j = 0; j < grid.numMatY; j++)
		{
			for (int i = 0; i < grid.numMatX; i++)
			{
				matIndex = i + j*grid.numMatX;

				centIndex = i + 1 + (j + 1)*grid.numX;
				leftIndex = i + (j + 1)*grid.numX;
				rightIndex = i + 1 + 1 + (j + 1)*grid.numX;
				bottomIndex = i + 1 + j*grid.numX;
				topIndex = i + 1 + (j + 1 + 1)*grid.numX;

				levelSet.unitNormal(i - 1, j, normalLeft);
				levelSet.unitNormal(i - 1, j, normalLeft);
				levelSet.unitNormal(i + 1, j, normalRight);
				levelSet.unitNormal(i, j, normalCenter);
				levelSet.unitNormal(i, j - 1, normalBottom);
				levelSet.unitNormal(i, j + 1, normalTop);

				if (levelSet.phi[leftIndex]>0 && levelSet.phi[centIndex] <= 0)
				{
					theta = abs(levelSet.phi[leftIndex]) / (abs(levelSet.phi[leftIndex]) + abs(levelSet.phi[centIndex]));
					aGamma = (jCondition1[leftIndex] * abs(levelSet.phi[centIndex]) + jCondition1[centIndex] * abs(levelSet.phi[leftIndex])) / (abs(levelSet.phi[leftIndex]) + abs(levelSet.phi[centIndex]));
					bGamma = (jCondition2[leftIndex] * normalLeft[0] * abs(levelSet.phi[centIndex]) + jCondition2[centIndex] * normalCenter[0] * abs(levelSet.phi[leftIndex])) / (abs(levelSet.phi[leftIndex]) + abs(levelSet.phi[centIndex]));
					tempBeta = beta[leftIndex] * beta[centIndex] * (abs(levelSet.phi[leftIndex]) + abs(levelSet.phi[centIndex])) / (beta[centIndex] * abs(levelSet.phi[leftIndex]) + beta[leftIndex] * abs(levelSet.phi[centIndex]));
					fL[matIndex] = tempBeta*aGamma / (grid.deltaX*grid.deltaX) - tempBeta*bGamma*theta / (beta[leftIndex] * grid.deltaX);
					//cout<<levelSet.phi[i]<<" "<<levelSet.phi[i+1]<<endl;
					//cout<<i<<" " <<i+1 <<endl;
					//cout << i<< " "<< fL[i]<<endl;
					//cout<<"theta"<<theta<<endl;
					//cout<<"aGamma : "<<aGamma<<endl;
					//cout<<"bGamma : "<<bGamma<<endl;
					//cout<<"tempBeta : "<<tempBeta<<endl;
					//cout<<endl;
				}
				else if (levelSet.phi[leftIndex] <= 0 && levelSet.phi[centIndex]>0)
				{
					theta = abs(levelSet.phi[leftIndex]) / (abs(levelSet.phi[leftIndex]) + abs(levelSet.phi[centIndex]));
					aGamma = (jCondition1[leftIndex] * abs(levelSet.phi[centIndex]) + jCondition1[centIndex] * abs(levelSet.phi[leftIndex])) / (abs(levelSet.phi[leftIndex]) + abs(levelSet.phi[centIndex]));
					bGamma = (jCondition2[leftIndex] * normalLeft[0] * abs(levelSet.phi[centIndex]) + jCondition2[centIndex] * normalCenter[0] * abs(levelSet.phi[leftIndex])) / (abs(levelSet.phi[leftIndex]) + abs(levelSet.phi[centIndex]));
					tempBeta = beta[leftIndex] * beta[centIndex] * (abs(levelSet.phi[leftIndex]) + abs(levelSet.phi[centIndex])) / (beta[centIndex] * abs(levelSet.phi[leftIndex]) + beta[leftIndex] * abs(levelSet.phi[centIndex]));
					fL[matIndex] = -tempBeta*aGamma / (grid.deltaX*grid.deltaX) + tempBeta*bGamma*theta / (beta[leftIndex] * grid.deltaX);
					//cout<<levelSet.phi[i]<<" "<<levelSet.phi[i+1]<<endl;
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

				if (levelSet.phi[centIndex] <= 0 && levelSet.phi[rightIndex]>0)
				{
					theta = abs(levelSet.phi[rightIndex]) / (abs(levelSet.phi[centIndex]) + abs(levelSet.phi[rightIndex]));
					aGamma = (jCondition1[centIndex] * abs(levelSet.phi[rightIndex]) + jCondition1[rightIndex] * abs(levelSet.phi[centIndex])) / (abs(levelSet.phi[centIndex]) + abs(levelSet.phi[rightIndex]));
					bGamma = (jCondition2[centIndex] * normalCenter[0] * abs(levelSet.phi[rightIndex]) + jCondition2[rightIndex] * normalRight[0] * abs(levelSet.phi[centIndex])) / (abs(levelSet.phi[centIndex]) + abs(levelSet.phi[rightIndex]));
					tempBeta = beta[centIndex] * beta[rightIndex] * (abs(levelSet.phi[centIndex]) + abs(levelSet.phi[rightIndex])) / (beta[rightIndex] * abs(levelSet.phi[centIndex]) + beta[centIndex] * abs(levelSet.phi[rightIndex]));
					fR[matIndex] = tempBeta*aGamma / (grid.deltaX*grid.deltaX) + tempBeta*bGamma*theta / (beta[rightIndex] * grid.deltaX);
					//cout<<levelSet.phi[i+1]<<" "<<levelSet.phi[i+1+1]<<endl;
					//cout<<i+1<<" " <<i+1+1 <<endl;
					//cout << i<< " "<< fR[i]<<endl;
					//cout<<"theta"<<theta<<endl;
					//cout<<"aGamma : "<<aGamma<<endl;
					//cout<<"bGamma : "<<bGamma<<endl;
					//cout<<"tempBeta : "<<tempBeta<<endl;
					//cout<<endl;

				}
				else if (levelSet.phi[centIndex]>0 && levelSet.phi[rightIndex] <= 0)
				{
					theta = abs(levelSet.phi[rightIndex]) / (abs(levelSet.phi[centIndex]) + abs(levelSet.phi[rightIndex]));
					aGamma = (jCondition1[centIndex] * abs(levelSet.phi[rightIndex]) + jCondition1[rightIndex] * abs(levelSet.phi[centIndex])) / (abs(levelSet.phi[centIndex]) + abs(levelSet.phi[rightIndex]));
					bGamma = (jCondition2[centIndex] * normalCenter[0] * abs(levelSet.phi[rightIndex]) + jCondition2[rightIndex] * normalRight[0] * abs(levelSet.phi[centIndex])) / (abs(levelSet.phi[centIndex]) + abs(levelSet.phi[rightIndex]));
					tempBeta = beta[centIndex] * beta[rightIndex] * (abs(levelSet.phi[centIndex]) + abs(levelSet.phi[rightIndex])) / (beta[rightIndex] * abs(levelSet.phi[centIndex]) + beta[centIndex] * abs(levelSet.phi[rightIndex]));
					fR[matIndex] = -tempBeta*aGamma / (grid.deltaX*grid.deltaX) - tempBeta*bGamma*theta / (beta[rightIndex] * grid.deltaX);
					//cout<<levelSet.phi[i+1]<<" "<<levelSet.phi[i+1+1]<<endl;
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

				if (levelSet.phi[bottomIndex]>0 && levelSet.phi[centIndex] <= 0)
				{
					theta = abs(levelSet.phi[bottomIndex]) / (abs(levelSet.phi[bottomIndex]) + abs(levelSet.phi[centIndex]));
					aGamma = (jCondition1[bottomIndex] * abs(levelSet.phi[centIndex]) + jCondition1[centIndex] * abs(levelSet.phi[bottomIndex])) / (abs(levelSet.phi[bottomIndex]) + abs(levelSet.phi[centIndex]));
					bGamma = (jCondition2[bottomIndex] * normalBottom[1] * abs(levelSet.phi[centIndex]) + jCondition2[centIndex] * normalCenter[1] * abs(levelSet.phi[bottomIndex])) / (abs(levelSet.phi[bottomIndex]) + abs(levelSet.phi[centIndex]));
					tempBeta = beta[bottomIndex] * beta[centIndex] * (abs(levelSet.phi[bottomIndex]) + abs(levelSet.phi[centIndex])) / (beta[centIndex] * abs(levelSet.phi[bottomIndex]) + beta[bottomIndex] * abs(levelSet.phi[centIndex]));
					fB[matIndex] = tempBeta*aGamma / (grid.deltaX*grid.deltaX) - tempBeta*bGamma*theta / (beta[bottomIndex] * grid.deltaX);
					//cout<<levelSet.phi[i]<<" "<<levelSet.phi[i+1]<<endl;
					//cout<<i<<" " <<i+1 <<endl;
					//cout << i<< " "<< fL[i]<<endl;
					//cout<<"theta"<<theta<<endl;
					//cout<<"aGamma : "<<aGamma<<endl;
					//cout<<"bGamma : "<<bGamma<<endl;
					//cout<<"tempBeta : "<<tempBeta<<endl;
					//cout<<endl;
				}
				else if (levelSet.phi[bottomIndex] <= 0 && levelSet.phi[centIndex]>0)
				{
					theta = abs(levelSet.phi[bottomIndex]) / (abs(levelSet.phi[bottomIndex]) + abs(levelSet.phi[centIndex]));
					aGamma = (jCondition1[bottomIndex] * abs(levelSet.phi[centIndex]) + jCondition1[centIndex] * abs(levelSet.phi[bottomIndex])) / (abs(levelSet.phi[bottomIndex]) + abs(levelSet.phi[centIndex]));
					bGamma = (jCondition2[bottomIndex] * normalBottom[1] * abs(levelSet.phi[centIndex]) + jCondition2[centIndex] * normalCenter[1] * abs(levelSet.phi[bottomIndex])) / (abs(levelSet.phi[bottomIndex]) + abs(levelSet.phi[centIndex]));
					tempBeta = beta[bottomIndex] * beta[centIndex] * (abs(levelSet.phi[bottomIndex]) + abs(levelSet.phi[centIndex])) / (beta[centIndex] * abs(levelSet.phi[bottomIndex]) + beta[bottomIndex] * abs(levelSet.phi[centIndex]));
					fB[matIndex] = -tempBeta*aGamma / (grid.deltaX*grid.deltaX) + tempBeta*bGamma*theta / (beta[bottomIndex] * grid.deltaX);
					//cout<<levelSet.phi[i]<<" "<<levelSet.phi[i+1]<<endl;
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

				if (levelSet.phi[centIndex] <= 0 && levelSet.phi[topIndex]>0)
				{
					theta = abs(levelSet.phi[topIndex]) / (abs(levelSet.phi[centIndex]) + abs(levelSet.phi[topIndex]));
					aGamma = (jCondition1[centIndex] * abs(levelSet.phi[topIndex]) + jCondition1[topIndex] * abs(levelSet.phi[centIndex])) / (abs(levelSet.phi[centIndex]) + abs(levelSet.phi[topIndex]));
					bGamma = (jCondition2[centIndex] * normalCenter[1] * abs(levelSet.phi[topIndex]) + jCondition2[topIndex] * normalTop[1] * abs(levelSet.phi[centIndex])) / (abs(levelSet.phi[centIndex]) + abs(levelSet.phi[topIndex]));
					tempBeta = beta[centIndex] * beta[topIndex] * (abs(levelSet.phi[centIndex]) + abs(levelSet.phi[topIndex])) / (beta[topIndex] * abs(levelSet.phi[centIndex]) + beta[centIndex] * abs(levelSet.phi[topIndex]));
					fT[matIndex] = tempBeta*aGamma / (grid.deltaX*grid.deltaX) + tempBeta*bGamma*theta / (beta[topIndex] * grid.deltaX);
					//cout<<levelSet.phi[i+1]<<" "<<levelSet.phi[i+1+1]<<endl;
					//cout<<i+1<<" " <<j+1 <<endl;
					//cout << fT[matIndex]<<endl;
					//cout<<"theta"<<theta<<endl;
					//cout<<"aGamma : "<<aGamma<<endl;
					//cout<<"bGamma : "<<bGamma<<endl;
					//cout<<"tempBeta : "<<tempBeta<<endl;
					//cout<<endl;

				}
				else if (levelSet.phi[centIndex]>0 && levelSet.phi[topIndex] <= 0)
				{
					theta = abs(levelSet.phi[topIndex]) / (abs(levelSet.phi[centIndex]) + abs(levelSet.phi[topIndex]));
					aGamma = (jCondition1[centIndex] * abs(levelSet.phi[topIndex]) + jCondition1[topIndex] * abs(levelSet.phi[centIndex])) / (abs(levelSet.phi[centIndex]) + abs(levelSet.phi[topIndex]));
					bGamma = (jCondition2[centIndex] * normalCenter[1] * abs(levelSet.phi[topIndex]) + jCondition2[topIndex] * normalTop[1] * abs(levelSet.phi[centIndex])) / (abs(levelSet.phi[centIndex]) + abs(levelSet.phi[topIndex]));
					tempBeta = beta[centIndex] * beta[topIndex] * (abs(levelSet.phi[centIndex]) + abs(levelSet.phi[topIndex])) / (beta[topIndex] * abs(levelSet.phi[centIndex]) + beta[centIndex] * abs(levelSet.phi[topIndex]));
					fT[matIndex] = -tempBeta*aGamma / (grid.deltaX*grid.deltaX) - tempBeta*bGamma*theta / (beta[topIndex] * grid.deltaX);
					//cout<<levelSet.phi[i+1]<<" "<<levelSet.phi[i+1+1]<<endl;
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

		//for (int j = 0; j < grid.numMatY; j++)
		//{
		//	for (int i = 0; i < grid.numMatX; i++)
		//	{
		//		ffl<<fL[i+j*grid.numMatX]<< " ";
		//		ffr<<fB[i+j*grid.numMatX]<< " ";
		//		ffb<<fB[i+j*grid.numMatX]<< " ";
		//		fft<<fT[i+j*grid.numMatX]<< " ";
		//		bbb<<-b[i+j*grid.numMatX]<< " ";
		//		fff<<f[i+j*grid.numMatX]<< " ";
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
	}
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

