#pragma once
#include "GridInfomation.h"

class LevelSet
{
public:
	GridInfo* grid;
	double* phi;
	double test;
	LevelSet(GridInfo& inputGrid);
	~LevelSet();

	int index(int i, int j);
	int index(int i, int j, int k);
private:

};

LevelSet::LevelSet(GridInfo& inputGrid)
{
	GridInfo* grid = &GridInfo(inputGrid);

	

	if (grid->dimension==1)
	{
		double* phi = new double[grid->numX];
	}
	if (grid->dimension==2)
	{
		double* phi = new double[grid->numX*grid->numY];
	}
	if (grid->dimension==3)
	{
		double* phi = new double[grid->numX*grid->numY*grid->numZ];
	}
	cout << grid->deltaX << "\n";
}

LevelSet::~LevelSet()
{
	delete grid;
	delete phi;

}

inline int LevelSet::index(int i, int j)
{
	return i + j*grid->numX;
}

inline int LevelSet::index(int i, int j, int k)
{
	return i + j*grid->numX + k*grid->numX*grid->numY;
}
