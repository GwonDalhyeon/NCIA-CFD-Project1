#pragma once
#include "GridInfomation.h"

class LevelSet
{
public:
	GridInfo grid;
	double* phi;
	double test;

	LevelSet(const GridInfo& inputGrid);
	~LevelSet();

	int index(int i, int j);
	int index(int i, int j, int k);
private:

};

LevelSet::LevelSet(const GridInfo& inputGrid)
{
	grid = inputGrid;
	//grid = GridInfo(inputGrid.X0, inputGrid.X1, inputGrid.numX);

	if (grid.dimension==1)
	{
		phi = new double[grid.numX];
	}
	if (grid.dimension==2)
	{
		phi = new double[grid.numX*grid.numY];
	}
	if (grid.dimension==3)
	{
		phi = new double[grid.numX*grid.numY*grid.numZ];
	}
	cout << grid.deltaX << "\n";
	test = 1;
}

LevelSet::~LevelSet()
{
	//delete grid;
	delete phi;

}

inline int LevelSet::index(int i, int j)
{
	return i + j*grid.numX;
}

inline int LevelSet::index(int i, int j, int k)
{
	return i + j*grid.numX + k*grid.numX*grid.numY;
}
