#pragma once
#include "GridInfomation.h"

class LevelSet
{
public:
	GridInfo grid;
	double* phi;

	LevelSet();
	LevelSet(const GridInfo& inputGrid);
	~LevelSet();

	int index(int i, int j);
	int index(int i, int j, int k);

	double unitNormal(int node);
	void unitNormal(int node1, int node2, double* normal);

	LevelSet& operator=(const LevelSet& originLevelSet)
	{
		grid = originLevelSet.grid;

		if (phi!=nullptr)
		{
			delete phi;
		}

		if (grid.dimension==1)
		{
			phi = new double[grid.numX];
			for (int i = 0; i < grid.numX; i++)
			{
				phi[i] = originLevelSet.phi[i];
			}
		}
		else if (grid.dimension==2)
		{
			phi = new double[grid.numX*grid.numY];
			for (int i = 0; i < grid.numX*grid.numY; i++)
			{
				phi[i] = originLevelSet.phi[i];
			}
		}
		else
		{
			phi = new double[grid.numX*grid.numY*grid.numZ];
			for (int i = 0; i < grid.numX*grid.numY*grid.numZ; i++)
			{
				phi[i] = originLevelSet.phi[i];
			}
		}


		return *this;
	}
private:

};

inline LevelSet::LevelSet()
{
}

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

inline double LevelSet::unitNormal(int node)
{
	int i = node;
	if (i<grid.numX - 1 && i>0)
	{
		return (phi[i + 1] - phi[i - 1]) / abs(phi[i + 1] - phi[i - 1]);
	}
	else if (i == 0)
	{
		return (phi[i + 1] - phi[i]) / abs(phi[i + 1] - phi[i]);
	}
	else if (i == grid.numX - 1)
	{
		return (phi[i] - phi[i - 1]) / abs(phi[i] - phi[i - 1]);
	}
	return 0;
}

inline void LevelSet::unitNormal(int node1, int node2, double* normal)
{
	int i = node1;
	int j = node2;

	if (i<grid.numX - 1 && i>0)
	{
		normal[0] = (phi[(i + 1) + j*grid.numX] - phi[i - 1 + j*grid.numX]) / (abs(phi[i + 1 + j*grid.numX] - phi[i - 1 + j*grid.numX]) + DBL_EPSILON);
	}
	if (i == 0)
	{
		normal[0] = (phi[(i + 1) + j*grid.numX] - phi[i + j*grid.numX]) / (abs(phi[i + 1 + j*grid.numX] - phi[i + j*grid.numX]) + DBL_EPSILON);
	}
	if (i == grid.numX - 1)
	{
		normal[0] = (phi[i + j*grid.numX] - phi[i - 1 + j*grid.numX]) / (abs(phi[i + j*grid.numX] - phi[i - 1 + j*grid.numX]) + DBL_EPSILON);
	}

	if (j<grid.numX - 1 && j>0)
	{
		normal[1] = (phi[i + (j + 1)*grid.numX] - phi[i + (j - 1)*grid.numX]) / (abs(phi[i + (j + 1)*grid.numX] - phi[i + (j - 1)*grid.numX]) + DBL_EPSILON);
	}
	if (j == 0)
	{
		normal[1] = (phi[i + (j + 1)*grid.numX] - phi[i + (j)*grid.numX]) / (abs(phi[i + (j + 1)*grid.numX] - phi[i + (j)*grid.numX]) + DBL_EPSILON);
	}
	if (j == grid.numX - 1)
	{
		normal[1] = (phi[i + (j)*grid.numX] - phi[i + (j - 1)*grid.numX]) / (abs(phi[i + (j)*grid.numX] - phi[i + (j - 1)*grid.numX]) + DBL_EPSILON);
	}
}
