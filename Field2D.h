#pragma once

#include "Grid2D.h"
#include "Array2D.h"

template<class TT>
class Field2D
{
public:
	Grid2D grid;
	Array2D<TT> dataArray;

	int iRes, jRes;
	int iStart, Start, iEnd, jEnd;
	double xMin, yMin, xMax, yMax;
	double xLength, yLength;
	double dx, dy;
	double twodx, twody;
	double dx2, dy2;
	double oneOverdx, oneOverdy;
	double oneOver2dx, oneOver2dy;
	double oneOverdx2, oneOverdy2;

	Field2D();
	~Field2D();

	Field2D(const Grid2D& ipGrid);
	Field2D(const double& ipXMin, const double& ipXmax, const int& ipiRes, const double& ipYMin, const double& ipYmax, const int& ipjRes);
	Field2D(const double & ipXMin, const double & ipXmax, const int & ipiStart, const int & ipiRes, const double & ipYMin, const double & ipYmax, const int & ipjStart, const int & ipjRes);

	void initialize(const Grid2D& ipGrid);
	void initialize(const double & ipXMin, const double & ipXmax, const int & ipiStart, const int & ipiRes, const double & ipYMin, const double & ipYmax, const int & ipjStart, const int & ipjRes);
	

private:

};

template<class TT>
Field2D<TT>::Field2D()
{
}

template<class TT>
Field2D<TT>::~Field2D()
{
}

template<class TT>
inline Field2D<TT>::Field2D(const Grid2D & ipGrid)
{
	grid = ipGrid;
	dataArray(grid);
}

template<class TT>
inline Field2D<TT>::Field2D(const double & ipXMin, const double & ipXmax, const int & ipiRes, const double & ipYMin, const double & ipYmax, const int & ipjRes)
{
	grid.initialize(ipXMin, ipXmax, 0, ipiRes, ipYMin, ipYmax, 0, ipjRes);
	initialize(grid);
	dataArray(grid);

}

template<class TT>
inline Field2D<TT>::Field2D(const double & ipXMin, const double & ipXmax, const int & ipiStart, const int & ipiRes, const double & ipYMin, const double & ipYmax, const int & ipjStart, const int & ipjRes)
{
	grid.initialize(ipXMin, ipXmax, ipiStart, ipiRes, ipYMin, ipYmax, ipjStart, ipjRes);
	initialize(grid);
	dataArray(grid);
}


template<class TT>
inline void Field2D<TT>::initialize(const Grid2D & ipGrid)
{
	initialize(ipGrid.xMin, ipGrid.xMax, ipGrid.iStart, ipGrid.iRes, ipGrid.yMin, ipGrid.yMax, ipGrid.iStart, ipGrid.iRes);
}

template<class TT>
inline void Field2D<TT>::initialize(const double & ipXMin, const double & ipXmax, const int & ipiStart, const int & ipiRes, const double & ipYMin, const double & ipYmax, const int & ipjStart, const int & ipjRes)
{
	iRes = ipiRes;
	jRes = ipjRes;
	iStart = ipiStart;
	jStart = ipjStart;
	iEnd = iStart + iRes - 1;
	jEnd = jStart + jRes - 1;
	xMin = ipXMin;
	yMin = ipYMin;
	xMax = ipXmax;
	yMax = ipYmax;
	xLength = xMax - xMin;
	yLength = yMax - yMin;
	dx = xLength / (double)iRes;
	dy = yLength / (double)jRes;
	twodx = 2.0 * dx;
	twody = 2.0*dy;
	dx2 = dx*dx;
	dy2 = dy*dy;
	oneOverdx = 1.0 / dx;
	oneOverdy = 1.0 / dy;
	oneOver2dx = 1.0 / twodx;
	oneOver2dy = 1.0 / twody;
	oneOverdx2 = oneOverdx*oneOverdx;
	oneOverdy2 = oneOverdy*oneOverdy;
}

