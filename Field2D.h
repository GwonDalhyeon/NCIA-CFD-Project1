#pragma once

#include "Grid2D.h"
#include "Array2D.h"
#include "VectorND.h"

template<class TT>
class Field2D
{
public:
	Grid2D grid;
	Array2D<TT> dataArray;

	int iRes, jRes;
	int iStart, jStart, iEnd, jEnd;
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
	
	inline TT& operator [](const int& i) const
	{
		assert(i >= 0 && i < iRes*jRes);
		return dataArray(i);
	}

	const int index(const Vector2D<int>& ipVector) const
	{
		assert(ipVector[0] >= iStart && ipVector[0] <= iEnd);
		assert(ipVector[1] >= jStart && ipVector[1] <= jEnd);
		return dataArray.index(ipVector);
	}

	const int index(const int& i, const int& j) const
	{
		assert(i >= iStart && i <= iEnd);
		assert(j >= jStart && j <= jEnd);
		return dataArray.index(i, j);
	}

	inline TT& operator ()(const int& i) const
	{
		assert(i >= iStart && i <= iEnd);
		return dataArray(i);
	}

	inline TT& operator ()(const Vector2D<int>& ipVector) const
	{
		assert(ipVector[0] >= iStart && ipVector[0] <= iEnd);
		assert(ipVector[1] >= jStart && ipVector[1] <= jEnd);
		return dataArray(ipVector[0],ipVector[1]);
	}

	inline TT& operator ()(const int& i, const int& j) const
	{
		assert(i >= iStart);
		assert(i <= iEnd);
		assert(j >= jStart);
		assert(j <= jEnd);
		return dataArray(i, j);
	}

	inline TT operator ()(const double& x, const double& y) const
	{
		assert(x >= xMin && x <= xMax);
		assert(y >= yMin && y <= yMax);
		//TT& a=x;// = interpolation(x, y);
		Vector2D<double> xy(x, y);
		Vector2D<int> cell = containedCell(x, y);
		TT value00, value10, value01, value11;
		value00 = dataArray(cell);
		value10 = dataArray(cell.i + 1, cell.j);
		value01 = dataArray(cell.i, cell.j + 1);
		value11 = dataArray(cell.i + 1, cell.j + 1);

		double distance00, distance10, distance01, distance11;
		distance00 = (grid(cell) - xy).magnitude();
		distance10 = (grid(cell.i + 1, cell.j) - xy).magnitude();
		distance01 = (grid(cell.i, cell.j + 1) - xy).magnitude();
		distance11 = (grid(cell.i + 1, cell.j + 1) - xy).magnitude();

		return ((dataArray(cell)*distance00 + dataArray(cell.i + 1, cell.j)*distance10 + dataArray(cell.i, cell.j + 1)*distance01 + dataArray(cell.i + 1, cell.j + 1)*distance11) / (distance00 + distance01 + distance10 + distance11));
		//return dataArray(1,1);
	}

	//inline TT& operator ()(const Vector2D <double>& ipVector) const
	//{
	//	assert(ipVector.x >= xMin && ipVector.x <= xMax);
	//	assert(ipVector.y >= yMin && ipVector.y <= yMax);

	//	return interpolation(ipVector.x, ipVector.y);
	//}

	inline TT interpolation(const double& x, const double& y) const;
	inline TT interpolation(const Vector2D<double>& ipVector) const;

	inline Vector2D<int> containedCell(const double& x, const double& y) const;

	inline Vector2D<double> gradient(const int& i, const int& j);

	inline TT minmod(const TT& constant1, const TT constant2) const;

	// Derivative
	inline TT dxxPhi(const int& i, const int& j) const;
	inline TT dxPhi(const int& i, const int& j) const;
	inline TT dxPlusPhi(const int& i, const int& j) const;
	inline TT dxMinusPhi(const int& i, const int& j) const;

	inline TT dyyPhi(const int& i, const int& j) const;
	inline TT dyPhi(const int& i, const int& j) const;
	inline TT dyPlusPhi(const int& i, const int& j) const;
	inline TT dyMinusPhi(const int& i, const int& j) const;

	inline TT dxyPhi(const int& i, const int& j) const;

	inline TT dxPlusPhiSubcell(const int& i, const int& j) const;
	inline TT dxMinusPhiSubcell(const int& i, const int& j) const;
	inline TT dyPlusPhiSubcell(const int& i, const int& j) const;
	inline TT dyMinusPhiSubcell(const int& i, const int& j) const;

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
	initialize(grid);
	dataArray=Array2D<TT>(grid);
}

template<class TT>
inline Field2D<TT>::Field2D(const double & ipXMin, const double & ipXmax, const int & ipiRes, const double & ipYMin, const double & ipYmax, const int & ipjRes)
{
	grid.initialize(ipXMin, ipXmax, 0, ipiRes, ipYMin, ipYmax, 0, ipjRes);
	initialize(grid);
	dataArray = Array2D<TT>(grid);
}

template<class TT>
inline Field2D<TT>::Field2D(const double & ipXMin, const double & ipXmax, const int & ipiStart, const int & ipiRes, const double & ipYMin, const double & ipYmax, const int & ipjStart, const int & ipjRes)
{
	grid.initialize(ipXMin, ipXmax, ipiStart, ipiRes, ipYMin, ipYmax, ipjStart, ipjRes);
	initialize(grid);
	dataArray = Array2D<TT>(grid);
}


template<class TT>
inline void Field2D<TT>::initialize(const Grid2D & ipGrid)
{
	grid = ipGrid;
	initialize(ipGrid.xMin, ipGrid.xMax, ipGrid.iStart, ipGrid.iRes, ipGrid.yMin, ipGrid.yMax, ipGrid.iStart, ipGrid.iRes);
	dataArray = Array2D<TT>(grid);
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
	dx = xLength / (double)(iRes - 1);
	dy = yLength / (double)(jRes - 1);
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

template<class TT>
inline TT Field2D<TT>::interpolation(const double & x, const double & y) const
{
	assert(x >= xMin && x <= xMax);
	assert(y >= yMin && y <= yMax);

	Vector2D<double> xy(x, y);
	Vector2D<int> cell = containedCell(x, y);

	double distance00, distance10, distance01, distance11;
	distance00 = (grid(cell) - xy).magnitude();
	distance10 = (grid(cell.i + 1, cell.j) - xy).magnitude();
	distance01 = (grid(cell.i, cell.j + 1) - xy).magnitude();
	distance11 = (grid(cell.i + 1, cell.j + 1) - xy).magnitude();

	return ((dataArray(cell)*distance00 + dataArray(cell.i + 1, cell.j)*distance10 + dataArray(cell.i, cell.j + 1)*distance01 + dataArray(cell.i + 1, cell.j + 1)*distance11) / (distance00 + distance01 + distance10 + distance11));
}

template<class TT>
inline TT Field2D<TT>::interpolation(const Vector2D<double>& ipVector) const
{
	return interpolation(ipVector.x, ipVector.y);
}

template<class TT>
inline Vector2D<int> Field2D<TT>::containedCell(const double & x, const double & y) const
{
	return Vector2D<int>(floor((x - xMin)*oneOverdx), floor((y - yMin)*oneOverdy));
}

template<class TT>
inline Vector2D<double> Field2D<TT>::gradient(const int & i, const int & j)
{
	return Vector2D<double>(dxPhi(i, j), dyPhi(i, j));
}

template<class TT>
inline TT Field2D<TT>::minmod(const TT & constant1, const TT constant2) const
{
	if (constant1*constant2<0)
	{
		return 0;
	}
	else if (abs(constant1) >= abs(constant2))
	{
		return constant2;
	}
	else
	{
		return constant1;
	}
	return TT();
}

template<class TT>
inline TT Field2D<TT>::dxxPhi(const int & i, const int & j) const
{
	assert(i >= grid.iStart && i <= grid.iEnd);
	assert(j >= grid.jStart && j <= grid.jEnd);

	if (i > grid.iStart && i < grid.iEnd)
	{
		return (dataArray(i + 1, j) - 2 * dataArray(i, j) + dataArray(i - 1, j))*grid.oneOver2dx;
	}
	else if (i == grid.iStart)
	{
		return (dataArray(i, j) - 2 * dataArray(i + 1, j) + dataArray(i + 2, j))*grid.oneOver2dx;
	}
	else
	{
		return (dataArray(i - 2, j) - 2 * dataArray(i - 1, j) + dataArray(i, j))*grid.oneOver2dx;
	}
}

template<class TT>
inline TT Field2D<TT>::dxPhi(const int & i, const int & j) const
{
	assert(i >= grid.iStart && i <= grid.iEnd);
	assert(j >= grid.jStart && j <= grid.jEnd);

	if (i > grid.iStart && i < grid.iEnd)
	{
		return (dataArray(i + 1, j) - dataArray(i - 1, j))*grid.oneOver2dx;// -grid.dx / 2.0*minmod(dxxPhi(i, j), dxxPhi(i + 1, j));
	}
	else if (i == grid.iStart)
	{
		return dxPlusPhi(i, j);
	}
	else
	{
		return dxMinusPhi(i, j);
	}
}

template<class TT>
inline TT Field2D<TT>::dxPlusPhi(const int & i, const int & j) const
{
	assert(i >= grid.iStart && i <= grid.iEnd);
	assert(j >= grid.jStart && j <= grid.jEnd);

	if (i < grid.iEnd)
	{
		return (dataArray(i + 1, j) - dataArray(i, j))*grid.oneOverdx;
	}
	else
	{
		return (dataArray(i, j) - dataArray(i - 1, j))*grid.oneOverdx;
	}
}

template<class TT>
inline TT Field2D<TT>::dxMinusPhi(const int & i, const int & j) const
{
	assert(i >= grid.iStart && i <= grid.iEnd);
	assert(j >= grid.jStart && j <= grid.jEnd);

	if (i > grid.iStart)
	{
		return (dataArray(i, j) - dataArray(i - 1, j))*grid.oneOverdx;
	}
	else
	{
		return (dataArray(i + 1, j) - dataArray(i, j))*grid.oneOverdx;
	}
}

template<class TT>
inline TT Field2D<TT>::dyyPhi(const int & i, const int & j) const
{
	assert(i >= grid.iStart && i <= grid.iEnd);
	assert(j >= grid.jStart && j <= grid.jEnd);

	if (j > grid.jStart && j < grid.jEnd)
	{
		return (dataArray(i, j + 1) - 2 * dataArray(i, j) + dataArray(i, j - 1))*grid.oneOver2dy;
	}
	else if (j == grid.jStart)
	{
		return (dataArray(i, j) - 2 * dataArray(i, j + 1) + dataArray(i, j + 2))*grid.oneOver2dy;
	}
	else
	{
		return (dataArray(i, j - 2) - 2 * dataArray(i, j - 1) + dataArray(i, j))*grid.oneOver2dy;
	}
}

template<class TT>
inline TT Field2D<TT>::dyPhi(const int & i, const int & j) const
{
	assert(i >= grid.iStart && i <= grid.iEnd);
	assert(j >= grid.jStart && j <= grid.jEnd);

	if (j > grid.jStart && j < grid.jEnd)
	{
		return (dataArray(i, j + 1) - dataArray(i, j - 1))*grid.oneOver2dy;// -grid.dy / 2.0*minmod(dyyPhi(i, j), dyyPhi(i, j + 1));
	}
	else if (j == grid.jStart)
	{
		return dyPlusPhi(i, j);
	}
	else
	{
		return dyMinusPhi(i, j);
	}
}

template<class TT>
inline TT Field2D<TT>::dyPlusPhi(const int & i, const int & j) const
{
	assert(i >= grid.iStart && i <= grid.iEnd);
	assert(j >= grid.jStart && j <= grid.jEnd);

	if (j < grid.jEnd)
	{
		return (dataArray(i, j + 1) - dataArray(i, j))*grid.oneOverdy;
	}
	else
	{
		return (dataArray(i, j) - dataArray(i, j - 1))*grid.oneOverdy;
	}
}

template<class TT>
inline TT Field2D<TT>::dyMinusPhi(const int & i, const int & j) const
{
	assert(i >= grid.iStart && i <= grid.iEnd);
	assert(j >= grid.jStart && j <= grid.jEnd);

	if (j > grid.jStart)
	{
		return (dataArray(i, j) - dataArray(i, j - 1))*grid.oneOverdy;
	}
	else
	{
		return (dataArray(i, j + 1) - dataArray(i, j))*grid.oneOverdy;
	}
}

template<class TT>
inline TT Field2D<TT>::dxyPhi(const int & i, const int & j) const
{
	assert(i >= grid.iStart && i <= grid.iEnd);
	assert(j >= grid.jStart && j <= grid.jEnd);

	if (j > grid.jStart && j < grid.jEnd)
	{
		return (dxPhi(i, j + 1) - dxPhi(i, j - 1))*grid.oneOver2dy;
	}
	else if (j == grid.iStart)
	{
		return (dxPhi(i, j + 1) - dxPhi(i, j))*grid.oneOverdy;
	}
	else
	{
		return (dxPhi(i, j) - dxPhi(i, j - 1))*grid.oneOverdy;
	}
}

template<class TT>
inline TT Field2D<TT>::dxPlusPhiSubcell(const int & i, const int & j) const
{
	return TT();
}

template<class TT>
inline TT Field2D<TT>::dxMinusPhiSubcell(const int & i, const int & j) const
{
	return TT();
}

template<class TT>
inline TT Field2D<TT>::dyPlusPhiSubcell(const int & i, const int & j) const
{
	return TT();
}

template<class TT>
inline TT Field2D<TT>::dyMinusPhiSubcell(const int & i, const int & j) const
{
	return TT();
}


