#pragma once
#include "Field2D.h"

class LevelSet2D
{
public:
	Grid2D grid;

	Field2D<double> phi;
	Field2D<Vector2D<double>> normal;
	Field2D<Vector2D<double>> unitNormal;
	Field2D<Vector2D<double>> tangential;
	Field2D<double> meanCurvature;


	LevelSet2D();
	~LevelSet2D();

	LevelSet2D(const Grid2D& ipGrid);

	const int index(const int& i, const int& j) const
	{
		assert(i >= phi.iStart && i <= phi.iEnd);
		assert(j >= phi.jStart && j <= phi.jEnd);
		return phi.index(i, j);
	}

	const int index(const Vector2D<int>& ipVector) const
	{
		assert(ipVector[0] >= phi.iStart && ipVector[0] <= phi.iEnd);
		assert(ipVector[1] >= phi.jStart && ipVector[1] <= phi.jEnd);
		return phi.index(ipVector);
	}

	inline double& operator ()(const int& i, const int& j) const
	{
		assert(i >= phi.iStart && i <= phi.iEnd);
		assert(j >= phi.jStart && j <= phi.jEnd);
		return phi(i, j);
	}

	inline double& operator ()(const Vector2D<int>& ipVector) const
	{
		assert(ipVector[0] >= phi.iStart && ipVector[0] <= phi.iEnd);
		assert(ipVector[1] >= phi.jStart && ipVector[1] <= phi.jEnd);
		return phi(ipVector);
	}

	inline double operator ()(const double& x, const double& y) const
	{
		assert(x >= grid.xMin && x <= grid.xMax);
		assert(y >= grid.yMin && y <= grid.yMax);

		return phi(x, y);
	}

	inline double operator ()(const Vector2D<double>& ipVector) const
	{
		assert(ipVector[0] >= grid.xMin && ipVector[0] <= grid.xMax);
		assert(ipVector[1] >= grid.yMin && ipVector[1] <= grid.yMax);
		return phi(ipVector.x, ipVector.y);
	}

	inline void computeNormal();
	inline Vector2D<double> computeNormal(const int& i, const int& j);
	inline Vector2D<double> computeNormal(const Vector2D<int> ipVector);

	inline void computeUnitNormal();
	inline Vector2D<double> computeUnitNormal(const int& i, const int& j);
	inline Vector2D<double> computeUnitNormal(const Vector2D<int> ipVector);

	inline void computeMeanCurvature();
	inline double computeMeanCurvature(const int& i, const int& j);
	inline double computeMeanCurvature(const Vector2D<int> ipVector);

	inline Vector2D<double> gradient(const int& i, const int& j);

	inline double interpolation(const double& x, const double& y);
	inline double interpolation(const Vector2D<double>& ipVector) ;

	// Derivative
	inline double dxxPhi(const int& i, const int& j);
	inline double dxPhi(const int& i, const int& j);
	inline double dxPlusPhi(const int& i, const int& j);
	inline double dxMinusPhi(const int& i, const int& j);

	inline double dyyPhi(const int& i, const int& j);
	inline double dyPhi(const int& i, const int& j);
	inline double dyPlusPhi(const int& i, const int& j);
	inline double dyMinusPhi(const int& i, const int& j);

	inline double dxyPhi(const int& i, const int& j);
	
private:

};

LevelSet2D::LevelSet2D()
{
}

LevelSet2D::~LevelSet2D()
{
}

inline LevelSet2D::LevelSet2D(const Grid2D & ipGrid)
{
	grid = ipGrid;
	phi = Field2D<double>(ipGrid);
	normal = Field2D<Vector2D<double>>(ipGrid);
	unitNormal = Field2D<Vector2D<double>>(ipGrid);
	tangential = Field2D<Vector2D<double>>(ipGrid);
	meanCurvature = Field2D<double>(ipGrid);
}


inline void LevelSet2D::computeNormal()
{
	for (int j = phi.jStart; j <= phi.jEnd; j++)
	{
		for (int i = phi.iStart; i <= phi.iEnd; i++)
		{
			normal.dataArray(index(i, j)) = computeNormal(i, j);
		}
	}
}

inline Vector2D<double> LevelSet2D::computeNormal(const int & i, const int & j)
{
	Vector2D<double> normal;

	if (j > phi.jStart && j < phi.jEnd)
	{

		if (i > phi.iStart && i < phi.iEnd)
		{
			normal.values[0] = (phi(i + 1, j) - phi(i - 1, j))*phi.oneOver2dx;
			normal.values[1] = (phi(i, j + 1) - phi(i, j - 1))*phi.oneOver2dy;
		}
		else if (i == phi.iStart)
		{
			normal.values[0] = (phi(i + 1, j) - phi(i, j))*phi.oneOverdx;
			normal.values[1] = (phi(i, j + 1) - phi(i, j - 1))*phi.oneOver2dy;
		}
		else if (i == phi.iEnd)
		{
			normal.values[0] = (phi(i, j) - phi(i - 1, j))*phi.oneOverdx;
			normal.values[1] = (phi(i, j + 1) - phi(i, j - 1))*phi.oneOver2dy;
		}
		else
		{
			cout << "Level set unitNormal error." << endl;
		}
	}
	else if (j == phi.jStart)
	{

		if (i > phi.iStart && i < phi.iEnd)
		{
			normal.values[0] = (phi(i + 1, j) - phi(i - 1, j))*phi.oneOver2dx;;
			normal.values[1] = (phi(i, j + 1) - phi(i, j))*phi.oneOverdy;
		}
		else if (i == phi.iStart)
		{
			normal.values[0] = (phi(i + 1, j) - phi(i, j))*phi.oneOverdx;
			normal.values[1] = (phi(i, j + 1) - phi(i, j))*phi.oneOverdy;
		}
		else if (i == phi.iEnd)
		{
			normal.values[0] = (phi(i, j) - phi(i - 1, j))*phi.oneOverdx;
			normal.values[1] = (phi(i, j + 1) - phi(i, j))*phi.oneOverdy;
		}
		else
		{
			cout << "Level set unitNormal error." << endl;
		}

	}
	else if (j == phi.jEnd)
	{

		if (i > phi.iStart && i < phi.iEnd)
		{
			normal.values[0] = (phi(i + 1, j) - phi(i - 1, j))*phi.oneOver2dx;;
			normal.values[1] = (phi(i, j) - phi(i, j - 1))*phi.oneOverdy;
		}
		else if (i == phi.iStart)
		{
			normal.values[0] = (phi(i + 1, j) - phi(i, j))*phi.oneOverdx;
			normal.values[1] = (phi(i, j) - phi(i, j - 1))*phi.oneOverdy;
		}
		else if (i == phi.iEnd)
		{
			normal.values[0] = (phi(i, j) - phi(i - 1, j))*phi.oneOverdx;
			normal.values[1] = (phi(i, j) - phi(i, j - 1))*phi.oneOverdy;
		}
		else
		{
			cout << "Level set unitNormal error." << endl;
		}

	}
	else
	{
		assert(j >= phi.jStart && j <= phi.jEnd);
	}


	return normal;
}

inline Vector2D<double> LevelSet2D::computeNormal(const Vector2D<int> ipVector)
{
	return computeNormal(ipVector[0], ipVector[1]);
}



inline void LevelSet2D::computeUnitNormal()
{
	for (int j = phi.jStart; j <= phi.jEnd; j++)
	{
		for (int i = phi.iStart; i <= phi.iEnd; i++)
		{
			unitNormal.dataArray(index(i, j)) = computeUnitNormal(i, j);
		}
	}
}

inline Vector2D<double> LevelSet2D::computeUnitNormal(const int & i, const int & j)
{
	Vector2D<double> normal;

	if (j > phi.jStart && j < phi.jEnd)
	{

		if (i > phi.iStart && i < phi.iEnd)
		{
			normal.values[0] = (phi(i + 1, j) - phi(i - 1, j)) / (abs(phi(i + 1, j) - phi(i - 1, j)) + DBL_EPSILON);
			normal.values[1] = (phi(i, j + 1) - phi(i, j - 1)) / (abs(phi(i, j + 1) - phi(i, j - 1)) + DBL_EPSILON);
		}
		else if (i == phi.iStart)
		{
			normal.values[0] = (phi(i + 1, j) - phi(i, j)) / (abs(phi(i + 1, j) - phi(i, j)) + DBL_EPSILON);
			normal.values[1] = (phi(i, j + 1) - phi(i, j - 1)) / (abs(phi(i, j + 1) - phi(i, j - 1)) + DBL_EPSILON);
		}
		else if (i == phi.iEnd)
		{
			normal.values[0] = (phi(i, j) - phi(i - 1, j)) / (abs(phi(i, j) - phi(i - 1, j)) + DBL_EPSILON);
			normal.values[1] = (phi(i, j + 1) - phi(i, j - 1)) / (abs(phi(i, j + 1) - phi(i, j - 1)) + DBL_EPSILON);
		}
		else
		{
			cout << "Level set unitNormal error." << endl;
		}
	}
	else if (j == phi.jStart)
	{

		if (i > phi.iStart && i < phi.iEnd)
		{
			normal.values[0] = (phi(i + 1, j) - phi(i - 1, j)) / (abs(phi(i + 1, j) - phi(i - 1, j)) + DBL_EPSILON);
			normal.values[1] = (phi(i, j + 1) - phi(i, j)) / (abs(phi(i, j + 1) - phi(i, j)) + DBL_EPSILON);
		}
		else if (i == phi.iStart)
		{
			normal.values[0] = (phi(i + 1, j) - phi(i, j)) / (abs(phi(i + 1, j) - phi(i, j)) + DBL_EPSILON);
			normal.values[1] = (phi(i, j + 1) - phi(i, j)) / (abs(phi(i, j + 1) - phi(i, j)) + DBL_EPSILON);
		}
		else if (i == phi.iEnd)
		{
			normal.values[0] = (phi(i, j) - phi(i - 1, j)) / (abs(phi(i, j) - phi(i - 1, j)) + DBL_EPSILON);
			normal.values[1] = (phi(i, j + 1) - phi(i, j)) / (abs(phi(i, j + 1) - phi(i, j)) + DBL_EPSILON);
		}
		else
		{
			cout << "Level set unitNormal error." << endl;
		}

	}
	else if (j == phi.jEnd)
	{

		if (i > phi.iStart && i < phi.iEnd)
		{
			normal.values[0] = (phi(i + 1, j) - phi(i - 1, j)) / (abs(phi(i + 1, j) - phi(i - 1, j)) + DBL_EPSILON);
			normal.values[1] = (phi(i, j) - phi(i, j - 1)) / (abs(phi(i, j) - phi(i, j - 1)) + DBL_EPSILON);
		}
		else if (i == phi.iStart)
		{
			normal.values[0] = (phi(i + 1, j) - phi(i, j)) / (abs(phi(i + 1, j) - phi(i, j)) + DBL_EPSILON);
			normal.values[1] = (phi(i, j) - phi(i, j - 1)) / (abs(phi(i, j) - phi(i, j - 1)) + DBL_EPSILON);
		}
		else if (i == phi.iEnd)
		{
			normal.values[0] = (phi(i, j) - phi(i - 1, j)) / (abs(phi(i, j) - phi(i - 1, j)) + DBL_EPSILON);
			normal.values[1] = (phi(i, j) - phi(i, j - 1)) / (abs(phi(i, j) - phi(i, j - 1)) + DBL_EPSILON);
		}
		else
		{
			cout << "Level set unitNormal error." << endl;
		}

	}
	else
	{
		assert(j >= phi.jStart && j <= phi.jEnd);
	}


	return normal;
}

inline Vector2D<double> LevelSet2D::computeUnitNormal(const Vector2D<int> ipVector)
{
	return computeUnitNormal(ipVector[0], ipVector[1]);
}

inline void LevelSet2D::computeMeanCurvature()
{
	for (int i = grid.iStart; i <= grid.iEnd; i++)
	{
		for (int j = grid.iStart; j <= grid.jEnd; j++)
		{
			meanCurvature(i, j) = computeMeanCurvature(i, j);
		}
	}
}

inline double LevelSet2D::computeMeanCurvature(const int & i, const int & j)
{
	return -(dxxPhi(i, j)*dyPhi(i, j)*dyPhi(i, j) - 2.0*dxyPhi(i, j)*dxPhi(i, j)*dyPhi(i, j) + dyyPhi(i, j)*dxPhi(i, j)*dxPhi(i, j)) / pow(dxPhi(i, j)*dxPhi(i, j) + dyPhi(i, j)*dyPhi(i, j), 3.0 / 2.0);
}

inline double LevelSet2D::computeMeanCurvature(const Vector2D<int> ipVector)
{
	return computeMeanCurvature(ipVector.i,ipVector.j);
}

inline Vector2D<double> LevelSet2D::gradient(const int & i, const int & j)
{
	return Vector2D<double>(dxPhi(i,j),dyPhi(i,j));
}

inline double LevelSet2D::interpolation(const double & x, const double & y)
{
	Vector2D<double> xy(x, y);
	Vector2D<int> cell = phi.containedCell(x, y);

	double distance00, distance10, distance01, distance11;
	distance00 = (grid(cell) - xy).magnitude();
	distance10 = (grid(cell.i + 1, cell.j) - xy).magnitude();
	distance01 = (grid(cell.i, cell.j + 1) - xy).magnitude();
	distance11 = (grid(cell.i + 1, cell.j + 1) - xy).magnitude();

	return ((phi(cell)*distance00 + phi(cell.i + 1, cell.j)*distance10 + phi(cell.i, cell.j + 1)*distance01 + phi(cell.i + 1, cell.j + 1)*distance11) / (distance00 + distance01 + distance10 + distance11));
}

inline double LevelSet2D::interpolation(const Vector2D<double>& ipVector)
{
	return interpolation(ipVector.x, ipVector.y);
}

inline double LevelSet2D::dxxPhi(const int & i, const int & j)
{
	assert(i >= grid.iStart && i <= grid.iEnd);
	assert(j >= grid.jStart && j <= grid.jEnd);

	if (i > grid.iStart && i < grid.iEnd)
	{
		return (phi(i + 1, j) - 2 * phi(i, j) + phi(i - 1, j))*grid.oneOver2dx;
	}
	else if (i == grid.iStart)
	{
		return (phi(i, j) - 2 * phi(i + 1, j) + phi(i + 2, j))*grid.oneOver2dx;
	}
	else
	{
		return (phi(i - 2, j) - 2 * phi(i - 1, j) + phi(i, j))*grid.oneOver2dx;
	}
}

inline double LevelSet2D::dxPhi(const int & i, const int & j)
{
	assert(i >= grid.iStart && i <= grid.iEnd);
	assert(j >= grid.jStart && j <= grid.jEnd);

	if (i > grid.iStart && i < grid.iEnd)
	{
		return (phi(i + 1, j) - phi(i - 1, j))*grid.oneOver2dx;
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

inline double LevelSet2D::dxPlusPhi(const int & i, const int & j)
{
	assert(i >= grid.iStart && i <= grid.iEnd);
	assert(j >= grid.jStart && j <= grid.jEnd);

	if (i < grid.iEnd)
	{
		return (phi(i + 1, j) - phi(i, j))*grid.oneOverdx;
	}
	else
	{
		return (phi(i, j) - phi(i - 1, j))*grid.oneOverdx;
	}
}

inline double LevelSet2D::dxMinusPhi(const int & i, const int & j)
{
	assert(i >= grid.iStart && i <= grid.iEnd);
	assert(j >= grid.jStart && j <= grid.jEnd);

	if (i > grid.iStart)
	{
		return (phi(i, j) - phi(i - 1, j))*grid.oneOverdx;
	}
	else
	{
		return (phi(i + 1, j) - phi(i, j))*grid.oneOverdx;
	}
}

inline double LevelSet2D::dyyPhi(const int & i, const int & j)
{
	assert(i >= grid.iStart && i <= grid.iEnd);
	assert(j >= grid.jStart && j <= grid.jEnd);

	if (j > grid.jStart && j < grid.jEnd)
	{
		return (phi(i, j + 1) - 2 * phi(i, j) + phi(i, j - 1))*grid.oneOver2dy;
	}
	else if (j == grid.jStart)
	{
		return (phi(i, j) - 2 * phi(i, j + 1) + phi(i, j + 2))*grid.oneOver2dy;
	}
	else
	{
		return (phi(i, j - 2) - 2 * phi(i, j - 1) + phi(i, j))*grid.oneOver2dy;
	}
}

inline double LevelSet2D::dyPhi(const int & i, const int & j)
{
	assert(i >= grid.iStart && i <= grid.iEnd);
	assert(j >= grid.jStart && j <= grid.jEnd);

	if (j > grid.jStart && j < grid.jEnd)
	{
		return (phi(i, j + 1) - phi(i, j - 1))*grid.oneOver2dy;
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

inline double LevelSet2D::dyPlusPhi(const int & i, const int & j)
{
	assert(i >= grid.iStart && i <= grid.iEnd);
	assert(j >= grid.jStart && j <= grid.jEnd);

	if (j < grid.jEnd)
	{
		return (phi(i, j + 1) - phi(i, j))*grid.oneOverdy;
	}
	else
	{
		return (phi(i, j) - phi(i, j - 1))*grid.oneOverdy;
	}
}

inline double LevelSet2D::dyMinusPhi(const int & i, const int & j)
{
	assert(i >= grid.iStart && i <= grid.iEnd);
	assert(j >= grid.jStart && j <= grid.jEnd);

	if (j > grid.jStart)
	{
		return (phi(i, j) - phi(i, j - 1))*grid.oneOverdy;
	}
	else
	{
		return (phi(i, j + 1) - phi(i, j))*grid.oneOverdy;
	}
}

inline double LevelSet2D::dxyPhi(const int & i, const int & j)
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
