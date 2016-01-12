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

	inline void computeNormal();
	inline Vector2D<double> computeNormal(const int& i, const int& j);
	inline Vector2D<double> computeNormal(const Vector2D<int> ipVector);

	inline void computeUnitNormal();
	inline Vector2D<double> computeUnitNormal(const int& i, const int& j);
	inline Vector2D<double> computeUnitNormal(const Vector2D<int> ipVector);
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
