#pragma once
#include <algorithm>

#include "AdvectionMethod2D.h"
#include "PoissonSolver.h"

template <class TT>
class SurfaceReconst
{
public:
	static double alpha; // = min(delta x, delta y) // This is a parameter for heaviside function and delta function.
	
	//SurfaceReconst() { cout << alpha << endl };
	SurfaceReconst();
	~SurfaceReconst();

	void distance2Data(const Field2D<TT>& ipData, Field2D<TT>& distance);
	static double heaviside(const double& ip);
	static double deltaFt(const double& ip);



private:

};

template <class TT>
SurfaceReconst<TT>::SurfaceReconst()
{
}

template <class TT>
SurfaceReconst<TT>::~SurfaceReconst()
{
}

template<class TT>
inline void SurfaceReconst<TT>::distance2Data(const Field2D<TT>& ipData, Field2D<TT>& distance)
{
	TT xMin, yMin;
	double h = min(distance.dx, distance.dy);

	for (int i = distance.iStart; i < distance.iEnd; i++)
	{
		for (int j = distance.jStart; j < distance.jEnd; j++)
		{
			if (i==distance.iStart)
			{
				xMin = min(distance(i, j), distance(i + 1, j));
			}
			else if (i==distance.iEnd)
			{
				xMin = min(distance(i-1, j), distance(i, j));
			}
			else
			{
				xMin = min(distance(i - 1, j), distance(i + 1, j));
			}

			if (j == distance.jStart)
			{
				yMin = min(distance(i, j), distance(i, j + 1));
			}
			else if (i == distance.iEnd)
			{
				yMin = min(distance(i, j - 1), distance(i, j));
			}
			else
			{
				yMin = min(distance(i, j - 1), distance(i, j + 1));
			}

			if (abs(xMin-yMin) >= h)
			{
				distance(i, j) = min(xMin, yMin) + h;
			}
			else
			{
				u(i, j) = (xMin + yMin + sqrt(2 * h*h - (xMin - yMin)*(xMin - yMin))) / 2.0;
			}
		}
	}

	for (int i = distance.iStart; i < distance.iEnd; i++)
	{
		for (int j = distance.jEnd; j < distance.jStart; j++)
		{
			if (i == distance.iStart)
			{
				xMin = min(distance(i, j), distance(i + 1, j));
			}
			else if (i == distance.iEnd)
			{
				xMin = min(distance(i - 1, j), distance(i, j));
			}
			else
			{
				xMin = min(distance(i - 1, j), distance(i + 1, j));
			}

			if (j == distance.jStart)
			{
				yMin = min(distance(i, j), distance(i, j + 1));
			}
			else if (i == distance.iEnd)
			{
				yMin = min(distance(i, j - 1), distance(i, j));
			}
			else
			{
				yMin = min(distance(i, j - 1), distance(i, j + 1));
			}

			if (abs(xMin - yMin) >= h)
			{
				distance(i, j) = min(xMin, yMin) + h;
			}
			else
			{
				u(i, j) = (xMin + yMin + sqrt(2 * h*h - (xMin - yMin)*(xMin - yMin))) / 2.0;
			}
		}
	}

	for (int i = distance.iEnd; i < distance.iStart; i++)
	{
		for (int j = distance.jStart; j < distance.jEnd; j++)
		{
			if (i == distance.iStart)
			{
				xMin = min(distance(i, j), distance(i + 1, j));
			}
			else if (i == distance.iEnd)
			{
				xMin = min(distance(i - 1, j), distance(i, j));
			}
			else
			{
				xMin = min(distance(i - 1, j), distance(i + 1, j));
			}

			if (j == distance.jStart)
			{
				yMin = min(distance(i, j), distance(i, j + 1));
			}
			else if (i == distance.iEnd)
			{
				yMin = min(distance(i, j - 1), distance(i, j));
			}
			else
			{
				yMin = min(distance(i, j - 1), distance(i, j + 1));
			}

			if (abs(xMin - yMin) >= h)
			{
				distance(i, j) = min(xMin, yMin) + h;
			}
			else
			{
				u(i, j) = (xMin + yMin + sqrt(2 * h*h - (xMin - yMin)*(xMin - yMin))) / 2.0;
			}
		}
	}

	for (int i = distance.iEnd; i < distance.iStart; i++)
	{
		for (int j = distance.jEnd; j < distance.jStart; j++)
		{
			if (i == distance.iStart)
			{
				xMin = min(distance(i, j), distance(i + 1, j));
			}
			else if (i == distance.iEnd)
			{
				xMin = min(distance(i - 1, j), distance(i, j));
			}
			else
			{
				xMin = min(distance(i - 1, j), distance(i + 1, j));
			}

			if (j == distance.jStart)
			{
				yMin = min(distance(i, j), distance(i, j + 1));
			}
			else if (i == distance.iEnd)
			{
				yMin = min(distance(i, j - 1), distance(i, j));
			}
			else
			{
				yMin = min(distance(i, j - 1), distance(i, j + 1));
			}

			if (abs(xMin - yMin) >= h)
			{
				distance(i, j) = min(xMin, yMin) + h;
			}
			else
			{
				u(i, j) = (xMin + yMin + sqrt(2 * h*h - (xMin - yMin)*(xMin - yMin))) / 2.0;
			}
		}
	}

	for (int i = distance.iStart; i < distance.iEnd; i++)
	{
		for (int j = distance.jStart; j < distance.jEnd; j++)
		{
			if (i == distance.iStart)
			{
				xMin = min(distance(i, j), distance(i + 1, j));
			}
			else if (i == distance.iEnd)
			{
				xMin = min(distance(i - 1, j), distance(i, j));
			}
			else
			{
				xMin = min(distance(i - 1, j), distance(i + 1, j));
			}

			if (j == distance.jStart)
			{
				yMin = min(distance(i, j), distance(i, j + 1));
			}
			else if (i == distance.iEnd)
			{
				yMin = min(distance(i, j - 1), distance(i, j));
			}
			else
			{
				yMin = min(distance(i, j - 1), distance(i, j + 1));
			}

			if (abs(xMin - yMin) >= h)
			{
				distance(i, j) = min(xMin, yMin) + h;
			}
			else
			{
				u(i, j) = (xMin + yMin + sqrt(2 * h*h - (xMin - yMin)*(xMin - yMin))) / 2.0;
			}
		}
	}
}

template<class TT>
inline double SurfaceReconst<TT>::heaviside(const double & ip)
{
	if (ip > alpha)
	{
		return 1.0;
	}
	else if (ip < -alpha)
	{
		return 0.0;
	}
	else
	{
		return (1 + ip / alpha + sin(PI*ip / alpha) / PI) / 2.0;
	}
}

template<class TT>
inline double SurfaceReconst<TT>::deltaFt(const double & ip)
{
	if (abs(ip)>alpha)
	{
		return 0.0;
	}
	else
	{
		return (1 + cos(PI*ip / alpha)) / (2 * alpha);
	}
}
