#pragma once
#include <assert.h>
#include <iostream>
#include "Grid2D.h"

template <class TT>
class Array2D
{
public:
	union
	{
		struct { int iStart, jStart, iEnd, jEnd; };
		struct { int ijStart[2], ijEnd[2]; };
	};

	int iRes, jRes;
	int ijRes;
	TT* values;

	Array2D();
	~Array2D();

	Array2D(const int& ipiRes);
	Array2D(const int& ipiRes, const int& ipjRes);
	Array2D(const int& ipiStart, const int& ipiEnd, const int& ipiRes);
	Array2D(const int& ipiStart, const int& ipiRes, const int& ipjStart, const int& ipjRes);
	Array2D(const Array2D<TT>& ipArray);

	Array2D(const Grid2D& ipGrid);

	void initialize(const int& iS, const int& iE, const int& iL, const int& jS, const int& jE, const int& jL, const int& ijL);
	void initialValues();

	const int index(const Vector2D<int>& ipVector) const
	{
		assert(ipVector[0] >= iStart && ipVector[0] <= iEnd);
		assert(ipVector[1] >= jStart && ipVector[1] <= jEnd);
		return (ipVector[0] - iStart) + (ipVector[1] - jStart)*iRes;
	}

	const int index(const int& i, const int& j) const
	{
		assert(i >= iStart && i <=iEnd);
		assert(j >= jStart && j <=jEnd);
		return (i - iStart) + (j - jStart)*iRes;
	}

	inline TT& operator [](const int& i) const
	{
		//assert(i >= iStart && i <=iEnd);
		return values[i];
	}
	inline TT& operator ()(const int& i)const
	{
		//assert(i >= iStart && i <=iEnd);
		return values[i];
	}

	//inline TT& operator [](const int& i, const int& j)
	//{
	//	assert(i >= iStart && i <=iEnd);
	//	assert(j >= jStart && j <=jEnd);
	//	return values[index(i,j)];
	//}

	inline TT& operator ()(const int& i, const int& j) const
	{
		assert(i >= iStart && i <=iEnd);
		assert(j >= jStart && j <=jEnd);

		return values[index(i, j)];
	}

	inline TT& operator ()(const Vector2D<int>& ipVector) const
	{
		assert(ipVector[0] >= iStart && ipVector[0] <= iEnd);
		assert(ipVector[1] >= jStart && ipVector[1] <= jEnd);

		return values[index(ipVector[0], ipVector[1])];
	}
	
	inline void operator =(const double& constant)
	{
		for (int i = 0; i < ijRes; i++)
		{
			values[i] = constant;
		}
	}

	inline void operator *=(const double& constant)
	{
		for (int i = 0; i < ijRes; i++)
		{
			values[i] *= constant;
		}
	}

	inline void operator +=(const double& constant)
	{
		for (int i = 0; i < ijRes; i++)
		{
			values[i] += constant;
		}
	}

	inline void operator -=(const double& constant)
	{
		for (int i = 0; i < ijRes; i++)
		{
			values[i] -= constant;
		}
	}

	inline void operator /=(const double& constant)
	{
		assert(constant != 0);
		
		for (int i = 0; i < ijRes; i++)
		{
			values[i] *= 1/constant;
		}
	}

	inline void operator =(const Array2D<TT>& ipArray)
	{
		if (values != nullptr)
		{
			delete[] values;
		}
		initialize(ipArray.iStart, ipArray.iEnd, ipArray.iRes, ipArray.jStart, ipArray.jEnd, ipArray.jRes, ipArray.ijRes);

		values = new TT[ijRes];

		for (int i = 0; i < ijRes; i++)
		{
			values[i] = ipArray.values[i];
		}
	}

	Array2D operator + (const Array2D<TT>& ipArray)
	{
		Array2D<TT> tempArray(ipArray.iStart, ipArray.iRes, ipArray.jStart, ipArray.jRes);
		
		for (int i = 0; i < tempArray.ijRes; i++)
		{
			tempArray.values[i] = value[i] + ipArray.values[i];
		}
		return tempArray;
	}

	Array2D operator - (const Array2D<TT>& ipArray)
	{
		Array2D<TT> tempArray(ipArray.iStart, ipArray.iRes, ipArray.jStart, ipArray.jRes);

		for (int i = 0; i < tempArray.ijRes; i++)
		{
			tempArray.values[i] = value[i] - ipArray.values[i];
		}
		return tempArray;
	}

	Array2D operator * (const Array2D<TT>& ipArray)
	{
		Array2D<TT> tempArray(ipArray.iStart, ipArray.iRes, ipArray.jStart, ipArray.jRes);

		for (int i = 0; i < tempArray.ijRes; i++)
		{
			tempArray.values[i] = value[i] * ipArray.values[i];
		}
		return tempArray;
	}

	Array2D operator / (const Array2D<TT>& ipArray)
	{
		Array2D<TT> tempArray(ipArray.iStart, ipArray.iRes, ipArray.jStart, ipArray.jRes);

		for (int i = 0; i < tempArray.ijRes; i++)
		{
			if (ipArray.values[i]!=0)
			{
				tempArray.values[i] = value[i] / ipArray.values[i];
			}
		}
		return tempArray;
	}

	Array2D operator + (const TT& constant)
	{
		Array2D<TT> tempArray = this;

		for (int i = 0; i < tempArray.ijRes; i++)
		{
			tempArray.values[i] = constant + ipArray.values[i];
		}
		return tempArray;
	}

	Array2D operator - (const TT& constant)
	{
		Array2D<TT> tempArray = this;

		for (int i = 0; i < tempArray.ijRes; i++)
		{
			tempArray.values[i] = ipArray.values[i] - constant;
		}
		return tempArray;
	}

	Array2D operator * (const TT& constant)
	{
		Array2D<TT> tempArray = this;

		for (int i = 0; i < tempArray.ijRes; i++)
		{
			tempArray.values[i] = ipArray.values[i] - constant;
		}
		return tempArray;
	}

	Array2D operator / (const TT& constant)
	{
		assert(constant != 0);
		Array2D<TT> tempArray = this;

		for (int i = 0; i < tempArray.ijRes; i++)
		{
			tempArray.values[i] = ipArray.values[i] - constant;
		}
		return tempArray;
	}

private:

};

template <class TT>
Array2D<TT>::Array2D()
{
	values = nullptr;
}
template <class TT>
Array2D<TT>::~Array2D()
{
	if (values != nullptr)
	{
		delete[] values;
	}
}

template<class TT>
inline Array2D<TT>::Array2D(const int & ipiRes)
{
	if (values != nullptr)
	{
		delete[] values;
	}
	initialize(0, ipiRes - 1, ipiRes, 0, 0, 0, ipiRes);

	assert(iRes > 0);

	values		= new TT[ijRes];

	initialValues();
}

template<class TT>
inline Array2D<TT>::Array2D(const int & ipiRes, const int& ipjRes)
{
	if (values!=nullptr)
	{
		delete[] values;
	}

	initialize(0, ipiRes - 1, ipiRes, 0, ipjRes - 1, ipjRes, ipiRes*ipjRes);

	assert(iRes > 0 && jRes > 0);

	values		= new TT[ijRes];

	initialValues();
}

template<class TT>
inline Array2D<TT>::Array2D(const int & ipiStart, const int & ipiEnd, const int & ipiRes)
{
	if (values != nullptr)
	{
		delete[] values;
	}

	initialize(ipiStart, ipiEnd, ipiRes, 0, 0, 0, ipiRes);

	assert(iRes > 0 && iEnd == iStart + iRes - 1);

	values = new TT[ijRes];

	initialValues();
}

template<class TT>
inline Array2D<TT>::Array2D(const int & ipiStart, const int & ipiRes, const int & ipjStart, const int & ipjRes)
{
	if (values != nullptr)
	{
		values = nullptr;
		//delete[] values;
	}
	initialize(ipiStart, ipiStart + ipiRes - 1, ipiRes, ipjStart, ipjStart + ipjRes - 1, ipjRes, ipiRes*ipjRes);

	assert(iRes > 0 && iEnd == iStart + iRes - 1);
	assert(jRes > 0 && jEnd == jStart + jRes - 1);
	values = new TT[ijRes];

	initialValues();
}

template<class TT>
inline Array2D<TT>::Array2D(const Array2D<TT>& ipArray)
{	
	if (values != nullptr)
	{
		delete[] values;
	}

	initialize(ipArray.iStart, ipArray.iEnd, ipArray.iRes, ipArray.jStart, ipArray.jEnd, ipArray.jRes, ipArray.ijRes);

	values = new TT[ijRes];

	for (int i = 0; i < ijRes; i++)
	{
		values[i] = ipArray.values[i];
	}
}

template<class TT>
inline Array2D<TT>::Array2D(const Grid2D & ipGrid)
{
	//if (values != nullptr)
	//{
	//	delete[] values;
	//}

	initialize(ipGrid.iStart, ipGrid.iEnd, ipGrid.iRes, ipGrid.jStart, ipGrid.jEnd, ipGrid.jRes, ipGrid.iRes*ipGrid.jRes);

	values = new TT[ijRes];

	for (int i = 0; i < ijRes; i++)
	{
		values[i] = 0;
	}
}

template<class TT>
inline void Array2D<TT>::initialize(const int & iS, const int & iE, const int & iL, const int & jS, const int & jE, const int & jL, const int& ijL)
{
	iStart = iS;
	iEnd = iE;
	iRes = iL;
	jStart = jS;
	jEnd = jE;
	jRes = jL;
	ijRes = ijL;
}

template<class TT>
inline void Array2D<TT>::initialValues()
{
	for (int i = 0; i < ijRes; i++)
	{
		values[i] = 0;
	}
}


