#pragma once
#include <assert.h>
#include <iostream>

template <class TT>
class Array2D
{
public:
	union
	{
		struct { int iStart, jStart, iEnd, jEnd; };
		struct { int ijStart[2], ijEnd[2]; };
	};

	int iLength, jLength;
	int ijLength;
	TT* values;

	Array2D();
	~Array2D();

	Array2D(const int& inputiLength);
	Array2D(const int& inputiLength, const int& inputjLength);
	Array2D(const int& inputiStart, const int& inputiEnd, const int& inputiLength);
	Array2D(const int& inputiStart, const int& inputiLength, const int& inputjStart, const int& inputjLength);
	Array2D(const Array2D<TT>& inputArray);
	void initialValues();

	inline int index(const int& i, const int& j)
	{
		assert(i >= iStart && i <=iEnd);
		assert(j >= jStart && j <=jEnd);
		return i + j*iLength;
	}

	inline TT& operator [](const int& i) const
	{
		assert(i >= iStart && i <=iEnd);
		return values[i];
	}
	inline TT& operator ()(const int& i)const
	{
		iStart;
		assert(i >= iStart && i <=iEnd);
		return values[i];
	}

	//inline TT& operator [](const int& i, const int& j)
	//{
	//	assert(i >= iStart && i <=iEnd);
	//	assert(j >= jStart && j <=jEnd);
	//	return values[index(i,j)];
	//}

	inline TT& operator ()(const int& i, const int& j) 
	{
		assert(i >= iStart && i <=iEnd);
		assert(j >= jStart && j <=jEnd);
		return values[index(i, j)];
	}
	
	inline void operator =(const double& constant)
	{
		for (int i = 0; i < ijLength; i++)
		{
			values[i] = constant;
		}
	}

	inline void operator *=(const double& constant)
	{
		for (int i = 0; i < ijLength; i++)
		{
			values[i] *= constant;
		}
	}

	inline void operator +=(const double& constant)
	{
		for (int i = 0; i < ijLength; i++)
		{
			values[i] += constant;
		}
	}

	inline void operator -=(const double& constant)
	{
		for (int i = 0; i < ijLength; i++)
		{
			values[i] -= constant;
		}
	}

	inline void operator /=(const double& constant)
	{
		assert(constant != 0);
		
		for (int i = 0; i < ijLength; i++)
		{
			values[i] *= 1/constant;
		}
	}

	inline void operator =(const Array2D<TT>& inputArray)
	{
		if (values != nullptr)
		{
			delete[] values;
		}
		
		iStart = inputArray.iStart;
		iEnd = inputArray.iEnd;
		iLength = inputArray.iLength;
		jStart = inputArray.jStart;
		jEnd = inputArray.jEnd;
		jLength = inputArray.jLength;
		ijLength = inputArray.ijLength;

		values = new TT[ijLength];

		for (int i = 0; i < ijLength; i++)
		{
			values[i] = inputArray.values[i];
		}
	}

	Array2D operator + (const Array2D<TT>& inputArray)
	{
		Array2D<TT> tempArray(inputArray.iStart, inputArray.iLength, inputArray.jStart, inputArray.jLength);
		
		for (int i = 0; i < tempArray.ijLength; i++)
		{
			tempArray.values[i] = value[i] + inputArray.values[i];
		}
		return tempArray;
	}

	Array2D operator - (const Array2D<TT>& inputArray)
	{
		Array2D<TT> tempArray(inputArray.iStart, inputArray.iLength, inputArray.jStart, inputArray.jLength);

		for (int i = 0; i < tempArray.ijLength; i++)
		{
			tempArray.values[i] = value[i] - inputArray.values[i];
		}
		return tempArray;
	}

	Array2D operator * (const Array2D<TT>& inputArray)
	{
		Array2D<TT> tempArray(inputArray.iStart, inputArray.iLength, inputArray.jStart, inputArray.jLength);

		for (int i = 0; i < tempArray.ijLength; i++)
		{
			tempArray.values[i] = value[i] * inputArray.values[i];
		}
		return tempArray;
	}

	Array2D operator / (const Array2D<TT>& inputArray)
	{
		Array2D<TT> tempArray(inputArray.iStart, inputArray.iLength, inputArray.jStart, inputArray.jLength);

		for (int i = 0; i < tempArray.ijLength; i++)
		{
			if (inputArray.values[i]!=0)
			{
				tempArray.values[i] = value[i] / inputArray.values[i];
			}
		}
		return tempArray;
	}

	Array2D operator + (const TT& constant)
	{
		Array2D<TT> tempArray = this;

		for (int i = 0; i < tempArray.ijLength; i++)
		{
			tempArray.values[i] = constant + inputArray.values[i];
		}
		return tempArray;
	}

	Array2D operator - (const TT& constant)
	{
		Array2D<TT> tempArray = this;

		for (int i = 0; i < tempArray.ijLength; i++)
		{
			tempArray.values[i] = inputArray.values[i] - constant;
		}
		return tempArray;
	}

	Array2D operator * (const TT& constant)
	{
		Array2D<TT> tempArray = this;

		for (int i = 0; i < tempArray.ijLength; i++)
		{
			tempArray.values[i] = inputArray.values[i] - constant;
		}
		return tempArray;
	}

	Array2D operator / (const TT& constant)
	{
		assert(constant != 0);
		Array2D<TT> tempArray = this;

		for (int i = 0; i < tempArray.ijLength; i++)
		{
			tempArray.values[i] = inputArray.values[i] - constant;
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
inline Array2D<TT>::Array2D(const int & inputiLength)
{
	if (values != nullptr)
	{
		delete[] values;
	}
	iStart		= 0;
	iEnd		= inputiLength - 1;
	jStart		= 0;
	jEnd		= 0;
	iLength		= inputiLength;
	jLength		= 0;
	ijLength	= iLength;

	assert(iLength > 0);

	values		= new TT[ijLength];

	initialValues();
}

template<class TT>
inline Array2D<TT>::Array2D(const int & inputiLength, const int& inputjLength)
{
	if (values!=nullptr)
	{
		delete[] values;
	}
	iStart		= 0;
	iEnd		= inputiLength-1;
	jStart		= 0;
	jEnd		= inputjLength-1;
	iLength		= inputiLength;
	jLength		= inputjLength;
	ijLength	= iLength*jLength;

	assert(iLength > 0 && jLength > 0);

	values		= new TT[ijLength];

	initialValues();
}

template<class TT>
inline Array2D<TT>::Array2D(const int & inputiStart, const int & inputiEnd, const int & inputiLength)
{
	if (values != nullptr)
	{
		delete[] values;
	}
	iStart		= inputiStart;
	iEnd		= inputiEnd;
	iLength		= inputiLength;
	jStart		= 0;
	jEnd		= 0;
	jLength		= 0;
	ijLength	= iLength;

	assert(iLength > 0 && iEnd == iStart + iLength - 1);

	values = new TT[ijLength];

	initialValues();
}

template<class TT>
inline Array2D<TT>::Array2D(const int & inputiStart, const int & inputiLength, const int & inputjStart, const int & inputjLength)
{
	if (values != nullptr)
	{
		delete[] values;
	}
	iStart		= inputiStart;
	iEnd		= inputiStart + inputiLength - 1;
	iLength		= inputiLength;
	jStart		= inputjStart;
	jEnd		= inputjStart + inputjLength - 1;
	jLength		= inputjLength;
	ijLength	= iLength*jLength;

	assert(iLength > 0 && iEnd == iStart + iLength - 1);
	assert(jLength > 0 && jEnd == jStart + jLength - 1);

	values = new TT[ijLength];

	initialValues();
}

template<class TT>
inline Array2D<TT>::Array2D(const Array2D<TT>& inputArray)
{	
	if (values != nullptr)
	{
		delete[] values;
	}

	iStart = inputArray.iStart;
	iEnd = inputArray.iEnd;
	iLength = inputArray.iLength;
	jStart = inputArray.jStart;
	jEnd = inputArray.jEnd;
	jLength = inputArray.jLength;
	ijLength = inputArray.ijLength;

	values = new TT[ijLength];

	for (int i = 0; i < ijLength; i++)
	{
		values[i] = inputArray.values[i];
	}
}

template<class TT>
inline void Array2D<TT>::initialValues()
{
	for (int i = 0; i < ijLength; i++)
	{
		values[i] = 0;
	}
}


