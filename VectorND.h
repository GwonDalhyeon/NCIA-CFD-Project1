#pragma once
#include <assert.h>
#include <iostream>

template <class TT>
class VectorND
{
public:

		TT* values;
		int iEnd;
		int iLength;

	VectorND();
	~VectorND();
	VectorND(const int& ipLength);
	VectorND(const VectorND<TT>& ipVector);
	VectorND(const TT* ipValues);

	inline TT operator [](const int& i) const
	{
		assert(i >= 0 || i < iLength);
		return values[i];
	}

	inline TT operator ()(const int& i) const
	{
		assert(i >= 0 || i < iLength);
		return values[i];
	}

	inline void operator = (const TT& constant)
	{
		for (int i = 0; i < iLength; i++)
		{
			values[i] = constant;
		}
	}

	inline void operator = (const VectorND<TT>& ipVector)
	{
		//assert(iLength == ipVector.iLength);
		iLength = ipVector.iLength;
		
		if (values != nullptr)
		{
			delete[] values;
		}
		
		values = new TT[iLength];

		for (int i = 0; i < iLength; i++)
		{
			values[i] =ipVector.values[i];
		}
	}

	inline void operator += (const VectorND<TT>& ipVector)
	{
		assert(iLength == ipVector.iLength);
		for (int i = 0; i < iLength; i++)
		{
			values[i] += ipVector.values[i];
		}
	}

	inline void operator -= (const VectorND<TT>& ipVector)
	{
		assert(iLength == ipVector.iLength);
		for (int i = 0; i < iLength; i++)
		{
			values[i] -= ipVector.values[i];
		}
	}

	inline void operator /= (const VectorND<TT>& ipVector)
	{
		assert(iLength == ipVector.iLength);
		for (int i = 0; i < iLength; i++)
		{
			values[i] /= ipVector.values[i];
		}
	}

	inline void operator *= (const VectorND<TT>& ipVector)
	{
		assert(iLength == ipVector.iLength);
		for (int i = 0; i < iLength; i++)
		{
			values[i] *= ipVector.values[i];
		}
	}

	inline void operator += (const TT& constant)
	{
		for (int i = 0; i < iLength; i++)
		{
			values[i] += constant;
		}
	}

	inline void operator -= (const TT& constant)
	{
		for (int i = 0; i < iLength; i++)
		{
			values[i] -= constant;
		}
	}

	inline void operator *= (const TT& constant)
	{
		for (int i = 0; i < iLength; i++)
		{
			values[i] *= constant;
		}
	}

	inline void operator /= (const TT& constant)
	{
		assert(constant>0 || constant<0)
		for (int i = 0; i < iLength; i++)
		{
			values[i] /= constant;
		}
	}

	VectorND operator + (const VectorND<TT>& ipVector)
	{
		assert(iLength == ipVector.iLength);

		VectorND<TT> returnVT(ipVector.iLength);

		for (int i = 0; i < iLength; i++)
		{
			returnVT.values[i] = values[i] + ipVector.values[i];
		}
		return returnVT;
	}

	VectorND operator - (const VectorND<TT>& ipVector)
	{
		assert(iLength == ipVector.iLength);

		VectorND<TT> returnVT(ipVector.iLength);

		for (int i = 0; i < iLength; i++)
		{
			returnVT.values[i] = values[i] - ipVector.values[i];
		}
		return returnVT;
	}

	VectorND operator * (const VectorND<TT>& ipVector)
	{
		assert(iLength == ipVector.iLength);

		VectorND<TT> returnVT(ipVector.iLength);

		for (int i = 0; i < iLength; i++)
		{
			returnVT.values[i] = values[i] * ipVector.values[i];
		}
		return returnVT;
	}

	VectorND operator / (const VectorND<TT>& ipVector)
	{
		assert(iLength == ipVector.iLength);

		VectorND<TT> returnVT(ipVector.iLength);

		for (int i = 0; i < iLength; i++)
		{
			assert(ipVector.values[i]>0 && ipVector.values[i]<0);
			returnVT.values[i] = values[i] / ipVector.values[i];
		}
		return returnVT;
	}

	VectorND operator + (const TT& constant)
	{
		VectorND<TT> returnVT(ipVector.iLength);

		for (int i = 0; i < iLength; i++)
		{
			returnVT.values[i] = ipVector.values[i] + constant;
		}
		return returnVT;
	}

	VectorND operator - (const TT& constant)
	{
		VectorND<TT> returnVT(ipVector.iLength);

		for (int i = 0; i < iLength; i++)
		{
			returnVT.values[i] = ipVector.values[i] - constant;
		}
		return returnVT;
	}

	VectorND operator * (const TT& constant)
	{
		VectorND<TT> returnVT(ipVector.iLength);

		for (int i = 0; i < iLength; i++)
		{
			returnVT.values[i] = ipVector.values[i] * constant;
		}
		return returnVT;
	}

	VectorND operator / (const TT& constant)
	{
		assert(constant != 0);
		VectorND<TT> returnVT(ipVector.iLength);

		for (int i = 0; i < iLength; i++)
		{
			returnVT.values[i] = ipVector.values[i] / constant;
		}
		return returnVT;
	}

	TT magnitude();

	void normalize();

private:

};

template <class TT>
VectorND<TT>::VectorND()
{
	values = nullptr;
}

template <class TT>
VectorND<TT>::~VectorND()
{
	if (iLength > 0 && values!=nullptr)
	{
		cout << iLength << endl;
		cout << values[0] << endl;
		delete[] values;
	}
	
}

template<class TT>
inline VectorND<TT>::VectorND(const int & ipLength)
{
	values = new TT[ipLength];
	iLength = ipLength;
	iEnd = ipLength - 1;
}

template<class TT>
inline VectorND<TT>::VectorND(const VectorND<TT>& ipVector)
{
	if (values != nullptr)
	{
		delete[] values;
	}

	assert(iLength > 0);

	values = new TT[iLength];

	for (int i = 0; i < iLength; i++)
	{
		values[i] = ipVector[i];
	}
}


template<class TT>
inline TT VectorND<TT>::magnitude()
{
	double mag2 = 0.0;

	for (int i = 0; i < iLength; i++)
	{
		mag2 = mag2 + values[i] * values[i];
	}
	return TT(sqrt(mag2));
}

template<class TT>
inline void VectorND<TT>::normalize()
{
	*this /= magnitude();
}



template<class TT>
inline VectorND<TT> operator + (const TT& constant, const VectorND<TT>& ipVector)
{
	VectorND<TT> returnVT(ipVector.iLength);

	for (int i = 0; i < iLength; i++)
	{
		returnVT.values[i] = constant + ipVector.values[i];
	}
	return returnVT;
}

template<class TT>
inline VectorND<TT> operator - (const TT& constant, const VectorND<TT>& ipVector)
{
	VectorND<TT> returnVT(ipVector.iLength);

	for (int i = 0; i < iLength; i++)
	{
		returnVT.values[i] = constant - ipVector.values[i];
	}
	return returnVT;
}

template<class TT>
inline VectorND<TT> operator *(const TT& constant, const VectorND<TT>& ipVector)
{
	VectorND<TT> returnVT(ipVector.iLength);

	for (int i = 0; i < iLength; i++)
	{
		returnVT.values[i] = constant * ipVector.values[i];
	}
	return returnVT;
}

template<class TT>
inline VectorND<TT> operator / (const TT& constant, const VectorND<TT>& ipVector)
{
	VectorND<TT> returnVT(ipVector.iLength);

	for (int i = 0; i < iLength; i++)
	{
		assert(ipVector.values[i] != 0);
		returnVT.values[i] = constant / ipVector.values[i];
	}
	return returnVT;
}

template<class TT>
inline TT dotProduct(const VectorND<TT>& ipVector1, const VectorND<TT>& ipVector2)
{
	assert(ipVector1.iLength == ipVector2.iLength);

	double dotPro = 0;

	for (int i = 0; i < ipVector1.iLength; i++)
	{
		dotPro = dotPro + ipVector1[i] * ipVector2[i];
	}

	return dotPro;
}

template<class TT>
inline std::ostream& operator << (std::ostream& output, const VectorND<TT>& ipVector)
{
	for (int i = 0; i < ipVector.iLength; i++)
	{
		output << i << " " << ipVector[i] << endl;;
	}
	return output;
}


template<class TT>
VectorND<TT> normalize(const VectorND<TT>& ipVector)
{
	VectorND<TT> returnVector=ipVector;
	returnVector.normalize();
	return returnVector;
}