#pragma once
#include <assert.h>
#include <iostream>

template <class TT>
class Vector2D
{
public:
	union
	{
		struct { TT x, y; };
		struct { TT i, j; };
		TT values[2];
	};

	Vector2D();
	~Vector2D();
	Vector2D(const TT& inputX, const TT& inputY);
	Vector2D(const Vector2D<TT>& inputVector);
	Vector2D(const TT inputValues[2]);

	inline TT operator [](const int& i)
	{
		assert(i == 0 || i == 1);
		return values[i];
	}
	
	inline TT operator ()(const int& i)
	{
		assert(i == 0 || i == 1);
		return values[i];
	}

	inline void operator = (const Vector2D<TT>& inputVector)
	{
		x = inputVector.x;
		y = inputVector.y;
	}

	inline void operator += (const Vector2D<TT>& inputVector)
	{
		x += inputVector.x;
		y += inputVector.y;
	}

	inline void operator -= (const Vector2D<TT>& inputVector)
	{
		x -= inputVector.x;
		y -= inputVector.y;
	}

	inline void operator /= (const Vector2D<TT>& inputVector)
	{
		x /= inputVector.x;
		y /= inputVector.y;
	}

	inline void operator *= (const Vector2D<TT>& inputVector)
	{
		x *= inputVector.x;
		y *= inputVector.y;
	}

	inline void operator += (const TT& value)
	{
		x += value;
		y += value;
	}

	inline void operator -= (const TT& value)
	{
		x -= value;
		y -= value;
	}

	inline void operator *= (const TT& value)
	{
		x *= value;
		y *= value;
	}

	inline void operator /= (const TT& value)
	{
		x /= value;
		y /= value;
	}

	Vector2D operator + (const Vector2D<TT>& inputVector)
	{
		return Vector2D(x + inputVector.x, y + inputVector.y);
	}

	Vector2D operator - (const Vector2D<TT>& inputVector)
	{
		return Vector2D(x - inputVector.x, y - inputVector.y);
	}
	
	Vector2D operator * (const Vector2D<TT>& inputVector)
	{
		return Vector2D(x *inputVector.x, y * inputVector.y);
	}
	
	Vector2D operator / (const Vector2D<TT>& inputVector)
	{
		assert(inputVector.x != 0 && inputVector.y != 0);
		return Vector2D(x / inputVector.x, y / inputVector.y);
	}
	
	Vector2D operator + (const TT& value)
	{
		return Vector2D(x + value, y + value);
	}
	
	Vector2D operator - (const TT& value)
	{
		return Vector2D(x  value, y - value);
	}
	
	Vector2D operator * (const TT& value)
	{
		return Vector2D(x * value, y * value);
	}
	
	Vector2D operator / (const TT& value)
	{
		assert(value != 0);
		return Vector2D(x / value, y / value);
	}

	TT magnitude();

	void normalize();

private:

};

template <class TT>
Vector2D<TT>::Vector2D()
{
}

template <class TT>
Vector2D<TT>::~Vector2D()
{
}

template<class TT>
inline Vector2D<TT>::Vector2D(const TT & inputX, const TT & inputY)
{
	x = inputX;
	y = inputY;
}

template<class TT>
inline Vector2D<TT>::Vector2D(const Vector2D<TT>& inputVector)
{
	x = inputVector.x;
	y = inputVector.y;
}

template<class TT>
inline Vector2D<TT>::Vector2D(const TT inputValues[2])
{
	x = inputValues[0];
	y = inputValues[1];
}

template<class TT>
inline TT Vector2D<TT>::magnitude()
{
	return TT(sqrt(x*x + y*y));
}

template<class TT>
inline void Vector2D<TT>::normalize()
{
	*this /= magnitude();
}



template<class TT>
inline Vector2D<TT> operator + (const TT& value, const Vector2D<TT>& inputVector)
{
	return Vector2D<TT>(value + inputVector.x, value + inputVector.y);
}

template<class TT>
inline Vector2D<TT> operator - (const TT& value, const Vector2D<TT>& inputVector)
{
	return Vector2D<TT>(value - inputVector.x, value - inputVector.y);
}

template<class TT>
inline Vector2D<TT> operator * (const TT& value, const Vector2D<TT>& inputVector)
{
	return Vector2D<TT>(value * inputVector.x, value * inputVector.y);
}

template<class TT>
inline Vector2D<TT> operator / (const TT& value, const Vector2D<TT>& inputVector)
{
	assert(inputVector.x != 0 && inputVector.y != 0);
	return Vector2D<TT>(value / inputVector.x, value / inputVector.y);
}

template<class TT>
inline TT dotProduct(const Vector2D<TT>& inputVector1, const Vector2D<TT>& inputVector2)
{
	return inputVector1.x*inputVector2.x + inputVector1.y*inputVector2.y;
}

template<class TT>
inline std::ostream& operator << (std::ostream& output, const Vector2D<TT>& v)
{
	return output << v.x << " " << v.y;
}


template<class TT>
Vector2D<TT> normalize(const Vector2D<TT>& inputVector)
{
	Vector2D<TT> returnVector(inputVector.x, inputVector.y);
	returnVector.normalize();
	return returnVector;
}