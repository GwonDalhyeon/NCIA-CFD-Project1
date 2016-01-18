#pragma once

#include "LevelSet2D.h"

template <class TT>
class AdvectionMethod2D
{
public:
	AdvectionMethod2D();
	~AdvectionMethod2D();

	static void advection()
	{
		cout << "static" << endl;
	}
private:

};

template <class TT>
AdvectionMethod2D<TT>::AdvectionMethod2D()
{
}

template <class TT>
AdvectionMethod2D<TT>::~AdvectionMethod2D()
{
}


