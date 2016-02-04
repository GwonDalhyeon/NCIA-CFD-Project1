#pragma once


//#ifndef EulerianFluidSolver_H
//#define EulerianFluidSolver_H
#include "AdvectionMethod2D.h"


class EulerianFluidSolver2D
{
public:
	Grid2D grid;

	LevelSet2D levelSet;

	Field2D<double> Pressure;

	Field2D<Vector2D<double>> velocity;

	double reynoldNum;
	double dt;

	EulerianFluidSolver2D();
	~EulerianFluidSolver2D();

private:

};



//#endif // !EulerianFluidSolver_H

EulerianFluidSolver2D::EulerianFluidSolver2D()
{
}

EulerianFluidSolver2D::~EulerianFluidSolver2D()
{
}
