#pragma once

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

EulerianFluidSolver2D::EulerianFluidSolver2D()
{
}

EulerianFluidSolver2D::~EulerianFluidSolver2D()
{
}