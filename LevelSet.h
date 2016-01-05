#pragma once
#include "GridInfomation.h"

class LevelSet
{
public:
	LevelSet();
	~LevelSet();

	GridInfo grid;
	double* phi;

private:

};

LevelSet::LevelSet()
{
}

LevelSet::~LevelSet()
{
}