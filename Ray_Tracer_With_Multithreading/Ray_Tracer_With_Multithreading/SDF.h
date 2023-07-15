#pragma once
#include "VectorMath.h"


class sdf
{
public:
	vec3 position;
	color colorOfSDF;
	virtual double DistanceFunction(vec3 point) { return 0; }
};

class SDFSphere :public sdf
{

public:
	double radius;
	
	SDFSphere(vec3 center,double r, color c)
	{
		position = center;
		radius = r;
		colorOfSDF = c;
	}

	double DistanceFunction(vec3 point) override
	{
		return distance(point, position)-radius;
	}
};