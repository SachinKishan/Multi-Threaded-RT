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
/*
 *union = min(a, b)
subtract = min(a, -b)
intersect = max(a, b)
*/

double SDFUnion(vec3 p, shared_ptr<sdf> sdf1, shared_ptr<sdf> sdf2)
{
	double d1 = sdf1->DistanceFunction(p);
	double d2 = sdf2->DistanceFunction(p);
	double dUnion = std::min(d1, d2);
	return dUnion;
}

double SDFIntersection(vec3 p, shared_ptr<sdf> sdf1, shared_ptr<sdf> sdf2)
{
	double d1 = sdf1->DistanceFunction(p);
	double d2 = sdf2->DistanceFunction(p);
	double dIntersection = std::max(d1, d2);
	return dIntersection;
}

double SDFSubtract(vec3 p, shared_ptr<sdf> sdf1, shared_ptr<sdf> sdf2)
{
	double d1 = sdf1->DistanceFunction(p);
	double d2 = sdf2->DistanceFunction(p);
	double dSubtract = std::max(d1, -d2);
	return dSubtract;
}




typedef double (*BooleanSDFOperation)(vec3, shared_ptr<sdf>, shared_ptr<sdf>);

class combinedSDF: public sdf
{
public:

	BooleanSDFOperation SDFOperation;
	shared_ptr<sdf> sdf1;
	shared_ptr<sdf> sdf2;

	combinedSDF(shared_ptr<sdf> _sdf1, shared_ptr<sdf> _sdf2, BooleanSDFOperation operations)
		: sdf1(_sdf1), sdf2(_sdf2), SDFOperation(operations)
	{
		colorOfSDF = (sdf1->colorOfSDF + sdf2->colorOfSDF) / 2;
	}

	double DistanceFunction(vec3 point) override
	{
		return SDFOperation(point, sdf1, sdf2);
	}

};