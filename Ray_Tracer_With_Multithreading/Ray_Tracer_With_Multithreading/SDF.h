#pragma once
#include "VectorMath.h"

class SDF_Material;

class sdf
{
public:
	vec3 position;
	color colorOfSDF;
	shared_ptr<SDF_Material> material;
	double distanceFound1;
	double distanceFound2;
	virtual double DistanceFunction(vec3 point) { return 0; }
};

class SDFSphere :public sdf
{

public:
	double radius;
	
	SDFSphere(vec3 center,double r, shared_ptr<SDF_Material> mat)
	{
		material = mat;
		position = center;
		radius = r;
	}

	double DistanceFunction(vec3 point) override
	{
		distanceFound1= distance(point, position) - radius;
		return distanceFound1;
	}
};

class SDFBox :public sdf
{

public:
	vec3 boxDimensions;

	SDFBox(vec3 boxPos, vec3 b, shared_ptr<SDF_Material> mat)
	{
		material = mat;
		position = boxPos;
		boxDimensions = b;

	}

	double DistanceFunction(vec3 point) override
	{
		const vec3 q = absolute(point - position) - boxDimensions;
		const vec3 max_q = vec3(std::max(q.x(), 0.0), std::max(q.y(), 0.0), std::max(q.z(), 0.0));
		vec3 a = max_q;
		vec3 b= std::min(std::max(q.x(), std::max(q.y(), q.z())), 0.0);
		return (a + b).length();
	}
};



float sminCubic(double a, double b, double k)
{
	double h = std::max(k - abs(a - b), 0.0) / k;
	return std::min(a, b) - h * h * h * k * (1.0 / 6.0);
}


double SDFUnion(vec3 p, shared_ptr<sdf> &sdf1, shared_ptr<sdf> &sdf2)
{
	sdf1->distanceFound1= sdf1->DistanceFunction(p);
	sdf2->distanceFound1= sdf2->DistanceFunction(p);
	double d1 = sdf1->DistanceFunction(p);
	double d2 = sdf2->DistanceFunction(p);


	
	double dUnion = sminCubic(d1, d2,0.1);
	return dUnion;
}

double SDFIntersection(vec3 p, shared_ptr<sdf> &sdf1, shared_ptr<sdf> &sdf2)
{
	double d1 = sdf1->distanceFound1 = sdf1->DistanceFunction(p);
	double d2 = sdf2->distanceFound1 = sdf2->DistanceFunction(p);
	double dIntersection = std::max(d1, d2);
	return dIntersection;
}

double SDFSubtract(vec3 p, shared_ptr<sdf> &sdf1, shared_ptr<sdf> &sdf2)
{
	double d1 = sdf1->distanceFound1 = sdf1->DistanceFunction(p);
	double d2 = sdf2->distanceFound1 = sdf2->DistanceFunction(p);
	double dSubtract = std::max(d1, -d2);
	return dSubtract;
}




typedef double (*BooleanSDFOperation)(vec3, shared_ptr<sdf>&, shared_ptr<sdf>&);

class combinedSDF: public sdf
{
public:

	BooleanSDFOperation SDFOperation;
	shared_ptr<sdf> sdf1;
	shared_ptr<sdf> sdf2;

	combinedSDF(shared_ptr<sdf> _sdf1, shared_ptr<sdf> _sdf2, BooleanSDFOperation operations, shared_ptr<SDF_Material> mat)
		: sdf1(_sdf1), sdf2(_sdf2), SDFOperation(operations)
	{
		std::move(material) = mat;
	}

	double DistanceFunction(vec3 point) override
	{
		distanceFound1 = sdf1->DistanceFunction(point);
		distanceFound2 = sdf2->DistanceFunction(point);
		
		return SDFOperation(point, sdf1, sdf2);
	}

};