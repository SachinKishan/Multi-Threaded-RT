#pragma once
#include "VectorMath.h"


class PointLight
{

	public:
		PointLight(float i, vec3 dir, color col):lightIntensity(i),lightPosition(dir),
	lightColor(col){}

public:
		float lightIntensity;//range from 0 to 1
		vec3 lightPosition;//position of source
		vec3 lightDirection;//position of source
		color lightColor;//color of light
};
