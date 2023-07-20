#pragma once
#include <vector>

#include "PointLight.h"
#include "SDF.h"

class SDF_List
{
public:
	//SDF_List
    SDF_List(color c) { defaultColor = c; }
    SDF_List(shared_ptr<sdf> object) { add(object); }

    void clear() { SDFObjects.clear(); }
    void add(shared_ptr<sdf> object) { SDFObjects.push_back(object); }

	std::vector<shared_ptr<sdf>> SDFObjects;
	std::vector<PointLight> lights;
    color defaultColor;

};
