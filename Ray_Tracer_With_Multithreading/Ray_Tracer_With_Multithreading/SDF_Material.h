#pragma once
#include <vector>

#include "PointLight.h"
#include "SDF.h"


class sdf;
class SDF_Material;
//base
/////
/// light "scatters"
/// ie light moves across implicit volumes, and the direction it moves depends on the material
///
///
struct SDF_hit_record {
    point3 p;
    vec3 normal;
    double t;
    bool front_face;

    shared_ptr<sdf> sdf_ptr;
};


class SDF_Material
{
public:
    virtual bool scatter(
        const ray& r_in, const SDF_hit_record& rec, color& attenuation, ray& scattered, std::vector<PointLight> &lights
    ) const
    {
        return false;
    }

    virtual color emitted(
        const ray& r_in, const SDF_hit_record& rec, const point3& p) const {
			return color(0, 0, 0);
    }
};


class SDF_material_emissive:public SDF_Material
{
public:

    color lightColor;
    double lightPower;
	SDF_material_emissive(double power,color a):lightPower(power),lightColor(a)
	{
	}

    bool scatter(const ray& r_in, const SDF_hit_record& rec, color& attenuation, ray& scattered, std::vector<PointLight> &lights) const override
	{
        
        return false;

	}
    color emitted(const ray& r_in, const SDF_hit_record& rec, const point3& p) const override
    {
        return lightColor*lightPower;

    }
};


//basic color


class SDF_material_basic_color :public SDF_Material
{
public:
    SDF_material_basic_color(color a, double diffuse, double specular, int sN) :albedo(a),diffuseCoeff(diffuse),specularCoeff(specular),specularN(sN)
    {}

    bool scatter(const ray& r_in, const SDF_hit_record& rec, color& attenuation, ray& scattered, std::vector<PointLight>& lights) const override
    {
       
        vec3 normal = rec.normal + vec3(0.001);

        for (PointLight l : lights)
        {
            vec3 lightVector = normalise(l.lightPosition);
            vec3 lightColor = l.lightColor;
            float lightPower = l.lightIntensity;

            double lightValue = diffuseCoeff * std::max(0.0, dot(lightVector, normal));
            vec3 h = unit_vector(r_in.direction() + lightVector);
            const float blinn = specularCoeff * std::pow((double)std::max(0.0, dot(normal, h)), (double)specularN);

            lightValue *= lightPower / dot(lightVector, lightVector);

            attenuation += blinn * lightColor + albedo * lightValue;
            
        }
        if(lights.empty())
        	attenuation = albedo;

        auto scatter_direction = normal+random_unit_vector();
        // Catch degenerate scatter direction
        if (scatter_direction.near_zero())
            scatter_direction = normal;

        scattered = ray(rec.p, scatter_direction);
        return true;
    }

    color albedo;
    double diffuseCoeff;
    double specularCoeff;
    int specularN;
};


double blend(double a, double b, double k)
{
    double h = std::max(k - abs(a - b), 0.0) / k;
    double m = h * h * h * 0.5;
    return (a < b) ? m : 1.0 - m;
}

double blendN(double a, double b, double k, double n)
{
    const double h = std::max(k - abs(a - b), 0.0) / k;
    const double m = pow(h, n) * 0.5;
    return (a < b) ? m : 1.0 - m;
}

class SDF_Material_Combined:public SDF_Material
{
public:
    shared_ptr<SDF_Material> material1;
    shared_ptr<SDF_Material> material2;
    SDF_Material_Combined(shared_ptr<SDF_Material>  m1, shared_ptr<SDF_Material>  m2) :material1(m1),material2(m2)
    {}

    bool scatter(const ray& r_in, const SDF_hit_record& rec, color& attenuation, ray& scattered, std::vector<PointLight>& lights) const override
    {
        color albedo1, albedo2;
        double m=0;
        ray s1, s2;
        material1->scatter(r_in, rec, albedo1, s1,lights);
        material2->scatter(r_in, rec, albedo2, s2,lights);

        //m = blend(rec.sdf_ptr->distanceFound1, rec.sdf_ptr->distanceFound2, 0.2);
        m = blendN(rec.sdf_ptr->distanceFound1, rec.sdf_ptr->distanceFound2, 0.5, 5);
        


        double d1 = abs(rec.sdf_ptr->distanceFound1);
        double d2 = abs(rec.sdf_ptr->distanceFound2);

        attenuation = m * albedo1 + (1.0 - m) * albedo2;

    	scattered = d1<d2 ? s1 : s2;
        return true;
    }
};
