#pragma once

#include "hittable.h";
#include "VectorMath.h"

class sphere:public hittable
{

public:
	point3 center;
	double radius;
    shared_ptr<material> mat_ptr;
public:
	sphere(){}
    sphere(point3 cen, double r, shared_ptr<material> m)
        : center(cen), radius(r), mat_ptr(m) {};

	virtual bool hit(const ray& r, double t_min, double t_max, hit_record& rec) const override;

};

bool sphere:: hit(const ray& r, double t_min, double t_max, hit_record& rec) const {
    vec3 oc = r.origin() - center;
    auto a = r.direction().length_squared();
    auto half_b = dot(oc, r.direction());
    auto c = oc.length_squared() - radius * radius;

    auto discriminant = half_b * half_b - a * c;
    if (discriminant < 0) return false;
    
    auto sqrtd = sqrt(discriminant);

    // Find the nearest root that lies in the acceptable range.

    auto root1 = (-half_b - sqrtd) / a;
    auto root2 = (-half_b + sqrtd) / a;
    if (root1 == root2)rec.hit_one_point = true;
    auto root = root1;
    if (root < t_min || t_max < root) {
        root = root2;
        if (root < t_min || t_max < root)
            return false;
    }

    rec.t = root;
    rec.p = r.at(rec.t);
    vec3 outward_normal = (rec.p - center) / radius;
    rec.set_face_normal(r, outward_normal);
    rec.mat_ptr = mat_ptr;

    return true;
}




