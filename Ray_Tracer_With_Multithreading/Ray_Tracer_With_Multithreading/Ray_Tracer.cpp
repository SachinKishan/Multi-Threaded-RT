// Ray_Tracer.cpp : This file contains the 'main' function. Program execution begins and ends there.
//




///Todo
//////All original RT Features
///Camera
///Objects
///Scene manipulation
///
#include<stdlib.h>
#include "VectorMath.h"
#include "Matrix.h"
#include <iostream>

#include "camera.h"
#include "hittable.h"
#include "hittable_list.h"
#include "lodepng.h"
#include "material.h"
#include "ray.h"
#include "sphere.h"

#pragma region PNG
void encodeOneStep(const char* filename, std::vector<unsigned char>& image, unsigned width, unsigned height) {
	//Encode the image
	unsigned error = lodepng::encode(filename, image, width, height);

	//if there's an error, display it
	if (error) std::cout << "encoder error " << error << ": " << lodepng_error_text(error) << std::endl;
}
#pragma endregion

color ray_color(const ray& r, const hittable_list& world, int depth, bool scattering=true) {
    hit_record rec;
    std::vector<bool> shouldLight;
    if (depth <= 0)
        return color(0, 0, 0);

    if (world.hit(r, 0.001, infinity, rec)) {
        ray scattered;
        color attenuation;
        color finalCol;
        color shadowing;
        for (PointLight l : world.lights)
        {
            bool blocked_by_all = true;
            hit_record shadowRec;
            ray shadow_ray(rec.p + (rec.normal * 1e-4), l.lightDirection);

            if (world.hit(shadow_ray, 0.001, infinity, shadowRec))
            {
                shouldLight.push_back(false);
            }
            else { shouldLight.push_back(true); blocked_by_all = false; }
            shadowing = blocked_by_all ? Black : White;

        }
        
        if (rec.mat_ptr->scatter(r, rec, attenuation, scattered,world.lights, shouldLight))
        {
            if (scattering)
                finalCol= attenuation * ray_color(scattered, world, depth - 1);
            else
                finalCol= attenuation;
        	return finalCol*shadowing;
        }
    	return color(0, 0, 0);
    }

    vec3 unit_direction = unit_vector(r.direction());
    auto t = 0.5 * (unit_direction.y() + 1.0);
    return color(0.5, 0.7, 1.0);
}


hittable_list random_scene() {
    hittable_list world;

    auto ground_material = make_shared<lambert>(color(0.5, 0.5, 0.5),0.5,0);
    world.add(make_shared<sphere>(point3(0, -1000, 0), 1000, ground_material));
    
    vec3 directionOfLight1(7, 9, 9);
    vec3 directionOfLight2(-9, 20, 9);
    color lightColor1(1, 1, 1);
    color lightColor2(1, 1, 1);
    float intensityOfLight1 = 0.1;
    float intensityOfLight2 = 0.1;

    PointLight l1(intensityOfLight1, directionOfLight1, lightColor1);
    PointLight l2(intensityOfLight2, directionOfLight2, lightColor2);
    world.lights.push_back(l1);
    //world.lights.push_back(l2);


	/*
    for (int a = -11; a < 11; a++) {
        for (int b = -11; b < 11; b++) {
            auto choose_mat = random_double();
            point3 center(a + 0.9 * random_double(), 0.2, b + 0.9 * random_double());

            if ((center - point3(4, 0.2, 0)).length() > 0.9) {
                shared_ptr<material> sphere_material;

                if (choose_mat < 0.8) {
                    // diffuse
                    auto albedo = color::random() * color::random();
                    sphere_material = make_shared<lambertian>(albedo);
                    world.add(make_shared<sphere>(center, 0.2, sphere_material));
                }
                else if (choose_mat < 0.95) {
                    // metal
                    auto albedo = color::random(0.5, 1);
                    auto fuzz = random_double(0, 0.5);
                    sphere_material = make_shared<metal>(albedo, fuzz);
                    world.add(make_shared<sphere>(center, 0.2, sphere_material));
                }
                else {
                    // glass
                    sphere_material = make_shared<dielectric>(1.5);
                    world.add(make_shared<sphere>(center, 0.2, sphere_material));
                }
            }
        }
    }*/

    //auto material1 = make_shared<dielectric>(1.5);
    //world.add(make_shared<sphere>(point3(0, 1, 0), 1, material1));

    //auto material2 = make_shared<lambertian>(color(0.4, 0.2, 0.1));
    //world.add(make_shared<sphere>(point3(-4, 1, 0), 1, material2));
    /*
    auto material4 = make_shared<lambert>(color(1, 0, 0),0.5,64);
    auto material5 = make_shared<lambert>(color(0, 0, 1),0.5,64);
    auto material6 = make_shared<lambert>(color(1, 1, 0),0.7,64);
    auto material3 = make_shared<metal>(color(0.7, 0.6, 0.5), 0.0);
    world.add(make_shared<sphere>(point3(4, 1, 0), 1, material4));
    world.add(make_shared<sphere>(point3(-4, 1, 0), 1, material5));
    world.add(make_shared<sphere>(point3(0, 1, 0), 1, material6));
    */

    double diff=0.35;
    double spec=32;
    auto john = make_shared<lambert>(Brown, diff, spec);
	auto mary = make_shared<lambert>(color(0.52, 0.8, 0.91), diff, spec);
	auto angel = make_shared<lambert>(color(0.7,0.7,0.7), diff, spec);
	auto jesus = make_shared<lambert>(White, diff, spec);
	auto king1 = make_shared<lambert>(Red, diff, spec);
	auto king2 = make_shared<lambert>(Green, diff , spec);
	auto king3 = make_shared<lambert>(Blue, diff, spec);

    //world.add(make_shared<sphere>(point3(-0.1, 2, 4.7), 2, jesus));
    world.add(make_shared<sphere>(point3(-4, 4, 0), 4, john));
    //world.add(make_shared<sphere>(point3(3.2, 3, 0.5), 3, mary));
    //world.add(make_shared<sphere>(point3(13, 4, 3.4 - 4 * 2.4), 4, king1));
    //world.add(make_shared<sphere>(point3(13+4*2, 4, 3.4-4*1.2), 4, king2));
    //world.add(make_shared<sphere>(point3(13+4*4, 4, 3.4), 4, king3));
    //world.add(make_shared<sphere>(point3(-16, 2.5, 1.3), 2.5, angel));



    return world;
}

int main()
{
	//generate image
	const char* filename = "out1.png";
	std::vector<unsigned char> image;
	const auto aspect_ratio = 1;
	const int image_width = 1024/2;
	const int image_height = static_cast<int>(image_width / aspect_ratio);
	int total = image_width * image_height;
    const int samples_per_pixel = 10;
    const int max_depth = 10;
	//image resizing
	image.resize(image_width * image_height * 4);

    // World

    auto world = random_scene();
    // Camera
    
    point3 lookfrom(6, 75, 150);
    //point3 lookfrom(1, 2, 18);
    point3 lookat(6, 1, 0);
    vec3 vup(0, 1, 0);
    auto dist_to_focus = 100.0;
    auto aperture = 0.1;

    camera cam(lookfrom, lookat, vup, 20, aspect_ratio, aperture, dist_to_focus);

	color defaultColor=Blue;
	for (unsigned y = image_height - 1; y >0; y--)
	{
		for (unsigned x = 0; x < image_width; x++) 
        {
            color col(0, 0, 0);
            for (int s = 0; s < samples_per_pixel; ++s) {
				auto u = (x + random_double()) / (image_width - 1);
				auto v = (y + random_double()) / (image_height - 1);

				ray ra = cam.get_ray_perspective(u, v);
				col += ray_color(ra, world, max_depth);

            }
            

#pragma region Image File Color Translation
			//conversion to unsigned char
            auto scale = 1.0 / samples_per_pixel;
            auto r = col.x();
            auto g = col.y();
            auto b = col.z();
            r = sqrt(scale * r);
            g = sqrt(scale * g);
            b = sqrt(scale * b);

			const unsigned char rf = static_cast<unsigned char>(std::min(1., r) * 255);
			const unsigned char gf = static_cast<unsigned char>(std::min(1., g) * 255);
			const unsigned char bf = static_cast<unsigned char>(std::min(1., b) * 255);

			//allot color
			int val = 4 * image_width * (image_height- y) + 4 * x;
			image[val + 0] = rf;
			image[val + 1] = gf;
			image[val + 2] = bf;
			image[val + 3] = 255;//alpha
#pragma endregion


        }
        if (y % 10 == 0)
        {
            system("cls");
            std::cout <<"Progress: " << int(100 * (double(image_height - y) / image_height)) << "%" << std::endl;
        }
	}
	encodeOneStep(filename, image, image_width, image_height);


	return 0;
}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
