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


#include <omp.h>
#include <thread>

#include "SDF.h"


std::vector<shared_ptr<sdf>> sdfWorld;


void explainOMP(std::string a)
{
    std::cout << "Thread number: " << omp_get_thread_num()<<" "<<a<<std::endl;
}


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
            for (PointLight l : world.lights)
            {
            	hit_record shadowRec;
                ray shadow_ray(rec.p + (rec.normal * 1e-4), l.lightDirection);

                if (world.hit(shadow_ray, 0.001, infinity, shadowRec))
                {
                    shouldLight.push_back(false);
                }
                else { shouldLight.push_back(true); }
            }
    
            if (rec.mat_ptr->scatter(r, rec, attenuation, scattered, world.lights, shouldLight))
            {
                if (scattering)
                    finalCol = attenuation * ray_color(scattered, world, depth - 1);
                else
                    finalCol = attenuation;
                return finalCol;
            }
            return color(0, 0, 0);
        }
    //vec3 unit_direction = unit_vector(r.direction());
    //auto t = 0.5 * (unit_direction.y() + 1.0);
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

    double diff=0.35;
    double spec=128;
    auto brown = make_shared<lambert>(Brown, diff, spec);

    world.add(make_shared<sphere>(point3(-4, 4, 0), 4, brown));



    return world;
}

auto world = random_scene();
void renderTile(int tno, int totaltiles, int image_width, int image_height, camera cam, std::vector<unsigned char> &image, int samples_per_pixel, int max_depth, int th, int tw)
{
    int t = tno;//tile number
    int ytno = t / totaltiles + 1;
    int xtno = t % totaltiles + 1;
    for (unsigned y = ytno * th ; y > (ytno - 1) * th; y--)
    {
        for (unsigned x = (xtno - 1) * tw; x < tw * xtno; x++)
        {
            color col(0, 0, 0);
            for (int s = 0; s < samples_per_pixel; ++s) {
                auto u = (x + random_double()) / (image_width - 1);
                auto v = (y + random_double()) / (image_width - 1);
                ray ra = cam.get_ray_perspective(u, v);
                col += ray_color(ra, world, max_depth);
            }
            //std::cout << std::endl<<col;
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
            int val = 4 * image_width * (image_height - y) + 4 * x;
            //allot color
            image[val + 0] = rf;
            image[val + 1] = gf;
            image[val + 2] = bf;
            image[val + 3] = 255;//alpha
#pragma endregion
        }
    }
}


double findSmallestSDFDistance(vec3 p, std::vector<shared_ptr<sdf>> sdfWorld, shared_ptr<sdf> &hitSDF)
{
  

    double minDistance = infinity;
    hitSDF = nullptr;  // Set the initial value of hitSDF to nullptr

    for (const auto& obj : sdfWorld)
    {
        double d = obj->DistanceFunction(p);
        if (d < minDistance)
        {
            minDistance = d;
            hitSDF = obj;
        }
    }

    return minDistance;
}



vec3 calcsdfNormal(vec3 p, std::vector<shared_ptr<sdf>> sdfWorld) {
    vec3 e = vec3(1.0, -1.0, 0) * 0.0005;

    vec3 a = vec3(e.x(), e.y(), e.y());
    vec3 b = vec3(e.y(), e.y(), e.x());
    vec3 c = vec3(e.y(), e.x(), e.y());
    vec3 d = vec3(e.x(), e.x(), e.x());
    shared_ptr<sdf> dummy;
    return normalise(
        a * findSmallestSDFDistance(p + a, sdfWorld,dummy) +
        b * findSmallestSDFDistance(p + b, sdfWorld,dummy) +
        c * findSmallestSDFDistance(p + c, sdfWorld,dummy) +
        d * findSmallestSDFDistance(p + d, sdfWorld,dummy));
}

void renderTileRayMarch(int tno, int totaltiles, int image_width, int image_height, camera cam, std::vector<unsigned char>& image, int samples_per_pixel, int max_depth, int th, int tw)
{
    int t = tno;//tile number
    int ytno = t / totaltiles + 1;
    int xtno = t % totaltiles + 1;
    for (unsigned y = ytno * th; y > (ytno - 1) * th; y--)
    {
        for (unsigned x = (xtno - 1) * tw; x < tw * xtno; x++)
        {
            color col(0, 0, 0);
            color raycol = Black;
            color defaultColor(0.3, 0.5, 0.7);

            for (int s = 0; s < samples_per_pixel; ++s) 
            {
                auto u = (x + random_double()) / (image_width - 1);
                auto v = (y + random_double()) / (image_width - 1);
                ray ra = cam.get_ray_perspective(u, v);
                float t = 1;
                float h = 0;
                shared_ptr<sdf> temp;
                for (int i = 0; i < 256; i++)
                {
                    vec3 pos = ra.orig + normalise(ra.dir) * t;
                	h = findSmallestSDFDistance(pos, sdfWorld, temp);
                    t += h;
                    if (h < 0.0001 || h>10) { /*std::cout << "found";*/ raycol = temp->colorOfSDF; break; }
                }
                if (h < 0.0001) 
                {
                    vec3 pos = ra.orig + normalise(ra.dir) * t;
                    vec3 normal = calcsdfNormal(pos, sdfWorld);
                    vec3 light = vec3(1, 2, 3);
                    float dif = clamp(dot(normal, normalise(light - pos)), 0., 1.);
                    dif *= 5. / dot(light - pos, light - pos);
                    col += dif*temp->colorOfSDF;
                }
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
            int val = 4 * image_width * (image_height - y) + 4 * x;
            //allot color
            image[val + 0] = rf;
            image[val + 1] = gf;
            image[val + 2] = bf;
            image[val + 3] = 255;//alpha
#pragma endregion
        }
    }
}





int main()
{
    //openmp set up
    double startTime = omp_get_wtime();
    omp_set_nested(1);
    omp_set_num_threads(8);

	//generate image
	
    const char* filename = "out1.png";
    std::vector<unsigned char> image;
    const auto aspect_ratio = 1;
    const int image_width = 256;
    const int image_height = static_cast<int>(image_width / aspect_ratio);
    int total = image_width * image_height;
    const int samples_per_pixel = 5;
    const int max_depth = 10;
    //image resizing
    image.resize(image_width * image_height * 4);

  

    // World
    shared_ptr<sdf> sphere2 = make_shared<SDFSphere>(vec3(1, 1, 0), 0.4, Green);
    shared_ptr<sdf> sphere = make_shared<SDFSphere>(vec3(0, 0, 0), 1, Blue);
    shared_ptr<sdf> sphere3 = make_shared<SDFSphere>(vec3(-1, 1, 0), 1, Red);

    //sdfWorld.push_back(sphere2);
    sdfWorld.push_back(sphere);
    //-sdfWorld.push_back(sphere3);


    //shared_ptr<sdf> cSDF = make_shared<combinedSDF>(sphere3, sphere, &SDFSubtract);


    //sdfWorld.push_back(cSDF);


    
    // Camera
    
    point3 lookfrom(0, 0, 10);
    point3 lookat(0, 0, 0);
    vec3 vup(0, 1, 0);
    auto dist_to_focus = 100.0;
    auto aperture = 0.1;

    camera cam(lookfrom, lookat, vup, 20, aspect_ratio, aperture, dist_to_focus);

	color defaultColor=Blue;

    
    
    int numberOfTiles = 128;
    int th = image_height / numberOfTiles;
    int tw = image_width / numberOfTiles;
    startTime = omp_get_wtime();

	#pragma omp parallel
    {
		#pragma omp single nowait
        {
            for (int i = 0; i < numberOfTiles*numberOfTiles; i++) 
            { 
				#pragma omp task
                {
                    renderTileRayMarch(i, numberOfTiles, image_width, image_height, cam, image, samples_per_pixel, max_depth, th,tw);
                }
				
            }
        }
	#pragma omp taskwait
    }
    std::cout <<"parallel: " << std::endl << omp_get_wtime() - startTime;
	//encodeOneStep("output_parallel", image, image_width, image_height);
    
    /*
    startTime = omp_get_wtime();
        for (unsigned y = image_height; y > 0; y--)
        {
            for (unsigned x = 0; x < image_width; x++)
            {
                color col(0, 0, 0);

                for (int s = 0; s < samples_per_pixel; ++s) {
                    auto u = (x + random_double()) / (image_width - 1);
                    auto v = (y + random_double()) / (image_width - 1);
                    ray ra = cam.get_ray_perspective(u, v);
                    float t = 1;
                    float h = 0;
                    shared_ptr<sdf> temp;
                    for (int i = 0; i < 256; i++)
                    {
                        vec3 pos = ra.orig + normalise(ra.dir) * t;
                        h = findSmallestSDFDistance(pos, sdfWorld, temp);
                        t += h;
                        if (h < 0.0001 || h>10) { break; }
                    }
                    if (h < 0.0001)
                    {
                        vec3 pos = ra.orig + normalise(ra.dir) * t;
                        vec3 normal = calcsdfNormal(pos, sdfWorld);
                        vec3 light = vec3(1, 2, 3);
                        float dif = clamp(dot(normal, normalise(light - pos)), 0., 1.);
                        dif *= 5. / dot(light - pos, light - pos);
                        col += dif * temp->colorOfSDF;
                    }

                }
                

                //std::cout << std::endl<<col;

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
                int val = 4 * image_width * (image_height - y) + 4 * x;
                //allot color
                //int val = 4 * image_width * (image_height- y) + 4 * x;
                image[val + 0] = rf;
                image[val + 1] = gf;
                image[val + 2] = bf;
                image[val + 3] = 255;//alpha
#pragma endregion
            }
            if (y % 10 == 0)
            {
                system("cls");
                std::cout << "Progress: " << int(100 * (double(image_height - y) / image_height)) << "%" << std::endl;
            }
        }
            std::cout << std::endl <<"Time: " << omp_get_wtime() - startTime;
            */
	encodeOneStep("ray_marched_output.png", image, image_width, image_height);
    
    

    
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
