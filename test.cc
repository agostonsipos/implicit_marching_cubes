#include "mc.h"

#include <cmath>
#include <chrono>

using namespace Geometry;

int main(int argc, char** argv)
{
	{
		std::function<double(Vector3D)> sphere = [](Vector3D p) {return p[0] * p[0] + p[1] * p[1] + p[2] * p[2] - 1.0; };
		std::array<Vector3D, 2> boundingBox1 = { Vector3D(-1,-1,-1), Geometry::Vector3D(1,1,1) };

		auto start = std::chrono::high_resolution_clock::now();

		TriMesh mesh = IMC::marching_cubes(sphere, 0, boundingBox1, { 200,200,200 }, true);

		auto end = std::chrono::high_resolution_clock::now();
		std::cerr << "Sphere (200x200) created in: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms" << std::endl;

		mesh.writeOBJ("sphere.obj");
	}
	{
		std::function<double(Vector3D)> something = [](Vector3D p) {return sin(p[0]) * cos(p[1]) + sin(p[1]) * cos(p[2]) + sin(p[2]) * cos(p[0]); };
		std::array<Vector3D, 2> boundingBox2 = { Vector3D(-10,-10,-10), Geometry::Vector3D(10,10,10) };

		auto start = std::chrono::high_resolution_clock::now();

		TriMesh mesh = IMC::marching_cubes(something, 0, boundingBox2, { 200,200,200 }, true);

		auto end = std::chrono::high_resolution_clock::now();
		std::cerr << "Gyroid (200x200) created in: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms" << std::endl;

		mesh.writeOBJ("gyroid.obj");
	}

	return 0;
}
