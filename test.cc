#include "mc.h"

#include <chrono>

using namespace Geometry;

int main(int argc, char** argv)
{
	auto start = std::chrono::high_resolution_clock::now();
	std::function<double(Vector3D)> sphere = [](Vector3D p) {return p[0] * p[0] + p[1] * p[1] + p[2] * p[2] - 1.0; };
	std::array<Vector3D, 2> boundingBox = { Vector3D(-1,-1,-1), Geometry::Vector3D(1,1,1) };

	TriMesh mesh = IMC::marching_cubes(sphere, 0, boundingBox, { 200,200,200 });

	auto end = std::chrono::high_resolution_clock::now();
	std::cerr << "Time elapsed: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms" << std::endl;

	mesh.writeOBJ("sphere.obj");

	return 0;
}