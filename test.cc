#include "mc.h"

int main(int argc, char** argv)
{
	std::function<double(Geometry::Vector3D)> sphere = [](Geometry::Vector3D p) {return p[0] * p[0] + p[1] * p[1] + p[2] * p[2] - 1.0; };
	std::array<Geometry::Vector3D, 2> boundingBox = { Geometry::Vector3D(-1,-1,-1), Geometry::Vector3D(1,1,1) };

	Geometry::TriMesh mesh = marching_cubes(sphere, 0, boundingBox, { 20,20,20 });
	mesh.writeOBJ("sphere.obj");
	return 0;
}