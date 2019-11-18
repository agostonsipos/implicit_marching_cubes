#include "mc.h"

using namespace Geometry;

int main(int argc, char** argv)
{
	std::function<double(Vector3D)> sphere = [](Vector3D p) {return p[0] * p[0] + p[1] * p[1] + p[2] * p[2] - 1.0; };
	std::array<Vector3D, 2> boundingBox = { Vector3D(-1,-1,-1), Geometry::Vector3D(1,1,1) };

	TriMesh mesh = IMC::marching_cubes(sphere, 0, boundingBox, { 100,100,100 });
	mesh.writeOBJ("sphere.obj");


	return 0;
}