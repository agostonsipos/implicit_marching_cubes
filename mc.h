#include "geometry.hh"

#include <functional>

// Constructs a triangle mesh approximating the isosurface of a given scalar function

namespace IMC { // implicit marching cubes

	Geometry::TriMesh marching_cubes(
		std::function<double(Geometry::Vector3D)> scalarFunc,
		double isolevel,
		std::array<Geometry::Vector3D, 2> boundingBox,
		std::array<int, 3> resolution
	);

}