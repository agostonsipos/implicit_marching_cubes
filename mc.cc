#include "mc.h"

#include <cmath>

/*
* Implementation based on http://paulbourke.net/geometry/polygonise/
* Geometric data structures from https://github.com/salvipeter/transfinite
*/
namespace IMC {
	namespace Tables {
		extern int edgePoints[12][2];
		extern int edgeTable[256];
		extern int triTable[256][16];
	}

	Geometry::Vector3D VertexInterp(double isolevel, Geometry::Vector3D p1, Geometry::Vector3D p2, double valp1, double valp2)
	{
		double mu;
		Geometry::Vector3D p;
		double eps = 1e-5;

		if (std::abs(isolevel - valp1) < eps)
			return p1;
		if (std::abs(isolevel - valp2) < eps)
			return p2;
		if (std::abs(valp1 - valp2) < eps)
			return (p1 + p2) / 2;
		mu = (isolevel - valp1) / (valp2 - valp1);
		p[0] = p1[0] + mu * (p2[0] - p1[0]);
		p[1] = p1[1] + mu * (p2[1] - p1[1]);
		p[2] = p1[2] + mu * (p2[2] - p1[2]);

		return p;
	}

	Geometry::TriMesh marching_cubes(std::function<double(Geometry::Vector3D)> scalarFunc, double isolevel, std::array<Geometry::Vector3D, 2> box, std::array<int, 3> res)
	{
		Geometry::TriMesh mesh;
		Geometry::PointVector points;
		int index = 0;
		double cellX = (box[1][0] - box[0][0]) / res[0];
		double cellY = (box[1][1] - box[0][1]) / res[1];
		double cellZ = (box[1][2] - box[0][2]) / res[2];
		for (int i = 0; i < res[0]; ++i)
			for (int j = 0; j < res[1]; ++j)
				for (int k = 0; k < res[2]; ++k)
				{
					std::array<Geometry::Vector3D, 8> corners;
					corners[0] = Geometry::Vector3D(
						box[0][0] + cellX * i,
						box[0][1] + cellY * j,
						box[0][2] + cellZ * k
					);
					corners[1] = corners[0] + Geometry::Vector3D(cellX, 0, 0);
					corners[2] = corners[0] + Geometry::Vector3D(cellX, cellY, 0);
					corners[3] = corners[0] + Geometry::Vector3D(0, cellY, 0);
					corners[4] = corners[0] + Geometry::Vector3D(0, 0, cellZ);
					corners[5] = corners[0] + Geometry::Vector3D(cellX, 0, cellZ);
					corners[6] = corners[0] + Geometry::Vector3D(cellX, cellY, cellZ);
					corners[7] = corners[0] + Geometry::Vector3D(0, cellY, cellZ);
					int cubeindex = 0;
					for (int c = 0; c < 8; ++c)
						if (scalarFunc(corners[c]) < isolevel) cubeindex |= (1 << c);

					std::array<Geometry::Vector3D, 12> vertlist;


					for (int i = 0; i < 12; ++i) {
						if(Tables::edgeTable[cubeindex] & 1 << i)
							vertlist[i] = VertexInterp(isolevel, 
								corners[Tables::edgePoints[i][0]], corners[Tables::edgePoints[i][1]], 
								scalarFunc(corners[Tables::edgePoints[i][0]]), scalarFunc(corners[Tables::edgePoints[i][1]]));
					}

					double tol = std::sqrt(cellX * cellX + cellY * cellY + cellZ * cellZ) / 1000.0;
					for (int e = 0; Tables::triTable[cubeindex][e] != -1; e += 3) {
						std::array<Geometry::Vector3D, 3> vertices = {
							vertlist[Tables::triTable[cubeindex][e]],
							vertlist[Tables::triTable[cubeindex][e + 1]],
							vertlist[Tables::triTable[cubeindex][e + 2]]
						};
						int indices[3] = { -1, -1, -1 };
						for (int k = 0; k < 3; ++k) {
							for (auto it = points.begin(); it != points.end(); ++it)
								if ((*it - vertices[k]).norm() < tol) {
									indices[k] = it - points.begin();
									break;
								}
							if (indices[k] == -1) {
								points.push_back(vertices[k]);
								indices[k] = index++;
							}
						}
						mesh.addTriangle(indices[0], indices[1], indices[2]);
					}
				}
		mesh.setPoints(points);
		return mesh;
	}
}