#include "mc.h"

/*
* Implementation based on http://paulbourke.net/geometry/polygonise/
* Geometric data structures from https://github.com/salvipeter/transfinite
*/
namespace IMC {
	namespace Tables {
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

					if (Tables::edgeTable[cubeindex] & 1)
						vertlist[0] =
						VertexInterp(isolevel, corners[0], corners[1], scalarFunc(corners[0]), scalarFunc(corners[1]));
					if (Tables::edgeTable[cubeindex] & 2)
						vertlist[1] =
						VertexInterp(isolevel, corners[1], corners[2], scalarFunc(corners[1]), scalarFunc(corners[2]));
					if (Tables::edgeTable[cubeindex] & 4)
						vertlist[2] =
						VertexInterp(isolevel, corners[2], corners[3], scalarFunc(corners[2]), scalarFunc(corners[3]));
					if (Tables::edgeTable[cubeindex] & 8)
						vertlist[3] =
						VertexInterp(isolevel, corners[3], corners[0], scalarFunc(corners[3]), scalarFunc(corners[0]));
					if (Tables::edgeTable[cubeindex] & 16)
						vertlist[4] =
						VertexInterp(isolevel, corners[4], corners[5], scalarFunc(corners[4]), scalarFunc(corners[5]));
					if (Tables::edgeTable[cubeindex] & 32)
						vertlist[5] =
						VertexInterp(isolevel, corners[5], corners[6], scalarFunc(corners[5]), scalarFunc(corners[6]));
					if (Tables::edgeTable[cubeindex] & 64)
						vertlist[6] =
						VertexInterp(isolevel, corners[6], corners[7], scalarFunc(corners[6]), scalarFunc(corners[7]));
					if (Tables::edgeTable[cubeindex] & 128)
						vertlist[7] =
						VertexInterp(isolevel, corners[7], corners[4], scalarFunc(corners[7]), scalarFunc(corners[4]));
					if (Tables::edgeTable[cubeindex] & 256)
						vertlist[8] =
						VertexInterp(isolevel, corners[0], corners[4], scalarFunc(corners[0]), scalarFunc(corners[4]));
					if (Tables::edgeTable[cubeindex] & 512)
						vertlist[9] =
						VertexInterp(isolevel, corners[1], corners[5], scalarFunc(corners[1]), scalarFunc(corners[5]));
					if (Tables::edgeTable[cubeindex] & 1024)
						vertlist[10] =
						VertexInterp(isolevel, corners[2], corners[6], scalarFunc(corners[2]), scalarFunc(corners[6]));
					if (Tables::edgeTable[cubeindex] & 2048)
						vertlist[11] =
						VertexInterp(isolevel, corners[3], corners[7], scalarFunc(corners[3]), scalarFunc(corners[7]));


					double tol = std::sqrt(cellX * cellX + cellY * cellY + cellZ * cellZ) / 1000.0;
					for (int e = 0; Tables::triTable[cubeindex][e] != -1; e += 3) {
						Geometry::Vector3D a = vertlist[Tables::triTable[cubeindex][e]];
						Geometry::Vector3D b = vertlist[Tables::triTable[cubeindex][e + 1]];
						Geometry::Vector3D c = vertlist[Tables::triTable[cubeindex][e + 2]];
						int indA = -1, indB = -1, indC = -1;
						for (auto it = points.begin(); it != points.end(); ++it)
							if ((*it - a).norm() < tol) {
								indA = it - points.begin();
								break;
							}
						if (indA == -1) {
							points.push_back(a);
							indA = index++;
						}
						for (auto it = points.begin(); it != points.end(); ++it)
							if ((*it - b).norm() < tol) {
								indB = it - points.begin();
								break;
							}
						if (indB == -1) {
							points.push_back(b);
							indB = index++;
						}
						for (auto it = points.begin(); it != points.end(); ++it)
							if ((*it - c).norm() < tol) {
								indC = it - points.begin();
								break;
							}
						if (indC == -1) {
							points.push_back(c);
							indC = index++;
						}
						mesh.addTriangle(indA, indB, indC);
					}
				}
		mesh.setPoints(points);
		return mesh;
	}
}