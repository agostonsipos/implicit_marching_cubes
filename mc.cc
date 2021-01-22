#include "mc.h"

#include <cmath>
#include <vector>

/*
* Implementation based on http://paulbourke.net/geometry/polygonise/
* Geometric data structures from https://github.com/salvipeter/libgeom
*/
namespace IMC {
	namespace Tables {
		extern std::array<std::array<int, 3>, 8> cornerPositions;
		extern std::array<std::array<int, 2>, 12> edgePoints;
		extern std::array<int, 256> edgeTable;
		extern std::array<std::array<int, 16>, 256> triTable;
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
		mu = (isolevel - valp1) / (valp2 - valp1);
		p = p1 + (p2 - p1) * mu;

		return p;
	}

	Geometry::Vector3D FunctionInterp(double isolevel, std::function<double(Geometry::Vector3D)> scalarFunc, Geometry::Vector3D p1, Geometry::Vector3D p2)
	{
		double eps = 1e-5;
		double valp1 = scalarFunc(p1);
		if (std::abs(isolevel - valp1) < eps)
			return p1;
		if (std::abs(isolevel - scalarFunc(p2)) < eps)
			return p2;
		auto p = (p1 + p2) / 2;
		int i = 0;
		double valp = scalarFunc(p);
		while (std::abs(valp) > eps && i < 10) {
			if (valp * valp1 < 0.0) {
				p2 = p;
			}
			else {
				p1 = p;
				valp1 = valp;
			}
			p = (p1 + p2) / 2;
			valp = scalarFunc(p);
			++i;
		}
		return p;
	}

	Geometry::TriMesh marching_cubes(std::function<double(Geometry::Vector3D)> scalarFunc, double isolevel, std::array<Geometry::Vector3D, 2> box, std::array<int, 3> res, bool findRoot)
	{
		Geometry::TriMesh mesh;
		Geometry::PointVector points;
		auto ind3d1 = [res](int x, int y, int z) { return x * res[1] * res[2] + y * res[2] + z; };
		auto ind3d2 = [res](int x, int y, int z) { return x * (res[1] + 1) * (res[2] + 1) + y * (res[2] + 1) + z; };

		std::vector<std::vector<int>> cells(res[0] * res[1] * res[2]);
		int index = 0;
		double cellX = (box[1][0] - box[0][0]) / res[0];
		double cellY = (box[1][1] - box[0][1]) / res[1];
		double cellZ = (box[1][2] - box[0][2]) / res[2];
		double tol = std::sqrt(cellX * cellX + cellY * cellY + cellZ * cellZ) / 1000.0;

		std::vector<Geometry::Vector3D> cellCorners((res[0] + 1) * (res[1] + 1) * (res[2] + 1));
		std::vector<double> cellValues((res[0] + 1) * (res[1] + 1) * (res[2] + 1));
		for (int i = 0; i <= res[0]; ++i) {
			for (int j = 0; j <= res[1]; ++j) {
				for (int k = 0; k <= res[2]; ++k) {
					cellCorners[ind3d2(i, j, k)] = Geometry::Vector3D(
						box[0][0] + cellX * i,
						box[0][1] + cellY * j,
						box[0][2] + cellZ * k
					);
					cellValues[ind3d2(i, j, k)] = scalarFunc(cellCorners[ind3d2(i, j, k)]);
				}
			}
		}

		for (int i = 0; i < res[0]; ++i)
			for (int j = 0; j < res[1]; ++j)
				for (int k = 0; k < res[2]; ++k)
				{
					std::array<Geometry::Vector3D, 8> corners;
					std::array<double, 8> values;
					unsigned char cubeindex = 0;
					for (int c = 0; c < 8; ++c) {
						int x = i + Tables::cornerPositions[c][0], y = j + Tables::cornerPositions[c][1], z = k + Tables::cornerPositions[c][2];
						corners[c] = cellCorners[ind3d2(x, y, z)];
						values[c] = cellValues[ind3d2(x, y, z)];
						if (values[c] < isolevel) cubeindex |= (1 << c);
					}
					std::array<Geometry::Vector3D, 12> vertlist;

					for (int e = 0; e < 12; ++e) {
						if (Tables::edgeTable[cubeindex] & (1 << e)) {
							if (findRoot) {
								vertlist[e] = FunctionInterp(isolevel, scalarFunc,
									corners[Tables::edgePoints[e][0]], corners[Tables::edgePoints[e][1]]);
							}
							else {
								vertlist[e] = VertexInterp(isolevel,
									corners[Tables::edgePoints[e][0]], corners[Tables::edgePoints[e][1]],
									values[Tables::edgePoints[e][0]], values[Tables::edgePoints[e][1]]);
							}
						}
					}

					for (int e = 0; Tables::triTable[cubeindex][e] != -1; e += 3) {
						std::array<Geometry::Vector3D, 3> vertices = {
							vertlist[Tables::triTable[cubeindex][e]],
							vertlist[Tables::triTable[cubeindex][e + 1]],
							vertlist[Tables::triTable[cubeindex][e + 2]]
						};
						int indices[3];
						for (int t = 0; t < 3; ++t) {
							bool foundPoint = false;
							for (int x = 0; x < 8 && !foundPoint; ++x) {
								if (i - ((x >> 2) & 1) < 0 || j - ((x >> 1) & 1) < 0 || k - (x & 1) < 0) continue;
								int cellInd = ind3d1(i - ((x >> 2) & 1), j - ((x >> 1) & 1), k - (x & 1));
								std::vector<int>& list = cells[cellInd];
								for (auto it = list.begin(); it != list.end(); ++it)
								{
									if ((points[*it] - vertices[t]).norm() < tol) {
										indices[t] = *it;
										cells[ind3d1(i, j, k)].push_back(indices[t]);
										foundPoint = true;
										break;
									}
								}
							}
							if (!foundPoint) {
								points.push_back(vertices[t]);
								indices[t] = index++;
								cells[ind3d1(i, j, k)].push_back(indices[t]);
							}
						}
						mesh.addTriangle(indices[0], indices[1], indices[2]);
					}
				}
		mesh.setPoints(points);
		return mesh;
	}
}
