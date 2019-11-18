#include "mc.h"

#include <cmath>
#include <map>
#include <vector>

/*
* Implementation based on http://paulbourke.net/geometry/polygonise/
* Geometric data structures from https://github.com/salvipeter/transfinite
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
		p[0] = p1[0] + mu * (p2[0] - p1[0]);
		p[1] = p1[1] + mu * (p2[1] - p1[1]);
		p[2] = p1[2] + mu * (p2[2] - p1[2]);

		return p;
	}

	Geometry::TriMesh marching_cubes(std::function<double(Geometry::Vector3D)> scalarFunc, double isolevel, std::array<Geometry::Vector3D, 2> box, std::array<int, 3> res)
	{
		Geometry::TriMesh mesh;
		Geometry::PointVector points;
		std::map<std::tuple<int,int,int>, std::vector<int>> cells;
		int index = 0;
		double cellX = (box[1][0] - box[0][0]) / res[0];
		double cellY = (box[1][1] - box[0][1]) / res[1];
		double cellZ = (box[1][2] - box[0][2]) / res[2];

		std::vector<std::vector<std::vector<Geometry::Vector3D>>> cellCorners(res[0]+1);
		std::vector<std::vector<std::vector<double>>> cellValues(res[0]+1);
		for (int i = 0; i <= res[0]; ++i) {
			cellCorners[i].resize(res[1] + 1);
			cellValues[i].resize(res[1] + 1);
			for (int j = 0; j <= res[1]; ++j) {
				cellCorners[i][j].resize(res[2] + 1);
				cellValues[i][j].resize(res[2] + 1);
				for (int k = 0; k <= res[2]; ++k) {
					cellCorners[i][j][k] = Geometry::Vector3D(
						box[0][0] + cellX * i,
						box[0][1] + cellY * j,
						box[0][2] + cellZ * k
					);
					cellValues[i][j][k] = scalarFunc(cellCorners[i][j][k]);
				}
			}
		}

		for (int i = 0; i < res[0]; ++i)
			for (int j = 0; j < res[1]; ++j)
				for (int k = 0; k < res[2]; ++k)
				{
					std::array<Geometry::Vector3D, 8> corners;
					std::array<double, 8> values;
					int cubeindex = 0;
					for (int c = 0; c < 8; ++c) {
						int x = i + Tables::cornerPositions[c][0], y = j + Tables::cornerPositions[c][1], z = k + Tables::cornerPositions[c][2];
						corners[c] = cellCorners[x][y][z];
						values[c] = cellValues[x][y][z];
						if (values[c] < isolevel) cubeindex |= (1 << c);
					}
					std::array<Geometry::Vector3D, 12> vertlist;

					for (int e = 0; e < 12; ++e) {
						if (Tables::edgeTable[cubeindex] & (1 << e)) {
							vertlist[e] = VertexInterp(isolevel,
								corners[Tables::edgePoints[e][0]], corners[Tables::edgePoints[e][1]],
								values[Tables::edgePoints[e][0]], values[Tables::edgePoints[e][1]]);
						}
					}

					double tol = std::sqrt(cellX * cellX + cellY * cellY + cellZ * cellZ) / 1000.0;
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
								std::vector<int>& list = cells[{i - ((x >> 2) & 1), j - ((x >> 1) & 1), k - (x & 1)}];
								for (auto it = list.begin(); it != list.end(); ++it)
								{
									if ((points[*it] - vertices[t]).norm() < tol) {
										indices[t] = *it;
										cells[{i,j,k}].push_back(indices[t]);
										foundPoint = true;
										break;
									}
								}
							}
							if (!foundPoint) {
								points.push_back(vertices[t]);
								indices[t] = index++;
								cells[{i,j,k}].push_back(indices[t]);
							}
						}
						mesh.addTriangle(indices[0], indices[1], indices[2]);
					}
				}
		mesh.setPoints(points);
		return mesh;
	}
}