#ifndef C3D_T_H
#define C3D_T_H
#include"Data_type.h"
#include "2Delaunay.h"
#include <vector>
using namespace std;

struct box_all_vertices
{
	int box_number;
	vector<int> vertice_index;
	int i, j, k;
	bool friend operator<(const box_all_vertices &C1, const box_all_vertices &C2)
	{
		if (C1.box_number<C2.box_number) return true;
		return false;
	}
};
struct box_one_vertice
{
	int box_number;
	int vertice_index;
	int i, j, k;
	bool friend operator<(const box_one_vertice &C1, const box_one_vertice &C2)
	{
		if (C1.box_number<C2.box_number) return true;
		if (C1.box_number == C2.box_number && C1.vertice_index<C2.vertice_index)return true;
		return false;
	}
};

struct box
{
	int box_number;
	int i, j, k;
};


class C3Delaunay_project
{
public:
	vector<C3DPoint> vertices;
	C3Delaunay_project(vector<C3DPoint>);
	vector<vector<int>> Kneighbourhood_generation();
	vector<C3DPoint> mesh_generation();
	



};

#endif