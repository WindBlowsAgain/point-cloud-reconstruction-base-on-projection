#ifndef DELAUNAY_H
#define DELAUNAY_H
#include"Data_type.h"
#include <vector>
using namespace std;

struct triangles_and_circumcircle
{
	double i, j, k, xc, yc, r2;//三角形三个点的索引，及外接圆圆心，半径的平方

}; 
struct supertriangle_vertices
{
	double x1, y1, x2, y2, x3, y3;
};
class CDelaunay
{
	//friend class COpengl1CView;
public:
	CDelaunay(vector<C2DPoint>);
	vector<C2DPoint> plane_vertices;
	supertriangle_vertices supertriangle();
	triangles_and_circumcircle circumcircle(int, int, int);
	vector<C3DPoint> mesh_generation();
};

#endif//!DELAUNAY_H