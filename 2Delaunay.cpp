#include "stdafx.h"
#include "2Delaunay.h"
#include <vector>
using namespace std;
CDelaunay::CDelaunay(vector<C2DPoint> temp_vertices)
{
	
	plane_vertices.insert(plane_vertices.begin(), temp_vertices.begin(), temp_vertices.end());
}
supertriangle_vertices CDelaunay::supertriangle()//构建超级三角形
{
	double 	xmin, xmax, ymin, ymax, xmid, ymid, dx, dy, dmax;
	int v_size=plane_vertices.size();
	xmin = plane_vertices[0].x;
	xmax = plane_vertices[v_size - 1].x;
	ymin = 0;
	ymax = 0;
	int i = 0;
	for (i = 0; i < v_size; i++)
	{
		if (plane_vertices[i].y >= ymax) ymax = plane_vertices[i].y;
		if (plane_vertices[i].y <= ymin) ymin = plane_vertices[i].y;
	}
	dx = xmax - xmin;
	dy = ymax - ymin;
	dmax = max(dx, dy);
	xmid = xmin + dx * 0.5;
	ymid = ymin + dy * 0.5;
	supertriangle_vertices temp;
	temp.x1 = xmid - 2 * dmax;
	temp.y1 = ymid - dmax;
	temp.x2 = xmid;
	temp.y2 = ymid + 2 * dmax;
	temp.x3 = xmid + 2 * dmax;
	temp.y3 = ymid - dmax;
	return temp;
}

triangles_and_circumcircle CDelaunay::circumcircle(int i, int j, int k)//构造外接圆
{
	
	   double x1, x2, x3, y1, y2, y3, abs_y1_y2, abs_y2_y3, xc, yc, m1, m2, mx1, mx2, my1, my2, dx, dy,r2;
	    x1 = plane_vertices[i].x;
		y1 = plane_vertices[i].y;
		x2 = plane_vertices[j].x;
		y2 = plane_vertices[j].y;
		x3 = plane_vertices[k].x;
		y3 = plane_vertices[k].y;
		
		double t1, t2, t3, temp, x, y;

		t1= x1*x1 + y1*y1;
		t2= x2*x2 + y2*y2;
		t3= x3*x3 + y3*y3;
		temp= x1*y2 + x2*y3 + x3*y1 - x1*y3 - x2*y1 - x3*y2;
		x= (t2*y3 + t1*y2 + t3*y1 - t2*y1 - t3*y2 - t1*y3) / temp / 2;
		y= (t3*x2 + t2*x1 + t1*x3 - t1*x2 - t2*x3 - t3*x1) / temp / 2;

		double a = sqrt((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2));
		double b = sqrt((x1 - x3)*(x1 - x3) + (y1 - y3)*(y1 - y3));
		double c = sqrt((x2 - x3)*(x2 - x3) + (y2 - y3)*(y2 - y3));
		double p = (a + b + c) / 2;
		double S = sqrt(p*(p - a)*(p - b)*(p - c));
		double radius = a*b*c / (4 * S);
		triangles_and_circumcircle c_temp;
		c_temp.i = i;
		c_temp.j = j;
		c_temp.k = k;
		c_temp.xc = x;
		c_temp.yc = y;
		c_temp.r2 = radius*radius;
		return c_temp;
}  

vector<C3DPoint> CDelaunay::mesh_generation()
{
	//将超级三角形存入顶点容器
	int n = plane_vertices.size();
	supertriangle_vertices temp; 
	temp=supertriangle();
	C2DPoint p;
	C2DPoint q;
	C2DPoint w;
	p.x = temp.x1;
	p.y = temp.y1;
	q.x = temp.x2;
	q.y = temp.y2;
	w.x = temp.x3;
	w.y = temp.y3;
	plane_vertices.push_back(p);
	plane_vertices.push_back(q);
	plane_vertices.push_back(w);
	
	vector<triangles_and_circumcircle> temp_triangles;//保存的是三角形的三个点及其外接圆的圆心，半径
	vector<triangles_and_circumcircle> triangles;
	vector<int> edges;
	triangles_and_circumcircle temp_t;
	temp_t = circumcircle( n, n + 1,  n + 2);
	temp_triangles.push_back(temp_t);
	for (int i = 0; i < n; i++)//递增地将每个顶点添加到网格中
	{
		
		for (int j = 0; j < temp_triangles.size(); j++)//遍历temp triangles中的每一个三角形
		{
			double dx, dy;
			dx = plane_vertices[i].x - temp_triangles[j].xc;
			//如果这个点在这个三角形的外圆的右边，那么这个三角形就不会再被检查了。从打开列表中删除它，将其添加到关闭列表中，并跳过。
			if (dx > 0.0 && dx * dx > temp_triangles[j].r2)
			{
				triangles.push_back(temp_triangles[j]);
				temp_triangles.erase(temp_triangles.begin() + j);
				j--;
				continue;
			}
			//如果这个点在圆周外，跳过这个三角形
			dy = plane_vertices[i].y - temp_triangles[j].yc;
			if (dx * dx + dy * dy >temp_triangles[j].r2)
				continue;
			//如果这个点在圆周内，则移除三角形并将它的边添加到边列表中
			edges.push_back(temp_triangles[j].i);
			edges.push_back(temp_triangles[j].j);
			edges.push_back(temp_triangles[j].j);
			edges.push_back(temp_triangles[j].k);
			edges.push_back(temp_triangles[j].k);
			edges.push_back(temp_triangles[j].i);
			temp_triangles.erase(temp_triangles.begin() + j);
			j--;//因为删除一个元素后，剩下的会前进一位
		}

		//删除重复边，注意并非只是将重复的部分去掉, 而是将发生重复的 edge 全部去除, 即不仅去除重复的部分, 原本的 edge 也要去掉).
		for (int i = 0; i < edges.size(); i+=2)
		{
			int a = edges[i];
			int b = edges[i+1];

			for (int j = i+2; j< edges.size();j+=2)
			{
				int m= edges[j];
				int n= edges[j+1];

				if ((a ==m && b ==n) || (a ==n && b ==m))
				{
					edges.erase(edges.begin() + i, edges.begin() + (i + 2));
					edges.erase(edges.begin() + (j-2), edges.begin() + j );
					i -= 2;//因为删除两个元素后，剩下的会前进两位
					break;
				}
			}
	   	}
		//为每条边添加一个新的三角形
		for (int j = 0; j < edges.size(); j+=2)
		{
			int a = edges[j];
			int b = edges[j + 1];
			temp_triangles.push_back(circumcircle(a, b, i));
		}
		vector<int>().swap(edges);
	}
	//删除与超三角形共享一个顶点的三角形，构建一个表示三角形的三联列表。
	for (int i = 0; i < temp_triangles.size(); i++)
	{
		triangles.push_back(temp_triangles[i]);
	}

	vector<C3DPoint>v_triangles;//保存三角形三个点的下标
	for (int j = 0; j < triangles.size(); j++)
	if (triangles[j].i < n && triangles[j].j < n && triangles[j].k < n)
	{
		C3DPoint p;
		p.x = triangles[j].i;
		p.y = triangles[j].j;
		p.z = triangles[j].k;
		v_triangles.push_back(p);
	}
	/*vector<C3DPoint> vertices_temp;
	for (int i = 0; i < v_triangles.size(); i++)
	{
		C3DPoint p;
		C3DPoint q;
		C3DPoint w;
		p.x = plane_vertices[v_triangles[i].x].x;
		p.y = plane_vertices[v_triangles[i].x].y;
		p.z = 0;
		q.x = plane_vertices[v_triangles[i].y].x;
		q.y = plane_vertices[v_triangles[i].y].y;
		q.z = 0;
		w.x = plane_vertices[v_triangles[i].z].x;
		w.y = plane_vertices[v_triangles[i].z].y;
		w.z = 0;
		vertices_temp.push_back(p);
		vertices_temp.push_back(q);
		vertices_temp.push_back(w);
	}
	return  vertices_temp;*/
	return  v_triangles;
}