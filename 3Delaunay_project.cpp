#include "stdafx.h"
#include "3Delaunay_project.h"
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <vector>
#include "windows.h"
using namespace Eigen;
using namespace std;
C3Delaunay_project::C3Delaunay_project(vector<C3DPoint> temp_vertices)
{
	vertices.insert(vertices.begin(), temp_vertices.begin(), temp_vertices.end());
}
vector<vector<int>>  C3Delaunay_project::Kneighbourhood_generation()
{
	vector<vector<int>> temp_Kneighbourhood;
	temp_Kneighbourhood.resize(vertices.size());
	vector<box> vertices_box;//对应索引点所在的包围盒的编号
	set<box_one_vertice>s_box_one_vertice;//保存包围盒编号和里面的一个点的索引
	vector<box_one_vertice>v_box_one_vertice;
	set<box_all_vertices>s_box_all_vertices; // 保存包围盒编号和里面的所有点的索引
	vertices_box.resize(vertices.size());
	/*******************************************************************************/
	//构造包围盒
	double xmin, xmax, ymin, ymax, zmin, zmax;
	double Lx, Ly, Lz;//大包围盒的边长
	double box_L; //小立方体的边长
	double nx, ny, nz;//小立方体在xyz方向上的个数；
	int  everybox_point_number = 20;
	int  k_point_number = 15;//K邻域数目
	xmin = 0;
	xmax = 0;
	ymin = 0;
	ymax = 0;
	zmin = 0;
	zmax = 0;
	for (int i = 0; i < vertices.size(); i++)
	{
		if (vertices[i].x < xmin) xmin = vertices[i].x;
		if (vertices[i].x > xmax) xmax = vertices[i].x;
		if (vertices[i].y < ymin) ymin = vertices[i].y;
		if (vertices[i].y > ymax) ymax = vertices[i].y;
		if (vertices[i].z < zmin) zmin = vertices[i].z;
		if (vertices[i].z > zmax) zmax = vertices[i].z;
	}
	
	Lx = xmax - xmin;
	Ly = ymax - ymin;
	Lz = zmax - zmin;
	double S = (Lx*Ly + Lx*Lz + Ly*Lz) * 2;
	double V = Lx*Ly*Lz;
	//box_L = sqrt((S*everybox_point_number) / vertices.size());
	//box_L = pow((V*everybox_point_number) / vertices.size(),1.0/3);
	//box_L = sqrt(S/ vertices.size());
	box_L = pow(V / vertices.size(), 1.0 / 3);

	xmin -= 1 * box_L;
	xmax += 1 * box_L;
	ymin -= 1* box_L;
	ymax += 1 * box_L;
	zmin -= 1 * box_L;
	zmax += 1 * box_L;
	nx = ceil((xmax - xmin) / box_L);
	ny = ceil((ymax - ymin) / box_L);
	nz = ceil((zmax - zmin) / box_L);
	//计算有多少个非空包围盒，为二次划分做准备
	set<int> s_box;
	for (int n = 0; n < vertices.size(); n++)
	{
		box_one_vertice temp_box_one_vertice;
		temp_box_one_vertice.i = ceil((vertices[n].x - xmin) / box_L) - 1;
		temp_box_one_vertice.j = ceil((vertices[n].y - ymin) / box_L) - 1;
		temp_box_one_vertice.k = ceil((vertices[n].z - zmin) / box_L) - 1;
		temp_box_one_vertice.box_number = temp_box_one_vertice.i*ny*nz + temp_box_one_vertice.j*nz + temp_box_one_vertice.k;
		s_box.insert(temp_box_one_vertice.box_number);
	}

	double p = vertices.size() / (s_box.size()*box_L*box_L*box_L);

	box_L = pow(everybox_point_number/p, 1.0 / 3);

	/******************************************/
	

	s_box.clear();
	for (int n = 0; n < vertices.size(); n++)
	{
		box_one_vertice temp_box_one_vertice;
		temp_box_one_vertice.i = ceil((vertices[n].x - xmin) / box_L) - 1;
		temp_box_one_vertice.j = ceil((vertices[n].y - ymin) / box_L) - 1;
		temp_box_one_vertice.k = ceil((vertices[n].z - zmin) / box_L) - 1;
		temp_box_one_vertice.box_number = temp_box_one_vertice.i*ny*nz + temp_box_one_vertice.j*nz + temp_box_one_vertice.k;
		s_box.insert(temp_box_one_vertice.box_number);
	}

	p = vertices.size() / (s_box.size()*box_L*box_L*box_L);
	box_L = pow(everybox_point_number / p, 1.0 / 3);
	

	/*******************************************/
	
	for (int n = 0; n < vertices.size(); n++)
	{
		box_one_vertice temp_box_one_vertice;
		temp_box_one_vertice.i = ceil((vertices[n].x - xmin) / box_L) - 1;
		temp_box_one_vertice.j = ceil((vertices[n].y - ymin) / box_L) - 1;
		temp_box_one_vertice.k = ceil((vertices[n].z - zmin) / box_L) - 1;
		temp_box_one_vertice.box_number = temp_box_one_vertice.i*ny*nz + temp_box_one_vertice.j*nz + temp_box_one_vertice.k;
		vertices_box[n].box_number = temp_box_one_vertice.box_number;
		vertices_box[n].i = temp_box_one_vertice.i;
		vertices_box[n].j = temp_box_one_vertice.j;
		vertices_box[n].k = temp_box_one_vertice.k;
		temp_box_one_vertice.vertice_index = n;
		s_box_one_vertice.insert(temp_box_one_vertice);
	}

	v_box_one_vertice.resize(s_box_one_vertice.size());
	copy(s_box_one_vertice.begin(), s_box_one_vertice.end(), v_box_one_vertice.begin());
	s_box_one_vertice.clear();


	box_all_vertices temp_box_all_vertices;
	temp_box_all_vertices.box_number = v_box_one_vertice[0].box_number;
	temp_box_all_vertices.i = v_box_one_vertice[0].i;
	temp_box_all_vertices.j = v_box_one_vertice[0].j;
	temp_box_all_vertices.k = v_box_one_vertice[0].k;
	for (int i = 0; i < v_box_one_vertice.size(); i++)
	{
		if (temp_box_all_vertices.box_number == v_box_one_vertice[i].box_number)
		{
			temp_box_all_vertices.vertice_index.push_back(v_box_one_vertice[i].vertice_index);
			if (i == (v_box_one_vertice.size() - 1))
				s_box_all_vertices.insert(temp_box_all_vertices);
		}

		else
		{
			s_box_all_vertices.insert(temp_box_all_vertices);
			//重置参数
			temp_box_all_vertices.box_number = v_box_one_vertice[i].box_number;
			temp_box_all_vertices.i = v_box_one_vertice[i].i;
			temp_box_all_vertices.j = v_box_one_vertice[i].j;
			temp_box_all_vertices.k = v_box_one_vertice[i].k;
			vector<int>().swap(temp_box_all_vertices.vertice_index);
			temp_box_all_vertices.vertice_index.push_back(v_box_one_vertice[i].vertice_index);
			if (i == (v_box_one_vertice.size() - 1))
				s_box_all_vertices.insert(temp_box_all_vertices);

		}
	}
	vector<box_one_vertice>().swap(v_box_one_vertice);

	/*******************************************************************************/
	//k邻域搜索
	

	long start_time = GetTickCount();//测试程序运行时间

	

	int error_vertice = 0;
	vector<int> test_searched_box_number;//记录搜索了多少个包围盒才得到K邻域
	vector<int> test_searched_vertices_number;//记录搜索了多少个点才得到K邻域
	test_searched_box_number.resize(vertices_box.size());
	test_searched_vertices_number.resize(vertices_box.size());
	for (int j = 0; j < vertices_box.size(); j++)
	{
		box_all_vertices  center_box;
		center_box.box_number = vertices_box[j].box_number;
		center_box.i = vertices_box[j].i;
		center_box.j = vertices_box[j].j;
		center_box.k = vertices_box[j].k;
		vector<int>().swap(center_box.vertice_index);

	
		set<box_all_vertices>::iterator it = s_box_all_vertices.find(center_box);
		
		

		vector<int> neighbourhood_vertices((*it).vertice_index);
		set<C2DPoint>s_distance_and_vertice;
		vector<C2DPoint>v_distance_and_vertice;

		
		
		for (int n = 0; n < neighbourhood_vertices.size(); n++)
		{
			
			C2DPoint temp_distance;
			temp_distance.x = sqrt((vertices[j].x - vertices[neighbourhood_vertices[n]].x)*(vertices[j].x - vertices[neighbourhood_vertices[n]].x) + (vertices[j].y - vertices[neighbourhood_vertices[n]].y)*(vertices[j].y - vertices[neighbourhood_vertices[n]].y) + (vertices[j].z - vertices[neighbourhood_vertices[n]].z)*(vertices[j].z - vertices[neighbourhood_vertices[n]].z));
			temp_distance.y = neighbourhood_vertices[n];
			
			s_distance_and_vertice.insert(temp_distance);
		}
		//求s_distance_and_vertice.insert【k_point_number】.x的值
		set<C2DPoint>::iterator it3 = s_distance_and_vertice.begin();
		if (s_distance_and_vertice.size() >= k_point_number + 1)
		{
			for (int i = 0; i <= k_point_number - 1; i++)
			{
				it3++;
			}
		}

		
		
		//判断是不是最优K邻域
		
		//计算该点到中心包围盒六个面的距离
		set<double> s_vertice_to_face_distance;
		s_vertice_to_face_distance.insert(vertices[j].x - (xmin + center_box.i*box_L) + 0.000000001*box_L);
		s_vertice_to_face_distance.insert((xmin + (center_box.i + 1)*box_L) - vertices[j].x + 0.000000002*box_L);
		s_vertice_to_face_distance.insert(vertices[j].y - (ymin + center_box.j*box_L) + 0.000000003*box_L);
		s_vertice_to_face_distance.insert((ymin + (center_box.j + 1)*box_L) - vertices[j].y + 0.000000004*box_L);
		s_vertice_to_face_distance.insert(vertices[j].z - (zmin + center_box.k*box_L) + 0.000000005*box_L);
		s_vertice_to_face_distance.insert((zmin + (center_box.k + 1)*box_L) - vertices[j].z + 0.000000006*box_L);
		if (s_distance_and_vertice.size() >= k_point_number + 1 && (*it3).x <= *(s_vertice_to_face_distance.begin()))
		{
			set<C2DPoint>::iterator it = s_distance_and_vertice.begin();
			it++;
			temp_Kneighbourhood[j].resize(k_point_number);
			for (int s = 0; s <k_point_number; s++)
			{
				temp_Kneighbourhood[j][s] = (*it).y;
				it++;
			}
			test_searched_box_number[j] = 1;
			test_searched_vertices_number[j] = s_distance_and_vertice.size();
		}
		//如果不是最优K邻域,就再扩充到多个立方体
		else
		{

			vector<int>().swap(neighbourhood_vertices);
			//s_distance_and_vertice.clear();
			//vector<C2DPoint>().swap(v_distance_and_vertice);
			//搜索其他26个立方体
			int a[27];
			//中间九个
			
			a[0] = center_box.box_number - ny*nz - nz;
			a[1] = center_box.box_number - ny*nz;
			a[2] = center_box.box_number - ny*nz + nz;
			a[3] = center_box.box_number - nz;
			a[4] = center_box.box_number;
			a[5] = center_box.box_number + nz;
			a[6] = center_box.box_number + ny*nz - nz;
			a[7] = center_box.box_number + ny*nz;
			a[8] = center_box.box_number + ny*nz + nz;
			//下面九个
			for (int k = 9; k <= 17; k++)
			{
				a[k] = a[k - 9] - 1;
			}
			//上面九个
			for (int k = 18; k <= 26; k++)
			{
				a[k] = a[k - 18] + 1;
			}

			//计算该点到中心包围盒六个面的距离
			double dx1, dx2, dy1, dy2, dz1, dz2;
			dx1 = vertices[j].x - (xmin + center_box.i*box_L);
			dx2 = (xmin + (center_box.i + 1)*box_L) - vertices[j].x;
			dy1 = vertices[j].y - (ymin + center_box.j*box_L);
			dy2 = (ymin + (center_box.j + 1)*box_L) - vertices[j].y;
			dz1 = vertices[j].z - (zmin + center_box.k*box_L);
			dz2 = (zmin + (center_box.k + 1)*box_L) - vertices[j].z;
			//计算该点到26个立方体的距离
			vector<C2DPoint> v_box_to_vertice_distance;
			v_box_to_vertice_distance.resize(27);
			for (int m = 0; m < 27; m++)
			{
				v_box_to_vertice_distance[m].y = a[m];
			}
			//中间九个
			
			v_box_to_vertice_distance[0].x = sqrt(dx1*dx1 + dy1*dy1);
			v_box_to_vertice_distance[1].x = dx1;
			v_box_to_vertice_distance[2].x = sqrt(dx1*dx1 + dy2*dy2);
			v_box_to_vertice_distance[3].x = dy1;
			v_box_to_vertice_distance[4].x = 0;
			v_box_to_vertice_distance[5].x = dy2;
			v_box_to_vertice_distance[6].x = sqrt(dx2*dx2 + dy1*dy1);
			v_box_to_vertice_distance[7].x = dx2;
			v_box_to_vertice_distance[8].x = sqrt(dx2*dx2 + dy2*dy2);
			//下面九个
			v_box_to_vertice_distance[9].x = sqrt(dx1*dx1 + dy1*dy1 + dz1*dz1);
			v_box_to_vertice_distance[10].x = sqrt(dx1*dx1 + dz1*dz1);
			v_box_to_vertice_distance[11].x = sqrt(dx1*dx1 + dy2*dy2 + dz1*dz1);
			v_box_to_vertice_distance[12].x = sqrt(dy1*dy1 + dz1*dz1);
			v_box_to_vertice_distance[13].x = dz1;
			v_box_to_vertice_distance[14].x = sqrt(dy2*dy2 + dz1*dz1);
			v_box_to_vertice_distance[15].x = sqrt(dx2*dx2 + dy1*dy1 + dz1*dz1);
			v_box_to_vertice_distance[16].x = sqrt(dx2*dx2 + dz1*dz1);
			v_box_to_vertice_distance[17].x = sqrt(dx2*dx2 + dy2*dy2 + dz1*dz1);
			//上面九个
			v_box_to_vertice_distance[18].x = sqrt(dx1*dx1 + dy1*dy1 + dz2*dz2);
			v_box_to_vertice_distance[19].x = sqrt(dx1*dx1 + dz2*dz2);
			v_box_to_vertice_distance[20].x = sqrt(dx1*dx1 + dy2*dy2 + dz2*dz2);
			v_box_to_vertice_distance[21].x = sqrt(dy1*dy1 + dz2*dz2);
			v_box_to_vertice_distance[22].x = dz2;
			v_box_to_vertice_distance[23].x = sqrt(dy2*dy2 + dz2*dz2);
			v_box_to_vertice_distance[24].x = sqrt(dx2*dx2 + dy1*dy1 + dz2*dz2);
			v_box_to_vertice_distance[25].x = sqrt(dx2*dx2 + dz2*dz2);
			v_box_to_vertice_distance[26].x = sqrt(dx2*dx2 + dy2*dy2 + dz2*dz2);

			
			set<C2DPoint> s_box_to_vertice_distance;
			s_box_to_vertice_distance.insert(v_box_to_vertice_distance.begin(), v_box_to_vertice_distance.end());
			copy(s_box_to_vertice_distance.begin(), s_box_to_vertice_distance.end(), v_box_to_vertice_distance.begin());
			//该点到中心包围盒六个面的距离
			vector<double> v_vertice_to_face_distance;
			v_vertice_to_face_distance.resize(s_vertice_to_face_distance.size()*2);
			copy(s_vertice_to_face_distance.begin(), s_vertice_to_face_distance.end(), v_vertice_to_face_distance.begin());
		
			for (int v = 0; v < 6; v++)
			{
				v_vertice_to_face_distance[v + 6] = v_vertice_to_face_distance[v] + box_L;
			}
	
			int t = 1;
			BOOL bool_set_erase = TRUE;
			for (int f = 1; f < 13;f++)
			{
				if (f == 12)
				{
					error_vertice++;
					break;

				}
				if (s_distance_and_vertice.size() < k_point_number + 1)
				{
					for (; t < 27 && v_box_to_vertice_distance[t].x < v_vertice_to_face_distance[f]; t++)
					{

						box_all_vertices  box;

						box.box_number = v_box_to_vertice_distance[t].y;
						set<box_all_vertices>::iterator it = s_box_all_vertices.find(box);

						if (it != s_box_all_vertices.end())
						{
							vector<int> v_box_inside_vertices((*it).vertice_index);

							for (int n = 0; n < v_box_inside_vertices.size(); n++)
							{

								C2DPoint temp_distance;
								temp_distance.x = sqrt((vertices[j].x - vertices[v_box_inside_vertices[n]].x)*(vertices[j].x - vertices[v_box_inside_vertices[n]].x) + (vertices[j].y - vertices[v_box_inside_vertices[n]].y)*(vertices[j].y - vertices[v_box_inside_vertices[n]].y) + (vertices[j].z - vertices[v_box_inside_vertices[n]].z)*(vertices[j].z - vertices[v_box_inside_vertices[n]].z));
								temp_distance.y = v_box_inside_vertices[n];

								s_distance_and_vertice.insert(temp_distance);
							}

						}
					}
				}
				else
				{
					if (bool_set_erase == TRUE)
					{
						set<C2DPoint>::iterator it3 = s_distance_and_vertice.begin();
						if (s_distance_and_vertice.size() >= k_point_number + 1)
						{

							for (int i = 0; i <= k_point_number; i++)
							{
								it3++;
							}
						}
						s_distance_and_vertice.erase(it3, s_distance_and_vertice.end());
						bool_set_erase = FALSE;
					}
					for (; t < 27 && v_box_to_vertice_distance[t].x < v_vertice_to_face_distance[f]; t++)
					{

						box_all_vertices  box;

						box.box_number = v_box_to_vertice_distance[t].y;
						set<box_all_vertices>::iterator it = s_box_all_vertices.find(box);

						if (it != s_box_all_vertices.end())
						{
							vector<int> v_box_inside_vertices((*it).vertice_index);

							for (int n = 0; n < v_box_inside_vertices.size(); n++)
							{

								C2DPoint temp_distance;
								temp_distance.x = sqrt((vertices[j].x - vertices[v_box_inside_vertices[n]].x)*(vertices[j].x - vertices[v_box_inside_vertices[n]].x) + (vertices[j].y - vertices[v_box_inside_vertices[n]].y)*(vertices[j].y - vertices[v_box_inside_vertices[n]].y) + (vertices[j].z - vertices[v_box_inside_vertices[n]].z)*(vertices[j].z - vertices[v_box_inside_vertices[n]].z));
								temp_distance.y = v_box_inside_vertices[n];
								set<C2DPoint>::iterator it = s_distance_and_vertice.end();
								it--;
								if (temp_distance.x < (*it).x)
								{
									s_distance_and_vertice.erase(it);
									s_distance_and_vertice.insert(temp_distance);
								}
								
							}

						}
					}

				}
				set<C2DPoint>::iterator it3 = s_distance_and_vertice.end();
				it3--;
				if (f <= 5)
				{

					if (!(s_distance_and_vertice.size() >= k_point_number + 1 && (*it3).x <= v_vertice_to_face_distance[f]) && t == 27)
					{
						error_vertice++;
						break;
					}

					if (s_distance_and_vertice.size() >= k_point_number + 1 && (*it3).x <= v_vertice_to_face_distance[f])
					{
						test_searched_box_number[j] = t;
						test_searched_vertices_number[j] = s_distance_and_vertice.size();
						break;
					}
				}
				else
				{
					if (!(s_distance_and_vertice.size() >= k_point_number + 1 && (*it3).x <= v_vertice_to_face_distance[6]) && t == 27)
					{
						error_vertice++;
						break;
					}

					if (s_distance_and_vertice.size() >= k_point_number + 1 && (*it3).x <= v_vertice_to_face_distance[6])
					{
						test_searched_box_number[j] = t;
						test_searched_vertices_number[j] = s_distance_and_vertice.size();
						break;
					}
				}
					
			} 
			temp_Kneighbourhood[j].resize(k_point_number);
			set<C2DPoint>::iterator it = s_distance_and_vertice.begin();
			it++;
			for (int s = 0; s <k_point_number; s++)
			{
				temp_Kneighbourhood[j][s] = (*it).y;
				it++;
			}

		}
	}

		long end_time = GetTickCount();
		long sum_time = end_time - start_time;//程序运行时间，单位ms
		
		
		int searched_box[28] = { 0 };
		int sum_searched_vertices = 0;
		for (int i = 0; i < test_searched_vertices_number.size(); i++)
		{
			sum_searched_vertices += test_searched_vertices_number[i];
		}
		double average_searched_vertices = (sum_searched_vertices + 0.1) / vertices.size();
		for (int d = 0; d < test_searched_box_number.size(); d++)
		{
			searched_box[test_searched_box_number[d]]++;
		}
		/*
		//验证结果准确性的测试程序
		long start_time2 = GetTickCount();//测试程序运行时间
		vector<vector<int>> test_Kneighbourhood;
		test_Kneighbourhood.resize(vertices.size());
		for (int j = 0; j < vertices.size(); j++)
		{
		set<C2DPoint> s_distance_and_vertice;
		for (int n = 0; n < vertices.size(); n++)
		{
		C2DPoint temp_distance;
		temp_distance.x = sqrt((vertices[j].x - vertices[n].x)*(vertices[j].x - vertices[n].x) + (vertices[j].y - vertices[n].y)*(vertices[j].y - vertices[n].y) + (vertices[j].z - vertices[n].z)*(vertices[j].z - vertices[n].z));
		temp_distance.y = n;
		s_distance_and_vertice.insert(temp_distance);
		}
		test_Kneighbourhood[j].resize(k_point_number);
		set<C2DPoint>::iterator it = s_distance_and_vertice.begin();
		it++;
		for (int i = 0; i < k_point_number; i++)
		{
		test_Kneighbourhood[j][i] = (*it).y;
		it++;
		}
		}

		long end_time2 = GetTickCount();
		long sum_time2 = end_time2 - start_time2;//程序运行时间，单位ms*/


		return  temp_Kneighbourhood;
	
}
vector<C3DPoint> C3Delaunay_project::mesh_generation()
{
	vector<C3DPoint> v_triangles;
	set<C3DPoint> s_triangles;
	vector<vector<int>> v_Kneighbourhood(Kneighbourhood_generation());//求K邻域
	for (int i = 0; i < v_Kneighbourhood.size(); i++)
	{
		//用最小二乘法拟合法向量
		double sum_x, sum_y, sum_z, average_x, average_y, average_z;
		sum_x = 0;
		sum_y = 0;
		sum_z = 0;
		vector<C3DPoint> temp_vertices;
		temp_vertices.resize(v_Kneighbourhood[i].size()+1);
		temp_vertices[0].x = vertices[i].x;
		temp_vertices[0].y = vertices[i].y;
		temp_vertices[0].z = vertices[i].z;
		sum_x += temp_vertices[0].x;
		sum_y += temp_vertices[0].y;
		sum_z += temp_vertices[0].z;
		for (int j = 0; j < v_Kneighbourhood[i].size(); j++)
		{
			temp_vertices[j + 1].x = vertices[v_Kneighbourhood[i][j]].x;
			temp_vertices[j + 1].y = vertices[v_Kneighbourhood[i][j]].y;
			temp_vertices[j + 1].z = vertices[v_Kneighbourhood[i][j]].z;
			sum_x += temp_vertices[j + 1].x;
			sum_y += temp_vertices[j + 1].y;
			sum_z += temp_vertices[j + 1].z;
		}
		average_x = sum_x / temp_vertices.size();
		average_y = sum_y / temp_vertices.size();
		average_z = sum_z / temp_vertices.size();

		
		double matrix[9] = { 0 };//定义矩阵A
		for (int k = 0; k < temp_vertices.size(); k++)
		{
			
			double deflection_x = temp_vertices[k].x - average_x;
			double deflection_y = temp_vertices[k].y - average_y;
			double deflection_z = temp_vertices[k].z - average_z;
			matrix[0] += deflection_x*deflection_x;
			matrix[1] += deflection_x*deflection_y;
			matrix[2] += deflection_x*deflection_z;
			matrix[3] += deflection_x*deflection_y;
			matrix[4] += deflection_y*deflection_y;
			matrix[5] += deflection_y*deflection_z;
			matrix[6] += deflection_x*deflection_z;
			matrix[7] += deflection_y*deflection_z;
			matrix[8] += deflection_z*deflection_z;
		}

		//用Eigen求特征值和特征向量
		Matrix3d A;
		A << matrix[0], matrix[1], matrix[2], matrix[3], matrix[4], matrix[5], matrix[6], matrix[7], matrix[8];
		EigenSolver<Matrix3d> es(A);
		Matrix3d Eigenvalue = es.pseudoEigenvalueMatrix();
		Matrix3d Eigenvectors = es.pseudoEigenvectors();
	
		double value[9];
		for (int i = 0; i < Eigenvalue.size(); i++)
			value[i] = *(Eigenvalue.data() + i);
		double vectors[9];
		for (int i = 0; i < Eigenvectors.size(); i++)
			vectors[i] = *(Eigenvectors.data() + i);

		//最小特征值对应的特征向量即为法向量
		int vectors_index=0;
		double min_value = value[0];
		if (value[4] < min_value)
		{
			min_value = value[4];
			vectors_index = 3;
		}
		if (value[8] < min_value)
		{
			min_value = value[8];
			vectors_index = 6;
		}
		//double normal_vector[3];//法向量
		//normal_vector[0] = vectors[vectors_index];
		//normal_vector[1] = vectors[vectors_index+1];
		//normal_vector[2] = vectors[vectors_index+2];
		double a, b, c;
		a = vectors[vectors_index];
		b = vectors[vectors_index + 1];
		c = vectors[vectors_index + 2];

		//將邻域旋转，使法向量与Z轴重合
		//求旋转矩阵
		double rotate_matrix[9];
		rotate_matrix[0] = sqrt(b*b + c*c);
		rotate_matrix[1] = (-a*b)/sqrt(b*b + c*c);
		rotate_matrix[2] = (-a*c) / sqrt(b*b + c*c);
		rotate_matrix[3] = 0;
		rotate_matrix[4] = c/ sqrt(b*b + c*c);
		rotate_matrix[5] = -b / sqrt(b*b + c*c);
		rotate_matrix[6] = a;
		rotate_matrix[7] = b;
		rotate_matrix[8] = c;

		//旋转K邻域里的所有点
		set<C3DPoint>s_sorted_plane_vertices_index;//保存转化后的点的xy和点的索引
	
		for (int m = 0; m < temp_vertices.size(); m++)
		{
			C3DPoint temp;
			double x, y, z;
			x = temp_vertices[m].x;
			y = temp_vertices[m].y;
			z = temp_vertices[m].z;
			temp.x = rotate_matrix[0] * x + rotate_matrix[1] * y + rotate_matrix[2] * z;
			temp.y = rotate_matrix[3] * x + rotate_matrix[4] * y + rotate_matrix[5] * z;
			if (m == 0)
			{
				temp.z = i;
			}
			else
			{
				temp.z = v_Kneighbourhood[i][m - 1];
			}
			s_sorted_plane_vertices_index.insert(temp);
		}

		vector<C3DPoint> v_sorted_plane_vertices_index;
		v_sorted_plane_vertices_index.resize(s_sorted_plane_vertices_index.size());
		copy(s_sorted_plane_vertices_index.begin(), s_sorted_plane_vertices_index.end(), v_sorted_plane_vertices_index.begin());
		vector<C2DPoint> v_plane_vertices;
		v_plane_vertices.resize(v_sorted_plane_vertices_index.size());

	


		vector<int>vertice_index;//保存v_plane_vertices中各点在vertices中的下标
		vertice_index.resize(v_sorted_plane_vertices_index.size());
		for (int n = 0; n < v_sorted_plane_vertices_index.size(); n++)
		{
			v_plane_vertices[n].x = v_sorted_plane_vertices_index[n].x;
			v_plane_vertices[n].y = v_sorted_plane_vertices_index[n].y;
			vertice_index [n]= v_sorted_plane_vertices_index[n].z;
		}

		//在平面上进行三角剖分
		CDelaunay delaunay(v_plane_vertices);
		vector<C3DPoint>v_plane_triangles(delaunay.mesh_generation());
		vector<C3DPoint>v_triangles_correct_index;//将v_plane_triangles中各点的下标转化为vertices中的
		v_triangles_correct_index.resize(v_plane_triangles.size());
		for (int t = 0; t < v_plane_triangles.size(); t++)
		{
			v_triangles_correct_index[t].x = vertice_index[v_plane_triangles[t].x];
			v_triangles_correct_index[t].y = vertice_index[v_plane_triangles[t].y];
			v_triangles_correct_index[t].z = vertice_index[v_plane_triangles[t].z];
		}

		//把含有中心点的三角形提取出来

		for (int t = 0; t < v_triangles_correct_index.size(); t++)
		{
			if (v_triangles_correct_index[t].x == i || v_triangles_correct_index[t].y == i || v_triangles_correct_index[t].z == i)
			{
				int a = v_triangles_correct_index[t].x;
				int b = v_triangles_correct_index[t].y;
				int c = v_triangles_correct_index[t].z;
				//将abc由小到大排序
				if (a>b)
				{
					int temp = a;
					a = b;
					b = temp;
				}
				if (a>c)
				{
					int temp = a;
					a = c;
					c = temp;
				}
				if (b>c)
				{
					int temp = b;
					b = c;
					c = temp;
				}
				
				C3DPoint temp;
				temp.x = a;
				temp.y = b;
				temp.z = c;
				s_triangles.insert(temp);
			}
		}
			

	}
	v_triangles.resize(s_triangles.size());
	copy(s_triangles.begin(), s_triangles.end(), v_triangles.begin());

	return v_triangles;
	 
}