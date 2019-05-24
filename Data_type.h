#ifndef DATA_H
#define DATA_H
#include <vector>
#include <set>
#include"Data_type.h"
#include<algorithm>

using namespace std;
#pragma once

struct C3DPoint
{
	double x, y, z;
	bool friend operator<(const C3DPoint &C1, const C3DPoint &C2)
	{
		if (C1.x<C2.x) return true;
		if (C1.x == C2.x && C1.y<C2.y)return true;
		if (C1.x == C2.x && C1.y == C2.y&C1.z<C2.z)return true;
		return false;
	}

};
struct C2DPoint
{
	double x, y;
	bool friend operator<(const C2DPoint &C1, const C2DPoint &C2)
	{
		if (C1.x<C2.x) return true;
		if (C1.x == C2.x && C1.y<C2.y)return true;
		return false;
	}
};
struct StlPoint_AND_OneFace    //����һ�������洢��.stl����ȡ�ĵ�����꼰�����ڵ�һ��������Ƭ������
{


	double x, y, z, f;
	bool friend operator<(const StlPoint_AND_OneFace  &S1, const StlPoint_AND_OneFace  &S2)
	{
		if (S1.x<S2.x) return true;
		if (S1.x == S2.x && S1.y<S2.y)return true;
		if (S1.x == S2.x && S1.y == S2.y&&S1.z<S2.z)return true;
		if (S1.x == S2.x && S1.y == S2.y&&S1.z == S2.z && S1.f<S2.f)return true;
		return false;
	}
};

struct StlPoint_AND_FaceNormals   //����һ�������洢��.stl����ȡ�Ĳ��ظ�������꼰���㷨����������
{


	double x, y, z, fx, fy, fz;
	bool friend operator<(const StlPoint_AND_FaceNormals  &S1, const StlPoint_AND_FaceNormals  &S2)
	{
		if (S1.x<S2.x) return true;
		if (S1.x == S2.x && S1.y<S2.y)return true;
		if (S1.x == S2.x && S1.y == S2.y&&S1.z<S2.z)return true;
		return false;
	}

};
struct tetrahedron
{
	double i, j, k, w, xc, yc, zc, r;//�������ĸ��������������������ģ��뾶
};
#endif