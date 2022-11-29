#pragma once

#ifndef PREDEFINE_H
#define PREDEFINE_H
#define varName(x) #x	// 宏，获取变量x的名字

using namespace std;
#include<string>
#include<vector>
#include <iostream>
#include <math.h>
#include <limits.h>
#include <cstdlib>
#include <iomanip>
#include <ilcplex/ilocplex.h>
#include <ilconcert/ilomodel.h>
#include <ctime>
#include <algorithm>
#include <random>

const double PI = 3.1415926535;
const double f = 2e9;		// 载波评率 4 Ghz
const double c = 299792458;			// 光速 m/s
const double eta_LOS = 1;	// shadow fading loss of LoS, 1 dB
const double N0 = -174;		// 背景噪声 −174 dBm/Hz    总功率 = -174 + 10log(1Mhz) = -114
static double N_I = 36;				// 干扰噪声，36dB  之后可能会更新
const double SINR_min = -7;	// 最小信噪比 10 dBm

union xu {
	double** d;
	int** i;
};
union yu {
	double* d;
	int* i;
};


class Point {
public:
	double X = 0;	// X坐标
	double Y = 0;	// Y坐标

	Point(double x, double y) { X = x; Y = y; }
	Point() {}
	~Point() {}
	double cal_distance(const Point& p);		//计算本对象与点p的距离 
	void set_location(double lx, double ly) { X = lx; Y = ly; }
	const int get_index();
	bool operator<(const Point& p) const;
};

class User : public Point
{
public:
	int BR = 0;		// 带宽需求
	int ID = -1;	// 用户id


	User(double lx, double ly, int id, int br);

	void print_user();
};

class Server : public Point
{
public:
	int BW = 0;
	int ID = -1;
	bool isSelected = false;
	double p = 33;			// 无人机功率 dBm		备选33 dBm
	double h = 300;			// 无人机飞行高度		300m


	Server(double lx, double ly, int id, int bw);

	void select() { isSelected = true; }	// 将该服务器标记为被选择
	void cancel() { isSelected = false; }	// 将该服务器标记为不被选择
	void print_server();
};

class Cover
{
public:
	int m = 0;
	int n = 0;
	vector<User> u;
	vector<Server> a;

	double** d = nullptr;	// 记录任意用户与服务器之间的距离
	double** L = nullptr;	// path loss
	double** SINR = nullptr;	// SINR
	double** Gp = nullptr;		// G*p = p - L


	/// <summary>
	/// 函数
	/// </summary>
	/// <param name="fname"></param>
	void initial(string fname);		// 根据fname文件中的数据初始化实例
	void read_file(string fname);	// 读取fname文件中的数据

	void cal_d();
	void cal_L();
	void cal_SINR();
	void cal_min_r();

	void print_all_usersServers();	// 输出所有用户与服务器
	template <class T>
	void print_array(T arr, string name);
	void print_all();
	template<class TT, class T>
	void print_cplex(TT x, T y);

	void LP();	// 线性规划
	void IP();	// 整数规划
};

class Result
{
public:
	int type = 0;		// x,y的类型 0为double，1为int
	xu x;
	yu y;
	int m;
	int n;

	Result(xu xx, yu yy, int t) { x = xx; y = yy; type = t; }
	Result() { ; }

	// 写将x，y结果中被选择的UAV，那个无人机覆盖那个用户输出到文件
};




#endif // PREDEFINE_H