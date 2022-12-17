#pragma once

#ifndef PREDEFINE_H
#define PREDEFINE_H
#define varName(x) #x	// 宏，获取变量x的名字

using namespace std;
#include<string>
#include<vector>
#include<map>
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
static double N_I = -52;				// 干扰噪声，36dB  之后可能会更新
const double SINR_min = 0;	// 最小信噪比 -7 dBm
const double power = 40;	// 无人机信号发射功率 33dBm
const double alpha = 0.5;
const double h = 100;		// 无人机飞行高度
const int E = 1e5;
// const double BW = 20 * 10e6;	// 无人机的带宽容量

union xu {
	double** dh;
	int** i;
};
union yu {
	double* dh;
	int* i;
};

typedef vector<int> Cluster;

class Point;
class User;
class Server;
class Cover;
class Result;

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
	double p = power;			// 无人机功率 dBm		备选33 dBm
	double h = 300;				// 无人机飞行高度		300m
	vector<int> served_user;

	Server(double lx, double ly, int id, int bw);

	void select() { isSelected = true; }	// 将该服务器标记为被选择
	void cancel() { isSelected = false; }	// 将该服务器标记为不被选择
	void print_server();
};



class Result
{
public:
	double ** x = nullptr;
	double* y = nullptr;
	int m = 0;
	int n = 0;
	int chosen_m = 0;

	Result() {}
	template<class TT, class T>
	Result(TT xx, T yy, int mm, int nn) {
		m = mm; n = nn;
		for (int i = 0; i < m; i++)
		{
			y[i] = yy[i];
			for (int j = 0; j < n; j++)
				x[i][j] = xx[i][j];
		}
	}
	Result(int mm, int nn) {
		m = mm;
		n = nn;
		x = new double* [m];
		y = new double[m];
		for (int i = 0; i < m; i++)
		{
			x[i] = new double[n];
			y[i] = 0;
			for (int j = 0; j < n; j++)
				x[i][j] = 0;
		}

	}

	template<class TT, class T>
	void init_(TT xx, T yy);
	// 写将x，y结果中被选择的UAV，那个无人机覆盖那个用户输出到文件
	void write_result_file(string fname);
	void write_cplex_file(string fname);

};

class Cover
{
public:
	int m = 0;
	int n = 0;
	vector<User> U;
	vector<int> all_u;			// 一个0~n-1的数组，记录所有用户id
	vector<Server> A;
	Result result;
	double ep_p = 1.0/6;		// p=33dBm,p'=34dBm
	double ep_SINR = -1;	// SINR_min=-7dB, SINR_min'=-8dB
	double ep = 0;
	double r = 0;				// 无人机在地面的覆盖半径




	double** NI = nullptr;
	double** dis_h = nullptr;	// 记录任意用户与服务器之间的距离，考虑高度h
	double** dis = nullptr;	// 记录任意用户与服务器之间的水平距离
	double** L = nullptr;	// path loss
	double** SINR = nullptr;	// SINR
	double** Gp = nullptr;		// G*p = p - L



	/// <summary>
	/// 函数
	/// </summary>
	/// <param name="fname"></param>
	~Cover() { }
	void initial(string fname);		// 根据fname文件中的数据初始化实例
	void read_file(string fname);	// 读取fname文件中的数据
	template<class TT, class T>
	void read_cplex_file(string fame, TT x, T y);

	void cal_d();			// 计算所有用户与无人机之间的距离
	void cal_L();			// 计算路径损耗
	void cal_SINR();		// 计算信噪比
	void cal_min_r();		// 计算无人机最小半径
	void cal_ep();			// 根据epsilon_p和epsilon_SINR计算epsilon
	template <class T>
	void cal_NI(T y);		// 根据被选中的无人机，计算所有的

	void print_all_usersServers();	// 输出所有用户与服务器
	template <class T>
	void print_array(T arr, string name);
	void print_all();
	template<class TT, class T>
	void print_cplex(TT x, T y);
	template<class TT, class T>
	void print_RB(TT x, T y, int is_y = 0);			// 打印当前状态下，y值不为零无人机对应的剩余容量
	void print_CC(map<int, Cluster>& CC);

	void LP(double** x0, double* y0);	// 线性规划
	void IP(int** x0, int* y0);	// 整数规划
	
	void GBTSR(int isLP = 1);	// Grid-based three-step rounding approximation algorithm
	void DSIS(double** x1, double* y1);						// Determining Superior and Inferior Servers (DSIS)
	void mergeII(double** x, double* y, double** x1, double* y1);
	map<int, Cluster> COS(double** x1, double* y1, double** x2, double* y2);						// Clusting of Servers(CoS)
	void SFS(double** x1, double** x2, double* y2, double** x3, double* y3, map<int, Cluster>& CC);												// Selecting the Final Servers(SFS)


	// template<class TT>
	// double sum(TT x, int len);				// 计算长度为len的x的和
	void construct_I_j(vector<int>& I_j,  double sum_I, double* y);
	void reroute(vector<int>& obj, vector<int>& src, int aim, double** x, double flow = 1); 	// 重引流
	template<class T>
	double get_sumBR(int i, T x);		// 计算无人机ai的已用容量
	template<class T>
	double get_RB(int i, T x);		// 计算无人机ai的剩余容量
	double get_BRmin(vector<int>& Ai);
	int get_max_xBR(int t, vector<int>& At, double RBt, double** x);		// 在无人机t当前服务的用户集合At中，找到最适合那个用户。这个用户将被重引流
	bool is_contain(vector<int>& At, vector<int>& Ai, int i);	// 判断At是否是Ai的子集，是返回true

	map<int, vector<int>> construct_GG(vector<int>& I, Point&cp, double L, double cl);		// 构造
	void merge_GG(map<int, vector<int>>& GG, double cl);

	map<int, vector<int>> get_ISorder(vector<int>& SS, vector<int>& II);
	map<int, vector<int>> get_SIorder(vector<int>& SS, vector<int>& II);
};






#endif // PREDEFINE_H