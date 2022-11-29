#include "predefine.h"
#include <iostream>
#include <cmath>

double Point::cal_distance(const Point& p)
{
	const double X1 = this->X;
	const double Y1 = this->Y;

	const double X2 = p.X;
	const double Y2 = p.Y;

	double dis = sqrt(pow(X1 - X2, 2) + pow(Y1 - Y2, 2));
	return dis;
}

User::User(double lx, double ly, int id, int br)
{
	X = lx;
	Y = ly;
	ID = id;
	BR = br;
}

void User::print_user()
{	
	cout << "u" << ID << ":(" << X << ", " << Y << "), " << BR;
}

Server::Server(double lx, double ly, int id, int bw)
{
	X = lx;
	Y = ly;
	ID = id;
	BW = bw;
}

void Server::print_server()
{
	cout << "a" << ID << ":(" << X << ", " << Y << "), " << BW;
}

void Cover::cal_d()
{
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++)
		{
			d[i][j] = sqrt(pow(a[i].X - u[j].X, 2) + pow(a[i].Y - u[j].Y, 2) + 300 * 300);
		}
	}
}

void Cover::cal_L()
{
	double eta = 20 * log10((4 * PI * f) / c) + eta_LOS;
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++)
			L[i][j] = 20 * log10(d[i][j]) + eta;
	}
}

void Cover::cal_SINR()
{
	// 
	double N = N0 + 10 * log10(a[0].BW * 10e6) + N_I;
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++)
			SINR[i][j] = a[i].p - L[i][j] - N;
	}
}

void Cover::cal_min_r()
{
	double p = 33;
	double BW = 20 * 10e6;
	double eta = 20 * log10((4 * PI * f) / c) + eta_LOS;
	double N = N0 + 10 * log10(BW) + 36;
	double zhishu = (p - SINR_min - eta - N) / 20;
	double r = pow(10, zhishu);

	cout << "N0 = " << N0 << '\n';
	cout << "eta = " << eta << '\n';
	cout << "zhishu = " << zhishu << '\n';
	cout << "r = " << r << '\n';
}

void Cover::initial(string fname)
{
	read_file(fname);

	d = new double* [m];
	L = new double* [m];
	SINR = new double* [m];
	Gp = new double* [m];

	for (int i = 0; i < m; i++)
	{
		d[i] = new double [n];
		L[i] = new double [n];
		SINR[i] = new double [n];
		Gp[i] = new double [n];

		for (int j = 0; j < n; j++)
		{
			d[i][j] = 0;
			L[i][j] = 0;
			SINR[i][j] = 0;
			Gp[i][j] = 0;
		}
	}

	cal_d();
	cal_L();
	cal_SINR();
}

void Cover::read_file(string fname)
{
	string data;
	ifstream infpoint;
	infpoint.open(fname);
	// 读取m, n
	getline(infpoint, data);
	istringstream istr1(data);
	istr1 >> m;
	istr1 >> n;

	for (int i = 0; i < m; i++) {
		double x;
		double y;
		int BW;
		getline(infpoint, data);
		istringstream istr2(data);
		istr2 >> x >> y >> BW;
		Server s(x, y, i, BW);
		a.push_back(s);
	}

	for (int j = 0; j < n; j++) {
		double x;
		double y;
		int BR;
		getline(infpoint, data);
		istringstream istr3(data);
		istr3 >> x >> y >> BR;
		User user(x, y, j, BR);
		u.push_back(user);
	}
}

void Cover::print_all_usersServers()
{
	cout << "UAV:\n";
	for (int i = 0; i < m; i++)
	{
		a[i].print_server();
		cout << '\t';
		if ((i + 1) % 5 == 0)
			cout << '\n';
	}
	cout << '\n';
	cout << "Users:\n";
	for (int i = 0; i < n; i++)
	{
		u[i].print_user();
		cout << '\t';
		if ((i + 1) % 5 == 0)
			cout << '\n';
	}
}

template <class T>
void Cover::print_array(T arr, string name)
{
	cout << name << ':' << endl;
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++)
			cout << arr[i][j] << " ";
		cout << '\n';
	}
	cout << '\n';
}

template<class TT, class T>
void Cover::print_cplex(TT x, T y)
{
	cout << "y:\n";
	for (int i = 0; i < m; i++)
	{
		cout << y[i] << ' ';
		if ((i + 1) % 20 == 0)
			cout << '\n';
	}
	cout << "\nx:\n";
	for (int i = 0; i < m; i++)
	{
		cout << "a" << i << ":\t";
		for (int j = 0; j < n; j++)
		{
			cout << x[i][j] << ' ';
			if ((j + 1) % 20 == 0)
				cout << "\n\t";
		}
	}
}

void Cover::print_all()
{
	print_all_usersServers();
	print_array(d, "d");
	print_array(L, "L");
	print_array(SINR, "SINR");
}

typedef IloArray <IloNumVarArray> IloNumVarArray2;
typedef IloArray <IloNumVarArray2> IloNumVarArray3;

ILOSTLBEGIN
void Cover::LP()
{
	clock_t start_, end_;
	start_ = clock();

	double** x0 = new double* [m];
	for (int i = 0; i < m; i++)
	{
		x0[i] = new double[n];
		for (int j = 0; j < n; j++)
			x0[i][j] = 0;
	}

	double* y0 = new double[n];
	
	IloEnv env_LP;
	double solution_value = 0;

	try {
		IloModel model(env_LP);

		IloNumVarArray2 x(env_LP, m);
		IloNumVarArray y(env_LP, m);
		for (IloInt i = 0; i < m; i++)
		{
			y[i] = IloNumVar(env_LP, 0, 1, ILOFLOAT);
			x[i] = IloNumVarArray(env_LP, n);
			for (IloInt j = 0; j < n; j++)
			{
				x[i][j] = IloNumVar(env_LP, 0, IloInfinity, ILOFLOAT);
				model.add(x[i][j]);
			}
		}

		// 约束1.1
		// 被u被a覆盖，a必须被选中
		// xij <= yi
		int ud = 0;
		IloExpr cons_x_le_y(env_LP);
		for (IloInt i = 0; i < m; i++)
		{
			for (IloInt j = 0; j < n; j++)
			{
				cons_x_le_y.clear();
				cons_x_le_y = x[i][j] - y[i];
				model.add(cons_x_le_y <= 0);
			}
		}

		// 约束1.2
		// 容量限制约束
		IloExpr capacity_cons(env_LP);
		IloExpr sumBR(env_LP);
		for (IloInt i = 0; i < m; i++)
		{
			capacity_cons.clear();
			for (IloInt j = 0; j < n; j++)
			{
				capacity_cons += x[i][j] * u[j].BR;
			}
			capacity_cons -= y[i] * a[i].BW;
			model.add(capacity_cons <= 0);
		}

		// 约束1.3
		// 每个用户都要被覆盖
		IloExpr cover_all_user_cons(env_LP);
		for (IloInt j = 0; j < n; j++)
		{
			cover_all_user_cons.clear();
			for (IloInt i = 0; i < m; i++)
			{
				cover_all_user_cons += x[i][j];
			}
			model.add(cover_all_user_cons == 1);
		}

		// 约束1.4
		// 通信质量约束
		IloExpr SINR_cons(env_LP);
		for (IloInt i = 0; i < m; i++)
			for (IloInt j = 0; j < n; j++)
				if (SINR[i][j] < SINR_min)
				{
					SINR_cons.clear();
					SINR_cons = x[i][j];
					model.add(SINR_cons == 0);
				}


		// 目标函数
		IloExpr obj(env_LP);
		for (IloInt i = 0; i < m; i++)
			obj += y[i];


		model.add(IloMinimize(env_LP, obj));

		// 创建求解对象
		IloCplex cplex(model);
		if (!cplex.solve()) {
			env_LP.error() << "Failed to optimize LP." << endl;
			throw(-1);
		}

		solution_value = cplex.getObjValue();
		printf("Solution value = %lf\n", solution_value);

		IloNumArray yy(env_LP);
		cplex.getValues(yy, y);
		for (IloInt i = 0; i < m; i++)
		{
			y0[i] = yy[i];

			IloNumArray xx(env_LP);
			cplex.getValues(xx, x[i]);
			for (int j = 0; j < n; j++) {
				x0[i][j] = xx[j];
			}
		}
		
		//IloNumArray vals(env);
		//env.out() << "Solution status = " << cplex.getStatus() << endl;
		//env.out() << "Solution value = " << cplex.getObjValue() << endl;
		//cplex.getValues(vals, x);
		//env.out() << "Value = " << vals << endl;
	}
	catch (IloException& e) { cerr << "Concert exception caught:" << e << endl; }
	catch (...) { cerr << "Unknuwn exception caught" << endl; }
	env_LP.end();


	print_cplex(x0, y0);
	xu x; yu y;
	x.d = x0;
	y.d = y0;
	Result r = Result(x, y, 0);
	r.m = m;
	r.n = n;
}

ILOSTLBEGIN
void Cover::IP()
{
	
	clock_t start_, end_;
	start_ = clock();

	int** x0 = new int* [m];
	for (int i = 0; i < m; i++)
	{
		x0[i] = new int[n];
		for (int j = 0; j < n; j++)
			x0[i][j] = 0;
	}

	int* y0 = new int[n];

	IloEnv env_IP;
	double solution_value = 0;

	try {
		IloModel model(env_IP);

		IloNumVarArray2 x(env_IP, m);
		IloNumVarArray y(env_IP, m);
		for (IloInt i = 0; i < m; i++)
		{
			y[i] = IloNumVar(env_IP, 0, 1, ILOINT);
			x[i] = IloNumVarArray(env_IP, n);
			for (IloInt j = 0; j < n; j++)
			{
				x[i][j] = IloNumVar(env_IP, 0, 1, ILOINT);
				model.add(x[i][j]);
			}
		}

		// 约束1.1
		// 被u被a覆盖，a必须被选中
		// xij <= yi
		int ud = 0;
		IloExpr cons_x_le_y(env_IP);
		for (IloInt i = 0; i < m; i++)
		{
			for (IloInt j = 0; j < n; j++)
			{
				cons_x_le_y.clear();
				cons_x_le_y = x[i][j] - y[i];
				model.add(cons_x_le_y <= 0);
			}
		}

		// 约束1.2
		// 容量限制约束
		IloExpr capacity_cons(env_IP);
		IloExpr sumBR(env_IP);
		for (IloInt i = 0; i < m; i++)
		{
			capacity_cons.clear();
			for (IloInt j = 0; j < n; j++)
			{
				capacity_cons += x[i][j] * u[j].BR;
			}
			capacity_cons -= y[i] * a[i].BW;
			model.add(capacity_cons <= 0);
		}

		// 约束1.3
		// 每个用户都要被覆盖
		IloExpr cover_all_user_cons(env_IP);
		for (IloInt j = 0; j < n; j++)
		{
			cover_all_user_cons.clear();
			for (IloInt i = 0; i < m; i++)
			{
				cover_all_user_cons += x[i][j];
			}
			model.add(cover_all_user_cons == 1);
		}

		// 约束1.4
		// 通信质量约束
		IloExpr SINR_cons(env_IP);
		for (IloInt i = 0; i < m; i++)
			for (IloInt j = 0; j < n; j++)
				if (SINR[i][j] < SINR_min)
				{
					SINR_cons.clear();
					SINR_cons = x[i][j];
					model.add(SINR_cons == 0);
				}

		// 目标函数
		IloExpr obj(env_IP);
		for (IloInt i = 0; i < m; i++)
			obj += y[i];

		model.add(IloMinimize(env_IP, obj));

		// 创建求解对象
		IloCplex cplex(model);
		if (!cplex.solve()) {
			env_IP.error() << "Failed to optimize LP." << endl;
			throw(-1);
		}

		solution_value = cplex.getObjValue();
		printf("Solution value = %lf\n", solution_value);

		IloNumArray yy(env_IP);
		cplex.getValues(yy, y);
		for (IloInt i = 0; i < m; i++)
		{
			y0[i] = yy[i];

			IloNumArray xx(env_IP);
			cplex.getValues(xx, x[i]);
			for (int j = 0; j < n; j++) {
				x0[i][j] = xx[j];
			}
		}

		//IloNumArray vals(env);
		//env.out() << "Solution status = " << cplex.getStatus() << endl;
		//env.out() << "Solution value = " << cplex.getObjValue() << endl;
		//cplex.getValues(vals, x);
		//env.out() << "Value = " << vals << endl;
	}
	catch (IloException& e) { cerr << "Concert exception caught:" << e << endl; }
	catch (...) { cerr << "Unknuwn exception caught" << endl; }
	env_IP.end();

	print_cplex(x0, y0);
}


