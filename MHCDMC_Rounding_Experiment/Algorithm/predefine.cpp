#include "predefine.h"
#include <iostream>
#include <cmath>

// 计算当前对象与参数p之间的距离
double Point::cal_distance(const Point& p)
{
	const double X1 = this->X;
	const double Y1 = this->Y;

	const double X2 = p.X;
	const double Y2 = p.Y;

	double d = sqrt(pow(X1 - X2, 2) + pow(Y1 - Y2, 2));
	return d;
}

// User类的构造函数，初始化坐标(x, y)，ID以及带宽需求BR
User::User(double lx, double ly, int id, int br)
{
	X = lx;
	Y = ly;
	ID = id;
	BR = br;
}

// 输出当前用户对象的信息
void User::print_user()
{	
	cout << "u" << ID << ":(" << X << ", " << Y << "), " << BR;
}

// Server类的构造函数，初始化坐标(x, y)，ID以及带宽容量BW
Server::Server(double lx, double ly, int id, int bw)
{
	X = lx;
	Y = ly;
	ID = id;
	BW = bw;
}

// 输出当前服务器对象的信息
void Server::print_server()
{
	cout << "a" << ID << ":(" << X << ", " << Y << "), " << BW;
}

// 计算所有用户与服务器之间的距离 dis_h, dis
void Cover::cal_d()
{
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++)
		{
			dis_h[i][j] = sqrt(pow(A[i].X - U[j].X, 2) + pow(A[i].Y - U[j].Y, 2) + h * h);
			dis[i][j] = sqrt(pow(A[i].X - U[j].X, 2) + pow(A[i].Y - U[j].Y, 2));
		}
	}
}

// 计算所有用户与服务器之间的路径损益 L
void Cover::cal_L()
{
	double eta = 20 * log10((4 * PI * f) / c) + eta_LOS;
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++)
			L[i][j] = 20 * log10(dis_h[i][j]) + eta;
	}
}

// 计算所有用户与服务器之间的信噪比 SINR
void Cover::cal_SINR()
{
	// 
	// double N = N0 + 10 * log10(A[0].BW * 10e6) + N_I;
	
	// cout << "N =  " << N << '\n';
	for (int j = 0; j < n; j++)
	{
		for (int i = 0; i < m; i++)
		{
			SINR[i][j] = 0;
			SINR[i][j] = A[i].p - L[i][j] - NI[i][j];
		}
	}
	
}

// 根据已知数据，计算无人机的最打覆盖半径
void Cover::cal_min_r()
{
	//double BW = A[0].BW * 10e6;
	double BW = 5 * E;
	double eta = 20 * log10((4 * PI * f) / c) + eta_LOS;
	double N = N0 + 10 * log10(BW / 10);
	double NN = 10 * log10(pow(10, N / 10) + pow(10, N_I / 10));
	double zhishu = (power - SINR_min - eta - NN) / 20;
	double R = pow(10, zhishu);
	r = sqrt(pow(R, 2) - pow(h, 2));

	cout << "N0 = " << N0 << '\n';
	cout << "eta = " << eta << '\n';
	cout << "zhishu = " << zhishu << '\n';
	cout << "r = " << r << '\n';
}

// 根据功率p与SINR_min的扩张系数，计算无人机半径的扩张系数 ep
void Cover::cal_ep()
{
	double zhishu = (power * ep_p + SINR_min - ep_SINR) / 20;

	ep = pow(10, zhishu) - 1;
	cout << "ep = " << ep << '\n';
}

// 根据文件fname来初始化一个Cover对象
void Cover::initial(string fname)
{
	read_file(fname);

	dis_h = new double* [m];
	dis = new double* [m];
	NI = new double* [m];
	L = new double* [m];
	SINR = new double* [m];
	Gp = new double* [m];

	for (int i = 0; i < m; i++)
	{
		dis_h[i] = new double [n];
		dis[i] = new double [n];
		NI[i] = new double [n];
		L[i] = new double [n];
		SINR[i] = new double [n];
		Gp[i] = new double [n];

		for (int j = 0; j < n; j++)
		{
			NI[i][j] = N_I;
			dis_h[i][j] = 0;
			dis[i][j] = 0;
			L[i][j] = 0;
			SINR[i][j] = 0;
			Gp[i][j] = 0;
		}
	}

	result = Result(m, n);
	for (int j = 0; j < n; j++)
		all_u.push_back(j);

	cout << "initial:\t";
	cout << "m = " << m << ", n = " << n << endl;
	cal_d();
	cal_L();
	cal_SINR();
	cal_min_r();
	cal_ep();

	cout << "initial:\t";
	cout << "m = " << m << ", n = " << n << endl;
}

// 从fname文件中读取m，n以及点坐标等信息
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
		A.push_back(s);
	}

	for (int j = 0; j < n; j++) {
		double x;
		double y;
		int BR;
		getline(infpoint, data);
		istringstream istr3(data);
		istr3 >> x >> y >> BR;
		User user(x, y, j, BR / 2.5);
		U.push_back(user);
	}
	
}

template<class TT, class T>
void Cover::read_cplex_file(string fname, TT x, T y)
{
	int chosen_m;
	cout << fname << '\n';
	string data;
	ifstream infpoint;
	infpoint.open(fname);
	// 读取m, n
	getline(infpoint, data);
	istringstream istr1(data);
	istr1 >> chosen_m;
	istr1 >> n;

	for (int ii = 0; ii < chosen_m; ii++) {
		int i = 0;
		getline(infpoint, data);
		istringstream istr2(data);
		cout << data << '\n';
		istr2 >> i;
		istr2 >> y[i];

		getline(infpoint, data);
		istringstream istr3(data);
		for (int j = 0; j < n; j++)
			istr3 >> x[i][j];
	}
	
	infpoint.close();
}

// 根据数组y来计算每个用户的干扰噪声
template<class T>
void Cover::cal_NI(T y)
{

	
	// cout << "N = " << N << '\n';
	vector<int> S;				// 被选中的服务器
	int* uBs = new int [n]; 		// 每个用户被离他最近的服务器服务

	cout << '\n';
	for (int i = 0; i < m; i++)
		if (y[i] > 0)
		{
			S.push_back(i);
			//cout << i << ' ';
		}

	//cout << '\n';
	

	for (int j = 0; j < n; j++)
	{
		double N = N0 + 10 * log10(U[j].BR * E);
		//cout << "N = " << N << '\n';
		for (int i = 0; i < m; i++)
		{
			double num = pow(10, N / 10.0);
			// cout << "num = " << num << '\n';
			//cout << "i  = " << i << ", j = " << j << "\n";
			for (auto a : S)
				if (a != i)
				{
					//cout << "\ta" << a << ": " << (A[a].p - L[a][j]) / 10.0 << '\n';
					//cout << "\tpow a" << a << ": " << pow(10, (A[a].p - L[a][j]) / 10.0) << '\n';
					num += pow(10, (A[a].p - L[a][j]) / 10.0);
				}

			//cout << "num = " << num << '\n';
			NI[i][j] = 10 * log10(num);
			//cout << "N[i][j] = " << NI[i][j] << '\n';
		}
	}
}



// 输出所有用户与服务器信息
void Cover::print_all_usersServers()
{
	cout << "UAV:\n";
	for (int i = 0; i < m; i++)
	{
		A[i].print_server();
		cout << '\t';
		if ((i + 1) % 5 == 0)
			cout << '\n';
	}
	cout << '\n';
	cout << "Users:\n";
	for (int i = 0; i < n; i++)
	{
		U[i].print_user();
		cout << '\t';
		if ((i + 1) % 5 == 0)
			cout << '\n';
	}
}

// 输出名为name的二维数组arr
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

// 将算法的x，y结果输出
template<class TT, class T>
void Cover::print_cplex(TT x, T y)
{
	double sum_y = 0;
	cout << "y:\n";
	for (int i = 0; i < m; i++)
	{
		if (y[i] != 0)
		{
			sum_y += y[i];
			cout << "y" << i << '=' << y[i] << '\t';
		}
	}
	cout << '\n';
	cout << "sum_y = " << sum_y << '\n';
	cout << "\nx:\n";
	for (int i = 0; i < m; i++)
	{
		if (y[i] != 0)
		{
			cout << "a" << i << ":\t";
			for (int j = 0; j < n; j++)
			{
				cout << setprecision(4) << x[i][j] << '\t';
				if ((j + 1) % 20 == 0)
					cout << "\n";
			}
			cout << '\n';
		}
	}
}

// 根据当前x,y状态，输出被选中服务器的剩余容量
template<class TT, class T>
void Cover::print_RB(TT x, T y, int is_y)
{
	// is_y			服务器容量是否乘上y这个系数		1表示乘，0表示不乘
	map<int, double> RBs;
	for (int i = 0; i < m; i++)
	{
		if (y[i] > 0)
		{
			double yy = y[i];
			if (is_y == 0)
				yy = 1;
			double RB = 0;
			for (int j = 0; j < n; j++)
				RB += x[i][j] * U[j].BR;
			RB = A[i].BW * yy - RB;
			RBs[i] = RB;
		}
	}
	cout << "RBs:\n";
	map<int, double>::iterator it = RBs.begin();
	while (it != RBs.end()) {
		cout << "a" << it->first << ": " << it->second << '\t';
		it++;
	}
	cout << '\n';
}

// 根据当前x状态，输出当前服务器i的已用容量
template<class T>
double Cover::get_sumBR(int i, T x)
{
	double sum = 0;
	for (int j = 0; j < n; j++)
		sum += x[i][j] * U[j].BR;
	return sum;
}

// 根据当前x状态，输出服务器i的可用容量
template<class T>
double Cover::get_RB(int i, T x)
{
	double BW = A[i].BW;
	double sum = get_sumBR(i, x);
	return BW - sum >= 0 ? BW - sum : 0;
}

//template<class TT>
//double Cover::sum(TT x, int len)
//{
//	double sum = 0;
//	for (int i = 0; i < len; i++)
//		sum += x[i];
//	return sum
//}

// 输出当前对象对应的所有数据，包括 用户与服务器信息， 距离dis_h，路径损益L，信噪比SINR 
void Cover::print_all()
{
	print_all_usersServers();
	print_array(dis_h, "d");
	print_array(L, "L");
	print_array(SINR, "SINR");
}

// 输出Cluster的集合 CC
void Cover::print_CC(map<int, Cluster>& CC)
{
	map<int, Cluster>::iterator it = CC.begin();
	cout << "CC:\n";
	while (it != CC.end())
	{
		cout << "\tC" << it->first << ": ";
		for (int j : it->second)
			cout << j << ' ';
		cout << '\n';
		it++;
	}
}


typedef IloArray <IloNumVarArray> IloNumVarArray2;
typedef IloArray <IloNumVarArray2> IloNumVarArray3;
ILOSTLBEGIN
// 线性规划
void Cover::LP(double** x0, double* y0)
{
	clock_t start_, end_;
	start_ = clock();

	
	// IloEnv env_LP;
	double solution_value = 0;
	IloEnv env_LP;
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
				//cons_x_le_y = x[i][j] * U[j].BR - y[i] * A[i].BW;
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
				capacity_cons += x[i][j] * U[j].BR;
			}
			capacity_cons -= y[i] * A[i].BW;
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
		for (IloInt i = 0; i < m; i++)
			for (IloInt j = 0; j < n; j++)
				if (SINR[i][j] < SINR_min)
				{
					IloRange SINR_con(env_LP, 0, x[i][j], 0);
					model.add(SINR_con);
				}


		// 目标函数
		IloExpr obj(env_LP);
		for (IloInt i = 0; i < m; i++)
			obj += y[i];


		model.add(IloMinimize(env_LP, obj));

		// 创建求解对象
		IloCplex cplex(model);
		cplex.setParam(IloCplex::NodeAlg, IloCplex::Barrier);
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
	


	
}

ILOSTLBEGIN
//整数规划
void Cover::IP(int** x0, int* y0)
{
	
	clock_t start_, end_;
	start_ = clock();


	// IloEnv env_IP;
	double solution_value = 0;
	IloEnv env_IP;
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
				// cons_x_le_y = x[i][j] * U[j].BR - y[i] * A[i].BW;
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
				capacity_cons += x[i][j] * U[j].BR;
			}
			capacity_cons -= y[i] * A[i].BW;
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
		obj.clear();
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

}

// 我们的算法入口
void Cover::GBTSR(int isLP)
{
	double** x;
	double** x1;
	double** x2;
	double** x3;
	int** x4;
	double* y;
	double* y1;
	double* y2; 
	double* y3;
	int* y4;
	x = new double* [m];
	x1 = new double* [m];
	x2 = new double* [m];
	x3 = new double* [m];
	x4 = new int* [m];
	y = new double[m];
	y1 = new double[m];
	y2 = new double[m];
	y3 = new double[m];
	y4 = new int[m];
	for (int i = 0; i < m; i++)
	{
		x[i] = new double[n];
		x1[i] = new double[n];
		x2[i] = new double[n];
		x3[i] = new double[n];
		x4[i] = new int[n];

		y[i] = 0;
		for (int j = 0; j < n; j++)
			x[i][j] = 0;
	}

	cout << "GBTSR:\t";
	cout << "m = " << m << ", n = " << n << endl;
	
	cout << "\n===========================LP=============================\n";
	if (isLP == 1)
	{
		LP(x, y);
		result.init_(x, y);
		result.write_cplex_file(exp_path + "Algorithm\\result\\test\\MultiItTest\\n200BW50\\UAV_CLSn200l200_r1LP.txt");
	}
	else
		read_cplex_file(exp_path + "Algorithm\\result\\test\\MultiItTest\\n200BW50\\UAV_CLSn200l200_r1LP.txt", x, y);
	//print_RB(x, y);
	//print_cplex(x, y);
	
	//result.write_result_file("D:\\Myschool\\graduate_school\\02Graduate\\Research\\My paper\\2_Papers\\005_MHCDMC_Rounding\\MHCDMC_Rounding_Experiment\\Algorithm\\result\\n50BW50\\rLP.txt");
	

	cout << "\n===========================MII=============================\n";

	/*for (int i = 0; i < m; i++)
	{
		y1[i] = y[i];
		for (int j = 0; j < n; j++)
			x1[i][j] = x[i][j];
	}*/
	mergeII(x, y, x1, y1);
	//print_cplex(x1, y1);
	result.init_(x1, y1);
	//result.write_result_file("D:\\Myschool\\graduate_school\\02Graduate\\Research\\My paper\\2_Papers\\005_MHCDMC_Rounding\\MHCDMC_Rounding_Experiment\\Algorithm\\result\\n50BW50\\rDSIS_MII.txt");

	cout << "\n===========================DSIS=============================\n";
	DSIS(x1, y1);
	//print_RB(x1, y1);
	print_cplex(x1, y1);
	result.init_(x1, y1);
	//result.write_result_file("D:\\Myschool\\graduate_school\\02Graduate\\Research\\My paper\\2_Papers\\005_MHCDMC_Rounding\\MHCDMC_Rounding_Experiment\\Algorithm\\result\\n50BW50\\rDSIS.txt");


	cout << "\n===========================COS=============================\n";
	map<int, Cluster> CC = COS(x1, y1, x2, y2);
	print_CC(CC);
	//print_RB(x2, y2);
	print_cplex(x2, y2);
	result.init_(x2, y2);
	//result.write_result_file("D:\\Myschool\\graduate_school\\02Graduate\\Research\\My paper\\2_Papers\\005_MHCDMC_Rounding\\MHCDMC_Rounding_Experiment\\Algorithm\\result\\n50BW50\\rCOS.txt");


	cout << "\n===========================SFS=============================\n";
	SFS(x1, x2, y2, x3, y3, CC);
	print_RB(x3, y3);
	print_cplex(x3, y3);
	result.init_(x3, y3);
	// result.write_result_file("D:\\Myschool\\graduate_school\\02Graduate\\Research\\My paper\\2_Papers\\005_MHCDMC_Rounding\\MHCDMC_Rounding_Experiment\\Algorithm\\result\\n50BW50\\rSFS.txt");

	cal_NI(y3);
	//print_array(NI, "NI");

	/*cout << "\n===========================IP=============================\n";
	IP(x4, y4);
	print_cplex(x4, y4);
	result.init_(x4, y4);
	result.write_result_file("D:\\Myschool\\graduate_school\\02Graduate\\Research\\My paper\\2_Papers\\005_MHCDMC_Rounding\\MHCDMC_Rounding_Experiment\\Algorithm\\result\\test\\MultiItTest\\n200BW50\\UAV_AVEn200l200_IP.txt");*/


	for (int i = 0; i < m; i++)
	{
		delete x[i];
		delete x1[i];
		delete x2[i];
		delete x3[i];
	}
	delete[] x;
	delete[] x1;
	delete[] x2;
	delete[] x3;
	delete[] y1;
	delete[] y2;
	delete[] y3;
}

// 算法第一步，确定inferior 和superior servers
void Cover::DSIS(double** x1, double* y1)
{
	vector<int> allA;			// 记录y非零的服务器id
	vector<int> ou = all_u;		// 按照服务用户的服务器数量降序排序的用户id
	vector<int> ucount;			// 每个用户被覆盖的次数
	vector<int> SS;				// superior server的集合
	vector<int> II;				// inferior server的集合
	vector<int> Mid;			// y值介于alpha与1之间的服务器的集合

	// 生成II，SS，Mid
	for (int i = 0; i < m; i++)
	{
		if (y1[i] > 0)
		{
			allA.push_back(i);
			if (y1[i] == 1)
				SS.push_back(i);
			else if (y1[i] <= alpha)
				II.push_back(i);
			else
				Mid.push_back(i);
		}
	}
	
	/*result.init_(x1, y1);
	result.write_result_file("D:\\Myschool\\graduate_school\\02Graduate\\Research\\My paper\\2_Papers\\005_MHCDMC_Rounding\\MHCDMC_Rounding_Experiment\\Algorithm\\result\\n50BW50\\rDSIS_MII.txt");*/

	// 生成ou
	for (int j = 0; j < n; j++)
	{
		ucount.push_back(0);
		for (auto i : allA) {
			if (x1[i][j] > 0)
				ucount[j]++;
		}
	}
	cout << "order user befor: \n";
	for (int j = 0; j < n; j++)
	{
		cout << ou[j] << ' ' << ucount[ou[j]] << ", ";
	}
	cout << '\n';
	// 根据用户被覆盖的次数对用户排序
	//for (int j = 0; j < n; j++)
	//{
	//	// max_id = uid 而不是j或者jj
	//	int max_id = j;		// 先假定被覆盖次数最多的id为ou[j]
	//	for (int jj = j + 1; jj < n; jj++)
	//	{
	//		int u = ou[jj];
	//		int mu = ou[max_id];
	//		if (ucount[u] > ucount[mu])
	//			max_id = jj;
	//	}
	//	int temp = ou[j];
	//	ou[j] = ou[max_id];
	//	ou[max_id] = temp;
	//}
	cout << "order user after: \n";
	for (int j = 0; j < n; j++)
	{
		cout << ou[j] << ' ' << ucount[ou[j]] << ", ";
	}
	cout << '\n';
	
	/*cout << "\nSS: ";
	for (auto i : SS)
		cout << i << ' ';
	cout << "\nII: ";
	for (auto i : II)
		cout << i << ' ';
	cout << "\nMid: ";
	for (auto i : Mid)
		cout << i << ' ';
	cout << '\n';*/

	double cl = r * ep / sqrt(2);
	// cl = 500; 
	cout << "r = " << r << '\n';
	cout << "ep = " << ep << '\n';
	cout << "cl = " << cl << '\n';
	for (auto j : ou)
	{
		double sum_I = 0;		// I_j 中服务器对应yi的和
		vector<int> I_j;		// 服务xij的inferior server的集合
		for (auto a : II)
			if (x1[a][j] > 0)
			{
				if (y1[a] != 1 and y1[a] != 0)
				{
					sum_I += y1[a];
					I_j.push_back(a);
				}
			}
		if (sum_I > alpha)
		{
			cout << "-----------------------------------------\n";
			cout << "chosen u: " << j << ' ';
			U[j].print_user(); cout << '\n';
			// 构造符合alpha<sum_{ai_\in Ij}{y_i}<= 2 * alpha
			// 在construct_I_j函数中，I_j被修改成符合上述约束的集合
			construct_I_j(I_j, sum_I, y1);

			cout << "I_j: ";
			for (auto i : I_j)
				cout << i << ' ';
			cout << '\n';

			// cp为被划分区域的中心点
			Point cp = Point(U[j].X, U[j].Y);

			
			map<int, vector<int>> GG = construct_GG(I_j, cp, 4 * r, cl);
			map<int, vector<int>>::iterator it = GG.begin();
			map<int, vector<int>>::iterator itEnd = GG.end();
			cout << "GG before:\n";
			while (it != itEnd)
			{
				vector<int> g = it->second;
				cout << "\tcell id: " << it->first;
				cout << "\n\t\tg: ";
				for (auto i : g)
				{
					cout << i << ' ';
					A[i].print_server();
					cout << '\t';
				}
				cout << '\n';
				it++;
			}
			// 将GG中，只有一个服务器的G合并到能被合并的G
			cl = 500;
			merge_GG(GG, cl);
			it = GG.begin();
			itEnd = GG.end();
			cout << "GG after:\n";
			while (it != itEnd)
			{
				vector<int> g = it->second;
				cout << "\tcell id: " << it->first;
				cout << "\n\t\tg: ";
				for (auto i : g)
				{
					cout << i << ' ';
					A[i].print_server();
					cout << '\t';
				}
				cout << '\n';
				it++;
			}

			it = GG.begin();
			itEnd = GG.end();
			while (it != itEnd)
			{

				// 对GG中的每一个group进行处理
				vector<int> g = it->second;
				int im = g.back();		// 将am定义为g中最后一个
				g.pop_back();			// 

				reroute(all_u, g, im, x1);	// 将所有用户all_u在群组g中的所有流重引流到im
				y1[im] = 1;
				for (auto i : g)
					y1[i] = 0;


				it++; // 迭代器迭代
			}
		}
	}
	
	for (auto i : Mid)
		y1[i] = 1;
}

// 算法第二步，聚类inferior服务器
map<int, Cluster> Cover::COS(double** x1, double* y1, double** x2, double* y2)
{
	for (int i = 0; i < m; i++)
	{
		y2[i] = y1[i];
		for (int j = 0; j < n; j++)
			x2[i][j] = x1[i][j];
	}
	vector<int> allu = all_u;
	vector<int> allA;			// 所有服务器id集合
	vector<int> SS;				// superior server的集合
	vector<int> II;				// inferior server的集合
	// inferior server 到各个 superior server距离升序排序，保存序号
	// ISorder[2] = {3, 5, 6, 1} 的含义是，id为2的服务器，由近及远的无人机id为3,5,6,1
	map<int, vector<int>> ISorder;
	map<int, vector<int>> SIorder;
	vector<int> OO;				// 大O集合
	map<int, Cluster> CC;
	for (int i = 0; i < m; i++)
	{
		if (y2[i] > 0)
		{
			allA.push_back(i);
			if (y2[i] == 1)
			{
				Cluster C;
				C.push_back(i);
				CC[i] = C;
				SS.push_back(i);
			}
			else if (y2[i] <= alpha)
				II.push_back(i);
		}
	}

	cout << "II: ";
	for (auto i : II)
		cout << i << " ";
	cout << '\n';
	ISorder = get_ISorder(SS, II);
	SIorder = get_SIorder(SS, II);

	cout << "SIorder:\n";
	map<int, vector<int>>::iterator it = SIorder.begin();
	while (it != SIorder.end())
	{
		cout << "\ta" << it->first << ": ";
		for (auto i : it->second)
			cout << i << " ";
		cout << '\n';
		it++;
	}

	int count = 0;
	cout << "COS:\n";
	while (II.size() != 0)
	{
		int ii_len = II.size();
		int* is_clus = new int[ii_len];
		for (int tt = 0; tt < II.size(); tt++)
			is_clus[tt] = -1;
		
		for (int tt = 0; tt < II.size(); tt++)
		{
			if (is_clus[tt] != -1)
				continue;
			int t = II[tt];
			Server at = A[t];
			double sumBRt = get_sumBR(t, x2);

			cout << "\tround" << count++ << ":\n";
			cout << "\t\tat" << t << ", sumBRt=" << sumBRt << '\t';
			for (auto i : SIorder[t])
			{
				Server ai = A[i];
				double RBi = get_RB(i, x2);
				double d = ai.cal_distance(at);
				if (RBi >= sumBRt and d <= 2 * r)
				{
					cout << "ai" << i << ", sumBRt=" << sumBRt << '\t';
					// 额外步骤，将当前x2[t][j]的状态保存至x1[t][j]中
					for (int j = 0; j < n; j++)
						x1[t][j] = x2[t][j];

					vector<int> src;
					src.push_back(t);
					reroute(allu, src, i, x2);
					CC[i].push_back(t);

					is_clus[tt] = t;
					// 在II中删除t
					//auto iter = std::remove(II.begin(), II.end(), t);
					//II.erase(iter, II.end());
					break;
				}
			}
			cout << '\n';
		}
		cout << "is_clus: ";
		
		for (int tt = 0; tt < ii_len; tt++)
		{
			if (is_clus[tt] != -1)
			{
				int t = is_clus[tt];
				auto iter = std::remove(II.begin(), II.end(), t);
				II.erase(iter, II.end());
			}
		}


		//for (auto i : SS)
		//{
		//	cout << "\tround" << count++ << ":\n";
		//	Server ai = A[i];
		//	double RBi = get_RB(i, x2);
		//	cout << "\t\tai" << i << ", RBi=" << RBi << ": ";
		//	for (auto t : ISorder[i])
		//	{
		//		vector<int>::iterator iter = find(II.begin(), II.end(), t);
		//		if (iter != II.end())
		//		{
		//			Server at = A[t];
		//			double sumBRt = get_sumBR(t, x2);
		//			double d = ai.cal_distance(at);
		//			if (RBi >= sumBRt and d <= 2 * r)
		//			{
		//				cout << "at" << t << ", sumBRt=" << sumBRt << '\t';
		//				// 额外步骤，将当前x2[t][j]的状态保存至x1[t][j]中
		//				for (int j = 0; j < n; j++)
		//					x1[t][j] = x2[t][j];
		//				vector<int> src;
		//				src.push_back(t);
		//				reroute(allu, src, i, x2);
		//				CC[i].push_back(t);
		//				// 在II中删除t
		//				auto iter = std::remove(II.begin(), II.end(), t);
		//				II.erase(iter, II.end());
		//			}
		//		}
		//		RBi = get_RB(i, x2);
		//	}
		//	cout << '\n';
		//}
		cout << "II: ";
		for (auto i : II)
			cout << i << " ";
		cout << '\n';

		if (II.size() == 0)
			break;
		else
		{
			// 保存II中所有无人机对应的Ai，Ai为ai服务的用户的id集合
			map<int, vector<int>> AA;
			// Ki保存II中所有无人机对应的ki
			vector<double> Ki;
			for (auto i : II)
			{
				// 对II中的无人机ai，构造器服务的用户的id集合Ai
				vector<int> Ai;
				for (int j = 0; j < n; j++)
					if (x2[i][j] > 0)Ai.push_back(j);
				AA[i] = Ai;

				// 计算ai对应的ki
				double ki = 0;
				for (auto j : Ai)
					ki += U[j].BR;
				if (A[0].BW < ki)
					ki = A[0].BW;
				Ki.push_back(ki);
			}

			// 在II中找到ki最小的无人机II[tt]及最小的ki值=Ki[tt]
			int tt = 0;
			for (int ii = 0; ii < II.size(); ii++)
				if (Ki[ii] > Ki[tt])
					tt = ii;
			// at = II[tt]
			int t = II[tt];
			double kt = Ki[tt];
			OO.push_back(t);
			// 在II中删除t
			auto iter = std::remove(II.begin(), II.end(), t);
			II.erase(iter, II.end());

			double BRmint = get_BRmin(AA[t]);
			double sumBRt = 0;
			for (auto j : AA[t])
				sumBRt += U[j].BR;
			// if kt == sumBRt <= BW
			if (kt == sumBRt and kt <= A[0].BW)
			{
				vector<int> src = allA;
				for (auto i : OO)
				{
					auto it = std::remove(src.begin(), src.end(), i);
					src.erase(it, src.end());
				}
				reroute(AA[t], src, t, x2);
			}
			// if kt == BW < sumBRt and kt >= 2*BRmint
			if (kt == A[0].BW and kt >= 2 * BRmint)
			{
				double RBt = get_RB(t, x2);
				int j = get_max_xBR(t, AA[t], RBt, x2);
				while (j != -1)
				{
					vector<int> uk;
					vector<int> src = allA;
					for (auto i : OO)
					{
						auto it = std::remove(src.begin(), src.end(), i);
						src.erase(it, src.end());
					}
					uk.push_back(j);
					reroute(uk, src, t, x2);
					RBt -= x2[t][j] * U[j].BR;
					j = get_max_xBR(t, AA[t], RBt, x2);
				}
			}
			// if kt == BW < sumBRt and kt >= 2*BRmint
			if (kt == A[0].BW and kt <= 2 * BRmint)
			{
				double RBt = get_RB(t, x2);
				int j = get_max_xBR(t, AA[t], RBt, x2);
				vector<int> uk;
				uk.push_back(j);
				reroute(uk, II, t, x2);
				RBt = get_RB(t, x2);

				double f = 0;
				for (auto i : OO)
					f += x2[i][j];
				if (RBt / U[j].BR < f)
					f = RBt / U[j].BR;

				vector<int> src;

				map<int, vector<int>>::iterator it = CC.begin();
				while (it != CC.end())
				{
					for (auto i : it->second)
						src.push_back(i);
					it++;
				}
				reroute(uk, src, t, x2, f);
				
				
			}
		}
	}

	for (auto i : OO) {
		y2[i] = 1;
		Cluster C;
		C.push_back(i);
		CC[i] = C;
	}

	return CC;
}

// 根据当前x，y状态，将一些可以被其它inferior服务器取代的服务器合并
void Cover::mergeII(double** x, double* y, double** x1, double* y1)
{
	// 初始化
	for (int i = 0; i < m; i++)
	{
		y1[i] = y[i];
		for (int j = 0; j < n; j++)
		{
			// cout << "i = " << i << ", j = " << j << '\t';
			x1[i][j] = x[i][j];
		}

	}
	vector<int> II;				// inferior server的集合
	map<int, vector<int>> AA;	// 保存各个被选中服务器服务的用户 
	for (int i = 0; i < m; i++)
	{
		if (y1[i] > 0)
		{
			// 生成II
			if (y1[i] < 1 and y1[i] > 0)
				II.push_back(i);

			// 生成被选中服务器服务的用户集合
			vector<int> Ai;
			for (int j = 0; j < n; j++)
				if (x1[i][j] > 0)
					Ai.push_back(j);
			AA[i] = Ai;
		}
	}

	cout << "AA:\n";
	map<int, vector<int>>::iterator it = AA.begin();
	while (it != AA.end())
	{
		cout << "\tA" << it->first << ": ";
		for (auto i : it->second)
			cout << i << ' ';
		cout << '\n';
		it++;
	}
	cout << "merge II befor: ";
	for (auto i : II)
		cout << i << ' ';
	cout << '\n';

	// 合并能够被取代的inferior server
	// is_merge为一个长度为II.size()的数组，记录对应下标服务器是否被merge
	// 0为“未被merge”		1为被merge
	int* is_merge = new int[II.size()];
	for (int i = 0; i < II.size(); i++)
		is_merge[i] = 0;

	for (int ii = 0; ii < II.size(); ii++)
	{
		// 若下标为ii的服务器被merge了，就下个循环
		if (is_merge[ii] == 1)
			continue;
		// i为当前ii对应的服务器id
		int i = II[ii];
		Server ai = A[i];
		// Ai为服务器i当前服务的用户
		vector<int> Ai = AA[i];
		for (int tt = 0; tt < II.size(); tt++)
		{
			// 合并下标为tt的服务器与ii的服务器
			// 若下标为tt的服务器被merge了，就下个循环
			if (is_merge[tt] == 1)
				continue;
			if (ii == tt)
				continue;
			int t = II[tt];
			Server at = A[t];
			vector<int> At = AA[t];

			if (ai.cal_distance(at) > 2 * r)
				continue;
			// 在合并的过程中，II中的服务器可能变得superior，因此需要此判断
			if (y1[i] + y1[t] > 1)
				continue;

			// 
			if (is_contain(At, Ai, i))
			{
				cout << "\tAt" << t << ": ";
				for (auto t0 : At)
					cout << t0 << ' ';
				cout << '\n' << "\tAi" << i << ": ";
				for (auto i0 : Ai)
					cout << i0 << ' ';
				cout << '\n';
				cout << '\n';


				y1[i] += y1[t];
				y1[t] = 0;

				for (auto j : At)
				{
					x1[i][j] += x1[t][j];
					x1[t][j] = 0;
				}

				
				is_merge[tt] = 1;

			}

		}
	}

	II.clear();
	for (int i = 0; i < m; i++)
	{
		if (y1[i] > 0)
		{
			// 生成II
			if (y1[i] <= alpha)
				II.push_back(i);

		}
	}

	cout << "merge II after: ";
	for (auto i : II)
		cout << i << ' ';
	cout << '\n';

	print_cplex(x1, y1);
}

// 算法第三步，在每个cluster中选择最终的服务器
void Cover::SFS(double** x1, double** x2, double* y2, double** x3, double* y3, map<int, Cluster>& CC)
{
	// 参数：
	// x1, y1		暂存某个请求被聚类到某个重球的cluster之前的x值
	// x2, y2		COS之后的x，y值
	// x3, y3		SFS之后的x，y值

	for (int i = 0; i < m; i++)
	{
		y3[i] = y2[i];
		for (int j = 0; j < n; j++)
			x3[i][j] = x2[i][j];
	}
	vector<int> allu = all_u;
	vector<int> allA;			// 所有服务器id集合
	vector<int> SS;				// superior server的集合
	vector<int> II;				// inferior server的集合
	for (int i = 0; i < m; i++)
	{
		if (y3[i] > 0)
		{
			allA.push_back(i);
			if (y3[i] == 1)
				SS.push_back(i);
			else if (y2[i] <= alpha)
				II.push_back(i);
		}
	}

	map<int, Cluster>::iterator it = CC.begin();

	while (it != CC.end())
	{
		Cluster Ch = it->second;
		if (Ch.size() == 1)
		{
			y3[Ch[0]] = 1;
		}
		else
		{
			// 将C中轻球对应inferior server对应的x值恢复到加入C之前的状态
			vector<int>C;
			C.push_back(Ch[0]);
			for (int ii = 1; ii < Ch.size(); ii++)
			{
				int i = Ch[ii];
				C.push_back(i);
				for (int j = 0; j < n; j++)
				{
					x3[i][j] = x1[i][j];
					x3[Ch[0]][j] -= x3[i][j];
					if (x3[Ch[0]][j] < 0)
						x3[Ch[0]][j] = 0;
				}
			}
			Server ah = A[Ch[0]];
			double cl = (r * ep) / 4;
			// cl = 200;
			map<int, vector<int>> GG = construct_GG(C, ah, 6 * r, cl);
			merge_GG(GG, cl);
			map<int, vector<int>>::iterator it = GG.begin();
			while (it != GG.end())
			{
				vector<int> G = it->second;
				int im = G.back();
				G.pop_back();
				reroute(allu, G, im, x3);

				y3[im] = 1;
				for (auto i : G)
					y3[i] = 0;
				it++;
			}

		}

		it++;
	}


}

// 构造符合alpha<sum_{ai_\in Ij}{y_i}<= 2 * alpha
// 在construct_I_j函数中，I_j被修改成符合上述约束的集合 
void Cover::construct_I_j(vector<int>& I_j, double sum_I, double* y)
{
	while (sum_I > 2 * alpha)
	{
		int i = I_j.back();
		sum_I -= y[i];
		I_j.pop_back();
	}
}

// 将obj中的用户的流，从src中reroute到aim服务器中，需要修改x；flow为留的大小，默认为1，若不为1则
void Cover::reroute(vector<int>& obj, vector<int>& src, int aim, double** x, double flow)
{
	// 将obj中的用户的flow这么多的流，从src重引流到aim
	// 
	if (flow == 1)
	{
		for (auto j : obj)
		{
			for (auto i : src)
			{
				x[aim][j] += x[i][j];
				x[i][j] = 0;
			}
		}
	}
	else {
		for (int j = 0; j < n; j++)
		{
			for (auto i : src)
			{
				if (flow != 0)
				{
					if (x[i][j] >= flow)
					{
						x[aim][j] += flow;
						x[i][j] -= flow;
						flow = 0;
					}
					else if (x[i][j] > 0)
					{
						x[aim][j] += x[i][j];
						flow -= x[i][j];
						x[i][j] = 0;
					}
				}
				else
					break;
			}
			if (flow == 0)
				break;
		}
	}
}

// 找到i服务的用户集合Ai中，需求最小的BR 
double Cover::get_BRmin(vector<int>& Ai)
{
	double min = A[0].BW;
	for (auto j : Ai)
		if (min > U[j].BR)
			min = U[j].BR;
	return min;
}

// 在无人机t当前服务的用户集合At中，找到最适合那个用户。这个用户将被重引流
int Cover::get_max_xBR(int t, vector<int>& At, double RBt, double** x)
{
	// 在无人机t当前服务的用户集合At中，找到最适合那个用户。这个用户将被重引流
	// t: 当前无人机id
	// At		无人机t当前服务的用户集合(x[t][j]>0)
	// RBt		无人机t当前的剩余容量
	// x
	//	返回值：max_u			若max_u==-1,这说明没有合适的用户
	int max_u = -1;		
	double max_xBR = 0;
	for (int j : At)
	{
		double xBRj = x[t][j] * U[j].BR;
		if (xBRj <= RBt)
		{
			if (xBRj > max_xBR)
				max_u = j;
		}
	}


	return max_u;
}

// 判断At是否是Ai的子集
bool Cover::is_contain(vector<int>& At, vector<int>& Ai, int i)
{
	Server ai = A[i];
	
	for (auto t : At)
	{
		/*if (find(Ai.begin(), Ai.end(), t) == Ai.end())
			return false;*/
		if (ai.cal_distance(U[t]) > r)
			return false;
	}

	return true;
}

// inferior server 到各个 superior server距离升序排序，保存序号
map<int, vector<int>> Cover::get_ISorder(vector<int>& SS, vector<int>& II)
{
	// inferior server 到各个 superior server距离升序排序，保存序号
	// ISorder[2] = {3, 5, 6, 1} 的含义是，id为2的服务器，由近及远的无人机id为3,5,6,1

	map<int, vector<int>> ISorder;

	// 初始化
	for (auto s : SS)
		ISorder[s] = II;


	map<int, vector<int>>::iterator it = ISorder.begin();
	/*cout << "ISorder:\n";
	it = ISorder.begin();
	while (it != ISorder.end())
	{
		cout << "a" << it->first << ": ";
		for (auto i : it->second)
			cout << i << ' ';
		cout << '\n';
		it++;
	}*/

	// cout << "过程：\n";
	it = ISorder.begin();
	while (it != ISorder.end())
	{
		int s = it->first;
		// cout << "a" << i << ": "; 
		Server as = A[s];
		for (int ii = 0; ii < it->second.size(); ii++)
		{
			int clo_i = ii;
			double clo_d = as.cal_distance(A[it->second[clo_i]]);
			for (int ii2 = ii + 1; ii2 < it->second.size(); ii2++)
			{
				double d = as.cal_distance(A[it->second[ii2]]);
				if (d < clo_d)
				{
					clo_i = ii2;
					clo_d = d;
				}
			}
			int temp = it->second[ii];
			it->second[ii] = it->second[clo_i];
			it->second[clo_i] = temp;
		}
		//for (int ss = 0; ss < it->second.size(); ss++)
		//{
		//	int clo_s = ss;
		//	double clo_d = ai.cal_distance(A[it->second[clo_s]]);
		//	// cout << '\t' << 'd' << it->second[clo_s] << ": " << clo_d << '\n';
		//	for (int ss2 = ss + 1; ss2 < it->second.size(); ss2++)
		//	{

		//		double d = ai.cal_distance(A[it->second[ss2]] );
		//		if (d < clo_d)
		//		{
		//			clo_s = ss2;
		//			clo_d = d;
		//		}
		//	}

		//	int temp = it->second[ss];
		//	it->second[ss] = it->second[clo_s];
		//	it->second[clo_s] = temp;
		//}
		//
		it++;
	}

	/*cout << "ISorder:\n";
	it = ISorder.begin();
	while (it != ISorder.end())
	{
		cout << "a" << it->first << ": ";
		for (auto i : it->second)
			cout << i << ' ';
		cout << '\n';
		it++;
	}*/
	return ISorder;
}

// superior server 到各个 inferior server距离升序排序，保存序号
map<int, vector<int>> Cover::get_SIorder(vector<int>& SS, vector<int>& II)
{
	// inferior server 到各个 superior server距离升序排序，保存序号
	// ISorder[2] = {3, 5, 6, 1} 的含义是，id为2的服务器，由近及远的无人机id为3,5,6,1

	map<int, vector<int>> SIorder;

	// 初始化
	for (auto i : II)
		SIorder[i] = SS;


	map<int, vector<int>>::iterator it = SIorder.begin();
	/*cout << "ISorder:\n";
	it = ISorder.begin();
	while (it != ISorder.end())
	{
		cout << "a" << it->first << ": ";
		for (auto i : it->second)
			cout << i << ' ';
		cout << '\n';
		it++;
	}*/

	// cout << "过程：\n";
	it = SIorder.begin();
	while (it != SIorder.end())
	{
		int i = it->first;
		// cout << "a" << i << ": "; 
		Server ai = A[i];
		for (int ss = 0; ss < it->second.size(); ss++)
		{
			int clo_s = ss;
			double clo_d = ai.cal_distance(A[it->second[clo_s]]);
			for (int ss2 = ss + 1; ss2 < it->second.size(); ss2++)
			{
				double d = ai.cal_distance(A[it->second[ss2]]);
				if (d < clo_d)
				{
					clo_s = ss2;
					clo_d = d;
				}
			}
			int temp = it->second[ss];
			it->second[ss] = it->second[clo_s];
			it->second[clo_s] = temp;
		}
		//for (int ss = 0; ss < it->second.size(); ss++)
		//{
		//	int clo_s = ss;
		//	double clo_d = ai.cal_distance(A[it->second[clo_s]]);
		//	// cout << '\t' << 'd' << it->second[clo_s] << ": " << clo_d << '\n';
		//	for (int ss2 = ss + 1; ss2 < it->second.size(); ss2++)
		//	{

		//		double d = ai.cal_distance(A[it->second[ss2]] );
		//		if (d < clo_d)
		//		{
		//			clo_s = ss2;
		//			clo_d = d;
		//		}
		//	}

		//	int temp = it->second[ss];
		//	it->second[ss] = it->second[clo_s];
		//	it->second[clo_s] = temp;
		//}
		//
		it++;
	}

	/*cout << "ISorder:\n";
	it = ISorder.begin();
	while (it != ISorder.end())
	{
		cout << "a" << it->first << ": ";
		for (auto i : it->second)
			cout << i << ' ';
		cout << '\n';
		it++;
	}*/
	return SIorder;
}

// 本函数将包含集合I的区域，该区域以点cp为中心，边长为L，划分为边长为cl的小方格
map<int, vector<int>> Cover::construct_GG(vector<int>& I, Point& cp, double L, double cl)
{
	// 本函数将包含集合I的区域，该区域以点cp为中心，边长为L，划分为边长为cl的小方格
	// 返回值为GG，即一个包含所有groups的集合。数据结构为map<int, vector<int>> (第一项为网格序号，第二项为group包含的服务器)
	// 参数列表:
	// I		将被划分的服务器集合
	// cp		I所处区域的中心点
	// L		I所处区域的边长
	// cl		网格边长
	map<int, vector<int>> GG;
	// 被划分区域左下角点坐标(X0, Y0)
	double X0 = cp.X - L / 2;
	double Y0 = cp.Y - L / 2;
	int lines = ceil(L / cl);


	for (auto i : I)
	{
		Server a = A[i];
		// ii, jj代表的是a所处的网格坐标（从下往上，从左往右，原点为X0，Y0）
		int ii = int((a.X - X0) / cl);
		int jj = int((a.Y - Y0) / cl);
		// key代表的是a所处的网格id（从下往上，从左往右，原点为X0，Y0）
		int key = ii + jj * lines;

		// 在GG中找以key为关键字的键值对，若没有，则GG.find(key)指向GG.end()
		if (GG.find(key) == GG.end())
		{
			// 未找到以key为键的键值对
			vector<int> g;
			g.push_back(i);
			GG[key] = g;
		}
		else
			GG[key].push_back(i);
	}

	// 合并GG中那些g只有一个服务器的g

	return GG;
}

// 将一些质包含一个服务器的G，合并到符合条件的G中
void Cover::merge_GG(map<int, vector<int>>& GG, double cl)
{
	map<int, vector<int>> GG_ = GG;			// GG副本，在后面GG中某个G合并新g之后，之后要跟他合并的g只能与G.back()比较距离

	map<int, vector<int>>::iterator it1 = GG.begin();
	map<int, vector<int>>::iterator itEnd = GG.end();

	vector<int> erased_key;	// 记录在GG中，被合并的长度为1的g的key值。在之后根据key值，在map中删除
	// 记录每个是否被合并，被合并的话值为被合并的key值
	int* is_mer = new int[GG.size()];
	for (int i = 0; i < GG.size(); i++)
		is_mer[i] = 0;

	int j = 0;
	while (it1 != itEnd)
	{
		cout << j << '\n';
		if (is_mer[j] == 1)
			break;
		// 外循环，找到长度为1的g
		vector<int> g = it1->second;
		// 记录当前g是否被合并，未被合并为0
		if (g.size() == 1)
		{
			// 只有当前的g的长度为1才能被合并，进行以下步骤
			// g对应的key
			int key = it1->first;
			// g中的唯一的server
			Server a1 = A[g[0]];

			// 从map开头为g寻找能够加入的G，it2指向it1将要并入的group
			map<int, vector<int>>::iterator it2 = GG.begin();
			int jj = 0;
			while (it2 != itEnd) {
				// 每次it2都要判断g是否被merge，上个迭代被merge后，这个迭代就会退出迭代
				if (is_mer[j] == 1)
					break;
				// 判断it2指向的group是否被merge，若被merge则扫描下一个
				if (is_mer[jj] == 1)
				{

				}
				else {
					if (it2 == it1)
					{
						// 若it2==it1，表明it2指向了当前g，则什么都不做
					}
					else
					{
						// 记录G对应的KEY值
						int KEY = it2->first;
						// G是g可能加入的对象
						vector<int> G = it2->second;
						// G_是G的原身，若他们相同，则说明G没有加入过g
						vector<int> G_ = GG_[KEY];
						// 若他们相同，则说明G没有加入过g，需要将g与G中每一个服务器比较距离
						if (G.size() == G_.size())
						{
							// 遍历G，ii是其中元素的下标
							for (int ii = 0; ii < G.size(); ii++)
							{
								// i为当前的服务器id
								int i = G[ii];
								// 若g与A[i]距离小于cl，就可以merge
								if (a1.cal_distance(A[i]) <= cl)
								{
									// 将g中服务器id g[0]保存到G中ii下标初，会覆盖i
									it2->second[ii] = g[0];
									// 在末尾加入i
									it2->second.push_back(i);
									// 在erased_key中加入当前被merge的g的键值
									erased_key.push_back(key);
									is_mer[j] = 1;
									break;
								}
							}
						}
						else
						{
							// 如果G中已经merge过服务器，那么就只能与G中的最后一个服务器比较距离
							cout << "before " << it2->first << ':';
							for (auto i : it2->second)
								cout << i << ' ';
							cout << '\n';
							int i = G.back();
							int ii = G.size() - 1;
							if (a1.cal_distance(A[i]) <= cl)
							{
								it2->second[ii] = g[0];
								it2->second.push_back(i);
								erased_key.push_back(key);
								is_mer[j] = 1;
							}
							cout << "after " << it2->first << ':';
							for (auto i : it2->second)
								cout << i << ' ';
							cout << '\n';
						}
					}
				}

				jj++;
				it2++;
			}

		}

		j++;
		it1++;
	}

	for (auto key : erased_key)
	{
		GG.erase(key);
	}
}

// 根据xx与yy的值，对对象初始化
template<class TT, class T>
void Result::init_(TT xx, T yy)
{
	for (int i = 0; i < m; i++)
	{
		y[i] = yy[i];
		if (yy[i] > 0)
		{
			
			for (int j = 0; j < n; j++)
			{
				x[i][j] = 0;
				x[i][j] = xx[i][j];
			}
		}
	}
}


// 将x，y解析成以下结果写入文件
// 第1行: chosen_m n |chosen_m: 被选择服务器数 n: 用户数
// 第2k+1行: aid ucont |k=0,1,...,chosen_m-1; aid: 服务器id; ucont: 
// 第2k+2行: u0 u1 ... u_ucont-1 |u0... : 被2k+1行中对应服务器服务的用户	
void Result::write_result_file(string fname)
{
	// 将x，y解析成以下结果写入文件
	// 第1行: chosen_m n |chosen_m: 被选择服务器数 n: 用户数
	// 第2k+1行: aid ucont |k=0,1,...,chosen_m-1; aid: 服务器id; ucont: 
	// 第2k+2行: u0 u1 ... u_ucont-1 |u0... : 被2k+1行中对应服务器服务的用户	
	vector<string> wcontent;
	wcontent.push_back("\n");	// 占位，第一行为"chosen_m n"，即被选择的服务器数和用户数
	chosen_m = 0;
	for (int i = 0; i < m; i++)
	{
		if (y[i] > 0)
		{
			chosen_m++;
			string con_u;	// 被i服务的用户id字符串
			int u_cont = 0;	// 被i服务的用户数
			for (int j = 0; j < n; j++)
			{
				if (x[i][j] > 0)
				{
					con_u += to_string(j) + ' ';
					u_cont++;
				}
			}
			con_u += '\n';
			wcontent.push_back(to_string(i) + ' ' + to_string(u_cont) + '\n');
			wcontent.push_back(con_u);
		}
	}
	wcontent[0] = to_string(chosen_m) + ' ' + to_string(n) + '\n';
	ofstream output;
	output.open(fname, ios::app);
	for (auto w : wcontent)
		output << w;
	output.close();
}

void Result::write_cplex_file(string fname)
{
	// 将x，y的值写入文件
	// 第一行：choosen_m n
	// 第二行：yid y[i]
	// 第三行：x[0] x[1] ... x[n]
	vector<string> wcontent;
	wcontent.push_back("\n");	// 占位，第一行为"chosen_m n"，即被选择的服务器数和用户数
	chosen_m = 0;
	for (int i = 0; i < m; i++)
	{
		if (y[i] > 0)
		{
			wcontent.push_back(to_string(i) + ' ' + to_string(y[i]) + '\n');
			chosen_m++;
			string x_values;	// 对应所有x[i]的值
			int u_cont = 0;		// 被i服务的用户数
			for (int j = 0; j < n; j++)
			{
				x_values += to_string(x[i][j]) + ' ';
			}
			x_values += '\n';
			wcontent.push_back(x_values);
		}
	}
	wcontent[0] = to_string(chosen_m) + ' ' + to_string(n) + '\n';
	ofstream output;
	output.open(fname, ios::app);
	for (auto w : wcontent)
		output << w;
	output.close();
}
