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
			d[i][j] = sqrt(pow(A[i].X - U[j].X, 2) + pow(A[i].Y - U[j].Y, 2) + 300 * 300);
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
	// double N = N0 + 10 * log10(A[0].BW * 10e6) + N_I;
	double N = N0 + 10 * log10(20 * 10e6) + N_I;
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++)
			SINR[i][j] = A[i].p - L[i][j] - N;
	}
}

void Cover::cal_min_r()
{
	//double BW = A[0].BW * 10e6;
	double BW = 20 * 10e6;
	double eta = 20 * log10((4 * PI * f) / c) + eta_LOS;
	double N = N0 + 10 * log10(BW) + 36;
	double zhishu = (power - SINR_min - eta - N) / 20;
	double R = pow(10, zhishu);
	r = sqrt(pow(R, 2) - pow(h, 2));

	cout << "N0 = " << N0 << '\n';
	cout << "eta = " << eta << '\n';
	cout << "zhishu = " << zhishu << '\n';
	cout << "r = " << r << '\n';
}

void Cover::cal_ep()
{
	double zhishu = (power * ep_p + SINR_min * ep_SINR) / 20;

	ep = pow(10, zhishu) - 1;
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

	result = Result(m, n);
	for (int j = 0; j < n; j++)
		all_u.push_back(j);

	cal_d();
	cal_L();
	cal_SINR();
	cal_min_r();
	cal_ep();
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
		Server s(x, y, i, 20);
		A.push_back(s);
	}

	for (int j = 0; j < n; j++) {
		double x;
		double y;
		int BR;
		getline(infpoint, data);
		istringstream istr3(data);
		istr3 >> x >> y >> BR;
		User user(x, y, j, BR);
		U.push_back(user);
	}
}

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
		if (y[i] != 0)
			cout << "y" << i << '=' << y[i] << '\t';
	}
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
		}
	}
}

template<class TT, class T>
void Cover::print_RB(TT x, T y)
{
	map<int, double> RBs;
	for (int i = 0; i < m; i++)
	{
		if (y[i] > 0)
		{
			double RB = 0;
			for (int j = 0; j < n; j++)
				RB += x[i][j] * U[j].BR;
			RB = A[i].BW * y[i] - RB;
			RBs[i] = RB;
		}
	}
	cout << "RBs:\n";
	map<int, double>::iterator it = RBs.begin();
	while (it != RBs.end()) {
		cout << "a" << it->first << ": " << it->second << '\t';
		it++;
	}
}

template<class T>
double Cover::get_sumBR(int i, T x)
{
	double sum = 0;
	for (int j = 0; j < n; j++)
		sum += x[i][j] * U[j].BR;
	return sum;
}

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
void Cover::LP(double** x0, double* y0)
{
	clock_t start_, end_;
	start_ = clock();

	
	// IloEnv env_LP;
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

	int* y0 = new int[m];

	// IloEnv env_IP;
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


	print_cplex(x0, y0);
	Cover::result.init(x0, y0, m, n);
}

void Cover::GBTSR()
{
	double** x;
	double** x1;
	double** x2;
	double** x3;
	double* y;
	double* y1;
	double* y2; 
	double* y3;
	x = new double* [m];
	x1 = new double* [m];
	x2 = new double* [m];
	x3 = new double* [m];
	y = new double[m];
	y1 = new double[m];
	y2 = new double[m];
	y3 = new double[m];
	for (int i = 0; i < m; i++)
	{
		x[i] = new double[n];
		x1[i] = new double[n];
		x2[i] = new double[n];
		x3[i] = new double[n];
	}


	LP(x, y);
	print_RB(x, y);
	print_cplex(x, y);
	result.int_(x, y);


	DSIS(x, y, x1, y1);
	print_RB(x1, y1);
	print_cplex(x1, y1);
	result.int_(x1, y1);


	map<int, Cluster> CC = COS(x1, y1, x2, y2);
	print_RB(x2, y2);
	print_cplex(x2, y2);
	result.int_(x2, y2);





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

void Cover::DSIS(double** x, double* y, double** x1, double* y1)
{
	for (int i = 0; i < m; i++)
	{
		y1[i] = y[i];
		for (int j = 0; j < n; j++)
			x1[i][j] = x[i][j];
	}
	vector<int> allA;			// 记录y非零的服务器id
	vector<int> ou = all_u;		// 按照服务用户的服务器数量降序排序的用户id
	vector<int> ucount;			// 每个用户被覆盖的次数
	vector<int> SS;				// superior server的集合
	vector<int> II;				// inferior server的集合
	vector<int> Mid;			// y值介于alpha与1之间的服务器的集合
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

	// 生成ou
	for (int j = 0; j < n; j++)
	{
		ucount.push_back(0);
		for (auto i : allA) {
			if (x[i][j] > 0)
				ucount[j]++;
		}
	}
	cout << "befor: \n";
	for (int j = 0; j < n; j++)
	{
		cout << ou[j] << ' ' << ucount[ou[j]] << ", ";
	}
	cout << '\n';
	for (int j = 0; j < n; j++)
	{
		// max_id = uid 而不是j或者jj
		int max_id = j;		// 先假定被覆盖次数最多的id为ou[j]
		for (int jj = j + 1; jj < n; jj++)
		{
			int u = ou[jj];
			int mu = ou[max_id];
			if (ucount[u] > ucount[mu])
				max_id = jj;
		}
		int temp = ou[j];
		ou[j] = ou[max_id];
		ou[max_id] = temp;
	}
	cout << "befor: \n";
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
	// cl = 1000; 
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
				cout << "cell id: " << it->first;
				cout << "\ng: ";
				for (auto i : g)
				{
					cout << i << ' ';
					A[i].print_server();
					cout << '\t';
				}
				cout << '\n';
				it++;
			}
			merge_GG(GG, cl);
			it = GG.begin();
			itEnd = GG.end();
			cout << "GG after:\n";
			while (it != itEnd)
			{
				vector<int> g = it->second;
				cout << "cell id: " << it->first;
				cout << "\ng: ";
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


	while (II.size() != 0)
	{
		for (auto i : SS)
		{
			Server ai = A[i];
			double RBi = get_RB(i, x2);
			for (auto t : II)
			{
				Server at = A[t];
				double sumBRt = get_sumBR(t, x2);
				if (RBi >= sumBRt)
				{
					vector<int> src;
					src.push_back(t);
					reroute(allu, src, i, x2);
					CC[i].push_back(t);
					// 在II中删除t
					auto iter = std::remove(II.begin(), II.end(), t);
					II.erase(iter, II.end());
				}
			}
		}

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

Result Cover::SFS(double** x2, double* y2, double** x3, double* y3, vector<Cluster>& CC)
{
	return Result();
}

void Cover::construct_I_j(vector<int>& I_j, double sum_I, double* y)
{
	while (sum_I > 2 * alpha)
	{
		int i = I_j.back();
		sum_I -= y[i];
		I_j.pop_back();
	}
}

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

double Cover::get_BRmin(vector<int>& Ai)
{
	double min = A[0].BW;
	for (auto j : Ai)
		if (min > U[j].BR)
			min = U[j].BR;
	return min;
}

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

void Cover::merge_GG(map<int, vector<int>>& GG, double cl)
{
	map<int, vector<int>> GG_ = GG;			// GG副本，在后面GG中某个G合并新g之后，之后要跟他合并的g只能与G.back()比较距离

	map<int, vector<int>>::iterator it1 = GG.begin();
	map<int, vector<int>>::iterator itEnd = GG.end();

	vector<int> erased_key;	// 记录在GG中，被合并的长度为1的g的key值。在之后根据key值，在map中删除


	while (it1 != itEnd)
	{
		// 外循环，找到长度为1的g
		vector<int> g = it1->second;
		// 记录当前g是否被合并，未被合并为0
		int is_merge = 0;
		if (g.size() == 1)
		{
			// 只有当前的g的长度为1才能被合并，进行以下步骤
			// g对应的key
			int key = it1->first;
			// g中的唯一的server
			Server a1 = A[g[0]];

			// 从map开头为g寻找能够加入的G
			map<int, vector<int>>::iterator it2 = GG.begin();
			while (it2 != itEnd) {
				if (it2 == it1)
				{
					// 若it2==it1，表明it2指向了当前g，则什么都不做
				}
				else
				{
					// 每次it2都要判断g是否被merge，上个迭代被merge后，这个迭代就会退出迭代
					if (is_merge == 0)
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
									is_merge = 1;
									break;
								}
							}
						}
						else
						{
							// 如果G中已经merge过服务器，那么就只能与G中的最后一个服务器比较距离
							int i = G.back();
							int ii = G.size() - 1;
							if (a1.cal_distance(A[i]) <= cl)
							{
								it2->second[ii] = g[0];
								it2->second.push_back(i);
								erased_key.push_back(key);
								is_merge = 1;
							}
						}

					}
					else
					{
						break;
					}
				}
				it2++;
			}

		}
		it1++;
	}

	for (auto key : erased_key)
	{
		GG.erase(key);
	}
}


void Result::int_(double** xx, double* yy)
{
	for (int i = 0; i < m; i++)
	{
		y[i] = 0;
		if (yy[i] > 0)
		{
			y[i] = 1;
			
			for (int j = 0; j < n; j++)
			{
				x[i][j] = 0;
				if (xx[i][j] > 0)
					x[i][j] = 1;
			}
		}
	}
}

void Result::write_result_file(string fname)
{
	// 将x，y解析成以下结果写入文件
	// 第1行: chosen_m n |chosen_m: 被选择服务器数 n: 用户数
	// 第2k+1行: aid ucont |k=0,1,...,chosen_m-1; aid: 服务器id; ucont: 
	// 第2k+2行: u0 u1 ... u_ucont-1 |u0... : 被2k+1行中对应服务器服务的用户	
	vector<string> wcontent;
	wcontent.push_back("\n");	// 占位，第一行为"chosen_m n"，即被选择的服务器数和用户数
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
