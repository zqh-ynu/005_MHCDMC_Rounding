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
	// ��ȡm, n
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

		// Լ��1.1
		// ��u��a���ǣ�a���뱻ѡ��
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

		// Լ��1.2
		// ��������Լ��
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

		// Լ��1.3
		// ÿ���û���Ҫ������
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

		// Լ��1.4
		// ͨ������Լ��
		for (IloInt i = 0; i < m; i++)
			for (IloInt j = 0; j < n; j++)
				if (SINR[i][j] < SINR_min)
				{
					IloRange SINR_con(env_LP, 0, x[i][j], 0);
					model.add(SINR_con);
				}


		// Ŀ�꺯��
		IloExpr obj(env_LP);
		for (IloInt i = 0; i < m; i++)
			obj += y[i];


		model.add(IloMinimize(env_LP, obj));

		// ����������
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

		// Լ��1.1
		// ��u��a���ǣ�a���뱻ѡ��
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

		// Լ��1.2
		// ��������Լ��
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

		// Լ��1.3
		// ÿ���û���Ҫ������
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

		// Լ��1.4
		// ͨ������Լ��
		IloExpr SINR_cons(env_IP);
		for (IloInt i = 0; i < m; i++)
			for (IloInt j = 0; j < n; j++)
				if (SINR[i][j] < SINR_min)
				{
					SINR_cons.clear();
					SINR_cons = x[i][j];
					model.add(SINR_cons == 0);
				}

		// Ŀ�꺯��
		IloExpr obj(env_IP);
		obj.clear();
		for (IloInt i = 0; i < m; i++)
			obj += y[i];

		model.add(IloMinimize(env_IP, obj));

		// ����������
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
	print_cplex(x, y);
	result.int_(x, y);


	DSIS(x, y, x1, y1);
	print_cplex(x1, y1);
	result.int_(x1, y1);
	/*vector<Cluster> CC = COS(x1, y1, x2, y2);
	Result resulth = SFS(x2, y2, x3, y3, CC);*/





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
	vector<int> allA;			// ��¼y����ķ�����id
	vector<int> ou = all_u;		// ���շ����û��ķ�������������������û�id
	vector<int> ucount;			// ÿ���û������ǵĴ���
	vector<int> SS;				// superior server�ļ���
	vector<int> II;				// inferior server�ļ���
	vector<int> Mid;			// yֵ����alpha��1֮��ķ������ļ���
	for (int i = 0; i < m; i++)
	{
		if (y[i] > 0)
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

	// ����ou
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
		// max_id = uid ������j����jj
		int max_id = j;		// �ȼٶ������Ǵ�������idΪou[j]
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
	cout << "r = " << r << '\n';
	cout << "ep = " << ep << '\n';
	cout << "cl = " << cl << '\n';
	for (auto j : ou)
	{
		double sum_I = 0;		// I_j �з�������Ӧyi�ĺ�
		vector<int> I_j;		// ����xij��inferior server�ļ���
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
			// �������alpha<sum_{ai_\in Ij}{y_i}<= 2 * alpha
			// ��construct_I_j�����У�I_j���޸ĳɷ�������Լ���ļ���
			construct_I_j(I_j, sum_I, y1);

			cout << "I_j: ";
			for (auto i : I_j)
				cout << i << ' ';
			cout << '\n';

			// cpΪ��������������ĵ�
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

				// ��GG�е�ÿһ��group���д���
				vector<int> g = it->second;
				int im = g.back();		// ��am����Ϊg�����һ��
				g.pop_back();			// 

				reroute(all_u, g, im, x1);	// �������û�all_u��Ⱥ��g�е���������������im
				y1[im] = 1;
				for (auto i : g)
					y1[i] = 0;


				it++; // ����������
			}
		}
	}
	
	for (auto i : Mid)
		y1[i] = 1;
}

vector<Cluster> Cover::COS(double** x1, double* y1, double** x2, double* y2)
{
	return vector<Cluster>();
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
	// ��obj�е��û���flow��ô���������src��������aim
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

map<int, vector<int>> Cover::construct_GG(vector<int>& I, Point& cp, double L, double cl)
{
	// ����������������I�����򣬸������Ե�cpΪ���ģ��߳�ΪL������Ϊ�߳�Ϊcl��С����
	// ����ֵΪGG����һ����������groups�ļ��ϡ����ݽṹΪmap<int, vector<int>> (��һ��Ϊ������ţ��ڶ���Ϊgroup�����ķ�����)
	// �����б�:
	// I		�������ֵķ���������
	// cp		I������������ĵ�
	// L		I��������ı߳�
	// cl		����߳�
	map<int, vector<int>> GG;
	// �������������½ǵ�����(X0, Y0)
	double X0 = cp.X - L / 2;
	double Y0 = cp.Y - L / 2;
	int lines = ceil(L / cl);


	for (auto i : I)
	{
		Server a = A[i];
		// ii, jj�������a�������������꣨�������ϣ��������ң�ԭ��ΪX0��Y0��
		int ii = int((a.X - X0) / cl);
		int jj = int((a.Y - Y0) / cl);
		// key�������a����������id���������ϣ��������ң�ԭ��ΪX0��Y0��
		int key = ii + jj * lines;

		// ��GG������keyΪ�ؼ��ֵļ�ֵ�ԣ���û�У���GG.find(key)ָ��GG.end()
		if (GG.find(key) == GG.end())
		{
			// δ�ҵ���keyΪ���ļ�ֵ��
			vector<int> g;
			g.push_back(i);
			GG[key] = g;
		}
		else
			GG[key].push_back(i);
	}

	// �ϲ�GG����Щgֻ��һ����������g

	return GG;
}

void Cover::merge_GG(map<int, vector<int>>& GG, double cl)
{
	map<int, vector<int>> GG_ = GG;			// GG�������ں���GG��ĳ��G�ϲ���g֮��֮��Ҫ�����ϲ���gֻ����G.back()�ȽϾ���

	map<int, vector<int>>::iterator it1 = GG.begin();
	map<int, vector<int>>::iterator itEnd = GG.end();

	vector<int> erased_key;	// ��¼��GG�У����ϲ��ĳ���Ϊ1��g��keyֵ����֮�����keyֵ����map��ɾ��


	while (it1 != itEnd)
	{
		// ��ѭ�����ҵ�����Ϊ1��g
		vector<int> g = it1->second;
		// ��¼��ǰg�Ƿ񱻺ϲ���δ���ϲ�Ϊ0
		int is_merge = 0;
		if (g.size() == 1)
		{
			// ֻ�е�ǰ��g�ĳ���Ϊ1���ܱ��ϲ����������²���
			// g��Ӧ��key
			int key = it1->first;
			// g�е�Ψһ��server
			Server a1 = A[g[0]];

			// ��map��ͷΪgѰ���ܹ������G
			map<int, vector<int>>::iterator it2 = GG.begin();
			while (it2 != itEnd) {
				if (it2 == it1)
				{
					// ��it2==it1������it2ָ���˵�ǰg����ʲô������
				}
				else
				{
					// ÿ��it2��Ҫ�ж�g�Ƿ�merge���ϸ�������merge����������ͻ��˳�����
					if (is_merge == 0)
					{
						// ��¼G��Ӧ��KEYֵ
						int KEY = it2->first;
						// G��g���ܼ���Ķ���
						vector<int> G = it2->second;
						// G_��G��ԭ����������ͬ����˵��Gû�м����g
						vector<int> G_ = GG_[KEY];
						// ��������ͬ����˵��Gû�м����g����Ҫ��g��G��ÿһ���������ȽϾ���
						if (G.size() == G_.size())
						{
							// ����G��ii������Ԫ�ص��±�
							for (int ii = 0; ii < G.size(); ii++)
							{
								// iΪ��ǰ�ķ�����id
								int i = G[ii];
								// ��g��A[i]����С��cl���Ϳ���merge
								if (a1.cal_distance(A[i]) <= cl)
								{
									// ��g�з�����id g[0]���浽G��ii�±�����Ḳ��i
									it2->second[ii] = g[0];
									// ��ĩβ����i
									it2->second.push_back(i);
									// ��erased_key�м��뵱ǰ��merge��g�ļ�ֵ
									erased_key.push_back(key);
									is_merge = 1;
									break;
								}
							}
						}
						else
						{
							// ���G���Ѿ�merge������������ô��ֻ����G�е����һ���������ȽϾ���
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
	// ��x��y���������½��д���ļ�
	// ��1��: chosen_m n |chosen_m: ��ѡ��������� n: �û���
	// ��2k+1��: aid ucont |k=0,1,...,chosen_m-1; aid: ������id; ucont: 
	// ��2k+2��: u0 u1 ... u_ucont-1 |u0... : ��2k+1���ж�Ӧ������������û�	
	vector<string> wcontent;
	wcontent.push_back("\n");	// ռλ����һ��Ϊ"chosen_m n"������ѡ��ķ����������û���
	for (int i = 0; i < m; i++)
	{
		if (y[i] > 0)
		{
			chosen_m++;
			string con_u;	// ��i������û�id�ַ���
			int u_cont = 0;	// ��i������û���
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
