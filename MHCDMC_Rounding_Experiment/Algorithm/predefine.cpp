#include "predefine.h"
#include <iostream>
#include <cmath>

// ���㵱ǰ���������p֮��ľ���
double Point::cal_distance(const Point& p)
{
	const double X1 = this->X;
	const double Y1 = this->Y;

	const double X2 = p.X;
	const double Y2 = p.Y;

	double d = sqrt(pow(X1 - X2, 2) + pow(Y1 - Y2, 2));
	return d;
}

// User��Ĺ��캯������ʼ������(x, y)��ID�Լ���������BR
User::User(double lx, double ly, int id, int br)
{
	X = lx;
	Y = ly;
	ID = id;
	BR = br;
}

// �����ǰ�û��������Ϣ
void User::print_user()
{	
	cout << "u" << ID << ":(" << X << ", " << Y << "), " << BR;
}

// Server��Ĺ��캯������ʼ������(x, y)��ID�Լ���������BW
Server::Server(double lx, double ly, int id, int bw)
{
	X = lx;
	Y = ly;
	ID = id;
	BW = bw;
}

// �����ǰ�������������Ϣ
void Server::print_server()
{
	cout << "a" << ID << ":(" << X << ", " << Y << "), " << BW;
}

// ���������û��������֮��ľ��� dis_h, dis
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

// ���������û��������֮���·������ L
void Cover::cal_L()
{
	double eta = 20 * log10((4 * PI * f) / c) + eta_LOS;
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++)
			L[i][j] = 20 * log10(dis_h[i][j]) + eta;
	}
}

// ���������û��������֮�������� SINR
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

// ������֪���ݣ��������˻�����򸲸ǰ뾶
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

// ���ݹ���p��SINR_min������ϵ�����������˻��뾶������ϵ�� ep
void Cover::cal_ep()
{
	double zhishu = (power * ep_p + SINR_min - ep_SINR) / 20;

	ep = pow(10, zhishu) - 1;
	cout << "ep = " << ep << '\n';
}

// �����ļ�fname����ʼ��һ��Cover����
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

// ��fname�ļ��ж�ȡm��n�Լ����������Ϣ
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
	// ��ȡm, n
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

// ��������y������ÿ���û��ĸ�������
template<class T>
void Cover::cal_NI(T y)
{

	
	// cout << "N = " << N << '\n';
	vector<int> S;				// ��ѡ�еķ�����
	int* uBs = new int [n]; 		// ÿ���û�����������ķ���������

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



// ��������û����������Ϣ
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

// �����Ϊname�Ķ�ά����arr
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

// ���㷨��x��y������
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

// ���ݵ�ǰx,y״̬�������ѡ�з�������ʣ������
template<class TT, class T>
void Cover::print_RB(TT x, T y, int is_y)
{
	// is_y			�����������Ƿ����y���ϵ��		1��ʾ�ˣ�0��ʾ����
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

// ���ݵ�ǰx״̬�������ǰ������i����������
template<class T>
double Cover::get_sumBR(int i, T x)
{
	double sum = 0;
	for (int j = 0; j < n; j++)
		sum += x[i][j] * U[j].BR;
	return sum;
}

// ���ݵ�ǰx״̬�����������i�Ŀ�������
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

// �����ǰ�����Ӧ���������ݣ����� �û����������Ϣ�� ����dis_h��·������L�������SINR 
void Cover::print_all()
{
	print_all_usersServers();
	print_array(dis_h, "d");
	print_array(L, "L");
	print_array(SINR, "SINR");
}

// ���Cluster�ļ��� CC
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
// ���Թ滮
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
				//cons_x_le_y = x[i][j] * U[j].BR - y[i] * A[i].BW;
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
//�����滮
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
				// cons_x_le_y = x[i][j] * U[j].BR - y[i] * A[i].BW;
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

}

// ���ǵ��㷨���
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

// �㷨��һ����ȷ��inferior ��superior servers
void Cover::DSIS(double** x1, double* y1)
{
	vector<int> allA;			// ��¼y����ķ�����id
	vector<int> ou = all_u;		// ���շ����û��ķ�������������������û�id
	vector<int> ucount;			// ÿ���û������ǵĴ���
	vector<int> SS;				// superior server�ļ���
	vector<int> II;				// inferior server�ļ���
	vector<int> Mid;			// yֵ����alpha��1֮��ķ������ļ���

	// ����II��SS��Mid
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

	// ����ou
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
	// �����û������ǵĴ������û�����
	//for (int j = 0; j < n; j++)
	//{
	//	// max_id = uid ������j����jj
	//	int max_id = j;		// �ȼٶ������Ǵ�������idΪou[j]
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
			// ��GG�У�ֻ��һ����������G�ϲ����ܱ��ϲ���G
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

// �㷨�ڶ���������inferior������
map<int, Cluster> Cover::COS(double** x1, double* y1, double** x2, double* y2)
{
	for (int i = 0; i < m; i++)
	{
		y2[i] = y1[i];
		for (int j = 0; j < n; j++)
			x2[i][j] = x1[i][j];
	}
	vector<int> allu = all_u;
	vector<int> allA;			// ���з�����id����
	vector<int> SS;				// superior server�ļ���
	vector<int> II;				// inferior server�ļ���
	// inferior server ������ superior server�����������򣬱������
	// ISorder[2] = {3, 5, 6, 1} �ĺ����ǣ�idΪ2�ķ��������ɽ���Զ�����˻�idΪ3,5,6,1
	map<int, vector<int>> ISorder;
	map<int, vector<int>> SIorder;
	vector<int> OO;				// ��O����
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
					// ���ⲽ�裬����ǰx2[t][j]��״̬������x1[t][j]��
					for (int j = 0; j < n; j++)
						x1[t][j] = x2[t][j];

					vector<int> src;
					src.push_back(t);
					reroute(allu, src, i, x2);
					CC[i].push_back(t);

					is_clus[tt] = t;
					// ��II��ɾ��t
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
		//				// ���ⲽ�裬����ǰx2[t][j]��״̬������x1[t][j]��
		//				for (int j = 0; j < n; j++)
		//					x1[t][j] = x2[t][j];
		//				vector<int> src;
		//				src.push_back(t);
		//				reroute(allu, src, i, x2);
		//				CC[i].push_back(t);
		//				// ��II��ɾ��t
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
			// ����II���������˻���Ӧ��Ai��AiΪai������û���id����
			map<int, vector<int>> AA;
			// Ki����II���������˻���Ӧ��ki
			vector<double> Ki;
			for (auto i : II)
			{
				// ��II�е����˻�ai��������������û���id����Ai
				vector<int> Ai;
				for (int j = 0; j < n; j++)
					if (x2[i][j] > 0)Ai.push_back(j);
				AA[i] = Ai;

				// ����ai��Ӧ��ki
				double ki = 0;
				for (auto j : Ai)
					ki += U[j].BR;
				if (A[0].BW < ki)
					ki = A[0].BW;
				Ki.push_back(ki);
			}

			// ��II���ҵ�ki��С�����˻�II[tt]����С��kiֵ=Ki[tt]
			int tt = 0;
			for (int ii = 0; ii < II.size(); ii++)
				if (Ki[ii] > Ki[tt])
					tt = ii;
			// at = II[tt]
			int t = II[tt];
			double kt = Ki[tt];
			OO.push_back(t);
			// ��II��ɾ��t
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

// ���ݵ�ǰx��y״̬����һЩ���Ա�����inferior������ȡ���ķ������ϲ�
void Cover::mergeII(double** x, double* y, double** x1, double* y1)
{
	// ��ʼ��
	for (int i = 0; i < m; i++)
	{
		y1[i] = y[i];
		for (int j = 0; j < n; j++)
		{
			// cout << "i = " << i << ", j = " << j << '\t';
			x1[i][j] = x[i][j];
		}

	}
	vector<int> II;				// inferior server�ļ���
	map<int, vector<int>> AA;	// ���������ѡ�з�����������û� 
	for (int i = 0; i < m; i++)
	{
		if (y1[i] > 0)
		{
			// ����II
			if (y1[i] < 1 and y1[i] > 0)
				II.push_back(i);

			// ���ɱ�ѡ�з�����������û�����
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

	// �ϲ��ܹ���ȡ����inferior server
	// is_mergeΪһ������ΪII.size()�����飬��¼��Ӧ�±�������Ƿ�merge
	// 0Ϊ��δ��merge��		1Ϊ��merge
	int* is_merge = new int[II.size()];
	for (int i = 0; i < II.size(); i++)
		is_merge[i] = 0;

	for (int ii = 0; ii < II.size(); ii++)
	{
		// ���±�Ϊii�ķ�������merge�ˣ����¸�ѭ��
		if (is_merge[ii] == 1)
			continue;
		// iΪ��ǰii��Ӧ�ķ�����id
		int i = II[ii];
		Server ai = A[i];
		// AiΪ������i��ǰ������û�
		vector<int> Ai = AA[i];
		for (int tt = 0; tt < II.size(); tt++)
		{
			// �ϲ��±�Ϊtt�ķ�������ii�ķ�����
			// ���±�Ϊtt�ķ�������merge�ˣ����¸�ѭ��
			if (is_merge[tt] == 1)
				continue;
			if (ii == tt)
				continue;
			int t = II[tt];
			Server at = A[t];
			vector<int> At = AA[t];

			if (ai.cal_distance(at) > 2 * r)
				continue;
			// �ںϲ��Ĺ����У�II�еķ��������ܱ��superior�������Ҫ���ж�
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
			// ����II
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

// �㷨����������ÿ��cluster��ѡ�����յķ�����
void Cover::SFS(double** x1, double** x2, double* y2, double** x3, double* y3, map<int, Cluster>& CC)
{
	// ������
	// x1, y1		�ݴ�ĳ�����󱻾��ൽĳ�������cluster֮ǰ��xֵ
	// x2, y2		COS֮���x��yֵ
	// x3, y3		SFS֮���x��yֵ

	for (int i = 0; i < m; i++)
	{
		y3[i] = y2[i];
		for (int j = 0; j < n; j++)
			x3[i][j] = x2[i][j];
	}
	vector<int> allu = all_u;
	vector<int> allA;			// ���з�����id����
	vector<int> SS;				// superior server�ļ���
	vector<int> II;				// inferior server�ļ���
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
			// ��C�������Ӧinferior server��Ӧ��xֵ�ָ�������C֮ǰ��״̬
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

// �������alpha<sum_{ai_\in Ij}{y_i}<= 2 * alpha
// ��construct_I_j�����У�I_j���޸ĳɷ�������Լ���ļ��� 
void Cover::construct_I_j(vector<int>& I_j, double sum_I, double* y)
{
	while (sum_I > 2 * alpha)
	{
		int i = I_j.back();
		sum_I -= y[i];
		I_j.pop_back();
	}
}

// ��obj�е��û���������src��reroute��aim�������У���Ҫ�޸�x��flowΪ���Ĵ�С��Ĭ��Ϊ1������Ϊ1��
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

// �ҵ�i������û�����Ai�У�������С��BR 
double Cover::get_BRmin(vector<int>& Ai)
{
	double min = A[0].BW;
	for (auto j : Ai)
		if (min > U[j].BR)
			min = U[j].BR;
	return min;
}

// �����˻�t��ǰ������û�����At�У��ҵ����ʺ��Ǹ��û�������û�����������
int Cover::get_max_xBR(int t, vector<int>& At, double RBt, double** x)
{
	// �����˻�t��ǰ������û�����At�У��ҵ����ʺ��Ǹ��û�������û�����������
	// t: ��ǰ���˻�id
	// At		���˻�t��ǰ������û�����(x[t][j]>0)
	// RBt		���˻�t��ǰ��ʣ������
	// x
	//	����ֵ��max_u			��max_u==-1,��˵��û�к��ʵ��û�
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

// �ж�At�Ƿ���Ai���Ӽ�
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

// inferior server ������ superior server�����������򣬱������
map<int, vector<int>> Cover::get_ISorder(vector<int>& SS, vector<int>& II)
{
	// inferior server ������ superior server�����������򣬱������
	// ISorder[2] = {3, 5, 6, 1} �ĺ����ǣ�idΪ2�ķ��������ɽ���Զ�����˻�idΪ3,5,6,1

	map<int, vector<int>> ISorder;

	// ��ʼ��
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

	// cout << "���̣�\n";
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

// superior server ������ inferior server�����������򣬱������
map<int, vector<int>> Cover::get_SIorder(vector<int>& SS, vector<int>& II)
{
	// inferior server ������ superior server�����������򣬱������
	// ISorder[2] = {3, 5, 6, 1} �ĺ����ǣ�idΪ2�ķ��������ɽ���Զ�����˻�idΪ3,5,6,1

	map<int, vector<int>> SIorder;

	// ��ʼ��
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

	// cout << "���̣�\n";
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

// ����������������I�����򣬸������Ե�cpΪ���ģ��߳�ΪL������Ϊ�߳�Ϊcl��С����
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

// ��һЩ�ʰ���һ����������G���ϲ�������������G��
void Cover::merge_GG(map<int, vector<int>>& GG, double cl)
{
	map<int, vector<int>> GG_ = GG;			// GG�������ں���GG��ĳ��G�ϲ���g֮��֮��Ҫ�����ϲ���gֻ����G.back()�ȽϾ���

	map<int, vector<int>>::iterator it1 = GG.begin();
	map<int, vector<int>>::iterator itEnd = GG.end();

	vector<int> erased_key;	// ��¼��GG�У����ϲ��ĳ���Ϊ1��g��keyֵ����֮�����keyֵ����map��ɾ��
	// ��¼ÿ���Ƿ񱻺ϲ������ϲ��Ļ�ֵΪ���ϲ���keyֵ
	int* is_mer = new int[GG.size()];
	for (int i = 0; i < GG.size(); i++)
		is_mer[i] = 0;

	int j = 0;
	while (it1 != itEnd)
	{
		cout << j << '\n';
		if (is_mer[j] == 1)
			break;
		// ��ѭ�����ҵ�����Ϊ1��g
		vector<int> g = it1->second;
		// ��¼��ǰg�Ƿ񱻺ϲ���δ���ϲ�Ϊ0
		if (g.size() == 1)
		{
			// ֻ�е�ǰ��g�ĳ���Ϊ1���ܱ��ϲ����������²���
			// g��Ӧ��key
			int key = it1->first;
			// g�е�Ψһ��server
			Server a1 = A[g[0]];

			// ��map��ͷΪgѰ���ܹ������G��it2ָ��it1��Ҫ�����group
			map<int, vector<int>>::iterator it2 = GG.begin();
			int jj = 0;
			while (it2 != itEnd) {
				// ÿ��it2��Ҫ�ж�g�Ƿ�merge���ϸ�������merge����������ͻ��˳�����
				if (is_mer[j] == 1)
					break;
				// �ж�it2ָ���group�Ƿ�merge������merge��ɨ����һ��
				if (is_mer[jj] == 1)
				{

				}
				else {
					if (it2 == it1)
					{
						// ��it2==it1������it2ָ���˵�ǰg����ʲô������
					}
					else
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
									is_mer[j] = 1;
									break;
								}
							}
						}
						else
						{
							// ���G���Ѿ�merge������������ô��ֻ����G�е����һ���������ȽϾ���
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

// ����xx��yy��ֵ���Զ����ʼ��
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


// ��x��y���������½��д���ļ�
// ��1��: chosen_m n |chosen_m: ��ѡ��������� n: �û���
// ��2k+1��: aid ucont |k=0,1,...,chosen_m-1; aid: ������id; ucont: 
// ��2k+2��: u0 u1 ... u_ucont-1 |u0... : ��2k+1���ж�Ӧ������������û�	
void Result::write_result_file(string fname)
{
	// ��x��y���������½��д���ļ�
	// ��1��: chosen_m n |chosen_m: ��ѡ��������� n: �û���
	// ��2k+1��: aid ucont |k=0,1,...,chosen_m-1; aid: ������id; ucont: 
	// ��2k+2��: u0 u1 ... u_ucont-1 |u0... : ��2k+1���ж�Ӧ������������û�	
	vector<string> wcontent;
	wcontent.push_back("\n");	// ռλ����һ��Ϊ"chosen_m n"������ѡ��ķ����������û���
	chosen_m = 0;
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

void Result::write_cplex_file(string fname)
{
	// ��x��y��ֵд���ļ�
	// ��һ�У�choosen_m n
	// �ڶ��У�yid y[i]
	// �����У�x[0] x[1] ... x[n]
	vector<string> wcontent;
	wcontent.push_back("\n");	// ռλ����һ��Ϊ"chosen_m n"������ѡ��ķ����������û���
	chosen_m = 0;
	for (int i = 0; i < m; i++)
	{
		if (y[i] > 0)
		{
			wcontent.push_back(to_string(i) + ' ' + to_string(y[i]) + '\n');
			chosen_m++;
			string x_values;	// ��Ӧ����x[i]��ֵ
			int u_cont = 0;		// ��i������û���
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
