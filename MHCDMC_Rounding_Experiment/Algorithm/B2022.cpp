// B2022.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//
#include <ilcplex/ilocplex.h>
#include <stdio.h>
#include <iostream>
#include "predefine.h"
using namespace std;
typedef IloArray <IloNumVarArray> IloNumVarArray2;
#include <iostream>

void one_example();


int main()
{
    one_example();
    /*Cover cover;
    cover.cal_min_r();*/
}


void one_example()
{
    string fname = "D:\\Myschool\\graduate_school\\02Graduate\\Research\\My paper\\2_Papers\\005_MHCDMC_Rounding\\MHCDMC_Rounding_Experiment\\Generate_Points\\data\\example\\oneInstance.txt";
    Cover cover;
    cover.initial(fname);
    //cover.print_all();
    // cout << "epsilon=" << cover.ep;
    cover.GBTSR();
    //cover.IP();
    //cover.result.write_result_file("D:\\Myschool\\graduate_school\\02Graduate\\Research\\My paper\\2_Papers\\005_MHCDMC_Rounding\\MHCDMC_Rounding_Experiment\\Algorithm\\result\\oneresultDSIS.txt");
}






// 运行程序: Ctrl + F5 或调试 >“开始执行(不调试)”菜单
// 调试程序: F5 或调试 >“开始调试”菜单

// 入门使用技巧: 
//   1. 使用解决方案资源管理器窗口添加/管理文件
//   2. 使用团队资源管理器窗口连接到源代码管理
//   3. 使用输出窗口查看生成输出和其他消息
//   4. 使用错误列表窗口查看错误
//   5. 转到“项目”>“添加新项”以创建新的代码文件，或转到“项目”>“添加现有项”以将现有代码文件添加到项目
//   6. 将来，若要再次打开此项目，请转到“文件”>“打开”>“项目”并选择 .sln 文件
//ILOSTLBEGIN
//void test1() {
//    IloEnv env;
//    try {
//        IloModel model(env);
//        IloNumVarArray vars(env);
//        vars.add(IloNumVar(env, 0.0, 40.0)); //0 <= x1 <= 40
//        vars.add(IloNumVar(env)); // 0 <= x2
//        vars.add(IloNumVar(env)); // 0 <= x3
//        //maximize x1 + 2 x2 + 3 x3
//        model.add(IloMaximize(env, vars[0] + 2 * vars[1] + 3 * vars[2]));
//        model.add(-vars[0] + vars[1] + vars[2] <= 20);
//        model.add(vars[0] - 3 * vars[1] + vars[2] <= 30);
//        IloCplex cplex(model);
//        if (!cplex.solve()) {
//            env.error() << "Failed to optimize LP." << endl;
//            throw(-1);
//        }
//        IloNumArray vals(env);
//        env.out() << "Solution status = " << cplex.getStatus() << endl;
//        env.out() << "Solution value = " << cplex.getObjValue() << endl;
//        cplex.getValues(vals, vars);
//        env.out() << "Value = " << vals << endl;
//    }
//    catch (IloException& e) { cerr << "Concert exception caught:" << e << endl; }
//    catch (...) { cerr << "Unknuwn exception caught" << endl; }
//    env.end();
//}