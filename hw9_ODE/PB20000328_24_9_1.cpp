#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>

using namespace std;

vector<double> RungeKutta(double f(double,double), double left, double right, double step, double f0);
vector<double> AdamsImplicit(double f(double,double), double left, double right, double step, double f0);
//vector<double> AdamsExplicit(vector<double> y, double f(double,double), double left, double right, double h)
double f(double,double);
double Fx(double);



int main()
{
    vector<int> step = {0,1,2,3};
    vector<double> error;
    double tmp;

    //四阶Runge-Kutta公式
    //误差计算
    for(int l : step) error.push_back(fabs(RungeKutta(f, 0, 1.5, 0.1/pow(2,l), 3).back() - Fx(1.5)));
    //输出
    cout << "四阶Runge-Kutta公式的误差和误差阶：" << endl;
    int j = 0;
    for(auto i : error) 
    {
        cout.precision(16);
        if(i==error.front()) cout << "h = "<< 0.1/pow(2,step[j]) <<", err = " << i << endl;
        else cout << "h = "<< 0.1/pow(2,step[j]) <<", err = " << i << ", ok = " << log(tmp/i)/log(2) <<endl;
        tmp = i;
        j++;
    }
    error.clear();

    //四阶隐式Adams公式
    //误差计算
    for(int l : step) error.push_back(fabs(AdamsImplicit(f, 0, 1.5, 0.1/pow(2,l), 3).back() - Fx(1.5)));
    //输出
    cout << "四阶隐式Adams公式的误差和误差阶：" << endl;
    j = 0;
    for(auto i : error) 
    {
        cout.precision(16);
        if(i==error.front()) cout << "h = "<< 0.1/pow(2,step[j]) <<", err = " << i << endl;
        else cout << "h = "<< 0.1/pow(2,step[j]) <<", err = " << i << ", ok = " << log(tmp/i)/log(2) <<endl;
        tmp = i;
        j++;
    }
    cout << "可以看出符合理论预期，整体截断误差为4阶。" << endl;

    system("pause");
}

double f(double x, double y)
{
    return -x*x*y*y;
}

double Fx(double x) 
{   
    return 3/(1+pow(x,3));
} 

vector<double> AdamsImplicit(double f(double,double), double left, double right, double h, double f0)
{
    vector<double> range;
    vector<double> y;
    double y_next,y_predicted;
    int j = 3;

    //定义计算范围
    for(int i = 0; i <= (int)((right-left)/h); i++) range.push_back(left+i*h);

    //R-K起3步
    y = RungeKutta(f, left, left+3*h, h, f0);
    //校准x4 
    auto i = range.begin() + 3;
    y_next = y[j-1] + h / 24 * (9 * f(*i,y[j]) + 19 * f(*(i-1),y[j-1]) - 5 * f(*(i-2),y[j-2]) + f(*(i-3),y[j-3]));
    y[j] = y_next;
    //校准后续数据
    
    for(; i!=range.end()-1; i++, j++)
    {
        //预估
        y_predicted = y[j] + h / 24 * (55 * f(*i,y[j]) - 59 * f(*(i-1),y[j-1]) + 37 * f(*(i-2),y[j-2]) - 9 * f(*(i-3),y[j-3]));
        //校准
        y_next = y[j] + h / 24 * (9 * f(*(i+1),y_predicted) + 19 * f(*i,y[j]) - 5 * f(*(i-1),y[j-1]) + f(*(i-2),y[j-2]));
        y.push_back(y_next);
    }
    return y;
    
}

vector<double> RungeKutta(double f(double,double), double left, double right, double h, double f0)
{
    vector<double> range;
    vector<double> y;
    double k1,k2,k3,k4,y_next;
    y.push_back(f0);
    //定义计算范围
    for(int i = 0; i <= (int)((right-left)/h); i++) range.push_back(left+i*h);

    //计算
    int j = 0;
    for(auto i = range.begin(); i!=range.end()-1; i++, j++)
    {
        k1 = f(*i,y[j]);
        k2 = f(*i+1/2.0*h, y[j]+1/2.0*h*k1);
        k3 = f(*i+1/2.0*h, y[j]+1/2.0*h*k2);
        k4 = f(*i+h, y[j]+h*k3);
        y_next = y[j] + h/6*(k1+2*k2+2*k3+k4);
        y.push_back(y_next);
    }
    return y;
}