#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <complex>
#include <iomanip>
#include <fstream>

using namespace std;

vector<complex<double>> FFT(double (*)(double), vector<complex<double>>);
double f(double);
complex<double> UnitRoot(double);

int main()
{
    //打开文件
    ofstream file("output.txt");

    vector<complex<double>> fn,gn;
    complex<double> tmp(0,0);
    int n[] = {128, 256}; //采样点个数
    //生成采样点
    for(int i = 0; i < 2; i++)
    {
        fn.clear();
        gn.clear();
        for(int j = 0; j < n[i]; j++)
        {
            tmp.real(f((double)j/n[i]));
            fn.push_back(tmp);
        }
        //计算
        gn = FFT(f, fn);

        //输出
        /* cout.precision(16);
        cout << "采样点为" << n[i] << "时：" << endl;
        for(int j = 0; j < n[i]; j++)
        {
            cout << fixed << "向量g的第" << j << "个分量, x_" << j << "=" << gn[j].real() << ", y_" << j << "=" << gn[j].imag() << endl;
        }
        cout << endl; */

        file.precision(16);
        file << "采样点为" << n[i] << "时：" << endl;
        for(int j = 0; j < n[i]; j++)
        {
            file << fixed << "向量g的第" << j << "个分量, x_" << j << "=" << gn[j].real() << ", y_" << j << "=" << gn[j].imag() << endl;
        }
        file << endl;
    }

    file.close();
    system("pause");
}

inline double f(double t)
{
    return 0.7*sin(2*M_PI * 2*t) + sin(2*M_PI * 5*t ) ;
}

vector<complex<double>> FFT(double (*f)(double), vector<complex<double>> fn)
{
    //递归边界
    if(fn.size() == 1) return fn;

    //生成分治数组
    vector<complex<double>> f0, f1;
    for(int i = 0; i < fn.size(); i++)
    {
        if(i%2 == 0) f0.push_back(fn[i]);
        else f1.push_back(fn[i]);
    }

    //分治
    vector<complex<double>> g0, g1;
    g0 = FFT(f, f0);
    g1 = FFT(f, f1);

    //计算
    vector<complex<double>> gn; //gn为最终输出
    double theta = 0;
    for(int i = 0; i < fn.size() / 2; i++)
    {
        gn.push_back((g0[i] + UnitRoot(theta) * g1[i]) / (complex<double>)2);
        theta -= 2*M_PI / fn.size(); 
    }
    for(int i = 0; i < fn.size() / 2; i++)
    {
        gn.push_back((g0[i] + UnitRoot(theta) * g1[i]) / (complex<double>)2);
        theta -= 2*M_PI / fn.size(); 
    }

    return gn;
}

inline complex<double> UnitRoot(double theta)
{
    complex<double> a(cos(theta), sin(theta));
    return a;
}