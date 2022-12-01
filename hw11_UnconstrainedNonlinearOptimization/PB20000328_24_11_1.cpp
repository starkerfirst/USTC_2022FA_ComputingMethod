#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <fstream>

using namespace std;

//定义2维向量
class Vector2D
{
public:
    double x;
    double y;
    Vector2D()
    {
        x = y = 0;
    }
    inline Vector2D operator+(Vector2D a)
    {
        Vector2D sum;
        sum.x = x + a.x;
        sum.y = y + a.y;
        return sum;
    }
    void operator= (const Vector2D& a)
    {
        x = a.x;
        y = a.y;
    }
    inline Vector2D operator- ()
    {
        Vector2D sum;
        sum.x = -x;
        sum.y = -y;
        return sum;
    }
    inline Vector2D operator*(double a)
    {
        Vector2D sum;
        sum.x = x * a;
        sum.y = y * a;
        return sum;
    }
};


void SteepestDescent(double (*)(double, double), Vector2D (*)(double, double));
void Newton(double (*)(double, double), Vector2D (*)(double, double), Vector2D (*)(double, double));
void Search(double (*)(double, double), Vector2D, Vector2D, double&);
double f(double, double);
Vector2D df(double, double);
Vector2D d2f(double, double);

int main()
{
    SteepestDescent(f, df);
    //cout << endl;
    //cout << endl; 
    Newton(f, df, d2f);
    system("pause");
}

inline double f(double x1, double x2)
{
    return 100 * pow(x2 - x1*x1, 2) + pow(1 - x1, 2);
}

inline Vector2D df(double x1, double x2)
{
    Vector2D grad;
    grad.x = (400 * x1 * (x1*x1 - x2) + 2 * (x1 - 1));
    grad.y = (200 * (x2 - x1*x1));
    return grad;
}

inline Vector2D d2f(double x1, double x2)
{
    //d2f矩阵的逆
    Vector2D pn;
    Vector2D gn = df(x1, x2);
    double d2f[2][2];
    double inv[2][2];

    d2f[0][0] = 1200 * x1 * x1 - 400 * x2 + 2;
    d2f[0][1] = -400 * x1;
    d2f[1][0] = -400 * x1;
    d2f[1][1] = 200;
    inv[0][0] = d2f[1][1] / (d2f[0][0]*d2f[1][1]-d2f[1][0]*d2f[0][1]);
    inv[0][1] = -d2f[0][1] / (d2f[0][0]*d2f[1][1]-d2f[1][0]*d2f[0][1]);
    inv[1][0] = -d2f[1][0] / (d2f[0][0]*d2f[1][1]-d2f[1][0]*d2f[0][1]);
    inv[1][1] = d2f[0][0] / (d2f[0][0]*d2f[1][1]-d2f[1][0]*d2f[0][1]);
    pn.x = -inv[0][0]*gn.x-inv[0][1]*gn.y;
    pn.y = -inv[1][0]*gn.x-inv[1][1]*gn.y;
    return pn;
}
void SteepestDescent(double (*f)(double, double), Vector2D (*Df)(double, double))
{
    ofstream file("output.txt");
    file << "最速下降法：" << endl;
    file.precision(16);
    Vector2D x;
    Vector2D pn=-Df(x.x, x.y),p;
    double lambda;
    int i = 0;
    while(sqrt(pow(pn.x,2)+pow(pn.y,2)) >= 1.0e-4 )
    {
        Search(f, x, pn, lambda);
        x = x + pn * lambda;
        //cout << i << "   "<< pn.x*p.x+pn.y*p.y <<endl;
        i++;
        file << fixed << "第" << i << "次迭代, f(x_" << i << ") = " << f(x.x,x.y) << ", x1 = " << x.x << ", x2 = " << x.y << endl;
        p = pn;
        pn = -Df(x.x, x.y); 
    }
}

void Newton(double (*f)(double, double), Vector2D (*df)(double, double), Vector2D (*d2f)(double, double))
{
    ofstream file("output.txt", ofstream::app);
    file << "牛顿法：" << endl;
    file.precision(16);
    Vector2D x;
    Vector2D pn = -d2f(x.x, x.y);
    Vector2D gn = df(x.x, x.y);
    double lambda;
    int i = 0;
    while(sqrt(pow(gn.x,2)+pow(gn.y,2)) >= 1.0e-4 )
    {  
        Search(f, x, pn, lambda);
        x = x + pn * lambda;
        //cout << i << "   "<< pn.x << "   "<< pn.y <<endl;
        i++;
        file << fixed << "第" << i << "次迭代, f(x_" << i << ") = " << f(x.x,x.y) << ", x1 = " << x.x << ", x2 = " << x.y << endl;
        pn = -d2f(x.x, x.y);
        gn = df(x.x, x.y);
    }
}

void Search(double (*)(double, double), Vector2D x, Vector2D pn, double& lambda)
{
    //进退法,初始步长设为0.0001
    double left, right;
    double h = 0.0001;
    double a = 0, b = 0, c = 0;
    int i=0;
    Vector2D x2, x1;
    while(1)
    {
        a = b + h;
        x2 = x + pn * a;
        x1 = x + pn * b;
        if(f(x2.x,x2.y) < f(x1.x, x1.y))
        {
            h = 2 * h;
            c = b;
            b = a;
            i++;
        }
        else 
        {
            if(i==0) 
            {
                h = -h;
                b = a;
                i++;
            }
            else
            {
                left = a>c? c :a;
                right = a<=c? c :a;
                break;
            }
        }
    }
    //黄金分割法，分割区间至1.0e-6
    double alpha,beta,fa,fb;
    alpha = left + 0.382 * (right - left);
    beta = left + 0.618 * (right - left);
    x2 = x + pn * alpha;
    x1 = x + pn * beta;
    fa = f(x2.x,x2.y);
    fb = f(x1.x,x1.y);
    while(right - left > 1.0e-6)
    {
        if(fa > fb)
        {
            left = alpha;
            alpha = beta;
            fa = fb;
            beta = left + 0.618 * (right - left); 
            x1 = x + pn * beta;
            fb = f(x1.x,x1.y);
        }
        else
        {
            right = beta;
            beta = alpha;
            fb = fa;
            alpha = left + 0.382 * (right - left); 
            x2 = x + pn * alpha;
            fa = f(x2.x,x2.y);
        }
    }

    //输出
    if(fa > fb) lambda = beta;
    else lambda = alpha;
}


