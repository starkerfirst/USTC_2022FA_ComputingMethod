#include<stdio.h>
#include<math.h>

long double simpson(double (*Fx)(double), long double left, long double right, int l)
{
    long double I = 0;
    long double x[4097];
    long double y[4097];
    //节点初始化
    for(int i=0; i<=pow(2,l); i++) 
        {
            x[i] = left + (right - left) * i / pow(2,l);
            y[i] = (long double)Fx(x[i]);
        }
    //边界
    I += 1.0 / 3 * (y[0] + y[(int)(pow(2,l))]); 
    //内部
    for(int i=1; i<(int)(pow(2,l)); i++)
    {
        if(i % 2 == 0)    I += 2.0 / 3 * y[i];
        else    I += 4.0 / 3 * y[i];
    }
    //乘上宽度
    I *= 4/pow(2,l);
    return I;
}

long double trapezoid(double (*Fx)(double), long double left, long double right, int l)
{
    long double I = 0;
    long double x[4097];
    long double y[4097];
    //节点初始化
    for(int i=0; i<=(int)pow(2,l); i++) 
        {
            x[i] = left + (right - left) * i / pow(2,l);
            y[i] = (long double)Fx(x[i]);
        }
    //边界
    I += 0.5 * (y[0] + y[(int)(pow(2,l))]); 
    //内部
    for(int i=1; i<(int)(pow(2,l)); i++)
    {
        I += y[i];
    }
    //乘上宽度
    I *= 4/pow(2,l);
    return I;
}

int main()
{
    long double Ir = cos(1) - cos(5);
    long double Is[12];
    long double It[12];
    long double es[12];
    long double et[12];
    long double os;
    long double ot;
    int l[12] = {1,2,3,4,5,6,7,8,9,10,11,12};
    
    for(int j=0; j<12; j++)
    {
        //积分计算
        Is[j] = simpson(sin, 1, 5, l[j]);
        It[j] = trapezoid(sin, 1, 5, l[j]);

        //误差
        es[j] = Ir - Is[j];
        et[j] = Ir - It[j];
        if(j != 0)
        {
            os = log((double)(es[j-1]/es[j]))/log(2);
            ot = log((double)(et[j-1]/et[j]))/log(2);
            printf("l=%d,simpson:%.15f,误差e=%.15e,误差阶o=%.15f\n   trapezoid:%.15f,误差e=%.15e,误差阶o=%.15f\n\n", l[j], (double)Is[j], (double)es[j], (double)os, (double)It[j], (double)et[j], (double)ot);
        }
        else printf("l=%d,simpson:%.15f,误差e=%.15e\n   trapezoid:%.15f,误差e=%.15e\n\n", l[j], (double)Is[j], (double)es[j], (double)It[j], (double)et[j]);
    }
    //输出
    
    system("pause");
}

