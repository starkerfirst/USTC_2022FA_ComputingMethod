#include<stdio.h>
#include<math.h>

double Fx(double x)
{
    double y = 1/(1+x*x);
    return y;
}

double Lagrange(double* x, int n, double input)
{
    double output = 0;
    double temp;
    
    for (int i = 0; i <= n; i++)
    {
        temp = 1;
        for(int j = 0; j <= n; j++)
            if(i != j)
                temp *= (input-x[j])/(x[i]-x[j]);
        output += temp*Fx(x[i]);
    }
    
    return output;
}

double ErrorCalculate(double* x, int n)
{
    double error = 0;
    double temp;
    double test[501];

    //数据准备
    for(int i=0; i<=500; i++)  test[i] = -5.0+10.0*i/500.0;

    //最大值查找
    for(int i=0; i<=500; i++)
    {
        temp = (Fx(test[i]) - Lagrange(x, n, test[i])) > 0 ? (Fx(test[i]) - Lagrange(x, n, test[i])) : (Lagrange(x, n, test[i]) - Fx(test[i]));
        if(temp > error)  error = temp;
    }
       
    return error;
}

int main()
{
    int n[4] = {5, 10, 20, 40}; //节点数目选取
    double x1[41];
    double x2[41];
    int j;

    for(int i=0; i<4; i++)
    {
        for(j=0; j<=n[i]; j++) //数据准备
        {
            x1[j] = -5.0+10.0*j/n[i];
            x2[j] = -5.0*cos((2*j+1)*M_PI/(2*n[i]+2));
        }

        printf("节点个数%d,第一组error=%.15lf,第二组error=%.15lf\n", n[i] + 1, ErrorCalculate(x1, n[i]), ErrorCalculate(x2, n[i])); //误差计算与输出
    }

    system("pause");
}



