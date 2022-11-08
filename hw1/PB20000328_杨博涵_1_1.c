#include<stdio.h>
#include<math.h>

int main()
{
    double x[7] = {0.0, 0.5, 1.0, sqrt(2), 10.0, 100.0, 300.0};
    for(int j = 0; j < 7; j ++)
    {
        double sum=0;
        for(int i=1; i<1e7; i++)  //要求截断误差小于1e-6,经过放缩计算以后取最大i值为1e7
        {
            sum += 1.0 / i / (i + x[j]);
        }
        printf("x=%.15f,y=%.15lf\n", x[j], sum);
    } 
}

