#include<stdio.h>
#include<math.h>

long double Fx(long double x)
{
    long double y = x * x * x / 3 - x;
    return y;
}

long double DFx(long double x)
{
    return x * x - 1;
}

long double newton(long double (*Fx)(long double), long double (*DFx)(long double), long double x0, long double epsilon)
{
    long double xi = x0;
    long double xj;
    long double error;
    long double error_back;
    int N = 0;
    do
    {
        xj = xi - Fx(xi) / DFx(xi);
        N++;
        error = (long double)fabs((double)(xi-xj));
        if(N != 0) printf("error[%d+1]/error[%d]^3 = %f\n", N, N, (double)(error/(long double)pow(error_back,3)));
    }
    while(error >= epsilon && (xi = xj) && (error_back = error));
    //printf("初值：x0=%.15f 根：x=%.15f 迭代步数：%d\n", (double)x0, (double)xj, N);
    return xj;
}

long double cutting(long double (*Fx)(long double), long double x0, long double x1, long double epsilon)
{
    long double xi1 = x0;
    long double xi2 = x1;
    long double xi3;
    long double diff;
    long double error;
    long double error_back;
    int N = 0;
    do
    {
        diff =(xi2 - xi1) / (Fx(xi2) - Fx(xi1));
        xi3 = xi2 - Fx(xi2) * diff;
        N++;
        error = (long double)fabs((double)(xi3-xi2));
        if(N != 0) printf("error[%d+1]/error[%d]^1.618 = %f\n", N, N, (double)(error/(long double)pow(error_back,(sqrt(5)+1)/2)));
    }
    while(error >= epsilon && (xi1 = xi2, 1) && (xi2 = xi3, 1) && (error_back = error));
    //printf("初值：x0=%.15f,x1=%.15f 根：x=%.15f 迭代步数：%d\n", (double)x0, (double)x1, (double)xi3, N);
    return xi3;
}

int main()
{
    long double newton_x0[] = {0.1, 0.2, 0.9, 9.0};
    long double cutting_x0[] = {-0.1, -0.2, -2.0, 0.9};
    long double cutting_x1[] = {0.1, 0.2, 0.9, 9.0};
    
    //newton
    printf("newton法：\n");
    for(int i=0; i<4; i++)
    {
        newton(Fx, DFx, newton_x0[i], 1.0e-8);        
        putchar('\n');
    }

    //cutting
    printf("弦截法：\n");
    for(int i=0; i<4; i++)
    {
        cutting(Fx, cutting_x0[i], cutting_x1[i], 1.0e-8);   
        putchar('\n');     
    }

}

