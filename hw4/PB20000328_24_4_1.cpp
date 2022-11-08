#include<iostream>
#include<cmath>
using namespace std;

long double simpson(double (*Fx)(double), long double left, long double right, int l)
{
    long double I = 0;
    long double x[4097];
    long double y[4097];
    for(int i=0; i<=pow(2,l); i++) 
        {
            x[i] = left + (right - left) * i / pow(2,l);
            y[i] = (long double)Fx(x[i]);
        }
    I += 1.0 / 3 * (y[0] + y[(int)(pow(2,l))]); 
    for(int i=1; i<(int)(pow(2,l)); i++)
    {
        if(i % 2 == 0)    I += 2.0 / 3 * y[i];
        else    I += 4.0 / 3 * y[i];
    }
    I *= 4/pow(2,l);
    return I;
}

long double trapezoid(double (*Fx)(double), long double left, long double right, int l)
{
    long double I = 0;
    long double x[4097];
    long double y[4097];
    for(int i=0; i<=(int)pow(2,l); i++) 
        {
            x[i] = left + (right - left) * i / pow(2,l);
            y[i] = (long double)Fx(x[i]);
        }
    I += 0.5 * (y[0] + y[(int)(pow(2,l))]); 
    for(int i=1; i<(int)(pow(2,l)); i++)
    {
        I += y[i];
    }
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
        Is[j] = simpson(sin, 1, 5, l[j]);
        It[j] = trapezoid(sin, 1, 5, l[j]);
        es[j] = Ir - Is[j];
        et[j] = Ir - It[j];
        if(j != 0)
        {
            os = log((double)(es[j-1]/es[j]))/log(2);
            ot = log((double)(es[j-1]/es[j]))/log(2);
            //printf("l=%d,simpson:%.15Lf,误差e=%.15Lf,误差阶o=%.15Lf\ntrapezoid:%.15Lf,误差e=%.15Lf,误差阶o=%.15Lf\n\n", l[j], Is[j], es[j], os, It[j], et[j], ot);
            //cout << "l=" << l[j] << ,simpson:
        }
        //else printf("l=%d,simpson:%.15Lf,误差e=%.15Lf\ntrapezoid:%.15Lf,误差e=%.15Lf\n\n", l[j], Is[j], es[j], It[j], et[j]);
    }
}

