#include<stdio.h>
#include<math.h>

int main()
{
    //节点初始化
    double x[] = {0.25, 0.50, 0.75, 1.00, 1.25, 1.50, 1.75, 2.00, 2.25, 2.50};
    double y[] = {1.284, 1.648, 2.117, 2.718, 3.427, 2.798, 3.534, 4.456, 5.465, 5.894};
    double X[2];
    double A[10][2];
    double A_tA[2][2];
    double A_ty[2][1];
    double Q=0;
    
    //A矩阵生成
    for(int i=0; i<10; i++) 
        {
            A[i][0] = sin(x[i]);
            A[i][1] = cos(x[i]);
        }

    //矩阵计算
    for(int i=0; i<2; i++) //A_tA
        for(int j=0; j<2; j++)
        {   
            double sum = 0;
            for(int k=0; k<10; k++)
            {
                sum += A[k][i] * A[k][j];
            }
            A_tA[i][j] = sum; 
        }

    for(int i=0; i<2; i++) //A_ty
        for(int j=0; j<1 ;j++)
        {   
            double sum = 0;
            for(int k=0; k<10; k++)
            {
                sum += A[k][i] * y[k];
            }
            A_ty[i][j] = sum; 
        }

    //Clamer法则运算线性方程组
    X[0] = (A_ty[0][0]*A_tA[1][1]-A_tA[0][1]*A_ty[1][0])/(A_tA[0][0]*A_tA[1][1]-A_tA[0][1]*A_tA[1][0]);
    X[1] = (A_tA[0][0]*A_ty[1][0]-A_tA[1][0]*A_ty[0][0])/(A_tA[0][0]*A_tA[1][1]-A_tA[0][1]*A_tA[1][0]);

    //均方误差
    for(int i=0; i<10; i++)
        Q += pow((X[0]*sin(x[i])+X[1]*cos(x[i])-y[i]),2);
    Q = Q / 10;

    //输出
    printf("a = %.15lf, b = %.15lf, 均方误差MSE = %.15lf\n", X[0], X[1], Q);
    system("pause");
}