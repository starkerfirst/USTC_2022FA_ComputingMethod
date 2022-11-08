#include<stdio.h>
#include<math.h>

double* Guass(double A[][9], double* b);

int main()
{
    //定义
    double A[9][9]={  \
        {31,-13,0,0,0,-10,0,0,0}, \
        {-13,35,-9,0,-11,0,0,0,0}, \
        {0,-9,31,-10,0,0,0,0,0}, \
        {0,0,-10,79,-30,0,0,0,-9}, \
        {0,0,0,-30,57,-7,0,-5,0}, \
        {0,0,0,0,-7,47,-30,0,0}, \
        {0,0,0,0,0,-30,41,0,0}, \
        {0,0,0,0,-5,0,0,27,-2}, \
        {0,0,0,-9,0,0,0,-2,29} 
    };
    double b[9]={-15,27,-23,0,-20,12,-7,7,10};
    double *x;
    x = Guass(A,b);
    for(int i=0; i<9; i++) printf("x%d=%.15f\n",i+1,x[i]);

    system("pause");
}

double* Guass(double A[][9], double* b)
{
    
    for(int i=0; i<9; i++)
        {
            //找到每列最大值
            int max = i;
            for(int j=i; j<9; j++)
            {
                if(fabs(A[j][i])>fabs(A[max][i]))  max = j;
            }
            if(max!=i) //交换两行数据
            {
                double tmp;
                for(int k=i; k<9; k++)
                    {
                        tmp = A[max][k];
                        A[max][k] = A[i][k];
                        A[i][k] = tmp;
                    }
                tmp = b[max];
                b[max] = b[i];
                b[i] = tmp;
            }

            //消列
            for(int j=i+1; j<9; j++)
            {
                b[j] = b[j] - A[j][i] / A[i][i] * b[i];
                for(int k=8; k>=i; k--) A[j][k] = A[j][k] - A[j][i] / A[i][i] * A[i][k];
            }
        }
    
    //回代
    for(int i=8; i>=0; i--)
        {
            for(int j=i-1; j>=0; j--)
                b[j] -= b[i] * A[j][i] / A[i][i];
            b[i] = b[i] / A[i][i];
        }

    return b;
}

