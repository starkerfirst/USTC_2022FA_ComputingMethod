#include<stdio.h>
#include<math.h>

#define epsilon pow(10,-4) //取控制精度为10^-4
#define N 5

void Jacobi(double A[][N]);
double* Guass(double A[][N-1], double* b);
double SquaredNondiagonal(double A[][N]);
void PrintA(double A[][N]);
void copy(double A[][N],double B[][N]); //B抄A

//定义
double A[N][N]={  \
    {3,2,5,4,6}, \
    {2,1,3,-7,8}, \
    {5,3,2,5,-4}, \
    {4,-7,5,1,3}, \
    {6,8,-4,3,8}
};


int main()
{
    /* double A[N][N]={  \
        {2,-1,3}, \
        {-1,5,0}, \
        {3,0,1}
    };
 */ 
    //运行
    double B[N][N];
    copy(A,B);
    //PrintA(B);
    Jacobi(B);
    
    system("pause");
}

void Jacobi(double B[N][N])
{
    double s,t,c,d;
    int p,q; 
    double C[N][N];
    while(SquaredNondiagonal(B) > epsilon)  
    {
        //判断p,q
        copy(B,C);
        double max=0;
        for(int i=0;i<N;i++)
            for(int j=i+1;j<N;j++)
                if(max<fabs(B[i][j]))  {max = fabs(B[i][j]); p=i; q=j;}

        //t生成
        s = (B[q][q]-B[p][p])/(2*B[p][q]);
        t = fabs(-s-sqrt(s*s+1)) > fabs(-s+sqrt(s*s+1)) ? -s+sqrt(s*s+1) : -s-sqrt(s*s+1);

        //c,d生成
        c = 1/sqrt(1+t*t);
        d = t*c;

        //迭代
        for(int i=0;i<N;i++)
            if(!(i==p||i==q)) 
            {
                B[i][p] = B[p][i] = c * C[p][i] - d * C[q][i];
                B[i][q] = B[q][i] = c * C[q][i] + d * C[p][i]; 
            }
        B[p][p] -= t*C[p][q];
        B[q][q] += t*C[p][q];
        B[p][q] = B[q][p] = 0;
        //PrintA(B);
        //putchar('\n');
    }
    
    //特征向量计算
    for(int i=0;i<N;i++)//第i个特征值
    {
        double S[N-1][N-1], b[N];
        //特征多项式矩阵生成
        for(int j=0;j<N-1;j++)
            for(int k=0;k<N-1;k++)
            {
                if(j==k) S[j][k] = A[j][k] - B[i][i];
                else S[j][k] = A[j][k];
            }
        //b生成
        for(int j=0;j<N-1;j++) b[j] = -A[j][N-1]; //设特征向量最后一位恒为1
        //Gauss
        Guass(S, b);
        //归一化
        b[N-1]=1;
        double sum=0;
        for(int j=0;j<N;j++) sum+=b[j]*b[j];
        for(int j=0;j<N;j++) b[j]/=sqrt(sum);
        //输出
        printf("r%d = %.15f, v = (%.15f,%.15f,%.15f,%.15f,%.15f)\n", i+1, B[i][i], b[0], b[1], b[2], b[3], b[4]);
        //for(int i=0;i<N;i++) printf("r%d = %.15f", i+1, A[i][i]);
    }
    
}

//计算非对角元平方和
double SquaredNondiagonal(double A[][N])
{
    int i,j;
    double sum=0;
    for(i=0;i<N;i++)
        for(j=0;j<N;j++)
            if(i!=j)  sum += A[i][j] * A[i][j];
    return sqrt(sum);
}


void PrintA(double A[][N])
{
    int i,j;
    for(i=0;i<N;i++)
        for(j=0;j<N;j++)
        {
            printf("%.5f ", A[i][j]);
            if(j==N-1) putchar('\n');
        }
}

void copy(double A[][N],double B[][N])
{
    for(int j=0;j<N;j++)
        for(int k=0;k<N;k++)
        {
            B[j][k] = A[j][k];
        }
}


#undef N
#define N 4
//高斯消元法
double* Guass(double A[][N], double* b)
{
    
    for(int i=0; i<N; i++)
        {
            //找到每列最大值
            int max = i;
            for(int j=i; j<N; j++)
            {
                if(fabs(A[j][i])>fabs(A[max][i]))  max = j;
            }
            if(max!=i) //交换两行数据
            {
                double tmp;
                for(int k=i; k<N; k++)
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
            for(int j=i+1; j<N; j++)
            {
                b[j] = b[j] - A[j][i] / A[i][i] * b[i];
                for(int k=N-1; k>=i; k--) A[j][k] = A[j][k] - A[j][i] / A[i][i] * A[i][k];
            }
        }
    
    //回代
    for(int i=N-1; i>=0; i--)
        {
            for(int j=i-1; j>=0; j--)
                b[j] -= b[i] * A[j][i] / A[i][i];
            b[i] = b[i] / A[i][i];
        }

    return b;
}

