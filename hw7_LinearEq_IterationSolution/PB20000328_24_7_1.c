#include<stdio.h>
#include<math.h>

#define epsilon pow(10,-7)

double* GuassSeidel(double A[][9], double* b);
double* SOR(double A[][9], double* b, double omega);
FILE * fp;

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
    double *y[99];
    fp = fopen("D://file_from_desktop//computing_method//hw7//result.txt","w");
    
    fprintf(fp, "Gauss-Seidel:\n");
    x = GuassSeidel(A,b);
    
    fprintf(fp, "SOR:\n");
    for(int i=1; i<100; i++) 
    {
        fprintf(fp, "omega = %d/50:\n",i);
        y[i] = SOR(A,b,i/50.0);
    }
    fclose(fp);
    system("pause");
}

double* GuassSeidel(double A[][9], double* b)
{
    //初始化为0
    static double x1[9] = {0}, x2[9] = {0};
    double norm;
    double sum;
    int count = 0;

    //迭代
    do
    {
        //x2->x1迁移
        for(int i=0; i<9; i++)  x1[i] = x2[i];

        //迭代格式
        for(int i=0; i<9; i++) 
        {
            sum = 0;
            for(int j=0; j<9; j++) sum += A[i][j] * x2[j];
            x2[i] = (b[i] - sum + A[i][i] * x2[i]) / A[i][i];
        }

        //计算范数
        double max = 0;
        for(int i=0; i<9; i++) if(max < fabs(x2[i] - x1[i])) max = fabs(x2[i] - x1[i]);
        norm = max;
        count += 1;
    }while(norm >= epsilon);
    fprintf(fp, "G-S迭代步数为：%d\n", count);
    for(int i=0; i<9; i++) fprintf(fp, "x%d=%.15f\n", i+1, x2[i]);
    return x2;
}

double* SOR(double A[][9], double* b, double omega)
{
    static double x1[9] = {0}, x2[9] = {0};
    double norm;
    double sum;
    int count = 0;

    //初始化
    for(int i=0; i<9; i++)  x1[i] = x2[i] = 0;
    
    //迭代
    do
    {
        //x2->x1迁移
        for(int i=0; i<9; i++)  x1[i] = x2[i];

        //迭代格式
        for(int i=0; i<9; i++) 
        {
            sum = 0;
            for(int j=0; j<9; j++) sum += A[i][j] * x2[j];
            x2[i] = (1 - omega) * x2[i] + omega * (b[i] - sum + A[i][i] * x2[i]) / A[i][i];
        }

        //计算范数
        double max = 0;
        for(int i=0; i<9; i++) if(max < fabs(x2[i] - x1[i])) max = fabs(x2[i] - x1[i]);
        norm = max;
        count += 1;
    }while(norm >= epsilon);
    fprintf(fp, "SOR迭代步数为：%d\n", count);
    for(int i=0; i<9; i++) fprintf(fp, "x%d=%.15f\n", i+1, x2[i]);
    return x2;
}

