// TRS4.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <iostream>
#include <math.h>
#include <vector>
#include <iomanip>
using namespace std;
namespace var9
{
    // info for first task
    double a = 0;
    double b = 2;
    double T = 3;
    double c = 1;
    double true_u(double x, double t)
    {
        return 1 + t * x - exp(x - t - 1);
    }
    double f(double x, double t)
    {
        return x + t;
    }
    double phi(double x)
    {
        return true_u(x, 0);
    }
    double ksi0(double t)
    {
        return true_u(a, t);
    }
    // info for second task
    double sa = 1; // square a
    double g(double x, double t) // it is f but renamed to g
    {
        return x * t;
    }
    double phi0(double x)
    {
        return x * x;
    }
    double phi1(double x)
    {
        return x;
    }
    double d2phi0(double x)
    {
        return 2;
    }
    double psi0(double t)
    {
        return  t * t;
    }
    double psi1(double t)
    {
        return 2 + t + pow(t, 3) / 6;
    }
    double nb = 1; // b for task 2
    double DalamberU(double x, double t)
    {
        return x*x+t*t+x*t+x*pow(t,3)/6;
    }
}
using namespace var9;
vector<double> Make_F(int n , int m , vector<double> X, vector<double> Time)
{
    vector<double> F((n - 1) * (m - 1), 0);
    for (int i = 0; i < m - 1; i++)
    {
        for (int j = 0; j < n - 1; j++)
        {
            F[i * (n - 1) + j] = f(X[j + 1], Time[i + 1]);
        }
    }
    return F;
}
vector<double> Make_G(int n, int m, vector<double> X, vector<double> Time)
{
    vector<double> G(n  * m, 0);
    for (int i = 0; i < m ; i++)
    {
        for (int j = 0; j < n ; j++)
        {
            G[i * n  + j] = g(X[j + 1], Time[i + 1]);
        }
    }
    return G;
}
void PrintMatrix(vector<vector<double>> Matrix)
{
    cout << fixed << std::setprecision(5);
    cout << "-------------------------------------------------------------" << endl;
    for (int i = 0; i < Matrix.size(); i++)
    {
        for (int j = 0; j < Matrix[i].size(); j++)
        {
            cout << setw(7) << Matrix[i][j] << "  ";
        }
        cout << endl;
    }
    cout << "-------------------------------------------------------------" << endl;
}
void PrintVector(vector<double> V)
{
    cout << fixed << std::setprecision(5);
    cout << "-------------------------------------------------------------" << endl;
    for (int i = 0; i < V.size(); i++)
    {
        cout << setw(7) << V[i] << "  " << endl;
    }
    cout << "-------------------------------------------------------------" << endl;
}
void PrintAllVectors(vector<vector<double>> V, int n, int m)
{
    cout << fixed << std::setprecision(5);
    cout << "-------------------------------------------------------------" << endl;
    for (int i = 0; i < V[0].size(); i++)
    {
        if ((i) % n == 0)
        {
            cout << "row num is:" << (i) / n << endl;
        }
        for (int n = 0; n < V.size(); n++)
        {
        cout << setw(7) << V[n][i] << "  ";
        }
        cout << endl;
    }
    cout << "-------------------------------------------------------------" << endl;
}
double MaxRazn(vector<double> a1, vector<double> a2)
{
    double raz = 0;
    for (int i = 0; i < a1.size(); i++)
    {
        if (abs(a1[i] - a2[i]) > raz)
            raz = abs(a1[i] - a2[i]);
    }
    return raz;
}

vector<double> Ex1(int n, int m) // сколько всего точек 
{
    double dx = (b - a) / (n - 1);
    double tau = T / (m - 1);
    vector<double> X(n, 0), Time(m, 0), otv((n - 1) * (m - 1), 0);
    for (int i = 0; i < n; i++)
        X[i] = a + i * dx;
    for (int i = 0; i < m; i++)
        Time[i] = i * tau;
    vector<double> F = Make_F(n, m, X, Time);
    double alpha = 1. - c * tau / dx , betta = c * tau / dx;
    // первая строчка значений 
    otv[0] = alpha * phi(X[1]) + betta * ksi0(X[0]) + tau * F[0];
    for (int i = 1; i <= n - 2; i++)
    {
        otv[i] = alpha * phi(X[i]) + betta * phi(X[i - 1]) + tau * F[i];
    }
    for (int i = 1; i < m - 1; i++)
    {
        otv[i * (n - 1)] = alpha * otv[(i - 1) * (n - 1)] + betta * ksi0(Time[i + 1]) 
            + tau * F[i * (n - 1)];
        for (int j = 1; j < n - 1; j++)
        {
            otv[i * (n - 1) + j] = alpha * otv[(i - 1) * (n - 1) + j ] + betta * otv[(i - 1) * (n - 1)  + j - 1]
                + tau * F[i * (n - 1) + j ];
        }
    }
    return otv;
}
vector<double> Ex2(int n, int m) // сколько всего точек 
{
    double dx = (b - a) / (n - 1);
    double tau = T / (m - 1);
    vector<double> X(n, 0), Time(m, 0), otv((n - 1) * (m - 1), 0);
    for (int i = 0; i < n; i++)
        X[i] = a + i * dx;
    for (int i = 0; i < m; i++)
        Time[i] = i * tau;
    vector<double> F = Make_F(n, m, X, Time);
    double alpha = 1. + c * tau / dx, betta = c * tau / dx;
    // первая строчка значений 
    otv[0] = ( betta * ksi0(Time[1]) + phi(X[1]) + tau * F[0] ) / alpha;
    for (int i = 1; i < n - 1; i++)
    {
        otv[i] = (betta * otv[i - 1] + phi(X[i + 1]) + tau * F[i] ) / alpha;
    }
    for (int i = 1; i < m - 1; i++)
    {
        otv[i * (n - 1)] = (  betta * ksi0(Time[i + 1]) + otv[(i - 1) * (n - 1)] + tau * F[i * ( n - 1)]  ) / alpha;
        for (int j = 1; j < n - 1; j++)
        {
            otv[i * (n - 1) + j] = ( betta * otv[i * (n - 1) + j - 1] + otv[i * (n - 1) + j - (n - 1)] + tau * F[i * (n - 1) + j]  ) / alpha;
        }
    }
    return otv;
}
vector<double> Ex3(int n, int m) // номер последней точки матрицы, но у вектора минус 1
{
    double h = (nb - a) / (n+1);
    double tau = T / (m);
    vector<double> X(n+2, 0), Time(m+1, 0), ans(n*m, 0); // мы не храним края и н у
    for (int i = 0; i <= n+1; i++)
        X[i] = a + i * h;
    for (int i = 0; i <= m; i++)
        Time[i] = i * tau;
    cout << "Tau is " << tau << "  " << "H is " << h << endl;
    vector<double> G = Make_G(n, m, X, Time);
    double p = sa * tau * tau / (h * h); // 
    // первая строчка значений 
    for (int i = 0; i <= n - 1; i++)
    {
       ans[i] = phi0(X[i+1])+tau*phi1(X[i+1])+tau*tau/2*( sa* d2phi0(X[i + 1])+ G[i]);
    }
    // вторая строчка
    ans[n] = p * (ans[1] - 2 * ans[0] + psi0(Time[1]) )+ tau * tau * G[n]
        + 2 * ans[0] - phi0(X[1]);
    for (int i = 1; i < n; i++)
    {
        ans[n+i] = p * (ans[i + 1] - 2 * ans[i] + ans[i - 1] ) + tau * tau * G[n+i]
            + 2 * ans[i] - phi0(X[i + 1]);
    }
    ans[2 * n - 1] = p * ((ans[n-1]+h*psi1(Time[1])) - 2 * ans[n-1] + ans[n-2]) + tau * tau * G[2 * n - 1]
        + 2 * ans[n-1] - phi0(X[n]);
    // остальные
    for (int i = 2; i < m; i++)
    {
        ans[i*n] = p * (ans[(i-1)*n+1] - 2 * ans[(i - 1) * n] +
            psi0(Time[i])) + tau * tau * G[i * n] + 2 * ans[(i - 1) * n] -
            ans[(i - 2) * n];
        for (int j = 1; j < n-1; j++)// цикл на последней итерации считает по краевому усл
        {
            ans[i*n + j] = p * (ans[(i - 1) * n + 1+j] - 2 * ans[(i - 1) * n+j] + ans[(i - 1) * n - 1 + j]) 
                + tau * tau * G[i * n + j]
                + 2 * ans[(i - 1) * n + j] - ans[(i - 2) * n + j];
        }
        ans[(i+1) * n-1] = p * ((ans[i*n - 1] + h * psi1(Time[i])) - 2 * ans[i * n-1] +
            ans[i*n-2]) + tau * tau * G[(i + 1) * n - 1] + 2 * ans[i * n-1] -
            ans[(i - 1) * n-1];
    }
    return ans;
}
vector<double> Ex4(int n, int m)
{
    const double sigma = 0.3;
    double h = (nb - a) / (n + 1);
    double tau = T / (m);
    vector<double> X(n + 2, 0), Time(m + 1, 0), ans(n * m, 0); // мы не храним края и н у
    for (int i = 0; i <= n + 1; i++)
    X[i] = a + i * h;
    for (int i = 0; i <= m; i++)
    Time[i] = i * tau;
    cout << "Tau is " << tau << "  " << "H is " << h << endl;
    vector<double> G = Make_G(n, m, X, Time);

    vector<double> A(n);
    vector<double> B(n);
    vector<double> C(n);
    vector<double> f(n);

    for (int i = 1; i < n - 1; i++)
    {
        A[i] = -double(a) * sigma / h / h;
        B[i] = 1. / tau / tau + 2 * a * sigma / h / h;
        C[i] = -double(a) * sigma / h / h;
    }

    A[0] = 0;
    B[0] = -1. / h;
    C[0] = 1. / h;
    A[n - 1] = 0;
    B[n - 1] = 1.;
    C[n - 1] = 0;


    for (int j = 2; j < *n; j++)
    {
        double t = j * tau;
        f[0] = ksi_0(t);
        for (int i = 1; i < *m - 1; i++)
        {
            double x = i * h;
            f[i] = 2 * u[j - 1][i] / tau / tau - u[j - 2][i] / tau / tau + a * (1 - 2 * sigma) / h / h * (u[j - 1][i + 1] - 2 * u[j - 1][i]
                + u[j - 1][i - 1]) + a * sigma / h / h * (u[j - 2][i + 1] - 2 * u[j - 2][i] + u[j - 2][i - 1]) + F(t, x);
        }
        f[*m - 1] = ksi_1(t);

        tridiagonal_matrix(*m, A, B, C, f, u[j]);

        u[j][0] = -ksi_0(t) * h + u[j][1];
        u[j][*m - 1] = ksi_1(t);
    }

    double diff, M = 0;
    for (int j = 0; j < *n; j++)
    {
        for (int i = 0; i < *m; i++)
        {
            diff = abs(u[j][i] - True_solve(j * tau, i * h));
            if (diff > M) {
                M = diff;

            }
            //cout << j * tau << "    " << i * h << "    " << M << endl;
            //cout << j * tau << "    " << i * h << "    " << u[j][i] << "   " << True_solve(j * tau, i * h) << endl;
        }
        //cout << endl;
        //if (j == 3) system("pause");
    }
    cout << "n: " << n << "  m: " << m << "  error:" << M << endl;
    return 0;
}
vector<double> vector_true_U(int n, int m)
{
    double dx = (b - a) / (n - 1);
    double tau = T / (m - 1);
    vector<double> X(n, 0), Time(m, 0), temp((n - 1) * (m - 1), 0);
    for (int i = 0; i < n; i++)
        X[i] = a + i * dx;
    for (int i = 0; i < m; i++)
        Time[i] = i * tau;
    //PrintVector(X);
    //PrintVector(Time);

    for (int i = 0; i < m - 1; i++)
    {
        for (int j = 0; j < n - 1; j++)
        {
            temp[i * (n - 1) + j] = true_u(X[j + 1], Time[i + 1]);
        }
    }
    return temp;
}
vector<double> vector_Dalamber_U(int n, int m)
{
    double dx = (nb - a) / (n + 1);
    double tau = T / (m);
    vector<double> X(n+2, 0), Time(m+1, 0), temp(n*m, 0);
    for (int i = 0; i < n+2; i++)
        X[i] = a + i * dx;
    for (int i = 0; i < m+1; i++)
        Time[i] = i * tau;
    //PrintVector(X);
    //PrintVector(Time);

    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n ; j++)
        {
            temp[i * n + j] = DalamberU(X[j + 1], Time[i + 1]);
        }
    }
    return temp;
}

int main()
{
    //PrintVector(Ex1(5, 10));
    //PrintVector(Ex2(5, 5));
    //PrintVector(vector_true_U(5, 10));
    /*vector<vector<double>> Temp;
    Temp.push_back(Ex1(10, 100));
    Temp.push_back(Ex2(10, 100));
    Temp.push_back(vector_true_U(10, 100));*/
    vector<vector<double>> Temp;
    /*Temp.push_back(Ex3(10, 100));
    Temp.push_back(vector_Dalamber_U(10, 100));
    PrintAllVectors(Temp,10,100);*/
    // не забывать про условие устойчивости для явной схемы tau < dx
    //cout << "ex1 max razn = " << MaxRazn(Ex1(100, 1000), vector_true_U(100, 1000))<< endl;
    //cout << "ex2 max razn = " << MaxRazn(Ex2(100, 1000), vector_true_U(100, 1000))<< endl;
    cout << "ex3 max razn = " << MaxRazn(Ex3(1500, 10000), vector_Dalamber_U(1500, 10000)) << endl;
    return 0;
}
