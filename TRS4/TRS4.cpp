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
    double sa = 4; // square a
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

int main()
{
    //PrintVector(Ex1(5, 10));
    //PrintVector(Ex2(5, 5));
    //PrintVector(vector_true_U(5, 10));
    // не забывать про условие устойчивости для явной схемы tau < dx
    cout << "ex1 max razn = " << MaxRazn(Ex1(100, 1000), vector_true_U(100, 1000))<< endl;
    cout << "ex2 max razn = " << MaxRazn(Ex2(100, 1000), vector_true_U(100, 1000))<< endl;
    return 0;
}

// Запуск программы: CTRL+F5 или меню "Отладка" > "Запуск без отладки"
// Отладка программы: F5 или меню "Отладка" > "Запустить отладку"

// Советы по началу работы 
//   1. В окне обозревателя решений можно добавлять файлы и управлять ими.
//   2. В окне Team Explorer можно подключиться к системе управления версиями.
//   3. В окне "Выходные данные" можно просматривать выходные данные сборки и другие сообщения.
//   4. В окне "Список ошибок" можно просматривать ошибки.
//   5. Последовательно выберите пункты меню "Проект" > "Добавить новый элемент", чтобы создать файлы кода, или "Проект" > "Добавить существующий элемент", чтобы добавить в проект существующие файлы кода.
//   6. Чтобы снова открыть этот проект позже, выберите пункты меню "Файл" > "Открыть" > "Проект" и выберите SLN-файл.
