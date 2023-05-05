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
    double psi0(double t)
    {
        return 16 * t * t;
    }
    double psi1(double t)
    {
        return  (1+t+ 16*t*t + pow(t,3) / 6);
    }
    double nb = 1; // b for task 2
    double DalamberU(double x, double t)
    {
        return x*x+16*t*t+x*t+x*pow(t,3)/6;
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
    vector<double> G((n - 1) * (m - 1), 0);
    for (int i = 0; i < m - 1; i++)
    {
        for (int j = 0; j < n - 1; j++)
        {
            G[i * (n - 1) + j] = g(X[j + 1], Time[i + 1]);
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
void PrintAllVectors(vector<vector<double>> V)
{
    cout << fixed << std::setprecision(5);
    cout << "-------------------------------------------------------------" << endl;
    for (int i = 0; i < V[0].size(); i++)
    {
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
vector<double> Ex3(int n, int m) // номер последней точки
{
    double h = (nb - a) / (n);
    double tau = T / (m);
    vector<double> X(n+1, 0), Time(m+1, 0), ans((n - 1) * (m - 1), 0); // мы не храним края и н у
    for (int i = 0; i <= n; i++)
        X[i] = a + i * h;
    for (int i = 0; i <= m; i++)
        Time[i] = i * tau;
    cout << "Tau is " << tau << "  " << "H is " << h << endl;
    vector<double> G = Make_G(n, m, X, Time);
    double p = sa * tau * tau / (h * h); // 
    // первая строчка значений 
    //ans[0] = p * (psi0(Time[2]) - 2 * psi0(Time[1]) + psi0(Time[0])) + tau * tau * G[0]/2 + psi0(Time[1]);
    for (int i = 0; i <= n - 2; i++)
    {
       ans[i] = p*(phi0(X[i+2])-2*phi0(X[i+1])+phi0(X[i]))+ tau*tau*G[i]+phi0(X[i+1]);
    }
    // вторая строчка
    ans[n-1] = p * (ans[1] - 2 * ans[0] + psi0(Time[1]) )+ tau * tau * G[n - 1]
        + 2 * ans[0] - phi0(X[1]);
    for (int i = 1; i < n - 2; i++)
    {
        ans[i + n - 1] = p * (ans[i + 1] - 2 * ans[i] + ans[i - 1] ) + tau * tau * G[i + n - 1]
            + 2 * ans[i] - phi0(X[i + 1]);
    }
    ans[2 * n - 3] = p * (psi1(Time[1]) - 2 * ans[2 * n - 3] + ans[2 * n - 4]) + tau * tau * G[2 * n - 3]
        + 2 * ans[2 * n - 3] - phi0(X[n - 2 + 1]);
    // остальные
    //for (int i = 2; i < m - 1; i++)
    //{
    //    /*ans[i*(n - 1)] = p * (ans[(i-1)*(n-1)+2] - 2 * ans[(i - 1) * (n - 1)+1] +
    //        phi0(X[i])) + tau * tau * G[i*(n - 1)] + 2 * ans[(i - 1) * (n - 1)+1] +
    //        ans[(i - 2) * (n - 1)];*/
    //    for (int j = 0; j < n - 1; j++)// цикл на последней итерации считает по краевому усл
    //    {
    //        if (j != n - 1)
    //        {
    //        
    //        ans[i * (n - 1) + j] = p * (ans[(i - 1) * (n - 1) + 2 + j] -
    //            2 * ans[(i - 1) * (n - 1) + 1 + j] + ans[(i - 1) * (n - 1) + j]) +
    //            tau * tau * G[i * (n - 1)] +
    //            2 * ans[(i - 1) * (n - 1) + 1 + j] + ans[(i - 2) * (n - 1) + j];
    //        }
    //        /*else
    //        {
    //        ans[i * (n - 1) + j] = p * (phi1(X[j]) -
    //            2 * ans[(i - 1) * (n - 1) + 1 + j] + ans[(i - 1) * (n - 1) + j]) +
    //            tau * tau * G[i * (n - 1)] +
    //            2 * ans[(i - 1) * (n - 1) + 1 + j] + ans[(i - 2) * (n - 1) + j];
    //        }*/
    //    }
    //}
    return ans;
}
//vector<double> Ex4(int n, int m)
//{
//    int n_big = int(l / h) + 1;
//    vector <double> X(n_big);
//    vector <double> T(n_big);
//    vector <vector<double>> U; // for U t is the first arg, x is second
//    for (int i = 0; i < n_big; i++)
//    {
//        X[i] = i * h;
//        T[i] = i * tau;
//    }
//    for (int m = 0; m < n_big; m++)
//    {
//        vector<double> temp(n_big, 0.);
//        U.push_back(temp);
//    }
//    for (int i = 0; i < n_big; i++)
//    {
//        U[0][i] = phi(X[i]);
//    }
//    for (int n = 1; n < n_big; n++)
//    {
//        vector <double> altha(n_big, 0), betta(n_big, 0);
//        altha[0] = 1; betta[0] = -h * psi0(T[n]);
//        for (int i = 1; i <= n_big - 2; i++)
//        {
//            double a = -tau / (h * h);
//            double b = 1 + 2 * tau / (h * h);
//            double c = -tau / (h * h);
//            double z = U[n - 1][i] + tau * f(T[n], X[i]);
//            altha[i] = -a / (b + c * altha[i - 1]);
//            betta[i] = (z - c * betta[i - 1]) / (b + c * altha[i - 1]);
//        }
//        U[n][n_big - 1] = ((3 * T[n] * T[n] / 2 + 2) * h + betta[n_big - 2]) / (1 + h - altha[n_big - 2]);
//        for (int i = n_big - 2; i >= 0; i--)
//        {
//            U[n][i] = altha[i] * U[n][i + 1] + betta[i];
//        }
//    }
//    PrintMatrix(U);
//    return U;
//}
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
    double dx = (nb - a) / (n - 1);
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
            temp[i * (n - 1) + j] = DalamberU(X[j + 1], Time[i + 1]);
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
    Temp.push_back(Ex3(15, 100));
    Temp.push_back(vector_Dalamber_U(15, 100));
    PrintAllVectors(Temp);
    // не забывать про условие устойчивости для явной схемы tau < dx
    //cout << "ex1 max razn = " << MaxRazn(Ex1(100, 1000), vector_true_U(100, 1000))<< endl;
    //cout << "ex2 max razn = " << MaxRazn(Ex2(100, 1000), vector_true_U(100, 1000))<< endl;
    cout << "ex3 max razn = " << MaxRazn(Ex3(15, 100), vector_Dalamber_U(15, 100)) << endl;
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
