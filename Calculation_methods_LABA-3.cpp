// Calculation_methods_LABA-3.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <iostream>
#include <fstream>
#include <vector>
#include <utility>
#include <cmath>
#include <numeric>
#include "QuadMatrix.h"
#include "Matrix.h"
#include "Norma.h"
#include "LinSolveAlgs.h"
#include "Diag3Matrix.h"
#include "IterativeLinSolveAlgs.h"
#include "Interpolation.h"
#include <string>


using namespace std;
ofstream out;


template<class T>
void printToFile(T(&func)(T), T a, T b, int n, const char* name_unif, const char* name_chebysh) {
    
    //Равномерная сетка
    vector<T> unifGrid = uniformGrid(a, b, n);

    auto calcValTable = calculateValueTable<T>(func, unifGrid);

    out.open(name_unif);
    for (int i = 0; i < calcValTable.x.size(); ++i) {
        out << calcValTable.x[i] << " ";
    }
    out << endl;

    for (int i = 0; i < calcValTable.x.size(); ++i) {
        out << calcValTable.f[i] << " ";
    }
    cout << endl;

    out.close();

    // Сетка Чебышева
    vector<double> chebyshGrid = chebyshevGrid(a, b, n);

    auto calcValTableChebysh = calculateValueTable<T>(func, chebyshGrid);

    out.open(name_chebysh);
    for (int i = 0; i < calcValTableChebysh.x.size(); ++i) {
        out << calcValTableChebysh.x[i] << " ";
    }
    out << endl;

    for (int i = 0; i < calcValTableChebysh.x.size(); ++i) {
        out << calcValTableChebysh.f[i] << " ";
    }
    cout << endl;

    out.close();
}

template<class T>
void printToFileSplain(T(&func)(T), T(&funcdiff)(T), T a, T b, int n, const char* name_unif, const char* name_chebysh) {

    //Равномерная сетка
    vector<T> unifGrid = uniformGrid(a, b, n);

    auto calcValTable = calculateValueTable<T>(func, unifGrid);
    auto calcValDiffTable = calculateValueTable<T>(funcdiff, unifGrid);

    auto splInterp = splaineInterpolation(calcValTable, calcValDiffTable, 1e-7);

    out.open(name_unif);
    for (int i = 0; i < splInterp.x.size(); ++i) {
        out << splInterp.x[i] << " ";
    }
    out << endl;

    for (int i = 0; i < splInterp.a.size(); ++i) {
        out << splInterp.a[i] << " ";
    }
    out << endl;

    for (int i = 0; i < splInterp.b.size(); ++i) {
        out << splInterp.b[i] << " ";
    }
    out << endl;

    for (int i = 0; i < splInterp.c.size(); ++i) {
        out << splInterp.c[i] << " ";
    }
    out << endl;

    for (int i = 0; i < splInterp.d.size(); ++i) {
        out << splInterp.d[i] << " ";
    }

    out.close();

    // Сетка Чебышева
    vector<double> chebyshGrid = chebyshevGrid(a, b, n);

    auto calcValTableChebysh = calculateValueTable<T>(func, chebyshGrid);
    auto calcValDiffTableChebysh = calculateValueTable<T>(funcdiff, chebyshGrid);

    auto splInterpChebysh = splaineInterpolation(calcValTableChebysh, calcValDiffTableChebysh, 1e-7);
        
    out.open(name_chebysh);
    for (int i = 0; i < splInterpChebysh.x.size(); ++i) {
        out << splInterpChebysh.x[i] << " ";
    }
    out << endl;

    for (int i = 0; i < splInterpChebysh.a.size(); ++i) {
        out << splInterpChebysh.a[i] << " ";
    }
    out << endl;

    for (int i = 0; i < splInterpChebysh.b.size(); ++i) {
        out << splInterpChebysh.b[i] << " ";
    }
    out << endl;

    for (int i = 0; i < splInterpChebysh.c.size(); ++i) {
        out << splInterpChebysh.c[i] << " ";
    }
    out << endl;

    for (int i = 0; i < splInterpChebysh.d.size(); ++i) {
        out << splInterpChebysh.d[i] << " ";
    }

    out.close();
}

int main(){

    setlocale(LC_ALL, "Russian");

    cout << "Полином Лагранжа:" << endl << endl;

    // Функция y = x^2, отрезок [-1,1]:
    auto sqrX = [](double x) { return x * x; };
    auto diff1sqrX = [](double x) { return 2 * x; };
    auto diffsqrX = [](double x) { return 2; };

    // Футкция y = 1/(1 + x^2), отрезок [-1, 1]:
    auto function2 = [](double x) { return 1 / (1 + x * x); };
    auto diff1function2 = [](double x) { return -((2 * x) / (1 + x * x) / (1 + x * x)); };
    auto difffunction2 = [](double x) { return (8 * x * x) / pow((1 + x * x), 3) - 2 / (1 + x * x) / (1 + x * x); };

    // Футкция y = 1 / atan(1 + 10 * x * x), отрезок [-3, 3]:
    auto function3 = [](double x) { return 1 / atan(1 + 10 * x * x); };
    auto diff1function3 = [](double x) { return -((20 * x) / ((1 + (1 + 10 * x * x) * (1 + 10 * x * x)) * pow(atan(1 + 10 * x * x), 2))); };
    auto difffunction3 = [](double x) { return (800 * pow(x, 2)) / (pow(1 + pow(1 + 10 * pow(x, 2), 2), 2) * pow(atan(1 + 10 * pow(x, 2)), 3)) +
        (800 * pow(x, 2) * (1 + 10 * pow(x, 2))) / (pow(1 + pow(1 + 10 * pow(x, 2), 2), 2) * pow(atan(1 + 10 * pow(x, 2)), 2)) -
        20 / ((1 + pow(1 + 10 * pow(x, 2), 2)) * pow(atan(1 + 10 * pow(x, 2)), 2)); };

    // Функция y = pow((4 * x * x * x + 2 * x * x - 4 * x + 2), sqrt(2)) + asin(1 / 5 + x - x * x), отрезок [-1, 1]:
    auto  function4 = [](double x) { return pow((4 * x * x * x + 2 * x * x - 4 * x + 2), sqrt(2)) + asin(1 / 5 + x - x * x); };
    auto  diff1function4 = [](double x) { 
        return (1 - 2 * x) / sqrt(1 - (-(x * x) + x + 1 / 5) * (-(x * x) + x + 1 / 5)) + sqrt(2) * (12 * x * x + 4 * x - 4) * pow((4 * x * x + 2 * x * x - 4 * x + 2), (sqrt(2) - 1)); };
    auto  difffunction4 = [](double x) {
        return sqrt(2) * (-1 + sqrt(2)) * pow(-4 + 4 * x + 12 * pow(x, 2), 2) * pow(2 - 4 * x + 2 * pow(x, 2) + 4 * pow(x, 3), -2 + sqrt(2)) +
            sqrt(2) * (4 + 24 * x) * pow(2 - 4 * x + 2 * pow(x, 2) + 4 * pow(x, 3), -1 + sqrt(2)) +
            (pow(1 - 2 * x, 2) * (0.2 + x - pow(x, 2))) / pow(1 - pow(0.2 + x - pow(x, 2), 2), 1.5) -
            2 / sqrt(1 - pow(0.2 + x - pow(x, 2), 2)); };


    // Функция y = x * x * x; отрезок [-1, 1]:
    auto cubX = [](double x) { return x * x * x; };
    auto diff1cubX = [](double x) { return 3 * x * x; };
    auto diffcubX = [](double x) { return 6 * x; };

    // Функция y = const; отрезок [-1, 1]:
    auto constX = [](double x) { return 3; };
    auto X = [](double x) { return 0; };



    vector<double> gridChebyshTest = chebyshevGrid(-1., 1., 3);

    auto testSqrX = calculateValueTable<double>(toFunction<double, decltype(&sqrX)>, gridChebyshTest);
    auto testdiffSqrX = calculateValueTable<double>(toFunction<double, decltype(&diffsqrX)>, gridChebyshTest);

    auto splain = splaineInterpolation(testSqrX, testdiffSqrX, 1e-7);
    cout << splain(0.92);

    //Для полинома Лагранжа
    printToFile<double>(toFunction<double, decltype(&sqrX)>, -1., 1., 3, "function1_unif_n=3.txt", "function1_chebysh_n=3.txt");
    printToFile<double>(toFunction<double, decltype(&sqrX)>, -1., 1., 10, "function1_unif_n=10.txt", "function1_chebysh_n=10.txt");
    printToFile<double>(toFunction<double, decltype(&sqrX)>, -1., 1., 50, "function1_unif_n=50.txt", "function1_chebysh_n=50.txt");

    printToFile<double>(toFunction<double, decltype(&function2)>, -1., 1., 3, "function2_unif_n=3.txt", "function2_chebysh_n=3.txt");
    printToFile<double>(toFunction<double, decltype(&function2)>, -1., 1., 10, "function2_unif_n=10.txt", "function2_chebysh_n=10.txt");
    printToFile<double>(toFunction<double, decltype(&function2)>, -1., 1., 50, "function2_unif_n=50.txt", "function2_chebysh_n=50.txt");

    printToFile<double>(toFunction<double, decltype(&function3)>, -3., 3., 3, "function3_unif_n=3.txt", "function3_chebysh_n=3.txt");
    printToFile<double>(toFunction<double, decltype(&function3)>, -3., 3., 10, "function3_unif_n=10.txt", "function3_chebysh_n=10.txt");
    printToFile<double>(toFunction<double, decltype(&function3)>, -3., 3., 50, "function3_unif_n=50.txt", "function3_chebysh_n=50.txt");

    printToFile<double>(toFunction<double, decltype(&function4)>, -1., 1., 3, "function4_unif_n=3.txt", "function4_chebysh_n=3.txt");
    printToFile<double>(toFunction<double, decltype(&function4)>, -1., 1., 10, "function4_unif_n=10.txt", "function4_chebysh_n=10.txt");
    printToFile<double>(toFunction<double, decltype(&function4)>, -1., 1., 50, "function4_unif_n=50.txt", "function4_chebysh_n=50.txt");

    printToFile<double>(toFunction<double, decltype(&cubX)>, -1., 1., 3, "cubX_unif_n=3.txt", "cubX_chebysh_n=3.txt");
    printToFile<double>(toFunction<double, decltype(&cubX)>, -1., 1., 10, "cubX_unif_n=10.txt", "cubX_chebysh_n=10.txt");
    printToFile<double>(toFunction<double, decltype(&cubX)>, -1., 1., 50, "cubX_unif_n=50.txt", "cubX_chebysh_n=50.txt");

    printToFile<double>(toFunction<double, decltype(&constX)>, -1., 1., 3, "const_unif_n=3.txt", "const_chebysh_n=3.txt");
    printToFile<double>(toFunction<double, decltype(&constX)>, -1., 1., 10, "const_unif_n=10.txt", "const_chebysh_n=10.txt");
    printToFile<double>(toFunction<double, decltype(&constX)>, -1., 1., 50, "const_unif_n=50.txt", "const_chebysh_n=50.txt");

    //Для сплайн интерполяции
    printToFileSplain(toFunction<double, decltype(&sqrX)>, toFunction<double, decltype(&diffsqrX)>, -1., 1., 3, "spl_function1_unif_n=3.txt", "spl_function1_chebysh_n=3.txt");
    printToFileSplain(toFunction<double, decltype(&sqrX)>, toFunction<double, decltype(&diffsqrX)>, -1., 1., 10, "spl_function1_unif_n=10.txt", "spl_function1_chebysh_n=10.txt");
    printToFileSplain(toFunction<double, decltype(&sqrX)>, toFunction<double, decltype(&diffsqrX)>, -1., 1., 50, "spl_function1_unif_n=50.txt", "spl_function1_chebysh_n=50.txt");

    printToFileSplain(toFunction<double, decltype(&function2)>, toFunction<double, decltype(&difffunction2)>, -1., 1., 3, "spl_function2_unif_n=3.txt", "spl_function2_chebysh_n=3.txt");
    printToFileSplain(toFunction<double, decltype(&function2)>, toFunction<double, decltype(&difffunction2)>, -1., 1., 10, "spl_function2_unif_n=10.txt", "spl_function2_chebysh_n=10.txt");
    printToFileSplain(toFunction<double, decltype(&function2)>, toFunction<double, decltype(&difffunction2)>, -1., 1., 50, "spl_function2_unif_n=50.txt", "spl_function2_chebysh_n=50.txt");

    printToFileSplain(toFunction<double, decltype(&function3)>, toFunction<double, decltype(&difffunction3)>, -3., 3., 3, "spl_function3_unif_n=3.txt", "spl_function3_chebysh_n=3.txt");
    printToFileSplain(toFunction<double, decltype(&function3)>, toFunction<double, decltype(&difffunction3)>, -3., 3., 10, "spl_function3_unif_n=10.txt", "spl_function3_chebysh_n=10.txt");
    printToFileSplain(toFunction<double, decltype(&function3)>, toFunction<double, decltype(&difffunction3)>, -3., 3., 50, "spl_function3_unif_n=50.txt", "spl_function3_chebysh_n=50.txt");

    printToFileSplain(toFunction<double, decltype(&function4)>, toFunction<double, decltype(&difffunction4)>, -1., 1., 3, "spl_function4_unif_n=3.txt", "spl_function4_chebysh_n=3.txt");
    printToFileSplain(toFunction<double, decltype(&function4)>, toFunction<double, decltype(&difffunction4)>, -1., 1., 10, "spl_function4_unif_n=10.txt", "spl_function4_chebysh_n=10.txt");
    printToFileSplain(toFunction<double, decltype(&function4)>, toFunction<double, decltype(&difffunction4)>, -1., 1., 50, "spl_function4_unif_n=50.txt", "spl_function4_chebysh_n=50.txt");

    printToFileSplain(toFunction<double, decltype(&cubX)>, toFunction<double, decltype(&diffcubX)>, -1., 1., 3, "spl_cubX_unif_n=3.txt", "spl_cubX_chebysh_n=3.txt");
    printToFileSplain(toFunction<double, decltype(&cubX)>, toFunction<double, decltype(&diffcubX)>, -1., 1., 10, "spl_cubX_unif_n=10.txt", "spl_cubX_chebysh_n=10.txt");
    printToFileSplain(toFunction<double, decltype(&cubX)>, toFunction<double, decltype(&diffcubX)>, -1., 1., 50, "spl_cubX_unif_n=50.txt", "spl_cubX_chebysh_n=50.txt");

    printToFileSplain(toFunction<double, decltype(&constX)>, toFunction<double, decltype(&X)>, -1., 1., 3, "spl_const_unif_n=3.txt", "spl_const_chebysh_n=3.txt");
    printToFileSplain(toFunction<double, decltype(&constX)>, toFunction<double, decltype(&X)>, -1., 1., 10, "spl_const_unif_n=10.txt", "spl_const_chebysh_n=10.txt");
    printToFileSplain(toFunction<double, decltype(&constX)>, toFunction<double, decltype(&X)>, -1., 1., 50, "spl_const_unif_n=50.txt", "spl_const_chebysh_n=50.txt");

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
