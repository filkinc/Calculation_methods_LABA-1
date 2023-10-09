#include <iostream>
#include <fstream>
#include <vector>
#include "QuadMatrix.h"
#include "LinSolveAlgs.h"
#include "Norma.h"

using namespace std;
ofstream ansFile;

template<class T>
void test() {
    
    vector<T> fileVector;

    
    ifstream matrixFile;
    matrixFile.open("DATA22.txt");

    while (!matrixFile.eof()) {
        T i;
        matrixFile >> i;
        fileVector.push_back(i);
    }
    //fileVector.pop_back();
    int k = 0;
    //n* (n + 1) = fileVector.size();
    size_t n = (-1 + sqrt(1 + 4 * fileVector.size())) / 2;

    QuadMatrix<T> matrix(n);
    vector<T> b(n);

    if (n * (n + 1) != fileVector.size()) {
        ansFile << "Матрица системы не является квадратной!" << endl;
        system("pause");
        exit(1);
    }

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (k == fileVector.size()) {
                break;
            }
            else {
                matrix(i, j) = fileVector[k];
                k++;
            }
        }
        b[i] = fileVector[k];
        ++k;
    }
    matrixFile.close();

    //ВЫВОД В КОНСОЛЬ

    /*cout << "Исходная матрица (A):";
    for (int i = 0; i < n; i++) {
        cout << endl;
        for (int j = 0; j < n; j++) {
            cout << matrix(i, j) << " ";
        }
    }
    cout << endl << endl;
    cout << "Свободный вектор (b):" << endl;
    for (int i = 0; i < n; ++i) {
        cout << b[i] << endl;
    }

    vector<double> res = gaussLinSolve(matrix, b);

    cout << endl;
    cout << "Результат полученный методом Гаусса (x):" << endl;
    for (int i = 0; i < n; ++i) {
        cout << res[i] << endl;
    }

    auto [Q, R] = qrDecomposition(matrix);
    vector<double> resQr = qrLinSolve(Q, R, b);

    cout << endl;
    cout << "Матрица (Q):";
    for (int i = 0; i < n; i++) {
        cout << endl;
        for (int j = 0; j < n; j++) {
            cout << Q(i, j) << " ";
        }
    }
    cout << endl << endl;

    cout << "Матрица (R):";
    for (int i = 0; i < n; i++) {
        cout << endl;
        for (int j = 0; j < n; j++) {
            cout << R(i, j) << " ";
        }
    }
    cout << endl << endl;

    cout << "Результат полученный QR-разложением (x):" << endl;
    for (int i = 0; i < n; ++i) {
        cout << resQr[i] << endl;
    }

    cout << endl;
    cout << "Сферическа норма исходной матрицы: " << norm_inf(matrix) << endl;
    cout  << "Октаэдральная норма исходной матрицы: " << norm_1(matrix) << endl;

    cout << endl;
    cout << "Обратная матрица: ";
    QuadMatrix<double> inversMatrix = matrix.inv();
    for (int i = 0; i < n; ++i) {
        cout << endl;
        for (int j = 0; j < n; ++j) {
            cout << inversMatrix(i, j) << ' ';
        }
    }

    cout << endl << endl;
    cout << "Число обусловленности: " << endl;
    double condMatrixInf = cond(matrix, norm_inf);
    double condMatrix_1 = cond(matrix, norm_1);
    cout << "При кубической норме: " << condMatrixInf << endl;
    cout << "При октаэдральной норме: " << condMatrix_1 << endl;*/


    //ВЫВОД В ФАЙЛ
    //ansFile << "Точность" << T << endl;

    ansFile << "Исходная матрица (A):";
    for (int i = 0; i < n; i++) {
        ansFile << endl;
        for (int j = 0; j < n; j++) {
            ansFile << matrix(i, j) << " ";
        }
    }
    ansFile << endl << endl;

    /*ansFile << "Свободный вектор (b):" << endl;
    for (int i = 0; i < n; ++i) {
        ansFile << b[i] << endl;
    }
    cout << endl;*/

    //vector<double> res = gaussLinSolve(matrix, b);
    auto [res, C] = gaussLinSolve(matrix, b);

    ansFile << "Прямой метод Гаусса: ";
    for (int i = 0; i < n; ++i) {
        ansFile << endl;
        for (int j = 0; j < n; ++j) {
            ansFile << C(i, j) << ' ';
        }
    }
    ansFile << endl;

    ansFile << endl;
    ansFile << "Результат полученный методом Гаусса (x):" << endl;
    for (int i = 0; i < n; ++i) {
        ansFile << res[i] << endl;
    }

    ansFile << endl << "Норма вектора невязки в методе Гаусса (||b - b1||): " << endl;
    ansFile << "При кубической норме: " << normDiscrepancyVectorGauss(matrix, b, norm_inf) << endl;
    ansFile << "При октаэдральной норме: " << normDiscrepancyVectorGauss(matrix, b, norm_1) << endl;

    auto [Q, R] = qrDecomposition(matrix);
    vector<T> resQr = qrLinSolve(Q, R, b);

    ansFile << endl;
    ansFile << "Матрица (Q):";
    for (int i = 0; i < n; i++) {
        ansFile << endl;
        for (int j = 0; j < n; j++) {
            ansFile << Q(i, j) << " ";
        }
    }
    ansFile << endl << endl;

    ansFile << "Матрица (R):";
    for (int i = 0; i < n; i++) {
        ansFile << endl;
        for (int j = 0; j < n; j++) {
            ansFile << R(i, j) << " ";
        }
    }
    ansFile << endl << endl;

    ansFile << "Результат полученный QR-разложением (x):" << endl;
    for (int i = 0; i < n; ++i) {
        ansFile << resQr[i] << endl;
    }

    ansFile << endl << "Норма вектора невязки в методе QR-разложения (||b - b1||): " << endl;
    ansFile << "При кубической норме: " << normDiscrepancyVectorQR(matrix, b, norm_inf) << endl;
    ansFile << "При октаэдральной норме: " << normDiscrepancyVectorQR(matrix, b, norm_1) << endl;

    ansFile << endl;
    ansFile << "Кубическая норма исходной матрицы: " << norm_inf(matrix) << endl;
    ansFile << "Октаэдральная норма исходной матрицы: " << norm_1(matrix) << endl;

    ansFile << endl;
    ansFile << "Обратная матрица: ";
    QuadMatrix<T> inversMatrix = matrix.inv();
    for (int i = 0; i < n; ++i) {
        ansFile << endl;
        for (int j = 0; j < n; ++j) {
            ansFile << inversMatrix(i, j) << ' ';
        }
    }

    ansFile << endl << endl << "Оценка числа обусловленности: " << endl;

    /*
    vector<double> db(n, 0.1);
    auto px = gaussLinSolve(matrix, sum(b, db));
    auto dx = diff(px, res);
    */

    ansFile << "При кубической норме: " << condEstimate(matrix, norm_inf) << endl;
    ansFile << "При октаэдральной норме: " << condEstimate(matrix, norm_1) << endl;

    ansFile << endl;
    ansFile << "Число обусловленности: " << endl;
    T condMatrixInf = cond(matrix, norm_inf);
    T condMatrix_1 = cond(matrix, norm_1);
    ansFile << "При кубической норме: " << condMatrixInf << endl;
    ansFile << "При октаэдральной норме: " << condMatrix_1 << endl;
    ansFile << endl << endl << endl;
    

    /*cout << endl;
    for (int i = 0; i < fileVector.size(); i++) {
        cout << fileVector[i] << " ";
    }*/
}


int main()
{
    setlocale(LC_ALL, "Russian");

    ansFile.open("AnswerFileDATA22.txt");
    ansFile << "Точность double:" << endl;
    test<double>();
    ansFile << "Точность float:" << endl;
    test<float>();
    ansFile.close();

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
