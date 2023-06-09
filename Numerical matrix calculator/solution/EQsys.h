#pragma once
#include <vector>
#include "SQmatrix.h"

#define LIMIT 200

#define EQSYS_FREETERMS_SIZE_MISSMATCH 100120
#define EQSYS_MATRIX_SIZE_MISSMATCH 100120
#define EQSYS_LUSOL_LUERROR 100201
#define EQSYS_SOLVE_WRONG_TYPE 100301
#define EQSYS_SOLVE_FPI_NOT_APPLYABLE 100402

class EQsys {
    //класс описывает систему уравнений в матричом виде
private:
    int varCnt;
    SQmatrix matr;
    std::vector<double> freeTerms;
    bool trace;

public:

    //пустая система
    EQsys();

    //незаполненная система
    EQsys(int s);

    //система, с заданной основной матрицей (задается через вектор)
    //система однородна
    EQsys(int s, std::vector<std::vector<double>>& v);

    //система, с заданной основной матрицей (задается через матрицу)
    //система однородна
    EQsys(int s, SQmatrix m);

    //система, с заданной основной матрицей (задается через вектор)
    //система неоднородна
    EQsys(int s, std::vector<std::vector<double>>& v, std::vector<double>& b);

    //система, с заданной основной матрицей (задается через матрицу)
    //система неоднородна
    EQsys(int s, SQmatrix m, std::vector<double>& b);

    //включена ли трассировка вычислений
    bool tracef();

    void traceSwitch();

    SQmatrix matrixf() const;

    std::vector<double> freeTermsf() const;

    //решить СЛАУ с помощью LU разложения
    std::vector<double> LUsolve(); 

    //решить СЛАУ методом прогонки
    std::vector<double> TMAsolve();

    //решить СЛАУ методом простых итераций
    std::pair<std::vector<double>, double> FPIsolve(double eps);

    //решить СЛАУ методом Гаусса-Зейделя
    std::pair<std::vector<double>, double> GSsolve(double eps);

    //автоматически выбирает метод решения
    std::vector<double> Solve();

};