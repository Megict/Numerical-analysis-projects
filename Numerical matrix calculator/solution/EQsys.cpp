#include <iostream>
#include <vector>
#include "VectorExten.h"
#include "SQmatrix.h"
#include "EQsys.h"

EQsys::EQsys() {
    //пустая система
    trace = 0;
    varCnt = 0;
    matr = SQmatrix(0);
    freeTerms = std::vector<double>(0, 0);
}

EQsys::EQsys(int s) {
    //незаполненная система
    trace = 0;
    varCnt = s;
    matr = SQmatrix(s);
    freeTerms = std::vector<double>(s, 0);
}

EQsys::EQsys(int s, std::vector<std::vector<double>>& v) {
    trace = 0;
    //система, с заданной основной матрицей (задается через вектор)
    //система однородна
    try {
        matr = SQmatrix(s, v);
    }
    catch (int err) {
        throw err;
    }

    varCnt = s;
    freeTerms = std::vector<double>(s, 0);
}

EQsys::EQsys(int s, SQmatrix m) {
    trace = 0;
    //система, с заданной основной матрицей (задается через матрицу)
    //система однородна

    if (m.sizef() != s) {
        throw EQSYS_MATRIX_SIZE_MISSMATCH;
    }

    varCnt = s;
    matr = m;
    freeTerms = std::vector<double>(s, 0);
}

EQsys::EQsys(int s, std::vector<std::vector<double>>& v, std::vector<double>& b) {
    trace = 0;
    //система, с заданной основной матрицей (задается через вектор)
    //система неоднородна
    varCnt = s;
    try {
        matr = SQmatrix(s, v);
    }
    catch (int err) {
        throw err;
    }

    if (b.size() != s) {
        throw EQSYS_FREETERMS_SIZE_MISSMATCH;
    }

    freeTerms = b;
}

EQsys::EQsys(int s, SQmatrix m, std::vector<double>& b) {
    trace = 0;
    //система, с заданной основной матрицей (задается через матрицу)
    //система неоднородна

    //как-то в дебажный вывод это надо отправить
    //std::cout << "inp matr| size: " << m.sizef() << "  type: " << m.typef() << "\n";

    if (m.sizef() != s) {
        if(m.typef() != "tridiagonalCompact")
            throw EQSYS_MATRIX_SIZE_MISSMATCH;
        else {
            if (m.sizef() != 3) {
                throw EQSYS_MATRIX_SIZE_MISSMATCH;
            }
        }
    }

    if (b.size() != s) {
        throw EQSYS_FREETERMS_SIZE_MISSMATCH;
    }

    varCnt = s;
    matr = m;
    freeTerms = b;
}

bool EQsys::tracef() {
    return trace;
}

void EQsys::traceSwitch() {
    trace = !trace;
}

SQmatrix EQsys::matrixf() const {
    return matr;
}

std::vector<double> EQsys::freeTermsf() const {
    return freeTerms;
}



std::vector<double> EQsys::LUsolve() {
    if (matr.typef() == "tridiagonalCompact") {
        throw EQSYS_SOLVE_WRONG_TYPE;
    }

    //перестановка, чтобы на ГД не было нулей
    /*std::pair<SQmatrix, std::pair<std::vector<double>,std::vector<int>>> arrMatr = matr.arrange(freeTerms);
    SQmatrix RSMmatr = arrMatr.first;
    std::vector<double> RSMfreeTerms = arrMatr.second.first;
    std::vector<int> resMap = arrMatr.second.second;*/

    std::vector<double> RSMfreeTerms = freeTerms;
    std::pair<std::pair<SQmatrix, SQmatrix>,std::vector<int>> preres = matr.LUdecomp(RSMfreeTerms);
    std::pair<SQmatrix, SQmatrix> LUres = preres.first;
    std::vector<int> resMap = preres.second;

    if (trace) {
        printf("LU decomp. results:\n");
        LUres.first.print();
        LUres.second.print();
    }



    //проверка правильности LU разложения
    if (!matr.LUcheck(LUres, 0.01, resMap)) {
        throw EQSYS_LUSOL_LUERROR;
    }

    std::vector<double> tmp(varCnt, 0);
    std::vector<double> res(varCnt, 0);

    for (int i = 0; i < varCnt; ++i) {
        double sumTmp = 0;
        for (int j = 0; j < i; ++j) {
            sumTmp += LUres.first[i][j] * tmp[j];
        }
        tmp[i] = RSMfreeTerms[i] - sumTmp;
    }

    for (int i = varCnt - 1; i >= 0; --i) {
        double sumTmp = 0;
        for (int j = i; j < varCnt; ++j) {
            sumTmp += LUres.second[i][j] * res[j];
        }
        res[i] = (tmp[i] - sumTmp) / (LUres.second[i][i]);
    }

    std::vector<double> realRes(varCnt);
    for (int i = varCnt - 1; i >= 0; --i) {
        //переставляем результат в правильном порядке
        realRes[resMap[i]] =  res[i];
    }

    return realRes;
}

std::vector<double> EQsys::TMAsolve() {
    //метод прогонки для компактной записи (3 вектора с коэфициентами на диагоналях)
    if (matr.typef() != "tridiagonalCompact") {
        throw EQSYS_SOLVE_WRONG_TYPE;
    }
    //первый вектор в пачке -- нижний, "a", второй -- центральный, "b", третий -- верхний, "c"

    std::vector<double> alpha(varCnt);
    std::vector<double> betha(varCnt);

    //отыскание коэфициентов
    alpha[0] = - matr[2][0] / matr[1][0];
    betha[0] = freeTerms[0] / matr[1][0];

    for (int i = 1; i < varCnt - 1; ++i) {
        alpha[i] = - matr[2][i] / (matr[1][i] + matr[0][i - 1] * alpha[i - 1]);
        betha[i] = (freeTerms[i] - matr[0][i - 1] * betha[i - 1]) / (matr[1][i] + matr[0][i - 1] * alpha[i - 1]);
        //первый вектор берется по коэфициентам i-1, т.к. из-за структуры матрицы он записан с увеличением коэфициента на 1
    }
    //пояснение к коэфициентам векторов
    /*
      1 1 0 0 0 0 ...
      1 2 2 0 0 0 ...
      0 2 3 3 0 0 ...
      ........... ...
    */

    std::vector<double> res(varCnt, 0);
    //значение последнего неизвестного
    res[varCnt - 1] = (freeTerms[varCnt - 1] - matr[0][varCnt - 2] * betha[varCnt - 2]) / (matr[1][varCnt - 1] + matr[0][varCnt - 2] * alpha[varCnt - 2]);
    //отычкание остальных неизвестных
    for (int i = varCnt - 2; i >= 0; --i) {
        res[i] = res[i + 1] * alpha[i] + betha[i];
    }

    return res;
}


std::pair<std::vector<double>, double> EQsys::FPIsolve(double eps) {
    if (matr.typef() == "tridiagonalCompact") {
        throw EQSYS_SOLVE_WRONG_TYPE;
    }
    //алгоритм простых итераций (fixed point iteration)
    if (!matr.JkbCompatable()) { //проверка сходимости
        throw EQSYS_SOLVE_FPI_NOT_APPLYABLE;
    }

    std::pair<SQmatrix, std::vector<double>> casted = matr.JkbCast(freeTerms); //выполнили преобразование
    std::vector<double> res = casted.second;
    std::vector<double> rp = res - std::vector<double>(res.size(), 1); //значния res на пердыдущем шаге.

    if (trace) {
        casted.first.print();
    }

    double alphaNorm = casted.first.norm();//лучше максимальная сумма строки
    double stepsEst = (log(eps) - VNorm(casted.second) + log(1 - alphaNorm)) / log(alphaNorm) - 1;
    if(trace) printf("steps estimated: %.0lf\n", stepsEst);

    int cnt = 0;
    while (VNorm(res - rp) > eps) {
        ++cnt;
        if (trace) {
            printf("step %d | X: ", cnt);
            for (int i = 0; i < varCnt; ++i) {
                printf("%lf ", res[i]);
            }
            printf("\n");
        }

        rp = res;
        res = casted.second + casted.first.mult(casted.first, res);
    }


    ++cnt;
    if (trace) {
        printf("step %d | X: ", cnt);
        for (int i = 0; i < varCnt; ++i) {
            printf("%lf ", res[i]);
        }
        printf("| final\n");


        printf("norm of alpha: %lf\n", alphaNorm);
        VecPrint(res - rp);
        printf("norm of vector: %lf\n", VNorm(res - rp));
    }

    double error = VNorm(res - rp) * alphaNorm / (1 - alphaNorm);
    return std::pair<std::vector<double>, double> (res,error);
}


//метод Гаусса-Зейделя
std::pair<std::vector<double>, double> EQsys::GSsolve(double eps) {
    if (matr.typef() == "tridiagonalCompact") {
        throw EQSYS_SOLVE_WRONG_TYPE;
    }
    if (!matr.JkbCompatable()) { //проверка сходимости
        throw EQSYS_SOLVE_FPI_NOT_APPLYABLE;
    }

    std::pair<SQmatrix, std::vector<double>> casted = matr.JkbCast(freeTerms); //выполнили преобразование
    std::pair<SQmatrix, SQmatrix> sliced = casted.first.slice();
    SQmatrix A = sliced.first.dif(SQmatrix(varCnt, "ident"), sliced.first).inverse();
    SQmatrix B = sliced.second;
    std::vector<double> res = casted.second;
    std::vector<double> rp = res - std::vector<double>(res.size(), 1); //значния res на пердыдущем шаге.
    SQmatrix Alpha = A.mult(A, B);
    std::vector<double> Betha = A.mult(A, casted.second);

    double alphaNorm = Alpha.norm();
    double stepsEst = (log(eps) - VNorm(Betha) + log(1 - alphaNorm)) / log(alphaNorm) - 1;
    if (trace) printf("steps estimated: %.0lf\n", stepsEst);

    int cnt = 0;
    while (VNorm(res - rp) > eps) {
        ++cnt;
        if (trace) {
            printf("step %d | X: ", cnt);
            for (int i = 0; i < varCnt; ++i) {
                printf("%lf ", res[i]);
            }
            printf("\n");
        }
        rp = res;
        res = Betha + casted.first.mult(Alpha, res);
    }

    ++cnt;
    if (trace) {
        printf("step %d | X: ", cnt);
        for (int i = 0; i < varCnt; ++i) {
            printf("%lf ", res[i]);
        }
        printf("| final\n");
    }

    double error = VNorm(res - rp) * B.norm() / (1 - Alpha.norm());
    return std::pair<std::vector<double>, double> (res,error);
}

//устаревший
std::vector<double> EQsys::Solve() {
    if (matr.typef() == "tridiagonalCompact") {
        return (*this).TMAsolve();
    }
    else {
        return (*this).LUsolve();
    }
}