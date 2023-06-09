#pragma once
#include <complex>
#include <vector>

#define LIMIT 200
#define EPS 0.001

#define ERR_MATR_IRR_SIZE 15
#define ERR_MATR_SIZE_MISSMATCH 12
#define ERR_MATR_TYPE_MISSMATCH 13
#define ERR_MATR_NOT_MATRIX 11
#define ERR_MATR_NULL_SIZE 10
#define ERR_MATR_MULT_SIZE_MISSMATCH 101
#define ERR_MATR_JKB_NOT_CASTABLE 301
#define ERR_MATR_NOT_INVERSIBLE 201
#define ERR_MATR_JRM_NO_CONVERGE 501
#define ERR_MATR_QR_NO_CONVERGE 502


#define ERR_VECTOR_MULT_SIZE_MISSMATCH 1012

class SQmatrix {
    //класс опиывает квадратную матрицу
private:
    bool trace;
    int size;
    std::string type;
    std::vector<std::vector<double>> matrix;
public:
    // онструкторы:
    SQmatrix();

    SQmatrix(int s);

    SQmatrix(int s, std::vector<std::vector<double>>& v);

    SQmatrix(int s, std::string type);

    SQmatrix(int s, std::vector<std::vector<double>>& v, std::string t);

    // омпоненты:

    bool tracef();

    void traceSwitch();

    const int sizef();

    const std::string typef();

    const std::vector<std::vector<double>> matrixf();

    //ѕровекрка:
    bool check();

    //визуализаци€:
    void print();

    //замена строк
    void swap(int lhs, int rhs);

    //координата максимального элемента
    std::pair<int, int> maxND();

    //посроить матрицу вращени€
    SQmatrix rotation(int i, int j, double phi);

    //транспонировать матрицу
    SQmatrix transpose();

    //арифметические операции:
    SQmatrix sum(const SQmatrix& lhs, const SQmatrix& rhs);

    SQmatrix dif(const SQmatrix& lhs, const SQmatrix& rhs);

    SQmatrix mult(const SQmatrix& lhs, const SQmatrix& rhs);

    SQmatrix mult(const SQmatrix& lhs, const double rhs);

    std::vector<double> mult(const SQmatrix& lhs, const std::vector<double>& rhs);

    //корень из суммы квадратов внедиагональных элементов (нужно только дл€ JRM)
    double nonDiagSqSm();

    //вернуть вектор диагональных элементов матрицы
    std::vector<double> diagonal();

    //метод вращений якоби (отыскание собственных значений и собственных векторов методом вращений
    std::pair < std::vector<double>, std::vector<std::vector<double>>> JRM(double acc);

    //корень из суммы квадратов поддиагональных элементов (дл€ метода QR)
    double underDiagSqSm();

    //отыскание собственных значений с помощью QR разложени€
    std::vector<std::complex<double>> QReigen(double prec);

    //помен€ть строки местами так, чтобы на главной диагонали не было нулевых элементов
    std::pair<SQmatrix, std::pair<std::vector<double>, std::vector<int>>> arrange(std::vector<double> ft) const;

    //разделить матрицу на нижнюю диагональную и верхнюю, чтобы L + U = A
    std::pair<SQmatrix, SQmatrix> slice();

    //LU разложение
    std::pair<std::pair<SQmatrix, SQmatrix>, std::vector<int>> LUdecomp(std::vector<double>& ft);

    bool LUcheck(const std::pair<SQmatrix, SQmatrix>& inp, double prec,std::vector<int> resMap);

    //QR разложение
    std::pair<SQmatrix, SQmatrix> QRdecomp();

    //приведение к специальному виду дл€ решеени€ методом простых итераций
    const std::pair<SQmatrix, std::vector<double>> JkbCast(const std::vector<double> freeTerms);

    //проверка применимости метода простых итераций
    bool JkbCompatable();

    //вычисление обратной матрицы методом ∆ордана-√аусса
    SQmatrix inverse();

    //вычисление определител€ методом √аусса
    double det();

    //вычисление нормы матрицы
    double norm();

    void swapf(int lhs, int rhs);

    //проверка правильности собственных значений и векторов
    bool eigenCheck(std::pair<std::vector<double>,std::vector<std::vector<double>>> eigen);

    std::vector<double>& operator[](int i);
};
