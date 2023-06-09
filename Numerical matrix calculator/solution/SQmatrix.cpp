#include <iostream>
#include <vector>
#include <algorithm>
#include "SQmatrix.h"
#include "VectorExten.h"
#include "equation.h"

#define EPS 0.001

SQmatrix::SQmatrix() {
    trace = false;
    //пустая матрица
    size = 0;
    type = "normal";
    matrix = std::vector<std::vector<double>>(0, std::vector<double>(0, 0));
}

SQmatrix::SQmatrix(int s) {
    trace = false;
    //незаполненная матрица

    if (s < 0) {
        //задан невозможный размер
        throw ERR_MATR_IRR_SIZE;
    }

    type = "normal";
    size = s;
    matrix = std::vector<std::vector<double>>(s, std::vector<double>(s, 0));
}

SQmatrix::SQmatrix(int s, std::vector<std::vector<double>>& v) {
    trace = false;
    //матрица с заданным заполнением

    size_t sz = v.size();
    if (sz > 0) {
        for (int i = 0; i < sz; ++i) {
            if (v[i].size() != sz) {
                //несоответствие размеров строк и столбцов
                throw ERR_MATR_NOT_MATRIX;
            }
        }
    }
    else {
        //матрица нулевой размерности
        throw ERR_MATR_NULL_SIZE;
    }

    if (s != sz) {
        //заявленная размерность не соответствует реальной
        throw ERR_MATR_SIZE_MISSMATCH;
    }

    size = s;
    matrix = v;
    type = "normal";
}

SQmatrix::SQmatrix(int s, std::string t) {
    trace = false;
    //матрица с заданным типом
    size = s;

    if (s < 0) {
        //задан невозможный размер
        throw ERR_MATR_IRR_SIZE;
    }

    matrix = std::vector<std::vector<double>>(s, std::vector<double>(s, 0));

    if (t == "ident") {
        for (int i = 0; i < size; i++) {
            matrix[i][i] = 1;
        }
    }

    type = "ident";
}

SQmatrix::SQmatrix(int s, std::vector<std::vector<double>>& v, std::string t) {
    trace = false;
    //матрица с заданным типом и заполнением
    size = s;

    //проверка входных данных
    if (t == "tridiagonalCompact") {
        //проверка входных данных
        size_t sz = v.size();
        if (sz == 3) {
            if (v[0].size() != s - 1 ||
                v[1].size() != s ||
                v[2].size() != s - 1) {
                throw ERR_MATR_NOT_MATRIX;
            }
        }
        else {
            //матрица нулевой размерности
            throw ERR_MATR_SIZE_MISSMATCH;
        }
    }
    else {
        size_t sz = v.size();
        if (sz > 0) {
            for (int i = 0; i < sz; ++i) {
                if (v[i].size() != sz) {
                    //несоответствие размеров строк и столбцов
                    throw ERR_MATR_NOT_MATRIX;
                }
            }
        }
        else {
            //матрица нулевой размерности
            throw ERR_MATR_NULL_SIZE;
        }

        if (s != sz) {
            //заявленная размерность не соответствует реальной
            throw ERR_MATR_SIZE_MISSMATCH;
        }
    }

    //трехдиагональная
    if (t == "tridiagonal") {
        //проверка тридиагональности матрицы
        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j) {
                //область, в которой должны быть только нулевые элескнты
                if (!(i == j || i + 1 == j || i - 1 == j)) {
                    if (abs(v[i][j] - 0) > 0.001) {
                        throw ERR_MATR_TYPE_MISSMATCH;
                    }
                }
            }
        }
        type = "tridiagonal";
    }


    if (t == "tridiagonalCompact") {
        matrix = v;
        type = "tridiagonalCompact";
    }
}

bool SQmatrix::check() { //добавить проверку типа
    //проверка корректности заполнения матрицы
    //возвращает true, если все ок

    if (matrix.size() != size) {
        return false;
    }

    for (int i = 0; i < size; ++i) {
        if (matrix[i].size() != size) {
            return false;
        }
    }

    return true;
}

bool SQmatrix::tracef() {
    return trace;
}

void SQmatrix::traceSwitch() {
    trace = !trace;
}

const std::string SQmatrix::typef() {
    return type;
}

const int SQmatrix::sizef() {
    if (matrix.size() == size) {
        return size;
    }
    else {
        if (type == "tridiagonalCompact") {
            return 3;
        }
        return -1;
    }
}

const std::vector<std::vector<double>> SQmatrix::matrixf() {
     return matrix;
}

void SQmatrix::print() {
    std::cout << "square matrix| size: " << (*this).sizef() << "  type: " << (*this).typef() << "\n";
    if (type == "tridiagonalCompact") {
        int sz = matrix[1].size();
        for (int i = 0; i < sz; ++i) {
            for (int j = 0; j < i - 1; ++j) {
                printf(" 0.00 ");
            }
            if (i != 0) {
                if (matrix[0][i - 1] >= 0) {
                    printf(" ");
                }
                printf("%.2lf ", matrix[0][i - 1]);
            }
            if (matrix[1][i] >= 0) {
                printf(" ");
            }
            printf("%.2lf ", matrix[1][i]);
            if (i != sz - 1) {
                if (matrix[2][i] >= 0) {
                    printf(" ");
                }
                printf("%.2lf ", matrix[2][i]);
            }
            for (int j = i + 2; j < sz; ++j) {
                printf(" 0.00 ");
            }
            printf("\n");
        }
    }
    else {
        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j) {
                if (matrix[i][j] >= 0) {
                    printf(" ");
                }
                printf("%.2lf ", matrix[i][j]);
            }
            printf("\n");
        }
    }
}

//замена сроки под номером lhs на строку под номером rhs
void SQmatrix::swap(int lhs, int rhs) {
    if (lhs < 0 || lhs >= size || rhs < 0 || rhs >= size || lhs == rhs) {
        return;
    }
    else {
        std::swap(matrix[lhs], matrix[rhs]);
        return;
    }
}

void SQmatrix::swapf(int lhs, int rhs) {
        return;
}

SQmatrix SQmatrix::sum(const SQmatrix& lhs, const SQmatrix& rhs) {
    //сложение двух матриц.
    // ошибка 1 - несответствие размеров матриц
    if (lhs.size != rhs.size) {
        throw 1;
    }

    SQmatrix res(lhs.size);
    for (int i = 0; i < res.size; ++i) {
        for (int j = 0; j < res.size; ++j) {
            res.matrix[i][j] = lhs.matrix[i][j] + rhs.matrix[i][j];
        }
    }

    return res;
}

SQmatrix SQmatrix::dif(const SQmatrix& lhs, const SQmatrix& rhs) {
    //сложение двух матриц.
    // ошибка 1 - несответствие размеров матриц
    if (lhs.size != rhs.size) {
        throw 1;
    }

    SQmatrix res(lhs.size);
    for (int i = 0; i < res.size; ++i) {
        for (int j = 0; j < res.size; ++j) {
            res.matrix[i][j] = lhs.matrix[i][j] - rhs.matrix[i][j];
        }
    }

    return res;
}

SQmatrix SQmatrix::mult(const SQmatrix& lhs, const SQmatrix& rhs) {
    // перемножение двух матриц.
    // ошибка 1 - несответствие размеров матриц
    if (lhs.size != rhs.size) {
        throw ERR_MATR_MULT_SIZE_MISSMATCH;
    }

    SQmatrix res(lhs.size);
    for (int i = 0; i < res.size; ++i) {
        for (int j = 0; j < res.size; ++j) {
            for (int k = 0; k < res.size; ++k) {
                res.matrix[i][j] += lhs.matrix[i][k] * rhs.matrix[k][j];
            }
        }
    }

    return res;
}

std::vector<double> SQmatrix::mult(const SQmatrix& lhs, const std::vector<double>& rhs) {
    // умножение матрицы на вектор
    // ошибка 1 - несответствие размеров
    if (lhs.size != rhs.size()) {
        throw ERR_MATR_MULT_SIZE_MISSMATCH;
    }

    std::vector<double> res(rhs.size(),0);
    for (int i = 0; i < res.size(); ++i) {
        for (int k = 0; k < res.size(); ++k) {
            res[i] += lhs.matrix[i][k] * rhs[k];
        }
    }

    return res;
}

SQmatrix SQmatrix::mult(const SQmatrix& lhs, const double rhs) {
    // умножение матрицы на вектор
    // ошибка 1 - несответствие размеров

    SQmatrix res(lhs.size);
    for (int i = 0; i < res.sizef(); ++i) {
        for (int k = 0; k < res.sizef(); ++k) {
            res[i][k] = lhs.matrix[i][k] * rhs;
        }
    }

    return res;
}

//разделить матрицу на нижнюю диагональную и верхнюю, чтобы L + U = A
std::pair<SQmatrix, SQmatrix> SQmatrix::slice() {
    SQmatrix L(size);
    SQmatrix U(size);

    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            if (j < i) {
                L[i][j] = matrix[i][j];
            }
            else {
                U[i][j] = matrix[i][j];
            }
        }
    }

    return std::pair<SQmatrix, SQmatrix>(L, U);
}

//координата максимального недиагонального элемента
std::pair<int, int> SQmatrix::maxND() {
    if (type == "tridiagonalCompact") { //может потом сделаю и для этого типа, но сейчас не нужно
        return std::pair<int, int>(-1, -1);
    }

    std::pair<int, int>  maxCoords(0, 1);
    double max = abs(matrix[0][1]);
    for (int i = 0; i < size; ++i) {
        for (int j = i + 1; j < size; ++j) {
            if (i != j && abs(matrix[i][j]) > max) {
                max = abs(matrix[i][j]);
                maxCoords = std::pair<int, int>(i, j);
            }
        }
    }
    return maxCoords;
}

//посроить матрицу вращения для данной матрицы
SQmatrix SQmatrix::rotation(int i, int j, double phi) {
    SQmatrix res(size, "ident");
    res[i][i] = cos(phi);
    res[i][j] = -sin(phi);
    res[j][i] = sin(phi);
    res[j][j] = cos(phi); 

    return res;
}

//транспонировать матрицу
SQmatrix SQmatrix::transpose() {
    SQmatrix res(size);
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            res[j][i] = matrix[i][j];
        }
    }

    return res;
}

//корень из суммы квадратов внедиагональных элементов (для метода вращений)
double SQmatrix::nonDiagSqSm() {
    double sm = 0;
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            if (i != j) {
                sm += matrix[i][j] * matrix[i][j];
            }
        }
    }

    return sqrt(sm/2); //в методе врепщений матрица симметричная, значение симметричной пары учитывается единожды
}

//вернуть вектор диагональных элементов матрицы
std::vector<double> SQmatrix::diagonal() {
    std::vector<double> diag(size);
    for (int i = 0; i < size; ++i) {
        diag[i] = matrix[i][i];
    }

    return diag;
}


//метод вращений Якоби (отыскание собственных значений и собственных векторов методом вращений
std::pair <std::vector<double>, std::vector<std::vector<double>>> SQmatrix::JRM(double acc) {
    if (type == "tridiagonalCompact") { //может потом сделаю и для этого типа, но сейчас не нужно
        throw ERR_MATR_TYPE_MISSMATCH;
    }
    //добавить проверку на симметричность

    std::pair<int, int>maxPos = maxND();
    double angle = (0.5) * atan(2 * matrix[maxPos.first][maxPos.second]/(matrix[maxPos.first][maxPos.first] - matrix[maxPos.second][maxPos.second]));
    SQmatrix rot = rotation(maxPos.first, maxPos.second, angle);
    //применение ротации к матрице
    SQmatrix cur = mult(rot.transpose(), *this);
    cur = mult(cur, rot);

    if (trace) {
        printf("max elm pos:\n");
        printf("%d %d\n", maxPos.first, maxPos.second);

        rot.print();
        cur.print();
    }
    SQmatrix eigen = rot;

    int i = 0;
    while (cur.nonDiagSqSm() > acc) {
        ++i;
        std::pair<int, int>maxPos = cur.maxND();

        if (trace) {
            printf("max elm pos:\n");
            printf("%d %d\n", maxPos.first, maxPos.second);
        }

        angle = (0.5) * atan(2 * cur[maxPos.first][maxPos.second] / (cur[maxPos.first][maxPos.first] - cur[maxPos.second][maxPos.second]));
        rot = rotation(maxPos.first, maxPos.second, angle);
        //применение ротации к матрице
        cur = mult(rot.transpose(), cur);
        cur = mult(cur, rot);
        eigen = mult(eigen, rot);

        if (trace) {
            rot.print();
            cur.print();
        }

        if (i > LIMIT) {
            throw ERR_MATR_JRM_NO_CONVERGE;
        }
    }

    SQmatrix eigenT = eigen.transpose();

    return std::pair<std::vector<double>, std::vector<std::vector<double>>>(cur.diagonal(), eigenT.matrixf());
}


//поменять строки местами так, чтобы на главной диагонали не было нулевых элементов
std::pair<SQmatrix, std::pair<std::vector<double>, std::vector<int>>> SQmatrix::arrange(std::vector<double> ft) const {
    if (type == "tridiagonalCompact") { //может потом сделаю и для этого типа, но сейчас не нужно
        throw ERR_MATR_TYPE_MISSMATCH;
    }

    std::vector<int> arrmap;
    for (int i = 0; i < size; ++i) {
        arrmap.push_back(i);
    }

    SQmatrix res = *this;
    //если на главной диагонали есть нулевые элементы, нужно заменить строки
    for (int i = 0; i < size; ++i) {
        int k = 0;
        //находим строку, в которой на том месте, которое нам нужно наибольший элемент
        double maxval = abs(res[k][i]);
        for (int j = 1; j < size; ++j) { 
            if (abs(res[j][i]) > maxval) {
                maxval = abs(res[j][i]);
                k = j;
            }
        }   

        if (maxval < EPS) {
            throw ERR_MATR_JKB_NOT_CASTABLE;
        }

        res.swap(k, i);
        VecSwap(arrmap, k, i);
        std::swap(ft[k], ft[i]);
    }

    return std::pair<SQmatrix, std::pair<std::vector<double>, std::vector<int>>>(res, std::pair<std::vector<double>,std::vector<int>>(ft,arrmap));
}

std::pair<std::pair<SQmatrix, SQmatrix>, std::vector<int>> SQmatrix::LUdecomp(std::vector<double>& ft) {
    if (type == "tridiagonalCompact") { //может потом сделаю и для этого типа, но сейчас не нужно
        throw ERR_MATR_TYPE_MISSMATCH;
    }
    SQmatrix L(size);
    SQmatrix U((*this));

    std::vector<int> arrmap;
    for (int i = 0; i < size; ++i) {
        arrmap.push_back(i);
    }

    for (int i = 0; i < size; ++i) {
        int k = 0;
        //находим строку, в которой на том месте, которое нам нужно наибольший элемент
        double maxval = abs(U[k][i]);
        for (int j = 1; j < size; ++j) {
            if (abs(U[j][i]) > maxval) {
                maxval = abs(U[j][i]);
                k = j;
            }
        }

        if (maxval < EPS) {
            throw ERR_MATR_JKB_NOT_CASTABLE;
        }

        U.swap(k, i);
        VecSwap(arrmap, k, i);
        std::swap(ft[k], ft[i]);
    }


    for (int i = 0; i < size; ++i) {
        //элементы на диагонали = 1
        L[i][i] = 1;
    }

    for (int j = 0; j < size - 1; ++j) {

        //Выполняется приведение матрицы U к верх.треуг. виду.
        //Коэфициенты, на которые для этого умножаются строки записываются в L.

        int k = 0;

        //находим строку, в которой на том месте, которое нам нужно наибольший элемент
        double maxval = abs(U[k][j]);
        for (int i = 1; i< size; ++i) {
            if (abs(U[j][i]) > maxval) {
                maxval = abs(U[j][i]);
                k = j;
            }
        }

        if (maxval < EPS) {
            throw ERR_MATR_JKB_NOT_CASTABLE;
        }

        U.swapf(k, j);

        for (int i = j + 1; i < size; ++i) {
            //вычисление коэфициента, на который надо умножить строку U[j], 
            //чтобы при сложении с U[i] эл-т U[i][j] обратился в 0
            double fst = U[i][j] / U[j][j];

            //запись коэфициента в L
            L[i][j] = fst;
            //выполнение сложения
            U[i] = U[i] - (U[j] * fst);
        }
    }

    return std::pair<std::pair<SQmatrix, SQmatrix>,std::vector<int>>(std::pair<SQmatrix, SQmatrix>(L, U), arrmap);
}

//на вход идет результат LU разложения, проверяется, равно ли произведение L и U изнач. матрицы
bool SQmatrix::LUcheck(const std::pair<SQmatrix, SQmatrix>& inp, double prec,std::vector<int> resMap){


    SQmatrix preRes = (*this).mult(inp.first, inp.second);

    SQmatrix multRes(size);
    for (int i = size - 1; i >= 0; --i) {
        //переставляем результат в правильном порядке
        multRes[resMap[i]] = preRes[i];
    }

    if (trace) {
        multRes.print();
    }

    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            if (multRes[i][j] - matrix[i][j] > prec) {
                return false;
            }
        }
    }
    //multRes.print();
    return true;
}

//умножение вектора не вектор (вертикального на горизонтальный)
SQmatrix multVH(std::vector<double> lhs, std::vector<double> rhs) {
    if (lhs.size() != rhs.size()) {
        throw ERR_VECTOR_MULT_SIZE_MISSMATCH;
    }
    int size = lhs.size();
    SQmatrix res(size);
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            res[i][j] = lhs[i] * rhs[j];
        }
    }
    return res;
}

//умножение вектора не вектор (горизонтального на вертикальный)
double multHV(std::vector<double> lhs, std::vector<double> rhs) {
    if (lhs.size() != rhs.size()) {
        throw ERR_VECTOR_MULT_SIZE_MISSMATCH;
    }
    int size = lhs.size();
    double res = 0;
    for (int i = 0; i < size; ++i) {
        res += lhs[i] * rhs[i];
    }
    return res;
}

//создание матрицы хаусхолдера на основе заданного вектора
SQmatrix hausholder(std::vector<double> vec) {
    int size = vec.size();
    SQmatrix res(size, "ident");


    return res.sum(res,  res.mult(multVH(vec, vec), - 2.0 / multHV(vec, vec)));
}

short sign(double a) {
    return a > 0 ? 1 : a < 0 ? -1 : 0;
}

//QR разложение матрицы
std::pair<SQmatrix, SQmatrix> SQmatrix::QRdecomp() {
    if (type == "tridiagonalCompact") { //может потом сделаю и для этого типа, но сейчас не нужно
        throw ERR_MATR_TYPE_MISSMATCH;
    }
    //Q - ортогональная мтарица, Q^-1 = QT
    SQmatrix Q(size);
    SQmatrix S = *this;
    if (trace) S.print();
    for (int i = 0; i < size - 1; ++i) {
        std::vector<double> vec(size, 0);  
        double sumfst = 0;

        for (int j = i; j < size; ++j) {
            sumfst += S.matrix[j][i]* S.matrix[j][i];
        }

        vec[i] = S.matrix[i][i] + sign(S.matrix[i][i])*sqrt(sumfst);

        for (int j = i + 1; j < size; ++j) {
            vec[j] = S.matrix[j][i];
        }

        SQmatrix H = hausholder(vec);
        if (i == 0) {
            Q = H;
        }
        else {
            Q = Q.mult(Q, H);
        }
        //H.print();
        S = S.mult(H,S);

        //S.print();
    }
    //Q.print();
    //S.print();

    //Q.mult(Q, S).print();

    return std::pair<SQmatrix, SQmatrix>(Q,S);
}

//корень из суммы квадратов поддиагональных элементов не последнего столбца (для метода QR)
double SQmatrix::underDiagSqSm() {
    double sm = 0;
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            if (i > j + 1) {
                sm += matrix[i][j] * matrix[i][j];
            }
        }
    }

    return sqrt(sm / 2); //в методе врепщений матрица симметричная, значение симметричной пары учитывается единожды
}

//отыскание собственных значений с помощью QR-разложения 
std::vector<std::complex<double>> SQmatrix::QReigen(double prec){
    if (type == "tridiagonalCompact") { //может потом сделаю и для этого типа, но сейчас не нужно
        throw ERR_MATR_TYPE_MISSMATCH;
    }
    std::pair<SQmatrix, SQmatrix> QR = QRdecomp();
    SQmatrix Q = QR.first,
             R = QR.second;
    SQmatrix RmQ = mult(R, Q);

    int i = 0;
    while (RmQ.underDiagSqSm() > prec) {
        ++i;
        QR = RmQ.QRdecomp();
        RmQ = mult(QR.second, QR.first);

        if(trace) RmQ.print();

        if (i > LIMIT) {
            throw ERR_MATR_QR_NO_CONVERGE;
        }
    }

    std::vector<std::complex<double>> eigenvalues(size);
    for (int i = 0; i < size; ++i) {
        if (i == size - 1) {
            eigenvalues[i] = RmQ[i][i];
            break;
        }

        if (pow(RmQ[i + 1][i],2) < prec) {
            eigenvalues[i] = RmQ[i][i];
        }
        else {//если найден блок, то надо вычислить два сз, соответствующих этому блоку и пропустить следующую ячейку
            std::pair<std::complex<double>, std::complex<double>> val = solveSQ(1, -(RmQ[i + 1][i + 1] + RmQ[i][i]),
                RmQ[i + 1][i + 1] * RmQ[i][i] - RmQ[i + 1][i] * RmQ[i][i + 1]);
            if (trace) std::cout << " cc " << val.first << " " << val.second << "\n";

            eigenvalues[i] = val.first;
            eigenvalues[i + 1] = val.second;

            ++i;
        }
    }

    return eigenvalues;
}




//принимает вектор свободных членов, возвращает приведеднный вид матрицы и вектора, позволяющий перейти к решению методом свободных итераций.
const std::pair<SQmatrix, std::vector<double>> SQmatrix::JkbCast(const std::vector<double> freeTerms) {
    if (type == "tridiagonalCompact") { //может потом сделаю и для этого типа, но сейчас не нужно
        throw ERR_MATR_TYPE_MISSMATCH;
    }
    SQmatrix matr = (*this);
    std::vector<double> ft = freeTerms; //работаем с копиями, чтобы не поломать ориг. матрицы.
    //если на главной диагонали есть нулевые элементы, нужно заменить строки
    for (int i = 0; i < size; ++i) {
        if (matr[i][i] == 0) {
            int j = 0;
            //находим строку, в которой на том месте, которое нам нужно не 0
            while (matr[j][i] == 0 /*в новой строке должен быть не ноль в нашем ГД элементе*/|| 
                   matr[i][j] == 0 /*в нашей строке не должно быть 0 в ГД элементе новой*/) {
                ++j; 
                if (j == matr.sizef()) {
                    //такой строки не нашли
                    throw ERR_MATR_JKB_NOT_CASTABLE;
                }
            }
            matr.swap(j, i);
            std::swap(ft[j], ft[i]);
        }
    }

    //само преобразование
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            if (i != j) {
                matr[i][j] = matr[i][j] / matr[i][i];
                matr[i][j] *= -1;
            }
        }
        ft[i] = ft[i] / matr[i][i];
        matr[i][i] = 0;
    }

    return std::pair<SQmatrix, std::vector<double>> (matr,ft);
}

//проверка диагонального преобладания
bool SQmatrix::JkbCompatable() {
    for (int i = 0; i < size; ++i) {
        double sumRow = 0;
        for (int j = 0; j < size; ++j) {
            if (j != i) {
                sumRow += abs(matrix[i][j]);
            }
        }

        if (abs(matrix[i][i]) <= sumRow) {
            return false;
        }
    }

    return true;
}


//вычисление обратной матрицы методом Жордана-Гаусса
SQmatrix SQmatrix::inverse() {
    if (type == "tridiagonalCompact") { //может потом сделаю и для этого типа, но сейчас не нужно
        throw ERR_MATR_TYPE_MISSMATCH;
    }
    if (det() == 0) {
        throw ERR_MATR_NOT_INVERSIBLE;
    }

    SQmatrix D(*this);
    //задание единичной матрицы
    SQmatrix A(size);
    for (int i = 0; i < size; ++i) {
        A[i][i] = 1;
    }

    for (int j = 0; j < size; ++j) {

        //Выполняется приведение матрицы D к верх.треуг. виду.

        //предварительно угловые элементы приводятся к 1.
        A[j] = A[j] / D[j][j];
        D[j] = D[j] / D[j][j];

        for (int i = j + 1; i < size; ++i) {
            //вычисление коэфициента, на который надо умножить строку 
            //чтобы при сложении первый эл-т обратился в 0
            double fst = D[i][j] / D[j][j];

            //выполнение сложения
            D[i] = D[i] - (D[j] * fst);
            A[i] = A[i] - (A[j] * fst);
        }

    }

    for (int j = size - 1; j >= 0; --j) {

        //Выполняется приведение матрицы D к диаг. виду

        for (int i = j - 1; i >= 0; --i) {

            //вычисление коэфициента, на который надо умножить строку 
            //чтобы при сложении первый эл-т обратился в 0
            double fst = D[i][j] / D[j][j];

            //выполнение сложения
            D[i] = D[i] - (D[j] * fst);
            A[i] = A[i] - (A[j] * fst);
        }

    }

    return A;
}


double SQmatrix::det() {
    if (type == "tridiagonalCompact") { //может потом сделаю и для этого типа, но сейчас не нужно
        throw ERR_MATR_TYPE_MISSMATCH;
    }
    SQmatrix D(*this);

    for (int j = 0; j < size; ++j) {

        //Выполняется приведение матрицы D к верх.треуг. виду.

        for (int i = j + 1; i < size; ++i) {
            //вычисление коэфициента, на который надо умножить строку 
            //чтобы при сложении первый эл-т обратился в 0
            if (D[j][j] == 0) {
                return 0;
            }

            double fst = D[i][j] / D[j][j];

            //выполнение сложения
            D[i] = D[i] - (D[j] * fst);
        }
    }

    double det = 1;
    //рассчет определителя, как произведения элементов на главной диагонали
    for (int i = 0; i < size; ++i) {
        det *= D[i][i];
    }
    return det;
}

double SQmatrix::norm() {
    double smmax = 0;
    for (std::vector <double> itr : matrix) {
        double sm = 0;
        for (double elm : itr) {
            sm += abs(elm);
        }
        if (sm > smmax) {
            smmax = sm;
        }
    }

    return smmax;
}

bool SQmatrix::eigenCheck(std::pair<std::vector<double>, std::vector<std::vector<double>>> eigen) {
    if (eigen.first.size() != eigen.first.size()) {
        throw ERR_MATR_SIZE_MISSMATCH;
    }
    if (eigen.first.size() != size) {
        throw ERR_MATR_SIZE_MISSMATCH;
    }

    for (int i = 0; i < size; ++i) {
        std::vector <double> resL = eigen.second[i] * eigen.first[i];
        std::vector <double> resR = mult(*this, eigen.second[i]);

        if (trace) {
            VecPrint(resL);
            VecPrint(resR);
            printf("\n");
        }

        if (!(resL == resR)) {
            return false;
        }
    }

    return true;
}



std::vector<double>& SQmatrix::operator[](int i) {
    return matrix[i];
}
