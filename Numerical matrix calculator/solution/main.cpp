//Гаврилов М.С. 8О-306
#include <iostream>
#include <fstream>
#include <vector>
#include <complex>
#include "VectorExten.h"
#include "SQmatrix.h"
#include "EQsys.h"
//вариант 7

int main() {


    if(true)
    while (true) {
        int n;
        printf("input number of variables (-1 for tridiagonal matr)\n");
        printf(" |> ");
        std::cin >> n;
        SQmatrix matr;
        if (n == -1) {
            printf("tdm| input number of variables\n");
            printf(" |> ");
            std::cin >> n;
            //записть матрицы в компактном виде
            std::vector<std::vector<double>> inp1 = std::vector<std::vector<double>>(3);
            inp1[0] = std::vector<double>(n - 1, 0);
            inp1[1] = std::vector<double>(n, 0);
            inp1[2] = std::vector<double>(n - 1, 0);
            std::vector<double> inp2 = std::vector<double>(n, 0);

            //заполнение матрицы
            printf("tdm| input matrix\n");
            for (int i = 0; i < n; ++i) {
                printf(" |> ");
                if (i == 0) {
                    std::cin >> inp1[1][i] >> inp1[2][i];
                }
                else
                    if (i == n - 1) {
                        std::cin >> inp1[0][i - 1] >> inp1[1][i];
                    }
                    else
                        std::cin >> inp1[0][i - 1] >> inp1[1][i] >> inp1[2][i];

            }
            //задание матрицы
            matr = SQmatrix(n, inp1, "tridiagonalCompact");
        }
        else {
            std::vector<std::vector<double>> inp1(n, std::vector<double>(n, 0));
            printf("input matrix\n");
            for (int i = 0; i < n; ++i) {
                printf(" |> ");
                for (int j = 0; j < n; ++j) {
                    std::cin >> inp1[i][j];
                }
            }
            //задание матрицы
            matr = SQmatrix(n, inp1);
        }
        int mode;
        while (true) {
            system("cls");
            printf("matrix:\n");
            matr.print();
            if (matr.tracef()) {
                printf(" (tracing matrix calcultions)\n");
            }
            printf("\n create system: 1\n calculate determinant: 2\n calculate inverse: 3\n calculate and check LU: 4\n calculate and check QR: 5\n \n find eigen with rotation method: 7\n find eigenvalues whith QR: 8\n switch trace mode: 9\n new matrix: 0\n\n");
            printf(" |> ");
            std::cin >> mode;
            if (mode == 1) {
                std::vector<double> inp2(n, 0);
                printf("input free terms\n");
                printf(" |> ");
                for (int i = 0; i < n; ++i) {
                    std::cin >> inp2[i];
                }
                EQsys system1(n, matr, inp2);
                
                while (true) {
                    system("cls");
                    printf("matrix of system:\n");
                    system1.matrixf().print();
                    if (matr.tracef()) {
                        printf(" (tracing matrix calcultions)\n");
                    }
                    printf("free terms of system:\n");
                    VecPrint(system1.freeTermsf());
                    if (system1.tracef()) {
                        printf(" (tracing system calcultions)\n");
                    }
                    printf("\n\n solve with LU: 1\n solve with tridiagonal matrix algorithm: 2\n solve with fixed-point iteration: 3\n solve with Gauss-Seidel: 4\n switch trace mode: 9\n return to matrices: 0\n\n");
                    printf(" |> ");
                    int mode;
                    std::cin >> mode;

                    if (mode == 1) {
                        try {
                            std::vector<double> res = system1.LUsolve();
                            printf("\nsolved with LU\n\nresult:\n\n");
                            for (int i = 0; i < n; ++i) {
                                printf("%.3lf\n", res[i]);
                            }
                        }
                        catch (int err) {
                            printf("error occured.\n error code: %d\n", err);
                        }
                        system("pause");
                    }
                    else if (mode == 2) {
                        try{
                            std::vector<double> res = system1.TMAsolve();
                            printf("\nsolved with TMA\n\nresult:\n\n");
                            for (int i = 0; i < n; ++i) {
                                printf("%.3lf\n", res[i]);
                            }
                        }
                        catch (int err) {
                            printf("error occured.\n error code: %d\n", err);
                        }
                        system("pause");
                    }
                    else if (mode == 3) {
                        try {
                            std::pair<std::vector<double>, double> res = system1.FPIsolve(0.01);
                            printf("\nsolved with FPI\n\nresult:\n\n");
                            for (int i = 0; i < n; ++i) {
                                printf("%.3lf\n", res.first[i]);
                            }
                            printf("error:\n%lf\n", res.second);
                        }
                        catch (int err) {
                            printf("error occured.\n error code: %d\n", err);
                        }
                        system("pause");
                    }
                    else if (mode == 4) {
                        try {
                            std::pair<std::vector<double>, double> res2 = system1.GSsolve(0.01);
                            printf("\nsolved with GS\n\nresult:\n\n");
                            for (int i = 0; i < n; ++i) {
                                printf("%.3lf\n", res2.first[i]);
                            }
                            printf("error:\n%lf\n", res2.second);
                        }
                        catch (int err) {
                            printf("error occured.\n error code: %d\n", err);
                        }
                        system("pause");
                    }
                    else if (mode == 9) {
                        system1.traceSwitch();
                    }
                    else if (mode == 0) {
                        break;
                    }
                    else {
                        continue;
                    }
                }


            }
            else if (mode == 2) {
                try {
                    double det = matr.det();
                    printf("determinant: %.2lf\n", det);
                }
                catch (int err) {
                    printf("error occured.\n error code: %d\n", err);
                }
                system("pause");
            }
            else if (mode == 3) {
                try {
                    SQmatrix res2 = matr.inverse();
                    printf("inverse matrix:\n");
                    res2.print();
                }
                catch (int err) {
                    printf("error occured.\n error code: %d\n", err);
                }
                system("pause");
            }
            else if (mode == 4) {
                try {
                    //вывод lu разложения
                    std::vector<double>uselessInput(matr.sizef(), 0);
                    std::pair<std::pair<SQmatrix, SQmatrix>, std::vector<int>> uselessResult = matr.LUdecomp(uselessInput);
                    std::pair<SQmatrix, SQmatrix> LUres = uselessResult.first;

                    printf("LU decomposition:\n");
                    LUres.first.print();
                    printf("\n");
                    LUres.second.print();
                    printf("\nLU decomposition multiplied:\n");
                    SQmatrix pr = LUres.first.mult(LUres.first, LUres.second);
                    pr.print();
                    bool ch = matr.LUcheck(LUres, 0.01, uselessResult.second);
                    if (ch) {
                        printf(" | autocheck: correct\n");
                    }
                    else {
                        printf(" | autocheck: wrong\n");
                    }
                }
                catch (int err) {
                    printf("error occured.\n error code: %d\n", err);
                }
                system("pause");
            }
            else if (mode == 5) {
                try {
                    //вывод йк разложения
                    std::pair<SQmatrix, SQmatrix> QRres = matr.QRdecomp();
                    printf("QR decomposition:\n");
                    QRres.first.print();
                    printf("\n");
                    QRres.second.print();
                    printf("\nQR decomposition multiplied:\n");
                    SQmatrix pr = QRres.first.mult(QRres.first, QRres.second);
                    pr.print();

                    std::vector<int> thisIsNotAnything(matr.sizef());
                    for (int i = 0; i < matr.sizef(); ++i) {
                        thisIsNotAnything[i] = i;
                    }
                    bool ch = matr.LUcheck(QRres, 0.01, thisIsNotAnything); //можно юзать и для QR, ибо он просто матрицы перемножает и сравнивает с оригом
                    if (ch) {
                        printf(" | autocheck: correct\n");
                    }
                    else {
                        printf(" | autocheck: wrong\n");
                    }
                }
                catch (int err) {
                    printf("error occured.\n error code: %d\n", err);
                }
                system("pause");
                }
            else if (mode == 6) {
                
            }
            else if (mode == 7) {
                //поиск собственных вектроров и числе с помощью метода вращений
                try {
                    std::pair < std::vector<double>, std::vector<std::vector<double>>> res = matr.JRM(0.01);

                    printf("eigenvalues:\n");
                    for (int i = 0; i < n; ++i) {
                        printf("%lf ", res.first[i]);
                    }
                    printf("\neigenvectors:\n");
                    for (int i = 0; i < n; ++i) {
                        printf(" | ");
                        for (int j = 0; j < n; ++j) {
                            printf("%lf ", res.second[i][j]);
                        }
                        printf("\n");
                    }
                    printf("\n");

                    if (matr.eigenCheck(res)) {
                        printf("eigen correct\n");
                    }
                    else {
                        printf("eigen wrong\n");
                    }
                }
                catch (int err) {
                    printf("error occured.\n error code: %d\n", err);
                }
                system("pause");
            }
            else if (mode == 8) {
                try {
                    //посик собственных векторов через QR разложение
                    std::vector<std::complex<double>> res = matr.QReigen(0.01);

                    printf("eigenvalues:\n");
                    for (int i = 0; i < n; ++i) {
                        std::cout << res[i] << " ";
                    }
                    printf("\n");
                }
                catch (int err) {
                    printf("error occured.\n error code: %d\n", err);
                }
                system("pause");
            }
            else if (mode == 9) {
                matr.traceSwitch();
            }
            else if (mode == 0) {
                system("cls");
                break;
            }
            else {
                continue;
            }
        }
    }



    

    int a;
    printf("test lab (n) ?\n");
    std::cin >> a;

    printf("input number of variables\n");
    int n;
    std::cin >> n;
    std::vector<std::vector<double>> inp1(n, std::vector<double>(n, 0));
    std::vector<double> inp2(n, 0);

    if (a == 1) { //решение уравнений с помощью LU разложения
        printf("input matrix\n");
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                std::cin >> inp1[i][j];
            }
        }
        printf("input free terms\n");
        for (int i = 0; i < n; ++i) {
            std::cin >> inp2[i];
        }
        try {
            //задание матрицы
            SQmatrix matr(n, inp1);
            //задание системы
            EQsys system1(n, matr, inp2);
            //решение системы
            std::vector<double> res = system1.LUsolve();
            printf("---------\n");
            printf("result:\n");
            for (int i = 0; i < n; ++i) {
                printf("%.3lf\n", res[i]);
            }
            printf("\n---------\n");

            //вывод lu разложения
            std::vector<double> uli = system1.freeTermsf();
            std::pair<std::pair<SQmatrix, SQmatrix>, std::vector<int>> uselessResult = system1.matrixf().LUdecomp(uli);
            std::pair<SQmatrix, SQmatrix> LUres = uselessResult.first;

            printf("LU decomposition:\n");
            printf("-----------------------------------------\n");
            LUres.first.print();
            printf("\n");
            LUres.second.print();
            printf("-----------------------------------------\n");
            SQmatrix pr = LUres.first.mult(LUres.first, LUres.second);
            pr.print();

            //отыскание обратной матрицы
            SQmatrix res2 = system1.matrixf().inverse();
            printf("-----------------------------------------\n");
            printf("inverse matrix:\n");
            res2.print();
            printf("-----------------------------------------\n");
            //отыскание определителя
            double det = system1.matrixf().det();
            printf("\ndeterminant: %.2lf\n", det);
        }
        catch (int err) {
            printf("error occured.\n error code: %d", err);
            return 1;
        }

        return 0;
    }
    if (a == 2) { //решение уравнения методом прогонки (на тридиагональной записи)
        inp1 = std::vector<std::vector<double>>(3);
        inp1[0] = std::vector<double>(n - 1, 0);
        inp1[1] = std::vector<double>(n, 0);
        inp1[2] = std::vector<double>(n - 1, 0);
        inp2 = std::vector<double>(n, 0);

        try {
            //заполнение матрицы
            printf("input matrix (compact tridiag)\n");
            for (int i = 0; i < n; ++i) {
                if (i == 0) {
                    std::cin >> inp1[1][i] >> inp1[2][i];
                }
                else
                    if (i == n - 1) {
                        std::cin >> inp1[0][i - 1] >> inp1[1][i];
                    }
                    else
                        std::cin >> inp1[0][i - 1] >> inp1[1][i] >> inp1[2][i];

            }
            //заполнение вектора свободных членов
            printf("input free terms\n");
            for (int i = 0; i < n; ++i) {
                std::cin >> inp2[i];
            }


            //задание матрицы
            SQmatrix matr(n, inp1, "tridiagonalCompact");
            printf("-----------------------------------------\n");
            matr.print();
            printf("-----------------------------------------\n");
            //задание системы
            EQsys system1(n, matr, inp2);
            //решение системы
            std::vector<double> res = system1.TMAsolve();

            printf("---------\n");
            printf("result:\n");
            for (int i = 0; i < n; ++i) {
                printf("%.3lf\n", res[i]);
            }
            printf("\n---------\n");
        }
        catch (int err) {
            printf("error occured.\n error code: %d", err);
            return 1;
        }

        return 0;
    }
    if (a == 3) {//метод простых итераций и Зейделя
        printf("input matrix\n");
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                std::cin >> inp1[i][j];
            }
        }
        printf("input free terms\n");
        for (int i = 0; i < n; ++i) {
            std::cin >> inp2[i];
        }
        try {
            //задание матрицы
            SQmatrix matr(n, inp1);
            //задание системы
            EQsys system1(n, matr, inp2);
            //решение системы
            std::pair<std::vector<double>, double> res = system1.FPIsolve(0.3);
            printf("---------\n");
            printf("result:\n");
            for (int i = 0; i < n; ++i) {
                printf("%.3lf\n", res.first[i]);
            }
            printf("\n---------\n");
            printf("error:\n%lf\n", res.second);

            //вывод jkb разложения
            std::pair<SQmatrix, std::vector<double>> cast = matr.JkbCast(inp2);
            printf("JkbCast:\n");
            printf("-----------------------------------------\n");
            cast.first.print();
            printf("-----------------------------------------\n");
            for (int i = 0; i < n; ++i) {
                printf("%lf\n", cast.second[i]);
            }
            printf("A * B:\n");
            printf("-----------------------------------------\n");

            printf("-----------------------------------------\n");
            printf("-----------------------------------------\n");

            std::pair<std::vector<double>, double> res2 = system1.GSsolve(0.01);
            printf("GS result:\n");
            for (int i = 0; i < n; ++i) {
                printf("%.3lf\n", res2.first[i]);
            }
            printf("\n---------\n");
            printf("error:\n%lf\n", res2.second);
        }
        catch (int err) {
            printf("error occured.\n error code: %d", err);
            return 1;
        }

    }
    if (a == 4) {//метод вращений


        printf("input matrix\n");
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                std::cin >> inp1[i][j];
            }
        }
        
        try {
            //задание матрицы
            SQmatrix matr(n, inp1);


            std::pair < std::vector<double>, std::vector<std::vector<double>>> res = matr.JRM(0.01);
            
            printf("eignvalues:\n");
            for (int i = 0; i < n; ++i) {
                printf("%lf ", res.first[i]);
            }
            printf("\neignvectors:\n");
            for (int i = 0; i < n; ++i) {
                printf("|> ");
                for (int j = 0; j < n; ++j) {
                    printf("%lf ", res.second[i][j]);
                }
                printf("\n");
            }


        }
        catch (int err) {
            printf("error occured.\n error code: %d", err);
            return 1;
        }

    }

    if (a == 5) {//QR разложения


        printf("input matrix\n");
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                std::cin >> inp1[i][j];
            }
        }

        try {
            //задание матрицы
            SQmatrix matr(n, inp1);

            std::vector<std::complex<double>> res = matr.QReigen(0.01);

            printf("eigenvalues:\n");
            for (int i = 0; i < n; ++i) {
                std::cout << res[i] << " ";
            }


        }
        catch (int err) {
            printf("error occured.\n error code: %d", err);
            return 1;
        }

    }
        
    return 0;
}
