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
    //����� ��������� ������� ��������� � �������� ����
private:
    int varCnt;
    SQmatrix matr;
    std::vector<double> freeTerms;
    bool trace;

public:

    //������ �������
    EQsys();

    //������������� �������
    EQsys(int s);

    //�������, � �������� �������� �������� (�������� ����� ������)
    //������� ���������
    EQsys(int s, std::vector<std::vector<double>>& v);

    //�������, � �������� �������� �������� (�������� ����� �������)
    //������� ���������
    EQsys(int s, SQmatrix m);

    //�������, � �������� �������� �������� (�������� ����� ������)
    //������� �����������
    EQsys(int s, std::vector<std::vector<double>>& v, std::vector<double>& b);

    //�������, � �������� �������� �������� (�������� ����� �������)
    //������� �����������
    EQsys(int s, SQmatrix m, std::vector<double>& b);

    //�������� �� ����������� ����������
    bool tracef();

    void traceSwitch();

    SQmatrix matrixf() const;

    std::vector<double> freeTermsf() const;

    //������ ���� � ������� LU ����������
    std::vector<double> LUsolve(); 

    //������ ���� ������� ��������
    std::vector<double> TMAsolve();

    //������ ���� ������� ������� ��������
    std::pair<std::vector<double>, double> FPIsolve(double eps);

    //������ ���� ������� ������-�������
    std::pair<std::vector<double>, double> GSsolve(double eps);

    //������������� �������� ����� �������
    std::vector<double> Solve();

};