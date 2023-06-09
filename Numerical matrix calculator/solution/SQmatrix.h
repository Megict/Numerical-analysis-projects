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
    //����� �������� ���������� �������
private:
    bool trace;
    int size;
    std::string type;
    std::vector<std::vector<double>> matrix;
public:
    //������������:
    SQmatrix();

    SQmatrix(int s);

    SQmatrix(int s, std::vector<std::vector<double>>& v);

    SQmatrix(int s, std::string type);

    SQmatrix(int s, std::vector<std::vector<double>>& v, std::string t);

    //����������:

    bool tracef();

    void traceSwitch();

    const int sizef();

    const std::string typef();

    const std::vector<std::vector<double>> matrixf();

    //���������:
    bool check();

    //������������:
    void print();

    //������ �����
    void swap(int lhs, int rhs);

    //���������� ������������� ��������
    std::pair<int, int> maxND();

    //�������� ������� ��������
    SQmatrix rotation(int i, int j, double phi);

    //��������������� �������
    SQmatrix transpose();

    //�������������� ��������:
    SQmatrix sum(const SQmatrix& lhs, const SQmatrix& rhs);

    SQmatrix dif(const SQmatrix& lhs, const SQmatrix& rhs);

    SQmatrix mult(const SQmatrix& lhs, const SQmatrix& rhs);

    SQmatrix mult(const SQmatrix& lhs, const double rhs);

    std::vector<double> mult(const SQmatrix& lhs, const std::vector<double>& rhs);

    //������ �� ����� ��������� ��������������� ��������� (����� ������ ��� JRM)
    double nonDiagSqSm();

    //������� ������ ������������ ��������� �������
    std::vector<double> diagonal();

    //����� �������� ����� (��������� ����������� �������� � ����������� �������� ������� ��������
    std::pair < std::vector<double>, std::vector<std::vector<double>>> JRM(double acc);

    //������ �� ����� ��������� ��������������� ��������� (��� ������ QR)
    double underDiagSqSm();

    //��������� ����������� �������� � ������� QR ����������
    std::vector<std::complex<double>> QReigen(double prec);

    //�������� ������ ������� ���, ����� �� ������� ��������� �� ���� ������� ���������
    std::pair<SQmatrix, std::pair<std::vector<double>, std::vector<int>>> arrange(std::vector<double> ft) const;

    //��������� ������� �� ������ ������������ � �������, ����� L + U = A
    std::pair<SQmatrix, SQmatrix> slice();

    //LU ����������
    std::pair<std::pair<SQmatrix, SQmatrix>, std::vector<int>> LUdecomp(std::vector<double>& ft);

    bool LUcheck(const std::pair<SQmatrix, SQmatrix>& inp, double prec,std::vector<int> resMap);

    //QR ����������
    std::pair<SQmatrix, SQmatrix> QRdecomp();

    //���������� � ������������ ���� ��� �������� ������� ������� ��������
    const std::pair<SQmatrix, std::vector<double>> JkbCast(const std::vector<double> freeTerms);

    //�������� ������������ ������ ������� ��������
    bool JkbCompatable();

    //���������� �������� ������� ������� �������-������
    SQmatrix inverse();

    //���������� ������������ ������� ������
    double det();

    //���������� ����� �������
    double norm();

    void swapf(int lhs, int rhs);

    //�������� ������������ ����������� �������� � ��������
    bool eigenCheck(std::pair<std::vector<double>,std::vector<std::vector<double>>> eigen);

    std::vector<double>& operator[](int i);
};
