#include <vector>
#include "VectorExten.h"
#define EPS 0.01


//���������� �������� � ���������
const std::vector<double> operator*(const double left, const std::vector<double>& right) { //��������� �� �����
    std::vector<double> res(0);
    for (double a : right) {
        res.push_back(left * a);
    }
    return res;
}

const std::vector<double> operator*(const std::vector<double>& left, const double right) { //��������� �� �����
    std::vector<double> res(0);
    for (double a : left) {
        res.push_back(right * a);
    }
    return res;
}

const std::vector<double> operator/(const std::vector<double>& left, const double right) { //������� �� �����
    std::vector<double> res(0);
    for (double a : left) {
        res.push_back(a / right);
    }
    return res;
}


const std::vector<double> operator+(const std::vector<double>& left, const std::vector<double>& right) { //�������� �������� ���������� �����
    if (left.size() != right.size()) {
        throw(1);
    }
    std::vector<double> res(left.size());
    for (size_t i = 0; i < left.size(); ++i) {
        res[i] = left[i] + right[i];
    }
    return res;
}

const std::vector<double> operator-(const std::vector<double>& left, const std::vector<double>& right) { //��������� �������� ���������� �����
    if (left.size() != right.size()) {
        throw(1);
    }
    std::vector<double> res(left.size());
    for (size_t i = 0; i < left.size(); ++i) {
        res[i] = left[i] - right[i];
    }
    return res;
}

const bool operator==(const std::vector<double>& left, const std::vector<double>& right) { //��������� �������� ���������� �����
    if (left.size() != right.size()) {
        throw(1);
    }

    for (size_t i = 0; i < left.size(); ++i) {

        if (abs(right[i] - left[i]) > EPS) {
            return false;
        }
    }
    return true;
}

void VecSwap(std::vector<int>& vec, int i, int j) {
    if (i >= vec.size() || j >= vec.size() || i < 0 || j < 0) {
        throw (1);
    }

    int k = vec[i];
    vec[i] = vec[j];
    vec[j] = k;
}

double VNorm(std::vector<double> inp) { //����� �������
    double max = 0;
    for (double a : inp) {
        if (abs(a) > max) {
            max = abs(a);
        }
    }
    return max;
}

void VecPrint(std::vector<double> inp) {
    for(double a : inp) {
        printf("%lf ", a);
    }
    printf("\n");
}