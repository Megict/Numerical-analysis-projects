#pragma once
#include <vector>

const std::vector<double> operator*(const double left, const std::vector<double>& right);
const std::vector<double> operator*(const std::vector<double>& left, const double right);
const std::vector<double> operator+(const std::vector<double>& left, const std::vector<double>& right);
const std::vector<double> operator-(const std::vector<double>& left, const std::vector<double>& right);
const std::vector<double> operator/(const std::vector<double>& left, const double right);
const bool operator==(const std::vector<double>& left, const std::vector<double>& right);
void VecSwap(std::vector<int>& vec, int i, int j);

double VNorm(std::vector<double> inp);
void VecPrint(std::vector<double> inp);