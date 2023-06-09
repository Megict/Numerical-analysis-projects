#include <iostream>
#include "equation.h"

std::pair<std::complex<double>, std::complex<double>> solveSQ(double a, double b, double c) {
    double D = b * b - 4 * a * c;
    return std::pair<std::complex<double>, std::complex<double>>
           ((std::complex<double>(-b) + sqrt(std::complex<double>(D))) / (2.0 * std::complex<double>(a)),
            (std::complex<double>(-b) - sqrt(std::complex<double>(D))) / (2.0 * std::complex<double>(a)));
}