#pragma once
#include <cmath>
#include <vector>
using namespace std;

double fx(double x);
double Integrals_rectangle(double h, const vector<pair<double, double>>& T);
double Integrals_trap(double h, double a, double b, const vector<pair<double, double>>& T);
double Integrals_simpson(double h, vector<pair<double, double>>& T);
vector<pair<double, double>> make_table(int N, const double a, const double b, double H);