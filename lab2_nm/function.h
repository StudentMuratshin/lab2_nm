#pragma once
#include <cmath>
#include <vector>
using namespace std;

double fx(double x);
double ft(double t);
double Integrals_rectangle(double h, const vector<pair<double, double>>& T);
double Integrals_trap(double h, double a, double b, const vector<pair<double, double>>& T);
double Integrals_simpson(double h, vector<pair<double, double>>& T);
vector<pair<double, double>> make_table(int N, const double a, const double b, double H);
double Legendre(double x, double i);
double Derivative(double x, double i);
double root(double i, double k, vector <double> &roots, double n);
double weights(double x, double n);
double I_tr(double xi, double xi_1, double h);
double I(double xi, double xi_1, double h);
double Eps(double xi, double xi_1, double h);
pair <double, int> opt_table(double a, double b);
double MC(vector<pair<double, double>>& T, double a, double b);