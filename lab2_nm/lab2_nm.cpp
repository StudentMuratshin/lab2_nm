#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include <numeric>
#include <cmath>
#include <Eigen/Dense>
#include "function.h"

constexpr auto M_PI = 3.14159265358979323846;
constexpr auto true_value = 2.69632737829770440658;

int main()
{
	
	const double a = 0, b = 3 * M_PI;
	const double H = abs(a - b) / 1e5;
	vector <vector<pair<double, double>>> table { make_table(1e5, 0, 3 * M_PI, H) };
	//vector<pair<double, double>> table1 = make_table(1e5, 0, 3 * M_PI, H);

	//прямоугольник
	std::cout << std::setprecision(15);
	cout << "rectangle integral: " << Integrals_rectangle(H, table[0]) << endl;
	// трапеция
	cout << "trapezoid integral: " << Integrals_trap(H, a, b, table[0]) << endl;

	// Тхомас Симпсон
	cout << "Simpson's integral: " << Integrals_simpson(H, table[0]) << endl;

	// Эйткин
	double k = 1 / 9.;
	table.push_back(make_table(1e5 / k, 0, 3 * M_PI, H * k));
	table.push_back(make_table(1e5 / (k * k), 0, 3 * M_PI, H * k * k));
	double B = (Integrals_rectangle(H, table[0]) - Integrals_rectangle(H * k, table[1])) /
		(Integrals_rectangle(H * k, table[1]) - Integrals_rectangle(H * k * k, table[2]));
	int p = round( - log(B) / log(k));
	cout << "p: " << p << endl;

	//метод Рунге
	double sigma = pow(k, p) / (pow(k, p) - 1);
	double Runge = sigma * Integrals_rectangle(H, table[0]) + (1 - sigma) * Integrals_rectangle(H * k, table[1]);
	cout << "Runge: " << Runge << endl;
	 
	// Ромберг
	const int q = table.size();
	Eigen::MatrixXd A1(q,q);
	Eigen::MatrixXd A2(q, q);
	Eigen::VectorXd vec_temp(q);
	for (int j = 0; j < q; j++)
	{
		for (int i = 0; i < q; i++)
		{
			if (j == 0) 
			{
				vec_temp[i] = Integrals_rectangle(H * pow(k, i), table[i]);
			}
			else
			{
				vec_temp[i] = pow(H * pow(k, i), p + j - 1);
			}
		}
		A1.col(j) = vec_temp;
	}
	vec_temp.setConstant(1);
	A2 = A1;
	A2.col(0) = vec_temp;
	double Romberg = A1.determinant() * (A2.inverse()).determinant();
	cout << "Romberg: " << Romberg << endl;
	cout << "abs error: " << abs(Romberg - true_value) << endl;

	//Legendre polynomials
	vector <double> roots;
	double n = 20;
	root(n, n, roots, n);
	double sum = 0;
	for (auto s: roots)
	{
		sum += weights(s, n) * ft(s);
	}
	cout << "Gaussian quadrature: " << (b - a) / 2. * sum << endl;

	//Optimal distribution
	cout << "Optimal distribution: " << opt_table(a, b).first << " n: " << opt_table(a, b).second << endl;

	//Монте_карло
	cout << "Monte Carlo: " << MC(table[2], a, b) << endl;
}