#include <iostream>
#include <iomanip>
#include <fstream>
#include "gnu.h"
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
	Gnuplot plot;

	
	const double a = 0, b = 3 * M_PI;
	int n = 1e5;
	const double H = abs(a - b) / n;
	vector <vector<pair<double, double>>> table { make_table(n, a, b, H) };

	//прямоугольник
	std::cout << std::setprecision(15);
	cout << "rectangle integral: " << Integrals_rectangle_with_gnu(H, table[0]) << endl;
	//plot("plot 'out.csv' using 1:2 w l lt 1 lw 1 title 'rectangle'"); // рисовать
	cout << "min n: " << min_n(a, b, "rectangle") << endl << endl;
	//system("pause"); cout << endl;

	// трапеция
	cout << "trapezoid integral: " << Integrals_trap_with_gnu(H, a, b, table[0]) << endl;
	//plot("plot 'out.csv' using 1:2 w l lt 1 lw 1 title 'trapezoi'"); // рисовать
	cout << "min n: " << min_n(a, b, "trapezoid") << endl << endl;
	//system("pause"); cout << endl;

	// Тхомас Симпсон
	cout << "Simpson's integral: " << Integrals_simpson(H, table[0]) << endl;
	//plot("plot 'out.csv' using 1:2 w l lt 1 lw 1 title 'simpson'"); // рисовать
	cout << "min n: " << min_n(a, b, "simpson") << endl << endl;
	//system("pause"); cout << endl;
	plot("exit gnuplot");

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
	cout << "Romberg: " << Romberg << endl << endl;

	//Legendre polynomials
	n = min_n(a, b, "Gaussian_quadrature");
	cout << "Gaussian quadrature: " << Gaussian_quadrature(a, b, n) << endl;
	cout << "min n: " << n << endl << endl;

	//Optimal distribution
	cout << "Optimal distribution: " << opt_table(a, b).first << endl; 
	cout << "min n: " << opt_table(a, b).second << endl << endl;

	//Монте_карло
	cout << "Monte Carlo: " << MC(table[2], a, b) << endl;
}