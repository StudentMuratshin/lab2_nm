#include <iostream>
#include <iomanip>
#include <fstream>
#include "gnu.h"
#include <vector>
#include <string>
#include <numeric>
#include <cmath>
#include "function.h"
#include <omp.h>


constexpr auto M_PI = 3.14159265358979323846;
constexpr auto true_value = 2.69632737829770440658;

int main()
{
	Gnuplot plot;
	double start = omp_get_wtime();
	
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
	//plot("exit gnuplot");

	// Эйткин
	double k = 1 / 9.;
	table.push_back(make_table(n / k, 0, 3 * M_PI, H * k));
	table.push_back(make_table(n / (k * k), 0, 3 * M_PI, H * k * k));
	int p = Aitken(k, table ,H);
	cout << "p: " << p << endl;

	//метод Рунге
	cout << "Runge: " << Runge(k, p, H, table) << endl;
	 
	// Ромберг
	cout << "Romberg: " << Romberg(table, H, k, p) << endl << endl;

	//Legendre polynomials
	n = min_n(a, b, "Gaussian_quadrature");
	cout << "Gaussian quadrature: " << Gaussian_quadrature(a, b, n) << endl;
	cout << "min n: " << n << endl << endl;

	//Optimal distribution
	pair <double, int> res = opt_table(a, b);
	cout << "Optimal distribution: " << res.first << endl; 
	cout << "min n: " << res.second << endl << endl;

	//Монте_карло
	n = 1000;
	MC_for_gnu(a, b, n);
	plot("plot 'out.csv' using 1:2 w l lt 1 lw 1 title 'Monte Carlo'"); // рисовать
	cout << "Monte Carlo: " << MC(a, b, n) << endl;
	system("pause"); cout << endl;

	double end = omp_get_wtime();
	cout << "Time: " << (end - start) << endl;
}