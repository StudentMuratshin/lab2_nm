#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <Eigen/Dense>
#include "function.h"

constexpr auto M_PI = 3.14159265358979323846;
using namespace std;

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
	//Eigen::MatrixXf A1(p + 1e5 - 2);
	int q = 3;
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
	cout << A1.determinant() * (A2.inverse()).determinant();
}