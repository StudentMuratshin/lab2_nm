#include "function.h"
#include <vector>
#include <fstream>
#include <cmath>
#include <iostream>
using namespace std;
constexpr auto M_PI = 3.14159265358979323846;
constexpr auto true_value = 2.69632737829770440658;
constexpr auto delta = 2e-3;

double fx(const double x)
{
	return sin(x) * log(x + 1);
}

double ft(double t)
{
	return sin(3 * M_PI / 2. * (t + 1)) * log(3 * M_PI / 2. * (t + 1) + 1);
}

double Integrals_rectangle_with_gnu(double h, const vector<pair<double, double>>& T)
{
	ofstream f_rectangle("out.csv");
	double sum = 0;
	for (int i = 1; i < T.size(); i++)
	{
		sum += fx((T[i - 1].first + T[i].first) / 2.);
		f_rectangle << i << " " << abs(true_value - h * sum) << endl; //в файл
	}
	f_rectangle.close();
	return h * sum;
}

double Integrals_rectangle(double h, const vector<pair<double, double>>& T)
{
	double sum = 0;
	for (int i = 1; i < T.size(); i++)
	{
		sum += fx((T[i - 1].first + T[i].first) / 2.);
	}
	return h * sum;
}

double Integrals_trap(double h, double a, double b, const vector<pair<double, double>>& T)
{
	double sum = 0;
	for (int i = 0; i < T.size(); i++)
	{
		sum += fx(T[i].first);
	}
	return h * sum + h * (fx(a) + fx(b)) / 2.;
}

double Integrals_trap_with_gnu(double h, double a, double b, const vector<pair<double, double>>& T)
{
	ofstream f_trap("out.csv");
	double sum = 0;
	for (int i = 0; i < T.size(); i++)
	{
		sum += fx(T[i].first);
		f_trap << i << " " << abs(true_value - h * sum + h * (fx(a) + fx(b)) / 2.) << endl; //в файл
	}
	f_trap.close();
	return h * sum + h * (fx(a) + fx(b)) / 2.;
}

double Integrals_simpson(double h, vector<pair<double, double>>& T)
{
	double sum4 = 0, sum2 = 0;
	for (int i = 0; i < T.size(); i ++)
	{
		if (i % 2 == 0) sum4 += fx(T[i].first);
		else sum2 += fx(T[i].first);
	}
	return h * (fx(T[0].first) + 4 * sum4 + 2 * sum2 + fx(T[T.size() - 1].first)) / 3.;
}

double Integrals_simpson_with_gnu(double h, vector<pair<double, double>>& T)
{
	ofstream f_simpson("out.csv");
	double sum4 = 0, sum2 = 0;
	for (int i = 0; i < T.size(); i++)
	{
		if (i % 2 == 0) sum4 += fx(T[i].first);
		else sum2 += fx(T[i].first);
		f_simpson << i << " " << abs(h * (fx(T[0].first) + 4 * sum4 + 2 * sum2 + fx(T[T.size() - 1].first)) / 3. - true_value) << endl; //в файл
	}
	f_simpson.close();
	return h * (fx(T[0].first) + 4 * sum4 + 2 * sum2 + fx(T[T.size() - 1].first)) / 3.;
}

vector<pair<double, double>> make_table(int N, const double a, const double b, double H)
{
	vector<pair<double, double>> T;
	for (int i = 0; i <= N; i++)
	{
		T.push_back({ i * H, fx(i * H) });
	}
	return T;
}

double Legendre(double x, double i)
{
	if (i == 0) return 1.;
	else if (i == 1) return x;
	else
	{
		return (2 * i - 1) / i * x * Legendre(x, i - 1) - (i - 1) / i * Legendre(x, i - 2);
	}
}

double Derivative(double x, double i)
{
	if (i == 1) return 1.;
	else if (i == 2) return 3*x;
	else return i / (1 - x * x) * (Legendre(x, i - 1) - x * Legendre(x, i));
}

double root(double i, double k, vector <double>& roots, double n)
{
	if (k == 0)
	{
		return cos(M_PI * (4 * i - 1) / (4 * n + 2));
	}
	else
	{
		double X = root(i, k - 1, roots, n);
		double T = X - Legendre(X, n) / Derivative(X, n);
		if (i == k) {
			roots.push_back(T);
			if (i!=1) root(i - 1, i - 1, roots, n);
		}
		return T;
	}

}

double weights(double x, double n)
{
	return 2 / ((1-x*x) * pow(Derivative(x,n),2));
}

double I_tr(double xi, double xi_1, double h)
{
	return h / 2. * (fx(xi_1) + fx(xi));
}

double I(double xi, double xi_1, double h)
{
	return h / 4. * (fx(xi_1) + 2 * fx((xi - h / 2.)) + fx(xi));
}

double Eps(double xi, double xi_1, double h)
{
	return abs(1 / 3. * (I(xi, xi_1, h) - I_tr(xi, xi_1,  h)));
}

pair <double, int> opt_table(double a, double b)
{
	int k = 1;
	bool key_out = false, key_in = false;
	double S = 0, x = a, eps, h = 2, I_trapezoid;

	while(!key_out)
	{
		if (key_in) key_out = true;

		while (true)
		{
			eps = Eps(x, x + h, h);
			I_trapezoid = I(x, x + h, h);
			if (eps > delta *  (b - a) / h) h /= 2;
			else break;
		}

		x += h;
		S += I_trapezoid;

		if (x + h > b)
		{
			h = b - x;
			key_in = true;
		}
		k++;
	}
	pair <double, int> r = { S,k };
	return r;
}

double MC(vector<pair<double, double>>& T, double a, double b)
{
	double sum = 0;
	for (auto s: T)
	{
		sum += s.second;
	}
	return (b - a) * sum / T.size();
}

int min_n(double a, double b, string name)
{
	double n = 1;
	double h = abs(a - b) / n;
	vector <pair<double, double>> T = make_table(n, a, b, h);

	if (name == "rectangle")
	{
		double I_ = Integrals_rectangle(h, T);
		while (abs(true_value - I_) >= delta)
		{
			n++;
			h = abs(a - b) / n;
			T = make_table(n, a, b, h);
			I_ = Integrals_rectangle(h, T);
		}
	}

	else if (name == "trapezoid")
	{
		double I_ = Integrals_trap(h, a, b, T);
		while (abs(true_value - I_) >= delta)
		{
			n++;
			h = abs(a - b) / n;
			T = make_table(n, a, b, h);
			I_ = Integrals_trap(h, a, b, T);
		}
	}

	else if (name == "simpson")
	{
		double I_ = Integrals_simpson(h, T);
		while (abs(true_value - I_) >= delta)
		{
			n++;
			h = abs(a - b) / n;
			T = make_table(n, a, b, h);
			I_ = Integrals_simpson(h, T);
		}
	}

	else if (name == "Gaussian_quadrature")
	{
		double I_ = Gaussian_quadrature(a, b, n);
		while (abs(true_value - I_) >= delta)
		{
			n++;
			h = abs(a - b) / n;
			T = make_table(n, a, b, h);
			I_ = Gaussian_quadrature(a, b, n);
		}
	}
	return n;
}

double Gaussian_quadrature(double a, double b, int n)
{
	vector <double> roots;
	root(n, n, roots, n);
	double sum = 0;
	for (auto s : roots)
	{
		sum += weights(s, n) * ft(s);
	}
	return (b - a) / 2. * sum;
}
