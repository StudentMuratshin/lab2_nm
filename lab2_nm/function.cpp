#include "function.h"
#include <vector>
using namespace std;

double fx(const double x)
{
	return sin(x) * log(x + 1);
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

double Integrals_simpson(double h, vector<pair<double, double>>& T)
{
	double sum4 = 0, sum2 = 0;
	for (int i = 0; i < T.size(); i += 2)
	{
		sum4 += fx(T[i].first);
	}


	for (int i = 1; i < T.size(); i += 2)
	{
		sum2 += fx(T[i].first);
	}


	return h * (fx(T[0].first) + 4 * sum4 + 2 * sum2 + fx(T[T.size() - 1].first)) / 3.;
}

vector<pair<double, double>> make_table(int N, const double a, const double b, double H)
{
	vector<pair<double, double>> T(N + 1);
	for (int i = 0; i <= N; i++)
	{
		T.push_back({ i * H, fx(i * H) });
	}
	return T;
}