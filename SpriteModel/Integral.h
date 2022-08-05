#pragma once
#include <cmath>
#include <functional>

double integral1(std::function<double(double)>  func, double a, double b, int n=1000)
{
	double h = (b - a) / n;
	double ans = 0;
	for (int i = 0; i < n; ++i)
	{
		ans += func(a + i* h) + 4 * func(a + i * h + h / 2)+ func(a + (1 + i) * h);
	}
	return h / 6 * ans;
}

double integral2(std::function<double(double, double)>  func, double a, double b, double c, double d, int n=1000, int m=1000)
{
	double h1 = (b - a) / n;
	double h2 = (d - c) / m;
	double ans = 0;
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < m; ++j)
		{
			ans += func(a + i * h1 + h1 / 2, c + j * h2  + h2 / 2);
		}
	}
	return h1 * h2 * ans;
}
