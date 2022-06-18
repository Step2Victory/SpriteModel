#pragma once
#include <iostream>
#include "Integral.h"
#include "SolverLE.h"
#include <vector>



class FiniteElements
{
	double a;
	double c;
	double d;
	double T;
	double h_r;
	double h_z;
	double h_t;
	int n;
	int m;
	int l;
	int N;
	std::vector<Matrix> x;
	std::function<std::vector<double>(double r, double z, double t)> j_ext;
	std::function<double(double r, double z, double t)> sigma;
	std::function<double(double r, double t)> upper_boundary;
	std::function<double(double r, double t)> lower_boundary;
	std::function<double(double z, double t)> right_boundary;
public:

	FiniteElements(double _a, double _c, double _d, double _T, int _n, int _m, int _l, 
	std::function<std::vector<double>(double r, double z, double t)> _j_ext,
	std::function<double(double r, double z, double t)> _sigma,
	std::function<double(double r, double t)> _upper_boundary, std::function<double(double r, double t)> _lower_boundary,
	std::function<double(double z, double t)> _right_boundary, std::function<double(double r, double z)> _init_conditions)
	{
		sigma = _sigma;
		j_ext = _j_ext;
		upper_boundary = _upper_boundary;
		lower_boundary = _lower_boundary;
		right_boundary = _right_boundary;
		a = _a;
		c = _c;
		d = _d;
		T = _T;
		l = _l;
		n = _n;
		m = _m;
		N = (n + 1) * (m + 1);
		x = std::vector<Matrix>(l, Matrix(N, 1));
		h_r = a / n;
		h_z = (d - c) / m;
		h_t = T / l;
		for (int s = 0; s < N; ++s)
		{
			int j_1 = s % (n + 1);
			int i_1 = s / (n + 1);
			x[0](s,0) = _init_conditions(j_1 * h_r, c + i_1 * h_z);
		}
	}



	void compute_matrix(Matrix& A, std::function<double(int i_1, int j_1, int i_2, int j_2, double r, double z)> func)
	{
		for (int s = 0; s < N; ++s)
		{
			for (int k = s; k < N; ++k)
			{
				int j_1 = s % (n + 1);
				int i_1 = s / (n + 1);
				int j_2 = k % (n + 1);
				int i_2 = k / (n + 1);
				std::function<double(double r, double z)> f = [&i_1, &j_1, &i_2, &j_2, &func](double r, double z)
				{
					return func(i_1, j_1, i_2, j_2, r, z);
				};
				if (i_1 == i_2 && j_1 == j_2)
				{
					A(s,k) = integral2(f, (j_1 - 1) * h_r, (j_1 + 1) * h_r, c + (i_1 - 1) * h_z, c + (i_1 + 1) * h_z, 2000, 2000);
				}
				else if(abs(i_1 - i_2) == 1 && j_1 == j_2)
				{
					A(s, k) = integral2(f, (j_1 - 1) * h_r, (j_1 + 1) * h_r, c + (std::max(i_1, i_2) - 1) * h_z, c + std::max(i_1, i_2) * h_z, 2000, 2000);
				}
				else if (i_1 == i_2 && abs(j_1 - j_2) == 1)
				{
					A(s, k) = integral2(f, (std::max(j_1, j_2) - 1) * h_r, (std::max(j_1, j_2)) * h_r, c + (i_1 - 1) * h_z, c + (i_1 + 1) * h_z, 2000, 2000);
				}
				else if (abs(i_1 - i_2) == 1 && abs(j_1 - j_2) == 1)
				{
					A(s, k) = integral2(f, (std::max(j_1, j_2) - 1) * h_r, (std::max(j_1, j_2)) * h_r, c + (std::max(i_1, i_2) - 1) * h_z, c + std::max(i_1, i_2) * h_z, 2000, 2000);
				}
				else
				{
					A(s, k) = 0;
				}
				A(k, s) = A(s, k);
			}
		}
	}

	void compute_r(Matrix &r, std::function<double(int i, int j, double r)> updown, std::function<double(int i, int j, double z)> leftright, std::function<double(int i, int j, double r, double z)> func)
	{

		for (int s = 0; s < N; ++s)
		{
			int j_1 = s % (n + 1);
			int i_1 = s / (n + 1);
			std::function<double(double r, double z)> f = [&i_1, &j_1, &func](double r, double z)
			{
				return func(i_1, j_1, r, z);
			};
			r(s, 0) = integral1([&i_1, &j_1, &updown](double r) {return updown(i_1, j_1, r); }, (j_1 - 1) * h_r, (j_1 + 1) * h_r, 2000) +
				integral1([&i_1, &j_1, &leftright](double z) {return leftright(i_1, j_1, z); }, c + (i_1 - 1) * h_z, c + (i_1 + 1) * h_z, 2000);// +
				//integral2(f, (j_1 - 1) * h_r, (j_1 + 1) * h_r, c + (i_1 - 1) * h_z, c + (i_1 + 1) * h_z, 100, 100);
		}
	}

	void compute()
	{
		Matrix A(N, N);
		std::function<double(int i_1, int j_1, int i_2, int j_2, double r, double z)> f = [&](int i_1, int j_1, int i_2, int j_2, double r, double z)
		{
			return r * v(i_1, j_1, r, z) * v(i_2, j_2, r, z);
		};
		compute_matrix(A, f);
		for (int i = 1; i < l; ++i)
		{
			Matrix B(N, N);
			Matrix r(N, 1);
			double t = i * h_t;


			std::function<double(int i_1, int j_1, int i_2, int j_2, double r, double z)> f_B = [&](int i_1, int j_1, int i_2, int j_2, double r, double z)
			{
				
				//return dvdz(i_1, j_1, r, z) * dvdz(i_2, j_2, r, z) * r / sigma(r, z, t) + drvdr(i_1, j_1, r, z) * drvdr(i_2, j_2, r, z) / (sigma(r, z, t) * r);
				return drvdr(i_1, j_1, r, z) * dvdr(i_2, j_2, r, z) / (sigma(r, z, t) * r);
			};
			compute_matrix(B, f_B);

			std::function<double(int i_1, int j_1, double r)> updown = [&](int i_1, int j_1, double r)
			{
				return (-upper_boundary(r, t) * v(i_1, j_1, r, d) * r + lower_boundary(r, t) * v(i_1, j_1, r, c)) * r ;
			};

			std::function<double(int i_1, int j_1, double z)> rightleft = [&](int i_1, int j_1, double z)
			{
				//return right_boundary(z, t) * v(i_1, j_1, a, z) * a;
				return right_boundary(z, t) * v(i_1, j_1, a, z);
			};

			std::function<double(int i_1, int j_1, double r, double z)> f_r = [&](int i_1, int j_1, double r, double z)
			{
				//return 1 / sigma(r, z, t) * j_ext(r, z, t)[1] * drvdr(i_1, j_1, r, z) - 1 / sigma(r, z, t) * j_ext(r, z, t)[0] * dvdz(i_1, j_1, r, z) * r;
				return 1.0 / sigma(r, z, t) * j_ext(r, z, t)[1] * dvdr(i_1, j_1, r, z);
			};

			compute_r(r, updown, rightleft, f_r);
			B.print();
			/*
			Matrix M = (1.0 / h_t) * A + (1.0 / 2) * B;
			
			Matrix b = r + (1.0 / h_t * A - (1.0 / 2) * B )* x[i - 1];
			*/
			Matrix M = B;

			Matrix b = r;
			x[i] = MinResidMethod(M, b, 10000, 10e-10);
		}
	}

	double v(int i, int j, double r, double z)
	{
		if (r < 0 || r > a || z < c || z > d  || r > (j + 1) * h_r || r < (j - 1) * h_r || z > c + (i + 1) * h_z || z < c + (i - 1) * h_z)
		{
			return 0;
		}

		return (r - (j + 1) * h_r) * (r - (j - 1) * h_r) / (h_r * h_r) * (z - c - (i + 1) * h_z) * (z - c - (i - 1) * h_z) / (h_z * h_z);
	}
	double dvdr(int i, int j, double r, double z)

	{
		if (r < 0 || r > a || z < c || z > d || r > (j + 1) * h_r || r < (j - 1) * h_r || z > c + (i + 1) * h_z || z < c + (i - 1) * h_z)
		{
			return 0;
		}
		return 2 * (r - j * h_r)/ (h_r * h_r) * (z - c - (i + 1) * h_z) * (z - c - (i - 1) * h_z) / (h_z * h_z);
	}
	double dvdz(int i, int j, double r, double z)
	{
		if (r < 0 || r > a || z < c || z > d || r > (j + 1) * h_r || r < (j - 1) * h_r || z > c + (i + 1) * h_z || z < c + (i - 1) * h_z)
		{
			return 0;
		}
		return (r - (j + 1) * h_r) * (r - (j - 1) * h_r) / (h_r * h_r) * 2 * (z - c - i * h_z) / (h_z * h_z);;
	}
	double drvdr(int i, int j, double r, double z)
	{
		if (r < 0 || r > a || z < c || z > d || r > (j + 1) * h_r || r < (j - 1) * h_r || z > c + (i + 1) * h_z || z < c + (i - 1) * h_z)
		{
			return 0;
		}
		return (3 * r *r - 4 * r * j * h_r + (j - 1) * (j + 1) * h_r * h_r) / (h_r * h_r) * (z - c - (i + 1) * h_z) * (z - c - (i - 1) * h_z) / (h_z * h_z);
	}

	void printanswer(double t)
	{
		int tau = int(t / h_t);
	
		for (int i = 0; i < (n + 1); ++i)
		{
			for (int j = 0; j < (m + 1); ++j)
			{
				std::cout << x[tau](i * (n + 1) + j) << ' ';
			}
			std::cout << std::endl;
		}
		std::cout << std::endl;
	}
};



