#pragma once
#include "SolverLE.h"
#include "Integral.h"
#include <functional>
#include <fstream>

class TestTask1
{
	double R;
	double T;
	double h;
	double dt;
	const double mu_0 = 0.00000125663706212;
	//const double mu_0 = 1;

	int n;
	int l;
	

	std::function< double(double r)> initial_conditions;
	std::function< double(double t)> E_z_R;
	std::function< double(double r, double t)> sigma;
	std::function< double(double r, double t)> j_ext;


	std::vector<Matrix> c;
	Matrix A;
	Matrix B;
	Matrix b;

public:
	TestTask1(double R_, double T_, double n_, double l_, std::function< double(double r)> initial_conditions_,
		std::function< double(double t)> E_z_R_,
		std::function< double(double r, double t)> sigma_,
		std::function< double(double r, double t)> j_ext_)
	{
		R = R_;
		T = T_;

		n = n_;
		l = l_;

		initial_conditions = initial_conditions_;
		E_z_R = E_z_R_;
		sigma = sigma_;
		j_ext = j_ext_;

		h = R / n;
		dt = T / l;
	
		c = std::vector<Matrix>(l, Matrix(n, 1));

		for (int i = 0; i < n; ++i)
		{
			
			double r = (i + 1) * h;
			c[0](i, 0) = initial_conditions(r);
		}
	}

	void compute_A()
	{
		A = Matrix(n, n);
		A(n - 1, n - 1) = integral1([&](double r) { return  (r - (n - 1) * h) * (r - (n - 1) * h) * r; }, (n - 1) * h, R) * mu_0 / (h * h);

		for (int i = 0; i < n - 1; ++i)
		{
			A(i, i) = (integral1([&](double r) { return  (r - i* h) * (r - i * h) * r; }, i * h, (i + 1) * h) 
				+ integral1([&](double r) { return  ((i + 2) * h - r) * ((i + 2) * h - r) * r; }, (i + 1) * h, (i + 2) * h))* mu_0 / (h * h);
			A(i, i + 1) = integral1([&](double r) { return  ((i + 2) * h - r) * (r - (i + 1) * h) * r; }, (i + 1) * h, (i + 2) * h)  * mu_0 / (h * h);
			A(i + 1, i) = A(i, i + 1);
		}
	}

	void compute_B(double t)
	{
		B = Matrix(n, n);
		B(0, 0) = integral1([&](double r) { return  4 * r  / (h * h * sigma(r, t)); }, 0, h) + integral1([&](double r) { return  (2 * h - 2 * r) * (2 * h - 2 * r) / (h * h * sigma(r, t) * r); }, h, 2 * h);
		B(n - 1, n - 1) = integral1([&](double r) { return  (2 * r - (n - 1) * h) * (2 * r - (n - 1) * h) / (h * h * sigma(r, t) * r); }, (n - 1) * h, R);

		for (int i = 0; i < n - 1; ++i)
		{
			
			if (i != 0)
				B(i, i) = (integral1([&](double r) { return  (2 * r - i * h) * (2 * r - i * h) / (sigma(r, t) * r); }, i * h, (i + 1) * h) + integral1([&](double r) { return  (2 * r - (i + 2) * h) * (2 * r - (i + 2) * h) / (sigma(r, t) * r); }, (i + 1) * h, (i + 2) * h)) / (h * h);
			B(i, i + 1) = integral1([&](double r) { return  (2 * r - (i + 1) * h) * ((i + 2) * h - 2 * r) / (sigma(r, t) * r); }, (i + 1) * h, (i + 2) * h) / (h * h);
			B(i + 1, i) = B(i, i + 1);
		}
	}

	void compute_b(double t)
	{
		b = Matrix(n, 1);
		for (int i = 0; i < n - 1; ++i)
		{
			b(i, 0) = (integral1([&](double r) { return  (2 * r - i * h) * j_ext(r, t)/ sigma(r, t); }, i * h, (i + 1) * h) + 
						integral1([&](double r) { return  ((i + 2) * h - 2 * r) * j_ext(r, t)/ sigma(r, t); }, (i + 1) * h, (i + 2) * h) ) / h;
		}
		b(n - 1, 0) = E_z_R(t) * R + integral1([&](double r) { return  (2 * r - (n - 1) * h) * j_ext(r, t) / sigma(r, t); }, (n - 1) * h, R) / h;
	}

	void compute(double alpha)
	{
		for (int i = 1; i < l; ++i)
		{
			double t = i * dt;
			compute_A();
			compute_B(t);
			compute_b(t);
			c[i] = MinResidMethod(1.0 / dt * A + alpha* B, b + (1.0 / dt * A - (1 - alpha) * B) * c[i - 1], 10000, 10e-5);
		}
	}

	void compute_const()
	{

		compute_B(0);
		compute_b(0);

		Matrix x = MinResidMethod(B, b, 10000, 10e-5);
		x.print();
		c[0].print();
	}

	void printanswer()
	{

		for (int i = 0; i < l; ++i)
		{
			//if (i % 2 == 0)
			{
				std::cout << 0 << " ";
				for (int j = 0; j < n; ++j)
				{
					std::cout << c[i](j, 0) << ' ';
				}
				std::cout << std::endl;
			}
		}
		std::cout << std::endl;
	}

	void write_answer_in_file(std::string a)
	{
		std::ofstream myfile;
		myfile.open(a);
		
		
		for (int i = 0; i < l; ++i)
		{
			//if (i % 2 == 0)
			{
				myfile << 0 << " ";
				for (int j = 0; j < n; ++j)
				{
					myfile << c[i](j, 0) << ' ';
				}
				myfile << std::endl;
			}
		}
		myfile << std::endl;
		myfile.close();
	}
};