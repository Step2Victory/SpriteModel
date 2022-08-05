#pragma once
#include "SolverLE.h"
#include "Integral.h"
#include <functional>
#include <fstream>

class TestTask2
{
	double R;
	double z1;
	double z2;
	double T;
	double h_r;
	double h_z;
	double dt;
	//const double mu_0 = 0.00000125663706212;
	const double mu_0 = 1;

	int n;
	int m;
	int l;


	std::function< double(double r, double z)> initial_conditions;
	std::function< double(double z, double t)> E_z_R;
	std::function< double(double r, double t)> E_r_h1;
	std::function< double(double r, double t)> E_r_h2;
	std::function< double(double r, double z, double t)> sigma;
	std::function< double(double r, double z, double t) > j_ext_z;
	std::function< double(double r, double z, double t) > j_ext_r;


	std::vector<Matrix> c;
	Matrix A;
	Matrix B;
	Matrix b;

public:
	TestTask2(double R_, double z1_, double z2_, double T_, double n_, double m_, double l_, std::function< double(double r, double z)> initial_conditions_,
		std::function< double(double z, double t)> E_z_R_, std::function< double(double r, double t)> E_r_h1_, std::function< double(double r, double t)> E_r_h2_,
		std::function< double(double r, double z, double t)> sigma_,
		std::function< double(double r, double z, double t)> j_ext_z_, std::function< double(double r, double z, double t) > j_ext_r_)
	{
		R = R_;
		z1 = z1_;
		z2 = z2_;
		T = T_;

		n = n_;
		m = m_;
		l = l_;

		initial_conditions = initial_conditions_;
		E_z_R = E_z_R_;
		sigma = sigma_;
		j_ext_z = j_ext_z_;
		j_ext_r = j_ext_r_;

		h_r = R / m;
		h_z = (z2 - z1) / n;
		dt = T / l;

		c = std::vector<Matrix>(l, Matrix((n + 1) * m, 1));

		for (int i = 0; i < n + 1; ++i)
		{
			for (int j = 0; j < m; ++j)
			{
				double r = (j + 1) * h_r;
				double z = i * h_z;
				c[0](index(i, j), 0) = initial_conditions(r, z);
			}
		}
	}

	int index(int i, int j)//i - узел по z, j - узел по r
	{
		return i * m + j;
	}

	std::vector<int> invert_index_func(int i)//для подстановки в функции
	{
		std::vector<int> a(2);
		a = { i / m, i % m + 1};
		return a;
	}
	std::vector<int> invert_index_matr(int i)//для индексирования матриц
	{
		std::vector<int> a(2);
		a = { i / m, i % m };
		return a;
	}
	void compute_A()
	{
		A = Matrix((n + 1) * m, (n + 1) * m);
		

		for (int i = 0; i < (n + 1) * m; ++i)
		{

			for (int j = i; j < (n + 1) * m; ++j)
			{

				std::vector<int> ind_grid_1 = invert_index_func(i);
				std::vector<int> ind_grid_2 = invert_index_func(j);


				std::function<double(double r, double z)> f = [&](double r, double z) {return v(ind_grid_1[0], ind_grid_1[1], r, z) * v(ind_grid_2[0], ind_grid_2[1], r, z) * r; };


				if (abs(ind_grid_1[0] - ind_grid_2[0]) < 2 && abs(ind_grid_1[1] - ind_grid_2[1]) < 2)
				{
					A(i, j) = mu_0 * integral2(f, (ind_grid_1[1] - 1) * h_r, (ind_grid_1[1] + 1) * h_r, z1 + (ind_grid_1[0] - 1) * h_z, z1 + (ind_grid_1[0] + 1) * h_z, 200, 200) / (h_r  * h_z) / (h_r * h_z);
				}
				else
				{
					A(i, j) = 0;
				}
				A(j, i) = A(i, j);
			}
		}
	}

	void compute_B(double t)
	{
		B = Matrix((n + 1) * m, (n + 1) * m);
		for (int i = 0; i < (n + 1) * m; ++i)
		{

			for (int j = i; j < (n + 1) * m; ++j)
			{

				std::vector<int> ind_grid_1 = invert_index_func(i);
				std::vector<int> ind_grid_2 = invert_index_func(j);


				std::function<double(double r, double z)> f = [&](double r, double z) {
				
					return dvdz(ind_grid_1[0], ind_grid_1[1], r, z) * dvdz(ind_grid_2[0], ind_grid_2[1], r, z) * r / sigma(r, z, t) 
																							+ drvdr(ind_grid_1[0], ind_grid_1[1], r, z) * drvdr(ind_grid_2[0], ind_grid_2[1], r, z) / (sigma(r,z,t) * r); };


				if (abs(ind_grid_1[0] - ind_grid_2[0]) < 2 && abs(ind_grid_1[1] - ind_grid_2[1]) < 2)
				{
					B(i, j) =  integral2(f, (ind_grid_1[1] - 1) * h_r, (ind_grid_1[1] + 1) * h_r, z1 + (ind_grid_1[0] - 1) * h_z, z1 + (ind_grid_1[0] + 1) * h_z, 1000, 1000) / (h_r * h_z) / (h_r  * h_z);
				}
				else
				{
					B(i, j) = 0;
				}
				B(j, i) = B(i, j);
			}
		}
	}
	

	void compute_b(double t)
	{
		b = Matrix((n + 1) * m, 1);
		for (int i = 0; i < (n + 1) * m; ++i)
		{

			std::vector<int> ind_grid = invert_index_func(i);
			std::function<double(double r, double z)> f2 = [&](double r, double z) {return drvdr(ind_grid[0], ind_grid[1], r, z) * j_ext_z(r, z, t) / sigma(r, z, t); };
			if (ind_grid[1] == m)
			{
				std::function<double(double z)> f1 = [&](double z) {return E_z_R(z, t) * v(ind_grid[0], ind_grid[1], R, z) * R; };
				

				b(i, 0) = integral1(f1, z1 + (ind_grid[0] - 1) * h_z, z1 + (ind_grid[0] + 1) * h_z, 1000) / (h_r * h_z) +
					integral2(f2, (ind_grid[1] - 1) * h_r, (ind_grid[1] + 1) * h_r, z1 + (ind_grid[0] - 1) * h_z, z1 + (ind_grid[0] + 1) * h_z, 1000, 1000) / (h_r  * h_z);
			}
			else
			{
				b(i, 0) = integral2(f2, (ind_grid[1] - 1) * h_r, (ind_grid[1] + 1) * h_r, z1 + (ind_grid[0] - 1) * h_z, z1 + (ind_grid[0] + 1) * h_z, 1000, 1000) / (h_r  * h_z);
			}
		}
	}

	void compute(double alpha)
	{
		for (int i = 1; i < l; ++i)
		{
			double t = i * dt;
			compute_A();
			compute_B(t);
			compute_b(t);
			c[i] = MinResidMethod(1.0 / dt * A + alpha * B, b + (1.0 / dt * A - (1 - alpha) * B) * c[i - 1], 10000, 10e-10);
		}
	}

	void compute_const()
	{
		compute_A();
		compute_B(0);
		compute_b(0);

		Matrix x = MinResidMethod(B, b, 10000, 10e-10);
		for (int i = 0; i < (n + 1) * m; ++i)
		{
			if (i % m == 0)
				std::cout << '\n' << 0 << ' ';
			std::cout << x(i, 0) << ' ';
			
		}
		
	}

	double v(int i, int j, double r, double z)
	{
		if (r < 0 || r > R || z < z1 || z > z2 || r > (j + 1) * h_r || r < (j - 1) * h_r || z > z1 + (i + 1) * h_z || z < z1 + (i - 1) * h_z)
		{
			return 0;
		}
		if (j * h_r <= r && r <= (j + 1) * h_r && z1 + i * h_z <= z && z <= z1 + (i + 1) * h_z)
		{
			return ((j + 1) * h_r - r) * (z1 + (i + 1) * h_z - z);
		}
		if ((j - 1) * h_r <= r && r <= j * h_r && z1 + i * h_z <= z && z <= z1 + (i + 1) * h_z)
		{
			return (r - (j - 1) * h_r) * (z1 + (i + 1) * h_z - z);
		}
		if ((j - 1) * h_r <= r && r <= j * h_r && z1 + (i - 1) * h_z <= z && z <= z1 + i * h_z)
		{
			return (r - (j - 1) * h_r) * (z - z1 - (i - 1) * h_z);
		}
		if (j * h_r <= r && r <= (j + 1) * h_r && z1 + (i - 1) * h_z <= z && z <= z1 + i * h_z)
		{
			return ((j + 1) * h_r - r) * (z - z1 - (i - 1) * h_z);
		}



	}
	double dvdr(int i, int j, double r, double z)

	{
		if (r < 0 || r > R || z < z1 || z > z2 || r > (j + 1) * h_r || r < (j - 1) * h_r || z > z1 + (i + 1) * h_z || z < z1 + (i - 1) * h_z)
		{
			return 0;
		}
		if (j * h_r <= r && r <= (j + 1) * h_r && z1 + i * h_z <= z && z <= z1 + (i + 1) * h_z)
		{
			return  -(z1 + (i + 1) * h_z - z);
		}
		if ((j - 1) * h_r <= r && r <= j * h_r && z1 + i * h_z <= z && z <= z1 + (i + 1) * h_z)
		{
			return (z1 + (i + 1) * h_z - z);
		}
		if ((j - 1) * h_r <= r && r <= j * h_r && z1 + (i - 1) * h_z <= z && z <= z1 + i * h_z)
		{
			return  (z - z1 - (i - 1) * h_z);
		}
		if (j * h_r <= r && r <= (j + 1) * h_r && z1 + (i - 1) * h_z <= z && z <= z1 + i * h_z)
		{
			return   -(z - z1 - (i - 1) * h_z);
		}
	}
	double dvdz(int i, int j, double r, double z)
	{
		if (r < 0 || r > R || z < z1 || z > z2 || r > (j + 1) * h_r || r < (j - 1) * h_r || z > z1 + (i + 1) * h_z || z < z1 + (i - 1) * h_z)
		{
			return 0;
		}
		if (j * h_r <= r && r <= (j + 1) * h_r && z1 + i * h_z <= z && z <= z1 + (i + 1) * h_z)
		{
			return -((j + 1) * h_r - r);
		}
		if ((j - 1) * h_r <= r && r <= j * h_r && z1 + i * h_z <= z && z <= z1 + (i + 1) * h_z)
		{
			return -(r - (j - 1) * h_r);
		}
		if ((j - 1) * h_r <= r && r <= j * h_r && z1 + (i - 1) * h_z <= z && z <= z1 + i * h_z)
		{
			return (r - (j - 1) * h_r);
		}
		if (j * h_r <= r && r <= (j + 1) * h_r && z1 + (i - 1) * h_z <= z && z <= z1 + i * h_z)
		{
			return ((j + 1) * h_r - r);
		}
	}
	double drvdr(int i, int j, double r, double z)
	{
		if (r < 0 || r > R || z < z1 || z > z2 || r > (j + 1) * h_r || r < (j - 1) * h_r || z > z1 + (i + 1) * h_z || z < z1 + (i - 1) * h_z)
		{
			return 0;
		}
		if (j * h_r <= r && r <= (j + 1) * h_r && z1 + i * h_z <= z && z <= z1 + (i + 1) * h_z)
		{
			return ((j + 1) * h_r - 2 * r) * (z1 + (i + 1) * h_z - z);
		}
		if ((j - 1) * h_r <= r && r <= j * h_r && z1 + i * h_z <= z && z <= z1 + (i + 1) * h_z)
		{
			return (2 * r - (j - 1) * h_r) * (z1 + (i + 1) * h_z - z);
		}
		if ((j - 1) * h_r <= r && r <= j * h_r && z1 + (i - 1) * h_z <= z && z <= z1 + i * h_z)
		{
			return (2 * r - (j - 1) * h_r) * (z - z1 - (i - 1) * h_z);
		}
		if (j * h_r <= r && r <= (j + 1) * h_r && z1 + (i - 1) * h_z <= z && z <= z1 + i * h_z)
		{
			return ((j + 1) * h_r - 2 * r) * (z - z1 - (i - 1) * h_z);
		}
	}

	void printanswer()
	{

		for (int t = 0; t < l; ++t)
		{
			for (int i = 0; i < (n + 1) * m; ++i)
			{
				if (i % m == 0)
					std::cout << '\n' << 0 << ' ';
				std::cout << c[t](i, 0) << ' ';

			}
		}
		std::cout << std::endl;
	}

	void write_answer_in_file(std::string a)
	{
		std::ofstream myfile;
		myfile.open(a);
		myfile << n << ' ' << m << ' ' << l << std::endl;

		for (int k = 0; k < l; ++k)
		{
			for (int i = 0; i < n + 1; ++i)
			{
				myfile << 0 << " ";
				for (int j = 0; j < m; ++j)
				{
					myfile << c[k](i * m + j, 0) << ' ';
				}
				myfile << std::endl;
			}
			myfile << std::endl;
		}
		myfile << std::endl;
		myfile.close();
	}
};