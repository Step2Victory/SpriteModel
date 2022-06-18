#pragma once
#include <vector>
#include <iostream>
class Matrix
{
	std::vector<std::vector<double>> res;
	int n;
	int m;
public:
	Matrix()
	{
		n = 0;
		m = 0;
	}
	Matrix(int n, int m)
	{
		res = std::vector<std::vector<double>>(n,
			std::vector<double>(m));
		this->n = n;
		this->m = m;
	}

	Matrix(const Matrix &a)
	{
		res = a.res;
		n = a.n;
		m = a.m;
	}
	Matrix(const std::vector<double>& arr)
	{
		res = std::vector<std::vector<double>> (arr.size(), 
			std::vector<double>(1));
		for (int i = 0; i < arr.size(); ++i)
		{
			res[i][0] = arr[i];
		}

		n = res.size();
		m = res[0].size();
	}

	Matrix(const std::vector<std::vector<double>> &arr)
	{
		res = arr;
		n = arr.size();
		m = arr[0].size();
	}


	double& operator()(int i, int j)
	{
		return res[i][j];
	}
	
	Matrix operator()(int i)
	{
		return Matrix(res[i]).T();
	}
	std::pair<int, int> shape()
	{
		return std::pair<int, int>(n, m);
	}

	Matrix T() const
	{
		std::vector<std::vector<double>> arr = 
			std::vector<std::vector<double>>(m,
				std::vector<double>(n));
		for (int i = 0; i < m; ++i)
		{
			for (int j = 0; j < n; ++j)
			{
				arr[i][j] = res[j][i];
			}
		}
		return Matrix(arr);
	}

	Matrix operator-() const
	{

		std::vector<std::vector<double>> arr =
			std::vector<std::vector<double>>(n,
				std::vector<double>(m));
		for (int i = 0; i < n; ++i)
		{
			for (int j = 0; j < m; ++j)
			{
				arr[i][j] = -res[i][j];
			}
		}
		return Matrix(arr);
	}

	Matrix operator+() const
	{
		return *this;
	}

	Matrix operator+(const Matrix & A) const
	{
		if (n != A.n || m != A.m)
		{
			throw "Wrong dimensions!!!";
		}
		std::vector<std::vector<double>> arr 
			= std::vector<std::vector<double>>(n, 
				std::vector<double>(m));
		for (int i = 0; i < n; ++i)
		{
			for (int j = 0; j < m; ++j)
			{
				arr[i][j] = res[i][j] + A.res[i][j];
			}
		}
		return Matrix(arr);
	}

	Matrix operator-(const Matrix& B) const
	{
		return *this + (-B);
	}
	Matrix operator*(const Matrix& A) const
	{
		if (A.m == 1 && A.n == 1)
			return *this * A.res[0][0];
		if (m == 1 && n == 1)
			return  A * res[0][0];
		if (m != A.n)
			throw "Wrong dimensions!!!";
		std::vector<std::vector<double>> arr =
			std::vector<std::vector<double>>(n, 
				std::vector<double>(A.m));
		for (int i = 0; i < n; ++i)
		{
			for (int j = 0; j < m; ++j)
			{
				for (int k = 0; k < A.m; ++k)
				{
					arr[i][k] += res[i][j] * A.res[j][k];
				}
				
			}
		}
		return Matrix(arr);
	}

	bool eq(const Matrix& A, double eps = 10e-5) const
	{
		if (n != A.n || m != A.m)
		{
			throw "Wrong dimensions!!!";
		}
		for (int i = 0; i < n; ++i)
		{
			for (int j = 0; j < m; ++j)
			{
				if (abs(A.res[i][j] - res[i][j]) > eps)
				{
					return false;
				}
			}
		}
		return true;
	}

	Matrix operator*(double a) const
	{
		std::vector<std::vector<double>> arr =
			std::vector<std::vector<double>>(n,
				std::vector<double>(m));
		for (int i = 0; i < n; ++i)
		{
			for (int j = 0; j < m; ++j)
			{
				arr[i][j] = a * res[i][j];
			}
		}
		return Matrix(arr);
	}

	operator double()
	{
		if (m != 1 || n != 1)
			throw "Wrong dimensions!!!";
		return res[0][0];
	}

	void print()
	{
		for (int i = 0; i < n; ++i)
		{
			for (int j = 0; j < m; ++j)
			{
				std::cout << res[i][j] << ' ';
			}
			std::cout << std::endl;
		}
		std::cout << std::endl;
	}
	
};

Matrix operator*(double a, const Matrix &A)
{
	return A * a;
}

