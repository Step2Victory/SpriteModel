#include "pch.h"
#include "../SpriteModel/Integral.h"
#include "../SpriteModel/SolverLE.h"
#include "../SpriteModel/FiniteElements.h"
#include "../SpriteModel/Matrix.h"

TEST(TestCaseName, TestName) {
  EXPECT_EQ(1, 1);
  EXPECT_TRUE(true);
}

TEST(Integral1, Test1) {
	double res = integral1([](double x) {return 1; }, 0, 1, 100);
	EXPECT_NEAR(res, 1, 10e-6);
}

TEST(Integral1, Test2) {
	double res = integral1([](double x) {return x; }, 0, 1, 100);
	EXPECT_NEAR(res, 0.5, 10e-6);
}

TEST(Integral1, Test3) {
	double res = integral1([](double x) {return sin(x); }, 0, 1, 100);
	EXPECT_NEAR(res, -cos(1) + 1, 10e-6);
}


TEST(Integral1, Test4) {
	double res = integral1([](double x) {return -(5 * x - 3) * (5 * x - 1); }, 0.2, 0.6, 100);
	EXPECT_NEAR(res, 4.0/15, 10e-6);
}

TEST(Integral2, Test1) {
	double res = integral2([](double x, double y) {return 1; }, 0, 1, 0, 1, 100, 100);
	EXPECT_NEAR(res, 1, 10e-6);
}

TEST(Integral2, Test2) {
	double res = integral2([](double x, double y) {return x * y; }, 0, 1, 1, 3, 100, 100);
	EXPECT_NEAR(res, 2, 10e-6);
}

TEST(Integral2, Test3) {
	double res = integral2([](double x, double y) {return y; }, 0, 1, 1, 3, 50, 100);
	EXPECT_NEAR(res, 4, 10e-6);
}
TEST(Integral2, Test4) {
	double res = integral2([](double x, double y) {return cos(x * y); }, 0, 1, 0, 2, 50, 100);
	EXPECT_NEAR(res, 1.60541297680269, 10e-3);
}


TEST(VectorArithm, Test1) {
	std::vector<std::vector<double>> A = { {1, 2} , {3, 4} };
	std::vector<double> b = { 5 , 6 };
	std::vector<double> cor = { 17, 39 };
	Matrix A_m(A);
	Matrix b_m(b);
	Matrix cor_m(cor);
	Matrix res = A_m * b_m;
	EXPECT_TRUE(cor_m.eq(res));
	
}

TEST(VectorArithm, Test2) {
	std::vector<std::vector<double>> A = { {1, 2, 1} , {3, 4, 2}, {1, 5, 2 } };
	std::vector<double> b = { 5 , 6 , 1};
	std::vector<double> cor = { 18, 41, 37 };
	Matrix A_m(A);
	Matrix b_m(b);
	Matrix cor_m(cor);
	Matrix res = A_m * b_m;
	EXPECT_TRUE(cor_m.eq(res));
}


TEST(VectorArithm, Test3) {
	std::vector<double> a = {1, 2, 1};
	std::vector<double> b = { 5 , 6 , 1 };
	Matrix a_m(a);
	Matrix b_m(b);
	Matrix res = a_m.T() * b_m;
	EXPECT_EQ(res(0,0), 18);
	
}

TEST(VectorArithm, Test4) {
	double a = 2;
	std::vector<double> b = { 5 , 6 , 1 };
	Matrix res = a * b;
	std::vector<double> cor = { 10, 12, 2 };
	Matrix b_m(b);
	Matrix cor_m(cor);
	EXPECT_TRUE(cor_m.eq(res));
}




TEST(MinResidMethod, Test1) {
	std::vector<std::vector<double>> A = { {1, 0} , {0, 1} };
	std::vector<double> b = { 1 , 1 };
	Matrix A_m(A);
	Matrix b_m(b);
	Matrix res = MinResidMethod(A, b, 100, 10e-3);
	Matrix a = res - b;
	EXPECT_NEAR( a.T() * a, 0, 10e-3);
}

TEST(MinResidMethod, Test2) {
	std::vector<std::vector<double>> A = { {2, 1} , {1, 2} };
	std::vector<double> b = { 3 , 4 };
	Matrix A_m(A);
	Matrix b_m(b);
	Matrix res = MinResidMethod(A, b, 100, 10e-3);
	Matrix a = res - Matrix(std::vector<double>{(double)2/3, (double)5/3});
	EXPECT_NEAR(a.T() * a, 0, 10e-3);
}

TEST(MinResidMethod, Test3) {
	std::vector<std::vector<double>> A = { {1, 2, 2} , {2, 1, 1}, {2, 1, 3} };
	std::vector<double> b = { 2 , 6, 8 };
	Matrix A_m(A);
	Matrix b_m(b);
	Matrix res = MinResidMethod(A, b, 1000, 10e-3);
	Matrix a = res - Matrix(std::vector<double>{(double)10 / 3, (double)-5 / 3, 1});
	EXPECT_NEAR(a.T() * a, 0, 10e-3);
}

/*
TEST(FiniteElements, Test1) {
	
	std::function<double(double r, double z, double t)> sigma = [](double r, double z, double t) {return 2; };
	std::function<std::vector<double>(double r, double z, double t)> j_ext = [](double r, double z, double t) {return std::vector<double>{ 0, 1 }; };

	std::function<double(double r, double t)> right_boundary = [](double r, double t) {return 1.0 / 2; };
	std::function<double(double r, double t)> up_boundary = [](double z, double t) {return 0; };
	std::function<double(double r, double t)> down_boundary = [](double z, double t) {return 0; };

	std::function<double(double r, double z)> init_conditions = [](double r, double z) {return r * r; };

	double a = 3;
	double c = 3;
	double d = 6;
	double T = 1;

	int n = 3;
	int m = 3;
	int l = 5;

	FiniteElements task(a, c, d, T, n, m, l, j_ext, sigma, up_boundary, down_boundary, right_boundary, init_conditions);
	Matrix A((n - 1) * (m - 1), (n - 1) * (m - 1));
	std::function<double(int i_1, int j_1, int i_2, int j_2, double r, double z)> f = [](int i_1, int j_1, int i_2, int j_2, double r, double z) {return 1; };
	task.compute_matrix(A, f);
	Matrix cor(std::vector<std::vector<double>>{ {4, 2, 2, 1}, { 2, 4, 1, 2 }, { 2, 1, 4,2 }, {1, 2, 2, 4}});
	EXPECT_TRUE(A.eq(cor));
 
}


TEST(FiniteElements, Test2) {

	std::function<double(double r, double z, double t)> sigma = [](double r, double z, double t) {return 2; };
	std::function<std::vector<double>(double r, double z, double t)> j_ext = [](double r, double z, double t) {return std::vector<double>{ 0, 1 }; };

	std::function<double(double r, double t)> right_boundary = [](double r, double t) {return 1.0 / 2; };
	std::function<double(double r, double t)> up_boundary = [](double z, double t) {return 0; };
	std::function<double(double r, double t)> down_boundary = [](double z, double t) {return 0; };

	std::function<double(double r, double z)> init_conditions = [](double r, double z) {return r * r; };

	double a = 3;
	double c = 3;
	double d = 6;
	double T = 1;

	int n = 3;
	int m = 3;
	int l = 5;

	FiniteElements task(a, c, d, T, n, m, l, j_ext, sigma, up_boundary, down_boundary, right_boundary, init_conditions);
	Matrix A((n - 1) * (m - 1), (n - 1) * (m - 1));
	std::function<double(int i_1, int j_1, int i_2, int j_2, double r, double z)> f = [](int i_1, int j_1, int i_2, int j_2, double r, double z) {return r; };
	task.compute_matrix(A, f);
	A.print();
	
	Matrix cor(std::vector<std::vector<double>>{ {4, 2, 2, 1}, { 2, 4, 1, 2 }, { 2, 1, 4,2 }, { 1, 2, 2, 4 }});
	EXPECT_TRUE(A.eq(cor));
	
}*/


