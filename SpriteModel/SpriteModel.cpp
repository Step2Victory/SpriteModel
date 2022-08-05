#include <iostream>
#include <functional>
#include "FiniteElements.h"
#include "TestTask1.h"
#include "TestTask2.h"

void testtask1()
{
	double E_z = 5;
	double sgm = 3;
	double j_z_ext = 1;
	std::function<double(double r, double t)> sigma = [&](double r, double t) {return sgm; };
	std::function<double(double r, double t)> j_ext = [&](double r, double t) {return  j_z_ext; };
	//std::function<double(double r)> init_conditions = [&](double r) {return r / 2 * (E_z * sgm + j_z_ext); };
	std::function<double(double r)> init_conditions = [&](double r) {return 8 * r  ; };

	std::function<double(double t)> E = [&](double t) {return E_z; };
	int n = 10;
	int l = 500;
	double T = 5;
	double R = 1;
	TestTask1 task(R, T, n, l, init_conditions, E, sigma, j_ext);
	//task.compute(1.0/2);
	//task.printanswer();
	//task.compute_const();
	//task.compute(1.0 / 2);
	//task.printanswer();
	task.compute_const();
	task.write_answer_in_file("1d.txt");
}

void testtask2()
{
	double E_z = 5;
	double sgm = 3;
	double j_z_ext = 1;
	std::function<double(double r, double z, double t)> sigma = [&](double r, double z, double t) {return sgm; };
	std::function<double(double r, double z, double t)> j_ext = [&](double r, double z, double t) {return  j_z_ext; };
	std::function<double(double r, double)> init_conditions = [&](double r, double z) {return 8 * r * r / 3 * (E_z * sgm + j_z_ext); };

	std::function<double(double z, double t)> E_z_R = [&](double z, double t) {return E_z; };
	std::function<double(double z, double t)> E_r_h1 = [&](double r, double t) {return 0; };
	std::function<double(double z, double t)> E_r_h2 = [&](double r, double t) {return 0; };

	std::function<double(double r, double z, double t)> j_ext_z = [&](double r, double z, double t) {return j_z_ext; };
	std::function<double(double r, double z, double t)> j_ext_r = [&](double r, double z, double t) {return 0; };

	int n = 10;
	int m = 10;
	int l = 40;
	double T = 4;
	double R = 1;
	double h1 = 0;
	double h2 = 1;
	TestTask2 task(R, h1, h2,  T, n, m, l, init_conditions, E_z_R, E_r_h1, E_r_h2, sigma, j_ext_z, j_ext_r);
	//task.compute(1.0/2);
	//task.printanswer();
	task.compute_const();
	task.compute(1.0 / 2);
	//task.printanswer();
	//task.compute_const();
	task.write_answer_in_file("2d.txt");
}

int main()
{
	/*
	double E_z = 5;
	double sgm = 3;
	double j_z_ext = 1;
	std::function<double(double r, double z, double t)> sigma = [&](double r, double z, double t) {return sgm; };
	std::function<std::vector<double>(double r, double z, double t)> j_ext = [&](double r, double z, double t) {return std::vector<double>{ 0, j_z_ext }; };

	std::function<double(double r, double t)> right_boundary = [&](double r, double t) {return E_z; };
	std::function<double(double r, double t)> up_boundary = [&](double z, double t) {return 0; };
	std::function<double(double r, double t)> down_boundary = [&](double z, double t) {return 0; };

	std::function<double(double r, double z)> init_conditions = [&](double r, double z) {return r / 2 * (E_z * sgm + j_z_ext); };

	double a = 1;
	double c = 0;
	double d = 1;
	double T = 1;

	int n = 10;
	int m = 10;
	int l = 2;
	
	FiniteElements task(a, c, d, T, n, m, l, j_ext, sigma, up_boundary, down_boundary, right_boundary, init_conditions);
	task.compute();
	task.printanswer(0);
	task.printanswer(0.8);
	*/
	testtask2();
	
}

