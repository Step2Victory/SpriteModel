#pragma once
#include "Matrix.h"



Matrix MinResidMethod(Matrix A, Matrix b, int max_iter, double eps)
{
	int n_iter = 0;
	Matrix x_last(b.shape().first, 1);
	Matrix x(b.shape().first , 1);
	while (n_iter < max_iter && abs((A * x_last - b).T() * (A * x_last - b)) >= eps)
	{
		
		Matrix r = A * x_last - b;
		double tau = (A * r).T() * r / ((A * r).T() * (A * r));
		x = x_last - tau * r;
		n_iter += 1;
		x_last = x;
	} 
	return x;
}

