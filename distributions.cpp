#include "stdafx.h"
#include <cstdlib>
#include <cmath>

#include "distributions.h"

namespace model
{
	double UniformDistrib(double m, double d)
	{
		double random;
		random = ((double)std::rand() / RAND_MAX)*(2*d) + m-d;
		return random;
	}

	double NormalDistrib(double m, double d)
	{
		double a = 0;
		for (int i = 0; i < 12; i++)
		{
			a += ((double)std::rand() / RAND_MAX);
		}
		a -= 6;
		return m + a*d;
	}

	double LogNormalDistrib(double m, double d)
	{
		return exp(NormalDistrib(m, d));
	}
}