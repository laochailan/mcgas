#include "random.h"
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <float.h>

double frand(void) {
	return random()/(double)RAND_MAX;
}

// polar method from https://en.wikipedia.org/wiki/Marsaglia_polar_method
double nfrand(void) {
	static bool hasSpare = false;
	static double spare;
#pragma omp threadprivate(spare)
#pragma omp threadprivate(hasSpare)

	if(hasSpare) {
		hasSpare = false;
		return spare;
	}

	hasSpare = true;
	static double u, v, s;
	do
	{
		u = (random() / ((double) RAND_MAX)) * 2.0 - 1.0;
		v = (random() / ((double) RAND_MAX)) * 2.0 - 1.0;
		s = u * u + v * v;
	}
	while( (s >= 1.0) || (s == 0.0) );

	s = sqrt(-2.0 * log(s) / s);
	spare = v * s;
	return u * s;
}
