#pragma once

#include <stdint.h>
#include <stdbool.h>
#include <stdbool.h>

#define GRAV -0. // gravitation
#define CE 0.03 // coupling energy
#define CR0 0.01 // coupling energy
#define R_CUTOFF 6*CR0

#define MASS 1.0

#define MC_SIGX CR0*2
#define MC_SIGP 1

enum {
	HIST_W = 30,
	HIST_H = 30,
	HIST_D = 40
};

typedef struct Particle{
	double x[2];
	double p[2];
} Particle;

typedef struct System {
	int N;
	double L;
	double T;

	double E;

	Particle *parts;
	int16_t *histgrid;
	bool hist_possible;
} System;
