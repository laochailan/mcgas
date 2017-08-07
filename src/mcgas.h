#pragma once

#include <stdint.h>
#include <stdbool.h>
#include <stdbool.h>

#define GRAV -0. // gravitation
#define CE 0.03 // coupling energy
#define CR0 0.01 // coupling energy
#define R_CUTOFF 2*CR0

#define MASS 1.0

#define MC_SIGX CR0*2
#define MC_SIGP 1

enum {
	HIST_W = 30,
	HIST_H = 30,
	HIST_D = 100
};

typedef struct Particle{
	double x[2];
	double p[2];

	int histidx;
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

void system_init(System *s, int N, double L, double T);
void system_free(System *s);
int hist_index(System *s, double x[2]);
void hist_impossible(const char *);
void system_set_L(System *s, double L);
void system_hist_parts(System *s);
