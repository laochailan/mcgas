#include "mcgas.h"
#include "random.h"
#include "render.h"
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <errno.h>

void system_init(System *s, int N, double L, double T) {
	s->T = T;
	s->E = 0;
	s->N = N;
	s->L = L;
	s->hist_possible = false;

	s->parts = calloc(N,sizeof(Particle));
	for(int i = 0; i < N; i++) {
		int sN = sqrt(N)+1;
		int ix = i/sN;
		int iy = i%sN;
		s->parts[i].x[0] = s->L/sN*ix;
		s->parts[i].x[1] = s->L/sN*iy;
		s->parts[i].p[0] = nfrand()*sqrt(MASS*s->T);
		s->parts[i].p[1] = nfrand()*sqrt(MASS*s->T);
	}
	s->histgrid = calloc(HIST_W*HIST_H*HIST_D,sizeof(int16_t));
	memset(s->histgrid,0,sizeof(int16_t)*HIST_W*HIST_H*HIST_D);
}

void system_free(System *s) {
	free(s->parts);
	free(s->histgrid);
}

int hist_index(System *s, double x[2]) {
	int ix = x[0]/s->L*HIST_W;
	int iy = x[1]/s->L*HIST_H;
	return HIST_D*(HIST_H*ix+iy);
}

void system_hist_parts(System *s) {
	s->hist_possible = true;
	for(int ix = 0; ix < HIST_W; ix++) {
		for(int iy = 0; iy < HIST_H; iy++) {
			for(int iz = 0; iz < HIST_D; iz++) {
				int idx = ix*HIST_H*HIST_D+iy*HIST_D+iz;
				if(s->histgrid[idx] == -1)
					break;
				s->histgrid[idx] = -1;
			}
		}
	}

	for(int i = 0; i < s->N; i++) {
		int idx = hist_index(s,s->parts[i].x);
		int nidx;
		for(nidx = idx; nidx < idx+HIST_D; nidx++) {
			if(s->histgrid[nidx] == -1) {
				s->histgrid[nidx] = i;
				break;
			}
		}
		if(nidx == idx + HIST_D) {
			s->hist_possible = false;
			break;
		}
	}
}

double potential(double x[2], double y[2]) {
	double d[2];
	d[0] = x[0] - y[0];
	d[1] = x[1] - y[1];
	double r2 = (d[0]*d[0]+d[1]*d[1])/CR0/CR0+1e-9;
	r2 = 1/(r2*r2*r2);
	return CE*(r2*r2-2*r2);
}

void system_mc_update(System *s) {
	int i = random()%s->N;
	double x[2],p[2],dp[2];
	x[0] = s->parts[i].x[0]+nfrand()*MC_SIGX;
	x[1] = s->parts[i].x[1]+nfrand()*MC_SIGX;

	if(x[0] < 0 || x[0] > s->L || x[1] < 0 || x[1] > s->L)
		return;
	p[0] = s->parts[i].p[0];
	p[1] = s->parts[i].p[1];
	dp[0] = nfrand()*MC_SIGP;
	dp[1] = nfrand()*MC_SIGP;

	double dE = (p[0]*dp[0]+p[1]*dp[1])/MASS + (dp[0]*dp[0]+dp[1]*dp[1])/MASS/2;
	dE += (x[1]-s->parts[i].x[1])*MASS*GRAV;

	if(s->hist_possible) {
		int rangex = R_CUTOFF/s->L*HIST_W+1;
		int rangey = R_CUTOFF/s->L*HIST_H+1;
		int idx = hist_index(s,x)/HIST_D;
		int ix0 = idx/HIST_H;
		int iy0 = idx%HIST_H;

		for(int ix = -rangex; ix <= rangex; ix++) {
			if(ix0+ix < 0 || ix0+ix >= HIST_W)
				continue;
			for(int iy = -rangey; iy <= rangey; iy++) {

				if(iy0+iy < 0 || iy0+iy >= HIST_H)
					continue;
				for(int iz = 0; iz < HIST_D; iz++) {
					int j = s->histgrid[HIST_D*(HIST_H*(ix0+ix)+iy0+iy)+iz];
					if(j == i)
						continue;
					if(j == -1)
						break;
					dE += potential(x,s->parts[j].x)-potential(s->parts[i].x,s->parts[j].x);
				}
			}
		}
	} else {
		for(int j = 0; j < s->N; j++) {
			if(j == i)
				continue;
			dE += potential(x,s->parts[j].x)-potential(s->parts[i].x,s->parts[j].x);
		}
	}

	if(dE < 0 || frand() < exp(-dE/s->T)) {
		s->parts[i].x[0] = x[0];
		s->parts[i].x[1] = x[1];
		s->parts[i].p[0] = p[0]+dp[0];
		s->parts[i].p[1] = p[1]+dp[1];
		s->E += dE;
	}
}

void system_mc_sweep(System *s) {
	system_hist_parts(s);
	if(!s->hist_possible)
		printf("WARNING: histogram not possible.\n");

	for(int i = 0; i < s->N; i++)
		system_mc_update(s);
}

void system_meanpos(System *s, double m[2]) {
	m[0]=0;
	m[1]=0;
	for(int i = 0; i < s->N; i++) {
		m[0]+=s->parts[i].x[0];
		m[1]+=s->parts[i].x[1];
	}
	m[0] /= s->N;
	m[1] /= s->N;
}

double system_energy(System *s) {
	double E=0, p[2];

	for(int i = 0; i < s->N; i++) {
		p[0] = s->parts[i].p[0];
		p[1] = s->parts[i].p[1];
		E += (p[0]*p[0]+p[1]*p[1])/2/MASS;
		E += GRAV*MASS*s->parts[i].x[1];

		if(s->hist_possible) {
			int rangex = R_CUTOFF/s->L*HIST_W+1;
			int rangey = R_CUTOFF/s->L*HIST_H+1;
			int idx = hist_index(s,s->parts[i].x)/HIST_D;
			int ix0 = idx/HIST_H;
			int iy0 = idx%HIST_H;

			for(int ix = -rangex; ix <= rangex; ix++) {
				if(ix0+ix < 0 || ix0+ix >= HIST_W)
					continue;
				for(int iy = -rangey; iy <= rangey; iy++) {

					if(iy0+iy < 0 || iy0+iy >= HIST_H)
						continue;
					for(int iz = 0; iz < HIST_D; iz++) {
						int j = s->histgrid[HIST_D*(HIST_H*(ix0+ix)+iy0+iy)+iz];
						if(j == i)
							continue;
						if(j == -1)
							break;
						E += potential(s->parts[i].x,s->parts[j].x);
					}
				}
			}
		} else {
			for(int j = 0; j < s->N; j++) {
				if(j == i)
					continue;
				E += potential(s->parts[i].x,s->parts[j].x);
			}
		}
		for(int j = 0; j < i; j++)
			E += potential(s->parts[i].x,s->parts[j].x);

	}
	return E;
}

int correlation_time(void) {
	System s;
	system_init(&s, 10000, 4,10);
	s.T = 30;
	double E0 = system_energy(&s)/s.N;
	for(int step = 0; step < 2000; step++) {
		system_mc_sweep(&s);
		if(system_energy(&s) > 0.95*3*E0) {
			system_free(&s);
			return step;
		}
	}
	printf("WARNING: Correlation time too long.\n");
	system_free(&s);
	return -1;
}

void measure_volume(int N, double L, double T0, double T1, int steps) {
	int equilibration_time = 1000;
	int measure_time = 5000;
	char filename[256];
	snprintf(filename,sizeof(filename),"data/N=%d_V=%f_T0=%f.csv",N,L*L,T0);
	FILE *out = fopen(filename,"a");
	if(out == 0) {
		perror("measure_volume");
		exit(1);
	}
	Renderer r;
	System s;
	system_init(&s, N, L,T0);
	for(int j = 0; j < equilibration_time; j++) // some more equilibration at the beginning
		system_mc_sweep(&s);
	for(int i = 0; i < steps; i++) {
		double T = T0+(T1-T0)/steps*i;
		s.T = T;
		for(int j = 0; j < equilibration_time; j++)
			system_mc_sweep(&s);
		for(int j = 0; j < measure_time; j++) {
			system_mc_sweep(&s);
			s.E = system_energy(&s);
			fprintf(out,"%.10e\t%.10e\t%.10e\n",L*L,T,s.E/N);
		}
		fflush(out);
		snprintf(filename,sizeof(filename),"pics/N=%d_V=%.2f_T=%.2f.webp",N,L*L,T);
		int rc = render_system(&r,&s,filename);
		if(rc != 0) {
			perror("render_system");
			exit(1);
		}
	}
	fclose(out);
	system_free(&s);
}

int main(int argc, char **argv) {
	int T_bigsteps = 13;
	int T_steps = 5;
	int V_steps = 129;
	if(argc != 3) {
		printf("Usage: %s Vstep[0-%d] Tstep[0-%d]\n",argv[0],V_steps-1,T_bigsteps-1);
		return 1;
	}

	int i = atoi(argv[1]);
	int j = atoi(argv[2]);
	if(i < 0 || j < 0 || i >= V_steps || j >= T_bigsteps) {
		printf("Usage: %s Vstep[0-%d] Tstep[0-%d]\n",argv[0],V_steps-1,T_bigsteps-1);
		return 1;
	}

	int N = 10000;
	double V0 = 8;
	double V1 = 0.05;
	double T0 = 30;
	double T1 = 0.01;

	double T0i = T0+(T1-T0)/T_bigsteps*j;
	double T1i = T0+(T1-T0)/T_bigsteps*(j+1);

	double V = V0 + (V1-V0)/(V_steps-1)*i;
	measure_volume(N, sqrt(V),T0i,T1i,T_steps);
	printf("V=%g complete\n",V);

	return 0;
}
