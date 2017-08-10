#include "mcgas.h"
#include "random.h"
#include "render.h"
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
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
		int sN = sqrt(N);
		int ix = i/sN;
		int iy = i%sN;
		s->parts[i].x[0] = s->L/sN*(ix+0.5+(frand()*0.5-0.25));
		s->parts[i].x[1] = s->L/sN*(iy+0.5+(frand()*0.5-0.25));
	}
	s->histgrid = calloc(HIST_W*HIST_H*HIST_D,sizeof(int16_t));
	memset(s->histgrid,0,sizeof(int16_t)*HIST_W*HIST_H*HIST_D);
	system_hist_parts(s);
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

void hist_impossible(const char *msg) {
	fprintf(stderr, "%s: ERROR: could not create histogram. Increase GRID_D\n",msg);
	exit(1);
}

void system_set_L(System *s, double L) {
	for(int i = 0; i < s->N; i++) {
		s->parts[i].x[0] *= L/s->L;
		s->parts[i].x[1] *= L/s->L;
	}
	s->L = L;
}
void system_hist_parts(System *s) {
	s->hist_possible = true;
	for(int ix = 0; ix < HIST_W; ix++) {
		for(int iy = 0; iy < HIST_H; iy++) {
			for(int iz = 0; iz < HIST_D; iz++) {
				int idx = ix*HIST_H*HIST_D+iy*HIST_D+iz;
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
				s->parts[i].histidx=nidx;
				break;
			}
		}
		if(nidx == idx + HIST_D) {
			hist_impossible("hist_parts");
		}
	}
}

/*
double potential(double x[2], double y[2]) {
	double d[2];
	d[0] = x[0] - y[0];
	d[1] = x[1] - y[1];
	double r2 = (d[0]*d[0]+d[1]*d[1])/CR0/CR0+1e-9;
	r2 = 1/(r2*r2*r2);
	return CE*(r2*r2-2*r2);
}*/

// https://spiral.imperial.ac.uk/bitstream/10044/1/262/1/edm2006jcp.pdf (1)
// computes both the potential and the virial interaction part for the pressure
void potential_rf(double x[2], double y[2], double *e, double *f) {
	double d[2];
	d[0] = x[0] - y[0];
	d[1] = x[1] - y[1];
	double r2 = (d[0]*d[0]+d[1]*d[1])/CR0/CR0+1e-9;
	r2 = 1/(r2*r2);

	double a = CE*r2*r2;
	double b = CE*r2;
	*e = a-2*b;
	*f = 8*(a-b);
}

void gather_particle_interaction(System *s, double x[2], int i, bool smaller_than_i, double *E, double *p) {
	int rangex = R_CUTOFF/s->L*HIST_W+1;
	int rangey = R_CUTOFF/s->L*HIST_H+1;
	int idx = hist_index(s,x)/HIST_D;
	int ix0 = idx/HIST_H;
	int iy0 = idx%HIST_H;

	*E = 0;
	*p = 0;
	for(int ix = -rangex; ix <= rangex; ix++) {
		if(ix0+ix < 0 || ix0+ix >= HIST_W)
			continue;
		for(int iy = -rangey; iy <= rangey; iy++) {

			if(iy0+iy < 0 || iy0+iy >= HIST_H)
				continue;
			for(int iz = 0; iz < HIST_D; iz++) {
				int j = s->histgrid[HIST_D*(HIST_H*(ix0+ix)+iy0+iy)+iz];
				if(j == -1)
					break;
				if((smaller_than_i && j>i) ||j == i || j < 0)
					continue;
				double de, f;
				potential_rf(x,s->parts[j].x, &de, &f);
				*E+=de;
				*p+=f;
			}
		}
	}
}

void system_mc_update(System *s) {
	int i = random()%s->N;
	double x[2];
	x[0] = s->parts[i].x[0]+nfrand()*MC_SIGX;
	x[1] = s->parts[i].x[1]+nfrand()*MC_SIGX;

	if(x[0] < 0 || x[0] > s->L || x[1] < 0 || x[1] > s->L)
		return;
	double dE = 0;
	dE += (x[1]-s->parts[i].x[1])*MASS*GRAV;

	double e, f;
	double dp = 0;
	gather_particle_interaction(s,x,i,false,&e,&f);
	dE+=e;
	dp+=f;
	gather_particle_interaction(s,s->parts[i].x,i,false,&e,&f);
	dE-=e;
	dp-=f;

	if(dE < 0 || frand() < exp(-dE/s->T)) {
		s->parts[i].x[0] = x[0];
		s->parts[i].x[1] = x[1];

		int nhistidx = hist_index(s,x);
		if(nhistidx/HIST_D != s->parts[i].histidx/HIST_D) {
			s->histgrid[s->parts[i].histidx] = -2;
			int n = 0;
			for(n = 0; n < HIST_D && s->histgrid[nhistidx+n] >= 0; n++);
			if(n == HIST_D)
				hist_impossible("update");
			s->histgrid[nhistidx+n] = i;
			s->parts[i].histidx = nhistidx+n;
		}
		s->E += dE;
		s->p += dp/2./s->L/s->L;
	}
}

void system_mc_sweep(System *s) {
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

void system_energypressure(System *s, double *sE, double *pressure) {
	double E=0, rf=0;

	for(int i = 0; i < s->N; i++) {
		E += GRAV*MASS*s->parts[i].x[1];

		double dE, drf;
		gather_particle_interaction(s,s->parts[i].x,i,true,&dE,&drf);
		E+=dE;
		rf+=drf;
	}

	*pressure = s->T*s->N/s->L/s->L + 1/2.*rf/s->L/s->L;
	*sE = E+s->T*s->N;
}

void measure(int N, double V0, double V1, int V_steps, double T0, double T1, int T_steps) {
	int equilibration_time = 1000+0*1000*(V1/N < 1.*CR0*CR0);
	int measure_time = 3000;
	char filename[256];
	srandom(time(0));
	snprintf(filename,sizeof(filename),"data/N=%d_V0=%f_T0=%f.csv",N,V0/N/CR0/CR0,T0);
	FILE *out = fopen(filename,"a");
	if(out == 0) {
		perror("measure_volume");
		exit(1);
	}
	Renderer r;
	System s;
	system_init(&s, N, sqrt(V0),T0);
	for(int j = 0; j < equilibration_time; j++) // some more equilibration at the beginning
		system_mc_sweep(&s);

	for(int k = 0; k < V_steps; k++) {
		double L = sqrt(V0+(V1-V0)/V_steps*k);
		system_set_L(&s,L);
		for(int i = 0; i < T_steps; i++) {
			double T = T0+(T1-T0)/T_steps*i;
			s.T = T;
			for(int j = 0; j < equilibration_time; j++)
				system_mc_sweep(&s);
			for(int j = 0; j < measure_time; j++) {
				system_mc_sweep(&s);
				if(j % 500 == 0) {
					double E, p;
					system_energypressure(&s, &E, &p);
					s.E=E;
					s.p=p;
				}

				// printf("%d: %g %g %g %g\n",j,E,s.E,p,s.p);
				fprintf(out,"%.8e\t%.8e\t%.8e\t%.8e\n",L*L/N/CR0/CR0,T,s.E/N,s.p*CR0*CR0);
			}
			fflush(out);
			snprintf(filename,sizeof(filename),"pics/N=%d_V=%.2f_T=%.2f.webp",N,L*L/N/CR0/CR0,T);
			int rc = render_system(&r,&s,filename);
			if(rc != 0) {
				perror("render_system");
				exit(1);
			}
		}
	}
	fclose(out);
	system_free(&s);
}
/*
int main_constant_V(int argc, char **argv) {
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
}*/

int main(int argc, char **argv) {
	int V_bigsteps = 6;
	int V_steps = 20;
	if(argc != 2) {
		printf("Usage: %s Vstep[0-%d]\n",argv[0],V_bigsteps-1);
		return 1;
	}

	int i = atoi(argv[1]);
	if(i < 0 || i >= V_bigsteps) {
		printf("Argument out of range\n");
		return 1;
	}

	int N = 10000;
	double V0 = 5;
	double V1 = 0.5;
	double T = 1.5;

	double V0i = V0+(V1-V0)/V_bigsteps*i;
	double V1i = V0+(V1-V0)/V_bigsteps*(i+1);

	measure(N, V0i,V1i,V_steps, T, T, 1);

	return 0;
}
