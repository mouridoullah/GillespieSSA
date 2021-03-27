#include <stdio.h>
#include <math.h>
#include <errno.h>
#include <time.h>
#include <stdlib.h>

#define NOMBRES_ESPECES 2
#define NOMBRES_DE_REACTIONS 4

/*--------------------------------------------------------*/
void ecrireDansData(FILE *f, double t, double *x, int N){
    int i;
    fprintf(f, "%f ", t);
    for (i = 0; i < N; i++) {
        fprintf(f, "%f ", x[i]);
    }
    fprintf(f, "\n");
}
/*-------------------------------------------------------------*/
double* calculDeH(double* x, int n) {
	double* h = malloc(n*sizeof(double));
    h[0] = 1.0;
    h[1] = x[0];
    h[2] = 0.5*x[0]*x[1]*(x[0] - 1);
    h[3] = x[0];

    return h;
}
/*-------------------------------------------------------------*/
double genererNombreAlea(){
	return (double)rand()/RAND_MAX;
}
/*-------------------------------------------------------------*/
double* calculDePropensity(double* h, double* c, int M) {
	double* a = malloc(M*sizeof(double));
	for (int i = 0; i < M; i++) {
		a[i] = h[i] * c[i];
	}
	return a;
}
/*-------------------------------------------------------------*/
double sommeDesA(double* a, int M) {
	double a0 = 0;
	for (int i = 0; i < M; i++) {
		a0 += a[i];
	}
	return a0;
}
/*-------------------------------------------------------------*/
int calculDMu(double* a, double r2, int M) {
	double r2a0 = r2 * sommeDesA(a, M);
	double sum = 0.0;
	
	for (int i = 0; i < M; i++) {
		sum += a[i];
		if (r2a0 <= sum) {
			return i;
		}
	}
	return 0;
}
/*-------------------------------------------------------------*/	
void miseAJourDesX(double* x, double v[][NOMBRES_ESPECES], int mu, int N) {
	for (int i = 0; i < N; i++) {
		x[i] += v[mu][i];
	}
}
/*-------------------------------------------------------------*/
void Gillespie(char *myfile, double* c, double v[][NOMBRES_ESPECES], double* x, int N, int M, double T) {
	FILE *data;
        data = fopen(myfile, "w");
	int mu = 0;
	double r1 = 0, r2 = 0, a0 = 0, t = 0;
	double* h;
	double* a;

	srand(time(NULL));

	while(t < T) {
		h = calculDeH(x, M);
		a = calculDePropensity(h, c, M);
		a0 = sommeDesA(a, M);

		r1 = genererNombreAlea();
		r2 = genererNombreAlea();

		if(a0 == 0 ) break;

		t += log(1.0 / r1) / a0;
		mu = calculDMu(a, r2, M);

		miseAJourDesX(x, v, mu, N);
		ecrireDansData(data, t, x, N);
	}
	free(a);
	free(h);
}
/*-------------------------------------------------------------*/
void RunSimulation(){
	double v[NOMBRES_DE_REACTIONS][NOMBRES_ESPECES];
	double x[NOMBRES_ESPECES]; double c[NOMBRES_DE_REACTIONS]; double TEMPS_MAX;

	FILE* f = NULL; 
	f = fopen("input.txt", "r");

	if (f != NULL){
        fscanf(f, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &TEMPS_MAX, 
        	     &x[0], &x[1], &c[0], &c[1], &c[2], &c[3], &v[0][0], &v[0][1], 
        	     &v[1][0], &v[1][1], &v[2][0], &v[2][1], &v[3][0], &v[3][1]);
        fclose(f);
    }

    Gillespie("data.txt", c, v, x, NOMBRES_ESPECES, NOMBRES_DE_REACTIONS, TEMPS_MAX);
}
/*-------------------------------------------------------------*/
int main(int argc, char const *argv[]){
	RunSimulation();
    return (0);
}
