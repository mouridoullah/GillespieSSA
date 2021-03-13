#include <stdio.h>
#include <math.h>
#include <errno.h>
#include <time.h>
#include <stdlib.h>

#define NOMBRES_ESPECES 3
#define NOMBRES_DE_REACTIONS 9

/*--------------------------------------------------------*/
void ecrireDansData(FILE *f, double t, double *x, int N){
    int i;
    fprintf(f, "%f ", t);
    for (i = 0; i < N; i++) {
        fprintf(f, "%f ", x[i]);
    }
    fprintf(f, "\n");
}
/*--------------------------------------------------------*/
void remplirTab(double t[], double n){
	for (int i = 0; i < n; ++i){
		printf("Donner la case numero %d : ", i+1);
		scanf("%lf", &t[i]);
	}
}
/*--------------------------------------------------------*/
void afficherTab(double t[], double n){
	int i;
	for ( i = 0, printf("["); i < n; ++i){
		printf("%f", t[i]);
		if(i < n-1) printf(", ");
	}
	printf("]\n");
}
/*--------------------------------------------------------*/
void remplirTab2D(double **t, int l, int c){
	for (int i = 0; i < l; ++i){
		printf("Veulliez remplir la %d la ligne : \n", i+1);
		remplirTab(t[i], c);
	}
}
/*--------------------------------------------------------*/
void afficherTab2D(double** t, int l, int c){
	for (int i = 0; i < l; ++i){
		afficherTab(t[i], c);
		//printf("\n");
	}
}
/*--------------------------------------------------------*/
void allouerTableau(double** t, int n){
	*t = (double*)malloc(n*sizeof(double));
}
/*-------------------------------------------------------------*/
double genererNombreAlea(){
	return (double)rand()/RAND_MAX;
}
/*-------------------------------------------------------------*/
double* calculDeH(double* x, int n) {
	double* h = malloc(n*sizeof(double));
    // h[0] = 1.0;
    // h[1] = x[0];
    // h[2] = 0.5*x[0]*x[1]*(x[0] - 1);
    // h[3] = x[0];
    h[0] = 1.0;
    h[1] = 1.0;
    h[2] = 1.0;
    h[3] = x[0];
    h[4] = x[1];
    h[5] = x[2];
    h[6] = 0.5*x[0]*x[1]*(x[0] - 1);
    h[7] = 0.5*x[1]*x[2]*(x[1] - 1);
    h[8] = 0.5*x[2]*x[0]*(x[2] - 1);

    return h;
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

		// printf("%f ", t);
	 //    for(int i = 0; i < N; ++i){
	 //    	printf("%f ", x[i]);
	 //    }
		// printf("\n");

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

int main(int argc, char const *argv[]){
	double T = 50.20;
	double c[] = {5000,5000,5000,5,5,5,0.00001,0.00001,0.00001};
	double x[] = {1000, 1000, 1000};
	// double v[NOMBRES_DE_REACTIONS][NOMBRES_ESPECES] = 	{ {1, 0, -1},
	// 											 	      {-1, 1, 0},
	// 											 	      {1, -1, 0},
	// 											 	      {-1, 0, 1},
	// 											 	      {1, 0, -1},
	// 											 	      {0, -1, 1},
	// 											 	      {-1, 1, 0},
	// 											 	      {1, 0, -1},
	// 											 	      {0, -1, 1} };

	double v[NOMBRES_DE_REACTIONS][NOMBRES_ESPECES] = 	{ {0,1,-1},
											 	     	  {-1,0,1},
											 	      	  {1,-1,0},
											 	      	  {1,-1,0},
											 	      	  {0,1,-1},
											 	      	  {-1,0,1},
											 	      	  {-1,0,1},
											 	      	  {1,-1,0},
											 	      	  {0,1,-1} };


	
	Gillespie("data.txt", c, v, x, NOMBRES_ESPECES, NOMBRES_DE_REACTIONS, T);


	/*--------------------------------------------------*/
	// double x[] = {500,2000};
	// double c[] = {5000.0, 50.0, 0.00005, 5.0};
	// afficherTab(x,NOMBRES_ESPECES);

	// double* h = calculDeH(x, NOMBRES_DE_REACTIONS);
	// afficherTab(h,NOMBRES_DE_REACTIONS);


	// double* a = calculDePropensity(h,c,NOMBRES_DE_REACTIONS);
	// afficherTab(a,NOMBRES_DE_REACTIONS);

	// double a0 = sommeDesA(a,NOMBRES_DE_REACTIONS);
	// printf("%f\n", a0);


	// srand(time(NULL));
	// double r2 = genererNombreAlea();
	// int mu = calculDMu(a, r2, NOMBRES_DE_REACTIONS);
	// printf("%d\n", mu);


	// // //afficherTab2D((double**)v, M, N);
	// // //afficherTab((double*)v[0], N);

	// miseAJourDesX(x,v,mu,NOMBRES_ESPECES);
	// afficherTab(x,NOMBRES_ESPECES);

	// free(h);
	// free(a);
	return 0;
}
