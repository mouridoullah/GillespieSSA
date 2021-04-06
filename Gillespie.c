#include <stdio.h>
#include <math.h>
#include <errno.h>
#include <time.h>
#include <stdlib.h>

#define NOMBRE_ESPECES 3
#define NOMBRE_DE_REACTIONS 9																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																															

/*----------------------------------------------------------------------------------*/
/*
Ecrit les données issue de la simulation dans un fichier pour la visualisation
*/
void ecrireDansData(FILE *f, double t, double *x, int N){
    int i;
    fprintf(f, "%f ", t);
    for (i = 0; i < N; i++) {
        fprintf(f, "%f ", x[i]);
    }
    fprintf(f, "\n");
}
/*-------------------------------------------------------------*/
/*
Nombre de combinaisons de réactifs moléculaires distinctes présente au temps t pour
réaction R.
*/
double* calculDeH(double* x, int M) {
	
	if (NOMBRE_ESPECES == 1){
		/*

		*  -> S1 => hu = 1
		S1 -> reactions products => hu = x1

		*/
		double* h = malloc(M*sizeof(double));
	    h[0] = 1.0;
	    h[1] = x[0];	

	    return h;	
	}else if(NOMBRE_ESPECES == 2){
		/*
	Formations
	*          -> S1                 => hu = 1
	*          -> 2S2                => hu = 1

	Destructions
	S1         -> reactions products => hu = X1
	S2         -> reactions products => hu = X2

	Bimoleculaire
	2S1        -> reactions products => hu = 0.5 * X1(X1 - 1)
	2S2        -> reactions products => hu = 0.5 * X2(X2 - 1)
	S1 + S2    -> reactions products => hu = X1 * X2

	Avec trois molecules
	S1 + 2S2   -> reactions products => hu = 0.5 * X1X2(X2 - 1)
	3S1        -> reactions products => hu = (1/6) * X1(X1 - 1)(X1 - 2)


	*/
		double* h = malloc(M*sizeof(double));
	    h[0] = 1.0;
	    h[1] = x[0];
	    h[2] = 0.5*x[0]*x[1]*(x[0] - 1);
	    h[3] = x[0];

	    return h;
	 } else if (NOMBRE_ESPECES == 3){
		double* h = malloc(M*sizeof(double));
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
}
/*-------------------------------------------------------------*/
/*
Génère un nombre aléatoire uniformément entre 0 et 1.
*/
double genererNombreAlea(){
	return (double)rand()/RAND_MAX;
}
/*-------------------------------------------------------------*/
/*
Calcule de propensions.
*/
double* calculDePropensity(double* h, double* c, int M) {
	double* a = malloc(M*sizeof(double));
	for (int i = 0; i < M; i++) {
		a[i] = h[i] * c[i];
	}
	return a;
}
/*-------------------------------------------------------------*/
/*
Somme les propensions. 
*/
double sommeDesA(double* a, int M) {
	double a0 = 0;
	for (int i = 0; i < M; i++) {
		a0 += a[i];
	}
	return a0;
}
/*-------------------------------------------------------------*/
/*
Calcule de l’indice de la prochaine réaction avec une transformation inverse.
*/
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
/*
Mets à jour le nombre de chaque espèces présent selon la matrice de stichométrique.
*/
void miseAJourDesX(double* x, int N, int v[][N], int mu) {
	for (int i = 0; i < N; i++) {
		x[i] += v[mu][i];
	}
}
/*-------------------------------------------------------------*/
/*
L’algorithme de Gillespie.
*/
void Gillespie(char *myfile, double* c, int N, int v[][N], double* x, int M, double T) {
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

		miseAJourDesX(x, N, v, mu);
		ecrireDansData(data, t, x, N);
	}
	free(a);
	free(h);
}
/*-------------------------------------------------------------*/
int main(){
	double TEMPS_MAX;
	double x[NOMBRE_ESPECES];
	double c[NOMBRE_DE_REACTIONS];
	int v[NOMBRE_DE_REACTIONS][NOMBRE_ESPECES];

	if (NOMBRE_ESPECES == 1){
		FILE* f = NULL; 
		f = fopen("input pour une espece.txt", "r");
		if (f != NULL){
	        fscanf(f, "%lf %lf %lf %lf %d %d", 
	        	&TEMPS_MAX, &x[0], &c[0], &c[1],
	        	&v[0][0], &v[1][0]);

	        fclose(f);
	    }
	}else if (NOMBRE_ESPECES == 2){
		FILE* f = NULL; 
		f = fopen("input pour deux especes.txt", "r");
		if (f != NULL){
	        fscanf(f, "%lf %lf %lf %lf %lf %lf %lf %d %d %d %d %d %d %d %d", 
	        	&TEMPS_MAX, &x[0], &x[1], &c[0], &c[1], &c[2], &c[3],
	        	&v[0][0], &v[0][1], &v[1][0], &v[1][1], &v[2][0], &v[2][1], &v[3][0], &v[3][1]);

	        fclose(f);
	    }
	}else if (NOMBRE_ESPECES == 3){
		FILE* f = NULL; 
		f = fopen("input pour trois especes.txt", "r");
		if (f != NULL){
	        fscanf(f, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d", 
        	&TEMPS_MAX, 
        	&x[0], &x[1], &x[2], 
        	&c[0], &c[1], &c[2], &c[3], &c[4], &c[5], &c[6], &c[7], &c[8],
        	&v[0][0], &v[0][1], &v[0][2],
        	&v[1][0], &v[1][1], &v[1][2],
        	&v[2][0], &v[2][1], &v[2][2],
        	&v[3][0], &v[3][1], &v[3][2],
        	&v[4][0], &v[4][1], &v[4][2],
        	&v[5][0], &v[5][1], &v[5][2],
        	&v[6][0], &v[6][1], &v[6][2],
        	&v[7][0], &v[7][1], &v[7][2],
        	&v[8][0], &v[8][1], &v[8][2] );

	        fclose(f);
	    }
	}

    Gillespie("data.txt", c, NOMBRE_ESPECES, v, x, NOMBRE_DE_REACTIONS, TEMPS_MAX);

    return 0;
	
}
	// printf("%lf\n", TEMPS_MAX);
	// for (int i = 0; i < NOMBRE_ESPECES; ++i){
	// 	printf("x[%d] = %lf ", i, x[i] );
	// }printf("\n");


	// for (int i = 0; i < NOMBRE_DE_REACTIONS; ++i){
	// 	printf("c[%d] = %lf ", i, c[i] );
	// }printf("\n");

	// double* h = calculDeH(x, NOMBRE_DE_REACTIONS);
	// for (int i = 0; i < NOMBRE_DE_REACTIONS; ++i){
	// 	printf("h[%d] = %lf ", i, h[i] );
	// }printf("\n");

	// double* a = calculDePropensity(h, c, NOMBRE_DE_REACTIONS);
	// for (int i = 0; i < NOMBRE_DE_REACTIONS; ++i){
	// 	printf("a[%d] = %lf ", i, a[i] );
	// }printf("\n");

	// for (int i = 0; i < NOMBRE_DE_REACTIONS; ++i){
	// 	for (int j = 0; j < NOMBRE_ESPECES; ++j){
	// 		printf("v[%d][%d] = %d ", i, j, v[i][j]);
	// 	}
	// }printf("\n");

	// srand(time(NULL));
	// double r2 = genererNombreAlea();
	// int mu = calculDMu(a, r2, NOMBRE_DE_REACTIONS);