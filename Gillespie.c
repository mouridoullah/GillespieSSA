#include <stdio.h>
#include <math.h>
#include <errno.h>
#include <time.h>
#include <stdlib.h>

																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																															

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
		double* h = malloc(M*sizeof(double));
	    h[0] = x[0]*x[1];
	    h[1] = 1.0;	

	    return h;
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

	ecrireDansData(data, t, x, N);
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

		for (int i = 0; i < N; ++i){
			if (x[i] <= 0){
				printf("Fin du programme\n");	
				exit(1);
			} 
		}


	}
	free(a);
	free(h);
}
/*-------------------------------------------------------------*/
int main(){
	double TEMPS_MAX = 1000.0;
	int NOMBRE_ESPECES, NOMBRE_DE_REACTIONS;
    FILE* out = NULL;
    out = fopen("test.txt", "r"); 
    if(out == NULL){
        printf("\nError cannot open file");
        exit(1);
    }

    fseek(out, 20, SEEK_SET);
    //printf("on est a la position: %ld\n", ftell(out));
    fscanf(out ,"%d", &NOMBRE_ESPECES);
    fseek(out, 44, SEEK_SET);
    //printf("on est a la position: %ld\n", ftell(out));
    fscanf(out, "%d", &NOMBRE_DE_REACTIONS);
    printf("Nombre de produits: %d\nNombre de reaction : %d\n", NOMBRE_ESPECES, NOMBRE_DE_REACTIONS);

    double x[NOMBRE_ESPECES], c[NOMBRE_DE_REACTIONS];
    int v[NOMBRE_DE_REACTIONS][NOMBRE_ESPECES];

    fseek(out, 85, SEEK_SET);
    //printf("on est a la position: %ld\n", ftell(out));
    fscanf(out,"%lf", &c[0]);
    printf("c1 = %lf\n", c[0]);

    fseek(out, 126, SEEK_SET);
    //printf("on est a la position: %ld\n", ftell(out));
    fscanf(out,"%lf", &c[1]);
    printf("c2 = %lf\n", c[1]);

    //char d;
    fseek(out, 175, SEEK_SET);
    fscanf(out,"%lf", &x[0]);
    printf("x1, nombre de molecules d'Oxygène: %lf\n", x[0]);

    fseek(out, 184, SEEK_SET);
    fscanf(out,"%lf", &x[1]);
    printf("x2, nombre de molecules d'Hydrogene: %lf\n", x[1]);

    fseek(out, 197, SEEK_SET);
    fscanf(out,"%lf", &x[2]);
    printf("x3, nombre de molecules d'Eau H2O: %lf\n", x[2]);

    fseek(out, 220, SEEK_SET);
    fscanf(out,"%d %d %d %d %d %d", &v[0][0], &v[0][1], &v[0][2], &v[1][0], &v[1][1], &v[1][2]);
    for (int i = 0; i < NOMBRE_DE_REACTIONS; ++i){
        for (int j = 0; j < NOMBRE_ESPECES; ++j){
            printf("%d ", v[i][j]);
        }
        printf("\n");
    }

    fclose(out);
/*-----------------------------------------------------------------------------------------------------------*/
	// double r1 = 0, r2 = 0, a0 = 0, t = 0;
	// double* h;
	// double* a;
	// int mu = 0;

	// srand(time(NULL));

	// for (int i = 0; i < NOMBRE_ESPECES; ++i){
	// 		printf("x[%d] = %lf ", i, x[i] );
	// 	}printf("\n");

	// while(t < TEMPS_MAX) {
	// 	h = calculDeH(x, NOMBRE_DE_REACTIONS);
	// 	for (int i = 0; i < NOMBRE_DE_REACTIONS; ++i){
	// 	printf("h[%d] = %lf ", i, h[i] );
	// 	}printf("\n");

	// 	a = calculDePropensity(h, c, NOMBRE_DE_REACTIONS);
	// 	for (int i = 0; i < NOMBRE_DE_REACTIONS; ++i){
	// 		printf("a[%d] = %lf ", i, a[i] );
	// 	}printf("\n");

	// 	a0 = sommeDesA(a, NOMBRE_DE_REACTIONS);

	// 	r1 = genererNombreAlea();
	// 	r2 = genererNombreAlea();

	// 	if(a0 == 0 ) break;

	// 	t += log(1.0 / r1) / a0;
	// 	mu = calculDMu(a, r2, NOMBRE_DE_REACTIONS);
	// 	printf("m = %d :\n", mu);

	// 	miseAJourDesX(x, NOMBRE_ESPECES, v, mu);

	// 	for (int i = 0; i < NOMBRE_ESPECES; ++i){
	// 		if (x[i] <= 0){
	// 			printf("Fin du programme\n");	
	// 			exit(1);
	// 		} 
	// 	}

	// 	for (int i = 0; i < NOMBRE_ESPECES; ++i){
	// 		printf("x[%d] = %lf ", i, x[i] );
	// 	}printf("\n");

	// }
	// free(a);
	// free(h);

 	//Gillespie("data.txt", c, NOMBRE_ESPECES, v, x, NOMBRE_DE_REACTIONS, TEMPS_MAX);
/*-----------------------------------------------------------------------------------------------------------*/
    printf("%lf\n", TEMPS_MAX);
	for (int i = 0; i < NOMBRE_ESPECES; ++i){
		printf("x[%d] = %lf ", i, x[i] );
	}printf("\n");


	for (int i = 0; i < NOMBRE_DE_REACTIONS; ++i){
		printf("c[%d] = %lf ", i, c[i] );
	}printf("\n");

	double* h = calculDeH(x, NOMBRE_DE_REACTIONS);
	for (int i = 0; i < NOMBRE_DE_REACTIONS; ++i){
		printf("h[%d] = %lf ", i, h[i] );
	}printf("\n");

	double* a = calculDePropensity(h, c, NOMBRE_DE_REACTIONS);
	for (int i = 0; i < NOMBRE_DE_REACTIONS; ++i){
		printf("a[%d] = %lf ", i, a[i] );
	}printf("\n");

	for (int i = 0; i < NOMBRE_DE_REACTIONS; ++i){
		for (int j = 0; j < NOMBRE_ESPECES; ++j){
			printf("v[%d][%d] = %d ", i, j, v[i][j]);
		}printf("\n");
	}printf("\n");

	srand(time(NULL));
	double r2 = genererNombreAlea();
	int mu;

	for (int i = 0; i < 10; ++i){
		mu = calculDMu(a, r2, NOMBRE_DE_REACTIONS);	
		printf("La reaction n° %d  a eu lieu\n", mu);
	}printf("\n");
/*-----------------------------------------------------------------------------------------------------------*/



    return 0;
}
