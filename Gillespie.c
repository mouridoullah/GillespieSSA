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
	// Formations
	// *          -> S1                 => hu = 1
	// *          -> 2S2                => hu = 1

	// Destructions
	// S1         -> reactions products => hu = X1
	// S2         -> reactions products => hu = X2

	// Bimoleculaire
	// 2S1        -> reactions products => hu = 0.5 * X1(X1 - 1)
	// 2S2        -> reactions products => hu = 0.5 * X2(X2 - 1)
	// S1 + S2    -> reactions products => hu = X1 * X2

	// Avec trois molecules
	// S1 + 2S2   -> reactions products => hu = 0.5 * X1X2(X2 - 1)
	// 3S1        -> reactions products => hu = (1/6) * X1(X1 - 1)(X1 - 2)
		double* h = malloc(M*sizeof(double));
	    h[0] = x[0]*x[3];
	    h[1] = 1.0;	
	    h[2] = x[18]*x[14];
	    h[3] = 1.0;
	    h[4] = x[15]*x[0];
	    h[5] = 1.0;
	    h[6] = x[20]*x[14];
	    h[7] = x[19]*x[14];
	    h[8] = 1.0;

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
				// On arrete le programme si l'une des especes est totalement consommé
				printf("Fin du programme\n");	
				//exit(1);
				return;
			} 
		}
	}
	free(a);
	free(h);
}
/*-------------------------------------------------------------*/
int main(){
	double start,stop;
	start = clock();


	double TEMPS_MAX = 2500.0;
	int NOMBRE_ESPECES, NOMBRE_DE_REACTIONS, pos = 20;

    FILE* out = NULL;
    out = fopen("input.txt", "r"); 
    if(out == NULL){
        printf("\nError cannot open file");
        exit(1);
    }

    fseek(out, pos, SEEK_SET);
    fscanf(out ,"%d", &NOMBRE_ESPECES);
    printf("Nombre de produits: %d\n", NOMBRE_ESPECES);

    pos += 25;
    fseek(out, pos, SEEK_SET);
    fscanf(out, "%d", &NOMBRE_DE_REACTIONS);
    printf("Nombre de reactions: %d\n", NOMBRE_DE_REACTIONS);

    printf("######################################################\n");

    double x[NOMBRE_ESPECES], c[NOMBRE_DE_REACTIONS];
    int v[NOMBRE_DE_REACTIONS][NOMBRE_ESPECES];

    pos += 53;
    fseek(out, pos, SEEK_SET);
    fscanf(out,"%lf", &x[0]);
    printf("x[0]: %f\n", x[0]);

    for (int i = 1; i < NOMBRE_ESPECES; ++i){
    	pos += 22;
    	fseek(out, pos, SEEK_SET);
    	fscanf(out,"%lf", &x[i]);
    	printf("x[%d]: %f\n", i, x[i]);
    }
    
    printf("\n######################################################\n\n");

    pos += 60;
    fseek(out, pos, SEEK_SET);
    fscanf(out,"%lf", &c[0]);
    printf("c[0]: %f\n", c[0]);

   for (int i = 1; i < NOMBRE_DE_REACTIONS; ++i){
   	pos += 74;
   	fseek(out, pos, SEEK_SET);
   	fscanf(out,"%lf", &c[i]);
   	printf("c[%d]: %f\n", i, c[i]);
   }    

    //printf("\n%ld\n", ftell(out));
   
    // char d;
    // for (pos = 1195; pos < 1413; ++pos){
    //   	fseek(out, pos, SEEK_SET);
    //  	fscanf(out, "%c", &d);
    //  	printf("%c", d);
    // }
    // // printf("\n%ld\n", ftell(out));

    printf("\n######################################################\n\n");

    pos += 220;
    fseek(out, pos, SEEK_SET);
    fscanf(out,"%d", &v[0][0]);
	printf("v[0][0]: %d\n", v[0][0]);
    for (int i = 0; i < NOMBRE_DE_REACTIONS; ++i){
    	for (int j = 1; j < NOMBRE_ESPECES; ++j){
    		pos += 8;
    	    fseek(out, pos, SEEK_SET);
   			fscanf(out,"%d", &v[i][j]);	
   			printf("v[%d][%d]: %d\n", i, j, v[i][j]);
    	}
    	printf("\n");
    	pos += 9;
        fseek(out, pos, SEEK_SET);
    	fscanf(out,"%d", &v[i][0]);
	 	if (i < NOMBRE_DE_REACTIONS - 1) printf("v[%d][0]: %d\n", i, v[i][0]);
    }

    fclose(out);

    printf("\n######################################################\n\n");

	printf("TEMPS_MAX: %lf\n", TEMPS_MAX);
	for (int i = 0; i < NOMBRE_ESPECES; ++i){
		printf("x[%d] = %lf\n", i, x[i] );
	}printf("\n");
    printf("\n######################################################\n\n");
	for (int i = 0; i < NOMBRE_DE_REACTIONS; ++i){
		printf("c[%d] = %lf\n", i, c[i] );
	}printf("\n");
    printf("\n######################################################\n\n");
	double* h = calculDeH(x, NOMBRE_DE_REACTIONS);
	for (int i = 0; i < NOMBRE_DE_REACTIONS; ++i){
		printf("h[%d] = %lf\n", i, h[i] );
	}printf("\n");
    printf("\n######################################################\n\n");
   	double* a = calculDePropensity(h, c, NOMBRE_DE_REACTIONS);
	for (int i = 0; i < NOMBRE_DE_REACTIONS; ++i){
		printf("a[%d] = %lf\n", i, a[i] );
	}printf("\n");

 //    printf("\n######################################################\n\n");


	Gillespie("data.txt", c, NOMBRE_ESPECES, v, x, NOMBRE_DE_REACTIONS, TEMPS_MAX);

	stop = clock();
	double k;
	k = (stop - start)/1000 ;
	printf("Temps d'exucution du programme: %f\n",k);


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

	// double TEMPS_MAX; int NOMBRE_DE_REACTIONS;
	// double x[NOMBRE_ESPECES];
	// double c[NOMBRE_DE_REACTIONS];
	// int v[NOMBRE_DE_REACTIONS][NOMBRE_ESPECES];