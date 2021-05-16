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
	    h[0] = x[0]*x[1];
	    h[1] = x[3]*x[4];	
	    h[2] = x[6]*x[4];
	    h[3] = x[9]*x[4];
	    h[4] = x[11];
	    h[5] = x[12]*x[13];
	    h[6] = x[5];
	    h[7] = x[14]*x[4];
	    h[8] = x[15]*x[16];
		h[9] = x[19]*x[4];
		h[10] = x[19]*x[4];
		h[11] = x[9]*x[4];
		h[12] = x[22]*x[4];
		h[13] = x[23]*x[4];
		h[14] = x[24]*x[4];
	    h[15] = x[25]*x[4];	
	    h[16] = x[26]*x[4];
	    h[17] = x[8]*x[5];
	    h[18] = x[28]*x[28];
	    h[19] = x[11]*x[4];
	    h[20] = x[32]*x[21];
	    h[21] = x[33]*x[34];
	    h[22] = x[37]*x[38];
		h[23] = x[38]*x[5];
		h[24] = x[19]*x[5];
		h[25] = x[35]*x[13];
		h[26] = x[42]*x[43];
		h[27] = x[45]*x[5];
		h[28] = x[32]*x[4];
	    h[29] = x[31]*x[32];	
	    h[30] = x[48]*x[6];
	    h[31] = x[50]*x[6];
	    h[32] = x[52]*x[6];
	    h[33] = x[48]*x[29];
	    h[34] = x[55]*x[56];
	    h[35] = x[58]*x[59];
	    h[36] = x[62]*x[6];

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
		//printf("\n######################################################\n");
		h = calculDeH(x, M);
		a = calculDePropensity(h, c, M);
		// for (int i = 0; i < M; ++i){
		// 	printf("a[%d] = %lf, h[%d] = %lf  \n", i, a[i], i, h[i] );
		// }
		a0 = sommeDesA(a, M);
		//printf("la somme Des propensions : %lf\n", a0);


		r1 = genererNombreAlea();
		r2 = genererNombreAlea();

		if(a0 == 0 ) break;

		t += log(1.0 / r1) / a0;

		mu = calculDMu(a, r2, M);
		printf("%d ", mu);


		miseAJourDesX(x, N, v, mu);
		ecrireDansData(data, t, x, N);

		for (int i = 0; i < N; ++i){
			if (x[i] <= 0){
				// On arrete le programme si l'une des especes est totalement consommé
				printf("\n\nFIN DU PROGRAMME\n");	
				//exit(1);
				return;
			} 
		}
	}
	printf("\n");
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

    //printf("######################################################\n");

    double x[NOMBRE_ESPECES], c[NOMBRE_DE_REACTIONS];
    int v[NOMBRE_DE_REACTIONS][NOMBRE_ESPECES];

    pos += 16;
    fseek(out, pos, SEEK_SET);
    fscanf(out,"%lf", &x[0]);
    //printf("x[0]: %f\n", x[0]);

    for (int i = 1; i < NOMBRE_ESPECES; ++i){
    	pos += 25;
    	fseek(out, pos, SEEK_SET);
    	fscanf(out,"%lf", &x[i]);
    	//printf("x[%d]: %f\n", i, x[i]);
    }
    
    //printf("\n######################################################\n\n");

    pos += 67;
    fseek(out, pos, SEEK_SET);
    fscanf(out,"%lf", &c[0]);
    //printf("c[0]: %f\n", c[0]);

   for (int i = 1; i < NOMBRE_DE_REACTIONS; ++i){
		pos += 78;
		fseek(out, pos, SEEK_SET);
		fscanf(out,"%lf", &c[i]);
		//printf("c[%d]: %f\n", i, c[i]);
   }    

    // printf("\n%ld\n", ftell(out));
   
    // char d;
    // for (pos = 4491; pos < 4784; ++pos){
    //   	fseek(out, pos, SEEK_SET);
    //  	fscanf(out, "%c", &d);
    //  	printf("%c", d);
    // }
    // printf("\n%ld\n", ftell(out));

    //printf("\n######################################################\n\n");

    pos += 297;
    fseek(out, pos, SEEK_SET);
    fscanf(out,"%d", &v[0][0]);
	//printf("v[0][0]: %d ", v[0][0]);
    for (int i = 0; i < NOMBRE_DE_REACTIONS; ++i){
    	for (int j = 1; j < NOMBRE_ESPECES; ++j){
    		pos += 4;
    	    fseek(out, pos, SEEK_SET);
   			fscanf(out,"%d", &v[i][j]);	
   			//printf("v[%d][%d]: %d ", i, j, v[i][j]);
    	}
    	//printf("\n");
    	pos += 6;
        fseek(out, pos, SEEK_SET);
    	fscanf(out,"%d", &v[i+1][0]);
	 	//if (i < NOMBRE_DE_REACTIONS - 1) printf("v[%d][0]: %d ", i, v[i][0]);
    }

    


    printf("\n------------------Les données de départ:--------------------\n");

	printf("TEMPS_MAX: %lf\n", TEMPS_MAX);
	for (int i = 0; i < NOMBRE_ESPECES; ++i){
		printf("x[%d] = %lf ", i, x[i] );
		if (i != 0 && i%10 == 0) printf("\n");
	}
    printf("\n\n----------------a, h c à t0:---------------------------------\n");
	// for (int i = 0; i < NOMBRE_DE_REACTIONS; ++i){
	//  	for (int j = 0; j < NOMBRE_ESPECES; ++j){
	//  		printf("v[%d][%d]: %d ", i, j, v[i][j]);
	//  	}printf("\n");
	// }printf("\n");
 //    printf("\n######################################################\n\n");
	double* h = calculDeH(x, NOMBRE_DE_REACTIONS);
	// for (int i = 0; i < NOMBRE_DE_REACTIONS; ++i){
	// 	printf("h[%d] = %lf\n", i, h[i] );
	// }printf("\n");
 //    printf("######################################################\n");
	// for (int i = 0; i < NOMBRE_DE_REACTIONS; ++i){
	// 	printf("c[%d] = %lf\n", i, c[i] );
	// } //printf("\n");

 //    printf("\n######################################################\n\n");
   	double* a = calculDePropensity(h, c, NOMBRE_DE_REACTIONS);
	for (int i = 0; i < NOMBRE_DE_REACTIONS; ++i){
		printf("a[%d] = %lf, h[%d] = %lf, c[%d] = %lf\n", i, a[i], i, h[i], i, c[i]);
	}
    printf("\n\n--------------------Numeros de reactions qui ont eu lieu:--------------\n");

 //    srand(time(NULL));
 //    double r2 = 0;
	// int mu = 0;
	// for (int i = 0; i < 105; ++i){
	// 	r2 = genererNombreAlea();
	// 	mu = calculDMu(a, r2, NOMBRE_DE_REACTIONS);	
	// 	printf("%d ", mu);
	// }printf("\n");

 //    printf("\n######################################################\n\n");

    fclose(out);

	Gillespie("data.txt", c, NOMBRE_ESPECES, v, x, NOMBRE_DE_REACTIONS, TEMPS_MAX);

	stop = clock();
	double k;
	k = (stop - start)/1000 ;
	printf("\nTemps d'exucution du programme: %f\n",k);



    return 0;
}





