#include <stdlib.h>
#include <stdio.h>

/*--------------------------------------------------------*/
void remplirTab(int t[], int n){
	for (int i = 0; i < n; ++i){
		printf("Donner la case numero %d : ", i+1);
		scanf("%d", &t[i]);
	}
}
/*--------------------------------------------------------*/
void afficherTab(int t[], int n){
	// for (int i = 0; i < n; ++i){
	// 	printf("%d\t ", t[i]);
	// }
	int i;
	for ( i = 0, printf("["); i < n; ++i){
		printf("%d", t[i]);
		if(i < n-1) printf(", ");
	}
	printf("]");
}
/*--------------------------------------------------------*/
void remplirTab2D(int **t, int l, int c){
	for (int i = 0; i < l; ++i){
		printf("Veulliez remplir la %d la ligne : \n", i+1);
		remplirTab(t[i], c);
	}
}
/*--------------------------------------------------------*/
void afficherTab2D(int** t, int l, int c){
	for (int i = 0; i < l; ++i){
		afficherTab(t[i], c);
		printf("\n");
	}
}
/*--------------------------------------------------------*/
void allouerTableau(int** t, int n){
	*t = (int*)malloc(n*sizeof(int));
}
/*--------------------------------------------------------*/
/*--------------------------------------------------------*/
/*--------------------------------------------------------*/
int main(int argc, char const *argv[]){

	int n = 2, m = 4;
	int** tab2D = (int**)malloc(n*sizeof(int*));
	for (int i = 0; i < n; ++i){
		allouerTableau( &(tab2D[i]), m);
		//remplirTab(tab2D[i], n);
	}
	remplirTab2D(tab2D,n,m);
	afficherTab2D(tab2D, 5, 5);
	free(tab2D);


	/*--------------------------------------------------------*/
	// int tab[5][4];
	// remplirTab2D(tab,3,3);
	// afficherTab2D(tab,3,3);
	/*--------------------------------------------------------*/
	return 0;
}