#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "metodos.h"
#include "operacoes.h"
#include "saida.h"
#include "entrada.h"

int main(){

    double **matrizA, *vetorB, *vetorX, *residuoPSO, *residuoSimulated;
    double tol, tempMax, tempMin, alfa, normaPSO, normaSimulated;
    int n, particulas, i, iteracoes;

    printf("Digite a dimensao da matriz A: ");
    scanf("%d%*c", &n);

    printf("\nDigite a quantidade de particulas para o PSO: ");
    scanf("%d%*c", &particulas);

    printf("\nDigite a tolerancia para o PSO: ");
    scanf("%lf%*c", &tol);

    printf("\nDigite a temperatura maxima: ");
    scanf("%lf%*c", &tempMax);

    printf("\nDigite a temperatura minima: ");
    scanf("%lf%*c", &tempMin);

    printf("\nDigite a variacao da temperatura: ");
    scanf("%lf%*c", &alfa);

    printf("\nDigite a quantidade de iteracoes: ");
    scanf("%d%*c", &iteracoes);
    
    srand((unsigned)time(NULL));

    matrizA = lerMatriz(n,n);
    vetorB  = lerVetor(n);

    vetorX           = criaVetor(n);
    residuoPSO       = criaVetor(n);
    residuoSimulated = criaVetor(n);

    PSO(matrizA, n, vetorB, vetorX, tol, particulas);

    printf("\nSolucao com PSO\n");

    imprimeVetor(vetorX, n);
    
    residuo(matrizA, vetorX, vetorB, residuoPSO, n);
    normaPSO = normaDois(residuoPSO, n);

    printf("\n");

    ///Simulated Annealing

    simulatedAnnealing(matrizA, vetorB, vetorX, alfa, tempMax, tempMin, iteracoes, n);

    printf("\nSolucao com Simulated Annealing\n");
    for(i = 0; i < n; i++)
    {
        printf("%lf\n", vetorX[i]);
    }


    residuo(matrizA, vetorX, vetorB, residuoSimulated, n);
    normaSimulated = normaDois(residuoSimulated, n);

    printf("\n\n");

    if(normaPSO < normaSimulated)
    {
        printf("O PSO foi mais eficiente nesse caso\n\n");
    }
    else
        printf("O Simulated Annealing foi mais eficiente nesse caso\n\n");
    
    printf("Norma 2 das solucoes:\n");
    printf("PSO: %lf\n", normaPSO);
    printf("Simulated Annealing: %lf\n", normaSimulated);
    
    getchar();
    
    liberaMatriz(matrizA, n);
    liberaVetor(vetorB);
    liberaVetor(vetorX);
    liberaVetor(residuoPSO);
    liberaVetor(residuoSimulated);

    return 0;
}