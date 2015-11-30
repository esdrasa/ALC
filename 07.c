#include <stdio.h>
#include "operacoes.h"
#include "metodos.h"
#include "saida.h"
#include "entrada.h"

int main()
{
    double **A, cLinha, cColuna, cNorma, cSassenfeld;
    int n;

    printf("Digite o tamanho do sistema: ");
    scanf("%d%*c", &n);

    A = lerMatriz(n, n);

    cLinha = criterioLinhas(A, n);
    cColuna = criterioColunas(A, n);
    cNorma = criterioNorma(A, n);
    cSassenfeld = criterioSassenfeld(A, n);

    printf("\nCriterio das Linhas: ");
    if(cLinha)
        printf("[X]");
    else
        printf("[ ]");
    printf(" %lf\n", cLinha);


    printf("Criterio das Colunas: ");
    if(cColuna)
        printf("[X]");
    else
        printf("[ ]");
    printf(" %lf\n", cColuna);


    printf("Criterio das Normas: ");
    if(cNorma)
        printf("[X]");
    else
        printf("[ ]");
    printf(" %lf\n", cNorma);


    printf("Criterio de Sassenfeld: ");
    if(cSassenfeld)
        printf("[X]");
    else
        printf("[ ]");
    printf(" %lf\n\n", cSassenfeld);

    if(cNorma || cColuna || cLinha)
        printf("O sistema pode ser resolvido usando o metodo de Jacobi.\n");
    

    if(cNorma || cColuna || cLinha || cSassenfeld)
        printf("O sistema pode ser resolvido usando o metodo de Gauss-Seidel.\n");        
    

    if(!(cSassenfeld || cNorma || cColuna || cLinha))
        printf("Nada podemos afirmar sobre a resolucao do sistema pelos metodos.\n");

    getchar();

    liberaMatriz(A, n);
    return 0;
}