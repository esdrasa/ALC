#include <stdio.h>
#include "operacoes.h"
#include "metodos.h"
#include "saida.h"
#include "entrada.h"

int main()
{   
    //declarando as variaveis a serem usadas
    double **A, *b, *x, *xAux, cLinha, cColuna, cNorma, cSassenfeld, tol;
    int n;
    unsigned long int iMax = 1000;

    //lendo as informacoes necessarias
    printf("Digite o tamanho do sistema: ");
    scanf("%d%*c", &n);

    //Alocando espaço para os vetores e matrizes
    x = criaVetor(n);

    xAux = criaVetor(n);

    A = lerMatriz(n, n);

    b = lerVetor(n);

    //executando os criterios
    cLinha = criterioLinhas(A, n);
    cColuna = criterioColunas(A, n);
    cNorma = criterioNorma(A, n);
    cSassenfeld = criterioSassenfeld(A, n);

    //prints para verem se os criterios foram atendidos
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

    //Verificadno se os metodos podem ser usados
    //a partir da garantia de convergencia
    //dos criterios

    if(cNorma || cColuna || cLinha)
    {
        printf("O sistema pode ser resolvido usando o metodo de Jacobi.\n");
        jacobi(A, b, x, 0, iMax, n);
        printf("A solucao do sistema pelo metodo de Jacobi e: \n");
        imprimeVetor(x, n);
    }
    

    if(cNorma || cColuna || cLinha || cSassenfeld)
    {
        printf("O sistema pode ser resolvido usando o metodo de Gauss-Seidel.\n");
        gaussSeidel(A, b, xAux, 0, iMax, n);
        printf("A solucao do sistema pelo metodo de Gauss-Seidel e: \n");
        imprimeVetor(xAux, n);
    }
    
    //Se os criterios nao garatem convergencia
    //nada podemos afirmar
    if(!(cSassenfeld || cNorma || cColuna || cLinha))
        printf("Nada podemos afirmar sobre a resolucao do sistema pelos metodos.\n");

    getchar();

    //liberando espaço alocado
    liberaVetor(x);
    liberaVetor(xAux);
    liberaVetor(b);
    liberaMatriz(A, n);

    return 0;
}