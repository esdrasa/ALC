#include <stdio.h>
#include "operacoes.h"
#include "metodos.h"
#include "saida.h"
#include "entrada.h"

int main()
{
    double **A, *coluna, produtoI, angle;
    int n, teste;

    printf("Digite o tamanho da matriz: ");
    scanf("%d%*c", &n);

    A = lerMatriz(n, n);

    if(teste = isVandermonde(A, n))
        printf("Eh uma matriz de Vandermonde.\n");
    else
        printf("Nao eh uma matriz de Vandermonde.\n");

    printf("Norma de Frobenius: %lf\n", normaF(A, n));
    printf("Norma Linha: %lf\n", normaLinha(A, n));
    printf("Norma Coluna: %lf\n", normaColuna(A, n));

    coluna = getColuna(A, 1, n);

    produtoI = produtoEscalar(A[0], coluna, n);

    angle = angulo(A[0], coluna, n);

    printf("Produto interno entre a primeira linha e a segunda coluna: %lf\n", produtoI);
    printf("Angulo em radianos entre a primeira linha e a segunda coluna: %lf\n", angle);

    if(teste)
        printf("Determinante: %lf\n", detVandermonde(A, n));

    getchar();

    return 0;
}