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
    coluna = lerVetor(n);

    if(teste = isVandermonde(A, n))
        printf("Eh uma matriz de Vandermonde.\n");
    else
        printf("Nao eh uma matriz de Vandermonde.\n");

    printf("Norma de Frobenius: %lf\n", normaF(A, n));
    printf("Norma Linha: %lf\n", normaLinha(A, n));
    printf("Norma Coluna: %lf\n", normaColuna(A, n));
    
    //Atribui os elementos da coluna 1 (segunda coluna) ao vetor "coluna"
    getColuna(A, coluna, 1, n);
    
    //Produto interno entre a primeira linha e a segunda coluna
    produtoI = produtoEscalar(A[0], coluna, n);
    
    //Ângulo entre a primeira e linha e a segunda coluna 
    angle = angulo(A[0], coluna, n);

    printf("Produto interno entre a primeira linha e a segunda coluna: %lf\n", produtoI);
    printf("Angulo em radianos entre a primeira linha e a segunda coluna: %lf\n", angle);
    
    //Se a matriz lida do arquivo for uma matriz de Vandermonde, será possível calcular seu determinante
    if(teste)
        printf("Determinante: %lf\n", detVandermonde(A, n));

    getchar();
    
    liberaMatriz(A, n);
    liberaVetor(coluna);

    return 0;
}