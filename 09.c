#include <stdio.h>
#include "entrada.h"
#include "saida.h"
#include "operacoes.h"
#include "metodos.h"

/*Resolver problemas de minimos quadrados usando A'Ax = A'b, onde A' seria A transposto*/

int main(void) {

    double **matriz, *vetorB, *vetorB2, **matrizTrans, **matrizMultA, *matrizMultB, *vetorX, **matrizL, **matrizU;

    int n, m;
    double x, resultado;

    printf("Digite a quantidade de linha(s): ");
    scanf("%d%*c",&n);

    printf("Digite a quantidade de coluna(s): ");
    scanf("%d%*c",&m);

    matriz = lerMatriz(n, m);
    vetorB = lerVetor(n);

    matrizL = criaMatrizI(n);
    matrizU = criaMatriz(n, n);

    vetorX = criaVetor(n);

    printf("\nResolvendo problemas de minimos quadrados da seguinte maneira:\n");
    printf("A'Ax = A'b.\n");
    printf("Obs: A' = A(transposto).\n");

    matrizTrans = criaMatriz(m, n);
    transposta(matriz, matrizTrans, n, m); // Transpor a matriz A

    matrizMultA = criaMatriz(m, m);
    multiplica(matrizTrans, m, n, matriz, n, m, matrizMultA); // Multiplica a matriz A' por A.

    matrizMultB = criaVetor(m);
    multiplicaVetor(matrizTrans, m, n, vetorB, matrizMultB, m); // Multiplica a matriz A' pelo vetor B

    /*Nessa etapa já tem o novo sistema Ax=b*/

    printf("\nMatriz A'A:\n");
    imprimeMatriz(matrizMultA, m, m); //Imprime o que para nós será a matriz A a partir de agora

    printf("\nVetor A'B.\n");
    imprimeVetor(matrizMultB, m); // Imprime o que para nós será o vetor B a partir de agora

    /*Agora só resta resolver o novo sistema*/

    printf("Utilizando o metodo LU.\n");


    lu(matrizMultA, matrizL, matrizU, m);

    printf("\nMatriz L:\n");
    imprimeMatriz(matrizL, m, m);

    printf("\nMatriz U:\n");
    imprimeMatriz(matrizU, m, m);

    //Explicação do método
    printf("Ax = b\n");
    printf("A = LU\n");
    printf("LUx = b\n");
    printf("Ux = y\n");
    printf("Ly = b    (resolucao por substituicao para frente).\n");
    printf("Ux = y    (resolucao por substituicao para tras).\n");

    vetorB2 = criaVetor(m);
    forwardSub(matrizL, matrizMultB, vetorB2, m);

    vetorX = criaVetor(m);
    backSub(matrizU, vetorB2, vetorX, m);

    printf("\nObtido o vetor solucao.\n");

    printf("\nVetor solucao(A e B, respectivamente):\n");
    imprimeVetor(vetorX, m); // Colocar para imprimir o A e o B separados

    printf("Digite um valor para X:\n");
    scanf("%lf%*c",&x);

    resultado = vetorX[0]*x + vetorX[1];
    printf("Resultado com o X solicitado: %lf",resultado);

    getchar();
    
    liberaMatriz(matriz, n);
    liberaMatriz(matrizL, n);
    liberaMatriz(matrizMultA, n);
    liberaMatriz(matrizTrans, n);
    liberaMatriz(matrizU, n);
    liberaVetor(vetorB);
    liberaVetor(vetorB2);
    liberaVetor(vetorX);
    liberaVetor(matrizMultB);

    return 0;
}