#include <stdio.h>
#include "entrada.h"
#include "saida.h"
#include "operacoes.h"
#include "metodos.h"

/*Resolver problemas de minimos quadrados usando A'Ax = A'b, onde A' seria A transposto*/

int main(void) {

   double **matriz, *vetorB, *vetorB2, **matrizTrans, **matrizMultA, *matrizMultB, *vetorX, **matrizL, **matrizU;

    int n, m;

   printf("Digite a quantidade de linha(s): ");
   scanf("%d",&n);

   printf("Digite a quantidade de coluna(s): ");
   scanf("%d",&m);

   matriz = lerMatriz(n, m);
   vetorB = lerVetor(n);

   matrizL = criaMatrizI(n);
   matrizU = criaMatriz(n, n);

   vetorX = criaVetor(n);

   printf("\nResolvendo problemas de minimos quadrados da seguinte maneira:\n");
   printf("A'Ax = A'b.\n");
   printf("Obs: A' = A(transposto).\n");

   matrizTrans = transposta(matriz, n, m); // Transpor a matriz A
   matrizMultA = multiplica(matrizTrans, m, n, matriz, n, m); // Multiplica a matriz A' por A.
   matrizMultB = multiplicaVetor(matrizTrans, m, n, vetorB, m); // Multiplica a matriz A' pelo vetor B

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


   vetorB2 = forwardSub(matrizL, matrizMultB, m);

   vetorX = backSub(matrizU, vetorB2, m);

   printf("\nObtido o vetor solucao.\n");

   printf("\nVetor solucao:\n");
   imprimeVetor(vetorX, m);

   system("pause");

    return 0;
}
