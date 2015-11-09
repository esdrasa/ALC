#include <stdio.h>
#include "entrada.h"
#include "saida.h"
#include "metodos.h"
#include "operacoes.h"
#define COLUNA1 0
#define COLUNA3 2

/* Descobrindo a matriz Laplaciana usando L = D - A
   D = matriz diagonal formada pelos graus dos vértices e
   A = matriz de Adjacência   */


int main(void) {

   double **matrizL, **matrizD, **matrizA, *vetorColuna1, *vetorColuna2;
   int n, i, j;

   printf("Digite a dimensao da matriz\n");
   scanf("%d%*c",&n);

   matrizD = criaMatriz(n, n);
   matrizL = criaMatriz(n, n);

   vetorColuna1 = criaVetor(n);
   vetorColuna2 = criaVetor(n);

   matrizA = lerMatriz(n, n); // Leu a matriz A

   matrizNula(matrizD, n, n); // Preencheu a matrizD com zeros

   encontraGraus(matrizA, n, n, matrizD); // Encontra grau de cada vértice e coloca na diagonal principal da matrizD

   printf("Matriz D\n\n");

   imprimeMatriz(matrizD, n, n);

   printf("\nMatriz Laplaciana de acordo com a matriz de adjacencia escolhida:\n\n");

   subtraiMatrizes(matrizD, matrizA, matrizL, n);// Faz L = D - A

   imprimeMatriz(matrizL, n, n);

   printf("\nNorma infinito da matriz L:\n");
   printf("%lf\n",normaLinha(matrizL, n));

   printf("\nDeterminante da matriz D:\n");
   printf("%lf\n",determinanteI(matrizD,n));

   getColuna(matrizL, vetorColuna1, COLUNA1, n); // Pegando a primeira coluna da matriz L

   getColuna(matrizL, vetorColuna2, COLUNA3, n); // Pegando a terceira coluna da matriz L

   printf("\nAngulo entre a coluna 1 e coluna 3 da matriz L\n");

   printf("%lf",angulo(vetorColuna1, vetorColuna2, n));

   getchar();

   return 0;
}