#include <stdio.h>
#include "operacoes.h"
#include "metodos.h"
#include "saida.h"
#include "entrada.h"

int main() {

  int n, i;
  double **A,
          *x,
          *q,
          tol = 0.000001,
          autovalor;

  //Lendo os dados necessários

  printf("Digite a dimensão da matriz\n");
  scanf("%d%*c", &n);
  A = lerMatriz(n, n);
  imprimeMatriz(A, n, n);

  printf("Digite o Vetor x\n");
  x = lerVetor(n);
  imprimeVetor(x, n);

  //Chamando o metodo da potencia.
  autovalor = potencia(n, tol, A, x);

  printf("\nO Autovalor Dominante e: %lf\n",autovalor );

  getchar();
  return 0;
}
