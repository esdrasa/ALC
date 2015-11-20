#include <stdio.h>
#include "operacoes.h"
#include "metodos.h"
#include "saida.h"
#include "entrada.h"

int main() {

  int n, i;
  double **A,
          *x,
          tol = 0.000001,
          autovalor;

  //Lendo os dados necessários

  printf("Digite a dimensão da matriz: \n");
  scanf("%d%*c", &n);
  printf("\n");

  A = lerMatriz(n, n);
  x = lerVetor(n);

  //Chamando o metodo da potencia.
  autovalor = potencia(n, tol, A, x);

  printf("\nO Autovalor Dominante e: %lf\n\n",autovalor );

   //liberando espaço alocado
  liberaVetor(x);
  liberaMatriz(A, n);
  
  return 0;
}