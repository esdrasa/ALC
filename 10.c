#include <stdio.h>
#include <math.h>
#include "operacoes.h"
#include "metodos.h"
#include "saida.h"
#include "entrada.h"

int main()
{
  double **A, **At, **AAt, **X, **Y, cofB, cofC, delta, rDelta, r1, r2, **U, **S, **V, **Vt, *u1, *u2, *v1, *v2, *v3;
  int n, i, j, k;
  printf("Decomposicao SVD para matrizes 2xn\n");
  printf("Digite o numero de colunas: ");
  scanf("%d%*c", &n);

  A = lerMatriz(2, n);
  At = criaMatriz(n, 2);
  AAt = criaMatriz(2, 2);
  S = criaMatriz(2, n);
  U = criaMatriz(2, 2);
  X = criaMatriz(2, 2);
  Y = criaMatriz(2, 2);
  V = criaMatriz(n, n);
  Vt = criaMatriz(n, n);
  u1 = criaVetor(2);
  u2 = criaVetor(2);
  v1 = criaVetor(n);
  v2 = criaVetor(n);
  v3 = criaVetor(n);
  
  matrizNula(S, 2, n);
  transposta(A, At, 2, n);
  multiplica(A, 2, n, At, n, 2, AAt);
  
  cofB = AAt[0][0]*(-1) + AAt[1][1]*(-1); 
  cofC = AAt[0][0]*AAt[1][1] - AAt[1][0]*AAt[0][1];
  delta = (cofB*cofB) - 4*cofC;
  rDelta = sqrt(delta);
  r1 = (-cofB + rDelta)/2; //Autovalor_1
  r2 = (-cofB - rDelta)/2; //Autovalor_2
  S[0][0] = sqrt(r1); //Valor Singular_1
  S[1][1] = sqrt(r2); //Valores Singular_2
  
  matrizNula(X, 2, 2); // 
  for (i = 0; i < 2; i++)// copia a matriz AAt para a matriz X
    for (j = 0; j < 2; j++)//
      X[i][j] = AAt[i][j];//

  matrizNula(Y, 2, 2);//
  for (i = 0; i < 2; i++)//copia a matriz AAt para a matriz Y
    for (j = 0; j < 2; j++)//
      Y[i][j] = AAt[i][j];//

  X[0][0] = r1 - AAt[0][0]; //(Lambida_1*I - AAt)
  X[1][1] = r1 - AAt[1][1]; //

  Y[0][0] = r2 - AAt[0][0]; //(Lambida_2*I - AAt)
  Y[1][1] = r2 - AAt[1][1]; //

  //Calculo de autovetores
  if (fabs(X[0][0]) < fabs(X[0][1]))
  {
     u1[0] = X[0][0]/X[0][0]*(-1);
     u1[1] = X[0][1]/X[0][0];
  }
  else
  if (fabs(X[0][0]) > fabs(X[0][1]))
  {
    u1[0] = X[0][0]/X[0][1];
    u1[1] = X[0][1]/X[0][1]*(-1);
  }
  else
  {
    u1[0] = X[0][0];
    u1[1] = X[0][1]*(-1);
  }

   if (fabs(Y[0][0]) < fabs(Y[0][1]))
  {
     u2[0] = Y[0][0]/Y[0][0]*(-1);
     u2[1] = Y[0][1]/Y[0][0];
  }
  else
  if (fabs(Y[0][0]) > fabs(Y[0][1]))
  {
    u2[0] = Y[0][0]/Y[0][1];
    u2[1] = Y[0][1]/Y[0][1]*(-1);
  }
  else
  { 
    u2[0] = Y[0][0];
    u2[1] = Y[0][1]*(-1);
  }

  if (u1[0] < 0 && u1[1] < 0)
  {
    u1[0] = u1[0]*(-1);
    u1[1] = u1[1]*(-1);
  }

  if (u2[0] < 0 && u2[1] < 0)
  {
    u2[0] = u2[0]*(-1);
    u2[1] = u2[1]*(-1);
  }
  matrizNula(U, 2, 2);
  U[0][0] = u2[0];
  U[0][1] = u1[0];
  U[1][0] = u2[1];
  U[1][1] = u1[1];

  matrizNula(V, 2, 2);

  printf("Matriz X:\n");
  imprimeMatriz(X, 2, 2);

  printf("Matriz Y:\n");
  imprimeMatriz(Y, 2, 2);

  printf("vetor X:\n");
  imprimeVetor(u1, 2);

  printf("vetor Y:\n");
  imprimeVetor(u2, 2);

  printf("Matriz ortogonal U:\n");
  imprimeMatriz(U, 2, 2);

  printf("Matriz diagonal S:\n");
  imprimeMatriz(S, 2, n);
// falta a ortogonalização

  printf("Matriz ortogonal V:\n");
  imprimeMatriz(V, n, n);


  getchar();
  return 0;




}