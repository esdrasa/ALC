#include <stdio.h>
#include <math.h>
#include "operacoes.h"
#include "metodos.h"
#include "saida.h"
#include "entrada.h"

int main()
{
  double **A, **At, **AAt, **AtA, **X, **Y, **X2, **Y2, **Uo, cofB, cofC, delta, rDelta, r1, r2;
  double **U, **S, **V, **Vt, **Av1, **Av2, *u1, *u2, *v1, *v2, *v3, vs1, vs2, *vetNulo, **UoS, **UoSVt;
  int n, i, j, k;
  printf("Decomposicao SVD para matrizes 2xn\n");
  printf("Digite o numero de colunas: ");
  scanf("%d%*c", &n);

  A = lerMatriz(2, n);
  At = criaMatriz(n, 2);
  AAt = criaMatriz(2, 2);
  AtA = criaMatriz(n, n);
  S = criaMatriz(2, n);
  U = criaMatriz(2, 2);
  Uo = criaMatriz(2, 2);
  UoS = criaMatriz(2, n);
  UoSVt = criaMatriz(2, n);  
  X = criaMatriz(2, 2);
  Y = criaMatriz(2, 2);
  X2 = criaMatriz(n, n);
  Y2 = criaMatriz(n, n);
  V = criaMatriz(n, n);
  Vt = criaMatriz(n, n);
  Av1 = criaMatriz(n, 2);
  Av2 = criaMatriz(n, 2);
  u1 = criaVetor(2);
  u2 = criaVetor(2);
  v1 = criaVetor(n);
  v2 = criaVetor(n);
  v3 = criaVetor(n);
  vetNulo = criaVetor(n);

  matrizNula(S, 2, n);                  // Zera a matriz S diagonal para ela ser completada com os valores singulares
  transposta(A, At, 2, n);             //
  multiplica(A, 2, n, At, n, 2, AAt); // Calculo das matriz a serem encontrados os autovalores
  multiplica(At, n, 2, A, 2, n, AtA);//
  
  cofB = AAt[0][0]*(-1) + AAt[1][1]*(-1);           //
  cofC = AAt[0][0]*AAt[1][1] - AAt[1][0]*AAt[0][1];// Encontrando os coeficientes do polinomio característico
  delta = (cofB*cofB) - 4*cofC;                   //
  rDelta = sqrt(delta);
  
  r1 = (-cofB + rDelta)/2; //Autovalor_1
  r2 = (-cofB - rDelta)/2; //Autovalor_2
  
  S[0][0] = sqrt(r1); 
  S[1][1] = sqrt(r2); 

  vs1 = S[0][0]; //Valores Singular_1
  vs2 = S[1][1]; //Valores Singular_2
  
  matrizNula(X, 2, 2);       // 
  for (i = 0; i < 2; i++)   // copia a matriz AAt para a matriz X
    for (j = 0; j < 2; j++)//
      X[i][j] = AAt[i][j];//

  matrizNula(Y, 2, 2);       //
  for (i = 0; i < 2; i++)   //copia a matriz AAt para a matriz Y
    for (j = 0; j < 2; j++)//
      Y[i][j] = AAt[i][j];//

  X[0][0] = r1 - AAt[0][0]; //(Lambida_1*I - AAt)
  X[1][1] = r1 - AAt[1][1];//

  Y[0][0] = r2 - AAt[0][0];  //(Lambida_2*I - AAt)
  Y[1][1] = r2 - AAt[1][1]; //

  //Calculo de autovetores da matriz U
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
  U[0][0] = u2[0];   //
  U[0][1] = u1[0];  // Formando a matriz U composta pelos autovetores ainda não ortonormais
  U[1][0] = u2[1]; //
  U[1][1] = u1[1];//

  schmidt(U, Uo, 2); // aplica schmidt para ortonormalizar os vetores da matriz U
 
  u1[0] = Uo[0][0];   //
  u2[0] = Uo[0][1];  // 
  u1[1] = Uo[1][0]; // substitiu os vetores
  u2[1] = Uo[1][1];//

  for (i = 0; i < n; i++)
    for (j = 0; j < 2; j++)
    {
      Av1[i][j] = At[i][j]*(1.0/vs1); // multiplica a matriz A por um escalar (1/valor singular_1)
    }

  multiplicaVetor(Av1, n, n, u1, v1, n); //[(1/vs1)*At*u1]

  for (i = 0; i < n; i++)
    for (j = 0; j < 2; j++)
    {
      Av2[i][j] = At[i][j]*(1.0/vs2); // multiplica a matriz A por um escalar (1/valor singular_2)
    }

  multiplicaVetor(Av2, n, 2, u2, v2, n);  //[(1/vs2)*At*u2]
  
  V = criaMatrizI(n);

  for (i = 0; i < n; i++)
    V[i][0] = v1[i];

  for (i = 0; i < n; i++)
    V[i][1] = v2[i];

  vetorNulo(vetNulo, n);

  PSO(AtA, n, vetNulo, v3, 0.000000000001, 800);  // encontra o terceiro alto vetor de AtA que é o vetor v3 ... vn esta dentro do for apenas para variar
  //  retirar PSO do ciclo aumenta a velocidade do algoritmo, mas deixa varias linhas da matriz Vt iguais

  for(i = 0; i < n; i++)
    for(j = 2; j < n; j++)  //
    {                      // Acrencenta ma matrix=z V o terceiro ou enésimo auto vetor de AtA
      V[i][j] = v3[i];                        
    }

  transposta(V, Vt, n, n);

  multiplica(Uo, 2, 2, S, 2, n, UoS);
  multiplica(UoS, 2, n, Vt, 2, n, UoSVt); //calculo da prova real SVD
 
  printf("\n"); 
  printf("Matriz ortogonal U:\n");
  imprimeMatriz(Uo, 2, 2);

  printf("Matriz diagonal S:\n");
  imprimeMatriz(S, 2, n);

  printf("Matriz ortogonal V:\n");
  imprimeMatriz(Vt, n, n);

  printf("Matriz SVD:\n");
  imprimeMatriz(UoSVt, 2, n);

  getchar();
  return 0;

}