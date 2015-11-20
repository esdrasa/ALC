#include <stdio.h>

void imprimeMatriz(double** matriz, int m, int n)
{
    int i, j;

    for(i = 0; i < m; i++)
    {
        for(j = 0; j < n; j++)
            printf("%lf\t", matriz[i][j]);
        printf("\n");
    }

    printf("\n");
}

void imprimeVetor(double* v, int n)
{
    int i;

    for(i = 0; i < n; i++)
        printf("%lf\n", v[i]);
    printf("\n");
}

void imprimePotencia(double *x, double sigma, int n, int k) {
  int i;

  if(k == 0)
  {
    printf("i\t Sigma \t\t\t Vetor q \n");
    printf("-----------------------------------------------------------------\n");
    }
  printf("%d\t ", k + 1);
  printf("%lf\t", sigma);
  for(i = 0; i < n; i++) {
    printf("%lf ", x[i]);
  }
  printf("\n");
}