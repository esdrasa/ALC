#include <stdio.h>

void imprimeMatriz(double** matriz, int n)
{
    int i, j;

    for(i = 0; i < n; i++)
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