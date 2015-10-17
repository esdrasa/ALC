#include "operacoes.h"

double* residuo(double** matriz, double* x, double* b, int n)
{
    return subtraiVetores(b, multiplicaVetor(matriz, x, n), n);
}

int lu(double** a, double** l, double** u, int n)
{
    int k, i, j;

    /* Os dois loops abaixo não são necessários caso a matriz "A" seja utilizada para calcular a matriz "U".
    A não utilização deles pode vir a ser útil no tempo de processamento, já que estes loops apenas atribuem os
    valores de A à U.*/
    for(i = 0; i < n; i++)
        for(j = 0; j < n; j++)
            u[i][j] = a[i][j];

    for(k = 1; k < n; k++)
    {
        if(u[k-1][k-1] == 0)
            return 0;
        for(i = k; i < n; i++)
        {
            l[i][k-1] = u[i][k-1] / u[k-1][k-1];

            for(j = k-1; j < n; j++)
            {
                u[i][j] = u[i][j] - l[i][k-1]*u[k-1][j];
            }
        }
    }

    if(u[n-1][n-1] == 0)
        return 0;

    return 1;
}
