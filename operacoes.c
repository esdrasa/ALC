#include <stdlib.h>
#include <math.h>

double** criaMatriz(int n)
{
    double **matriz;
    int i;

    matriz = (double**) malloc(sizeof(double*) * n);

    for(i = 0; i < n; i++)
        matriz[i] = (double*) malloc(sizeof(double) * n);

    return matriz;
}

double** criaMatrizI(int n)
{
    double **matriz;
    int i, j;

    matriz = criaMatriz(n);

    for(i = 0; i < n; i++)
        for(j = 0; j < n; j++)
        {
            if(i == j)
            {
                matriz[i][j] = 1;
            }
            else
            {
                    matriz[i][j] = 0;
            }
        }

    return matriz;
}

double* criaVetor(int n)
{
    double *v;

    v = (double*) malloc(sizeof(double) * n);

    return v;
}

void transpor(double** matriz, int n)
{
    int i, j;
    double aux;

    for(i = 0; i < n; i++)
        for(j = 0; j < n; j++)
        {
            aux = matriz[i][j];
            matriz[i][j] = matriz[j][i];
            matriz[j][i] = aux;
        }
}

double* backSub(double** matriz, double* b, int n)
{
    double *x;
    int i, j;

    x = criaVetor(n);

    for(i = n - 1; i >= 0; i--)
    {
        x[i] = b[i];
        for(j = n - 1; j >= i + 1; j--)
            x[i] -= matriz[i][j] * x[j];
        x[i] /= matriz[i][i];
    }

    return x;
}

double* forwardSub(double** matriz, double* b, int n)
{
    double *x;
    int i, j;

    x = criaVetor(n);

    for(i = 0; i < n; i++)
    {
        x[i] = b[i];
        for(j = 0; j < i; j++)
            x[i] -= matriz[i][j] * x[j];
        x[i] /= matriz[i][i];
    }

    return x;
}

double** multiplica(double** a, double** b, int n)
{
    int i, j, k;
    double **result;

    result = criaMatriz(n);

    for(i = 0; i < n; i++)
        for(j = 0; j < n; j++)
        {
            result[i][j] = 0;
            for(k = 0; k < n; k++)
                result[i][j] += a[i][k] * b[k][j];
        }

    return result;
}

double* multiplicaVetor(double** a, double* v, int n)
{
    int i, j;
    double *result;

    result = criaVetor(n);

    for(i = 0; i < n; i++)
    {
        result[i] = 0;
        for(j = 0; j < n; j++)
            result[i] += a[i][j] * v[j];
    }

    return result;
}

double produtoEscalar(double* a, double* b, int n)
{
    int i;
    double produto = 0;

    for(i = 0; i < n; i++)
        produto += a[i] * b[i];

    return produto;
}

double* subtraiVetores(double* a, double* b, int n)
{
    int i;
    double *result;

    result = criaVetor(n);

    for(i = 0; i < n; i++)
        result[i] = a[i] - b[i];

    return result;
}

double* somaVetores(double* a, double* b, int n)
{
    int i;
    double *result;

    result = criaVetor(n);

    for(i = 0; i < n; i++)
        result[i] = a[i] + b[i];

    return result;
}

double** somaMatrizes(double** a, double** b, int n)
{
    int i, j;
    double **result;

    result = criaMatriz(n);

    for(i = 0; i < n; i++)
        for(j = 0; j < n; j++)
            result[i][j] = a[i][j] + b[i][j];

    return result;
}

double** subtraiMatrizes(double** a, double** b, int n)
{
    int i, j;
    double **result;

    result = criaMatriz(n);

    for(i = 0; i < n; i++)
        for(j = 0; j < n; j++)
            result[i][j] = a[i][j] - b[i][j];

    return result;
}

double normaF(double** matriz, int n)
{
    int i, j;
    double norma = 0;

    for(i = 0; i < n; i++)
        for(j = 0; j < n; j++)
            norma += matriz[i][j] * matriz[i][j];

    norma = sqrt(norma);

    return norma;
}

double normaLinha(double** matriz, int n)
{
    int i, j;
    double norma = -1, soma = 0;

    for(i = 0; i < n; i++)
    {
        for(j = 0; j < n; j++)
            soma += fabs(matriz[i][j]);

        if(soma > norma)
            norma = soma;

        soma = 0;
    }

    return norma;
}

double normaColuna(double** matriz, int n)
{
    int i, j;
    double norma = -1, soma = 0;

    for(i = 0; i < n; i++)
    {
        for(j = 0; j < n; j++)
            soma += fabs(matriz[j][i]);

        if(soma > norma)
            norma = soma;

        soma = 0;
    }

    return norma;
}