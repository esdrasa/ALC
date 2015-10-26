#include <stdlib.h>
#include <math.h>

double** criaMatriz(int m, int n)
{
    double **matriz;
    int i;

    matriz = (double**) malloc(sizeof(double*) * m);

    for(i = 0; i < n; i++)
        matriz[i] = (double*) malloc(sizeof(double) * n);

    return matriz;
}

double** criaMatrizI(int n)
{
    double **matriz;
    int i, j;

    matriz = criaMatriz(n, n);

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

double** transposta(double** matriz, int m, int n)
{
    int i, j;
    double **matrizT;

    matrizT = criaMatriz(n, m);

    for(i = 0; i < m; i++)
        for(j = 0; j < n; j++)
            matrizT[j][i] = matriz[i][j];

    return matrizT;
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

double** multiplica(double** a, int m, int n, double** b, int p, int q)
{
    int i, j, k;
    double **result;

    result = criaMatriz(m, q);

    for(i = 0; i < m; i++)
        for(j = 0; j < q; j++)
        {
            result[i][j] = 0;
            for(k = 0; k < n; k++)
                result[i][j] += a[i][k] * b[k][j];
        }

    return result;
}

double* multiplicaVetor(double** a, int m, int n, double* v, int s)
{
    int i, j;
    double *result;

    result = criaVetor(m);

    for(i = 0; i < m; i++)
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

    result = criaMatriz(n, n);

    for(i = 0; i < n; i++)
        for(j = 0; j < n; j++)
            result[i][j] = a[i][j] + b[i][j];

    return result;
}

double** subtraiMatrizes(double** a, double** b, int n)
{
    int i, j;
    double **result;

    result = criaMatriz(n, n);

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

double normaDois(double* vetor, int n)
{
    int i;
    double norma = 0;

    for(i = 0; i < n; i++)
        norma += vetor[i] * vetor[i];

    norma = sqrt(norma);

    return norma;
}

double angulo(double* v1, double* v2, int n)
{
    double t;

    t = acos(produtoEscalar(v1, v2, n) / (normaDois(v1, n) * normaDois(v2, n)));

    return t;
}

int isVandermonde(double** a, int n)
{
    int i, j;

    for(i = 0; i < n; i++)
        for(j = 0; j < n; j++)
            if(pow(a[i][1], j) != a[i][j])
                return 0;

    return 1;
}

double detVandermonde(double** a, int n)
{
    int i, j;
    double d = 1;

    for(j = 0; j < n-1; j++)
        for(i = j+1; i < n; i++)
            d *= (a[j][1] - a[i][1]);

    return d;
}

double* getColuna(double** matriz, int coluna, int n)
{
    double* vetor;
    int i;

    vetor = criaVetor(n);

    for(i = 0; i < n; i++)
        vetor[i] = matriz[i][coluna];

    return vetor;
}

double normaInfinito(double* vetor, int n)
{
    int i;
    double maior = fabs(vetor[0]), elemento;

    for(i = 1; i < n; i++)
    {
        elemento = fabs(vetor[i]);
        if(elemento > maior)
            maior = elemento;
    }

    return maior;
}

double erroRelativo(double x1, double x0)
{
    double erro;

    erro = fabs(x1 - x0) / fabs(x1);

    return erro;
}