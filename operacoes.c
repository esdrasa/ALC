#include <stdlib.h>
#include <math.h>

double** criaMatriz(int m, int n)
{
    double **matriz;
    int i;

    matriz = (double**) malloc(sizeof(double*) * m);

    for(i = 0; i < m; i++)
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

void transposta(double** matriz, double** t, int m, int n)
{
    int i, j;

    for(i = 0; i < m; i++)
        for(j = 0; j < n; j++)
            t[j][i] = matriz[i][j];
}

void backSub(double** matriz, double* b, double* x, int n)
{
    int i, j;

    for(i = n - 1; i >= 0; i--)
    {
        x[i] = b[i];
        for(j = n - 1; j >= i + 1; j--)
            x[i] -= matriz[i][j] * x[j];
        x[i] /= matriz[i][i];
    }
}

void forwardSub(double** matriz, double* b, double* x, int n)
{
    int i, j;

    for(i = 0; i < n; i++)
    {
        x[i] = b[i];
        for(j = 0; j < i; j++)
            x[i] -= matriz[i][j] * x[j];
        x[i] /= matriz[i][i];
    }
}

void multiplica(double** a, int m, int n, double** b, int p, int q, double** result)
{
    int i, j, k;

    for(i = 0; i < m; i++)
        for(j = 0; j < q; j++)
        {
            result[i][j] = 0;
            for(k = 0; k < n; k++)
                result[i][j] += a[i][k] * b[k][j];
        }
}

void multiplicaVetor(double** a, int m, int n, double* v, double* result, int s)
{
    int i, j;

    for(i = 0; i < m; i++)
    {
        result[i] = 0;
        for(j = 0; j < n; j++)
            result[i] += a[i][j] * v[j];
    }

}

double produtoEscalar(double* a, double* b, int n)
{
    int i;
    double produto = 0;

    for(i = 0; i < n; i++)
        produto += a[i] * b[i];

    return produto;
}

void subtraiVetores(double* a, double* b, double* result, int n)
{
    int i;

    for(i = 0; i < n; i++)
        result[i] = a[i] - b[i];
}

void somaVetores(double* a, double* b, double* result, int n)
{
    int i;

    for(i = 0; i < n; i++)
        result[i] = a[i] + b[i];
}

void somaMatrizes(double** a, double** b, double** result, int n)
{
    int i, j;

    for(i = 0; i < n; i++)
        for(j = 0; j < n; j++)
            result[i][j] = a[i][j] + b[i][j];
}

void subtraiMatrizes(double** a, double** b, double** result, int n)
{
    int i, j;

    for(i = 0; i < n; i++)
        for(j = 0; j < n; j++)
            result[i][j] = a[i][j] - b[i][j];
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

void getColuna(double** matriz, double* vColuna, int coluna, int n)
{
    int i;

    for(i = 0; i < n; i++)
        vColuna[i] = matriz[i][coluna];
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

void liberaMatriz(double** a, int nLinhas)
{
    int i;
    
    for(i = 0; i < nLinhas; i++)
	free(a[i]);
    
    free(a);
}

void liberaVetor(double* a)
{
    free(a);
}

void derivada(double* y, double* d, int n)
{
    int i;
    
    for(i = 0; i < n; i++)
	d[i] = y[i] * (n-i);
}

double ordenada(double* polinomio, double abscissa, int n)
{
    int i;
    double y = 0;
    
    for(i = 0; i < n; i++)
	y += polinomio[i] * pow(abscissa, n-i);
    
    y += polinomio[n];
    
    return y;
}

void encontraGraus(double **matrizA, int n, int m, double **matrizD)
{
    int i, j;

    for(i=0; i<n; i++)
    {
        for(j=0; j<m; j++)
        {
            if(matrizA[i][j] > 0)
            {
                matrizD[i][i]++;
            }
        }
    }
}

void **matrizNula(double **matriz, int n, int m)
{
    int i,j;
    for(i=0; i<n; i++)
    {
        for(j=0; j<m; j++)
        {
            matriz[i][j] = 0;
        }
    }
}