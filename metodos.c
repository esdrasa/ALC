#include "operacoes.h
#include <math.h>

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

int jacobi(double** A, double* b, double* x, double tolerancia, unsigned long int iMax, int n)
{
    int i, j;
    double erro = tolerancia + 1, *xant, soma;
    unsigned long int k = 0;

    xant = criaVetor(n);

    for(i = 0; i < n; i++)
    {
        /*
        Se algum elemento da diagonal principal for zero, não haverá solução,
        a menos que as linha da matriz sejam permutadas.
        */
        if(A[i][i] == 0)
            return 0;

        //primeiro vetor para solução do método iterativo
        x[i] = b[i] / A[i][i];
    }

    while(k < iMax && erro >= tolerancia)
    {
        if(erro == 0 && tolerancia == 0)
            return 1;

        erro = 0;

        for(i = 0; i < n; i++)
            xant[i] = x[i];

        for(i = 0; i < n; i++)
        {
            soma = 0;

            for(j = 0; j < n; j++)
                if(j != i)
                    soma += A[i][j] * xant[j];

            x[i] = (b[i] - soma) / A[i][i];

            if(erroRelativo(x[i], xant[i]) > erro)
                erro = erroRelativo(x[i], xant[i]);
        }

        k++;
    }

    /*
    Se erro for maior que a tolerancia no final das iterações, significa que
    não ocorreu convergência para o número máximo de iterações definidas.
    */
    if(erro > tolerancia)
        return 0;

    return 1;
}

int gaussSeidel(double** A, double* b, double* x, double tolerancia, unsigned long int iMax, int n)
{
    int i, j;
    double erro = tolerancia + 1, *xant, soma;
    unsigned long int k = 0;

    xant = criaVetor(n);

    for(i = 0; i < n; i++)
    {
        /*
        Se algum elemento da diagonal principal for zero, não haverá solução,
        a menos que as linha da matriz sejam permutadas.
        */
        if(A[i][i] == 0)
            return 0;

        //primeiro vetor para solução do método iterativo
        x[i] = b[i] / A[i][i];
    }

    while(k < iMax && erro >= tolerancia)
    {
        if(erro == 0 && tolerancia == 0)
            return 1;

        erro = 0;

        for(i = 0; i < n; i++)
            xant[i] = x[i];

        for(i = 0; i < n; i++)
        {
            soma = 0;

            for(j = 0; j < n; j++)
                if(j != i)
                    soma += A[i][j] * x[j];

            x[i] = (b[i] - soma) / A[i][i];

            if(erroRelativo(x[i], xant[i]) > erro)
                erro = erroRelativo(x[i], xant[i]);
        }

        k++;
    }

    /*
    Se erro for maior que a tolerancia no final das iterações, significa que
    não ocorreu convergência para o número máximo de iterações definidas.
    */
    if(erro > tolerancia)
        return 0;

    return 1;
}

double criterioLinhas(double** A, int n)
{
    double soma;
    int i, j;

    for(i = 0; i < n; i++)
    {
        soma = 0;
        for(j = 0; j < n; j++)
            if(i != j)
                soma += fabs(A[i][j]);

        if(soma >= fabs(A[i][i]))
            return 0;
    }

    return soma;
}

double criterioColunas(double** A, int n)
{
    double soma;
    int i, j;

    for(j = 0; j < n; j++)
    {
        soma = 0;
        for(i = 0; i < n; i++)
            if(i != j)
                soma += (fabs(A[i][j]) / fabs(A[i][i]));

        if(soma >= 1)
            return 0;
    }

    return soma;
}

double criterioNorma(double** A, int n)
{
    double norma;

    norma = normaF(A,n);

    if(norma >= 1)
        return 0;

    return norma;
}

double criterioSassenfeld(double** A, int n)
{
    double *b, maior;
    int i, j;

    b = criaVetor(n);

    for(i = 0; i < n; i++)
    {
        b[i] = 0;
        for(j = 0; j < n; j++)
        {
            if(j < i)
            {
                b[i] += fabs(A[i][j]) * b[j];
            }
            else if(j > i)
            {
                b[i] += fabs(A[i][j]);
            }
        }
        b[i] /= A[i][i];
    }

    maior = normaInfinito(b, n);

    if(maior >= 1)
        return 0;

    return maior;
}