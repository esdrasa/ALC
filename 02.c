#include <stdio.h>
#include "operacoes.h"
#include "metodos.h"
#include "saida.h"
#include "entrada.h"

int main()
{
    double **A, *b, **R, **Rt, *y, *x1, *x2, tol, w;
    int n, teste, teste2;
    unsigned long int iMax;

    printf("Digite o tamanho da matriz: ");
    scanf("%d%*c", &n);

    A = lerMatriz(n, n);

    b = lerVetor(n);

    R = criaMatriz(n, n);

    teste = cholesky(A, R, n);

    if(!teste)
    {
        printf("Nao foi possivel determinar o fator de Cholesky\n");
    }
    else
    {
        printf("Fator de Cholesky encontrado:\n");
        imprimeMatriz(R, n, n);

        Rt = transposta(R, n, n);

        y = forwardSub(Rt, b, n);

        x1 = backSub(R, y, n);

        printf("Solucao do sistema obtida com o fator de Cholesky:\n");

        imprimeVetor(x1, n);
    }

    printf("Resolucao do sistema com o metodo SOR:\n");

    printf("\nDigite a tolerancia para o erro relativo: ");
    scanf("%lf%*c", &tol);

    printf("Digite um valor para w entre 0 e 2: ");
    scanf("%lf%*c", &w);

    printf("Digite o numero maximo de iteracoes permitidas: ");
    scanf("%lu%*c", &iMax);

    x2 = criaVetor(n);

    if(teste2 = SOR(A, b, x2, tol, w, iMax, n))
    {
        printf("\nSolucao do sistema com o metodo SOR:\n");
        imprimeVetor(x2, n);
    }
    else
    {
        printf("\nNao foi possivel resolver o sistema pelo metodo SOR.\n");
    }

    if(teste && teste2)
    {
        printf("\nDiferenca entre os resultados obtidos com os metodos SOR e Cholesky:\n");

        imprimeVetor(subtraiVetores(x1, x2, n), n);
    }

    getchar();

    return 0;
}