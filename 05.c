#include <stdio.h>
#include "operacoes.h"
#include "metodos.h"
#include "saida.h"
#include "entrada.h"

/*
Questão 05 não completa
*/
int main()
{
    double **a, *b, **l, **u, *x, *b2;
    int n;

    printf("Digite o tamanho da matriz: ");
    scanf("%d%*c", &n);

    a = lerMatriz(n, n);
    b = lerVetor(n);

    u = criaMatriz(n, n);

    l = criaMatrizI(n);

    imprimeMatriz(a, n, n);

    imprimeVetor(b, n);

    lu(a, l, u, n);
    
    b2 = criaVetor(n);
    forwardSub(l, b, b2, n);
    
    x = criaVetor(n);
    backSub(u, b2, x, n);

    imprimeMatriz(u, n, n);

    imprimeMatriz(l, n, n);

    imprimeVetor(b2, n);

    imprimeVetor(x, n);

    getchar();

    return 0;
}