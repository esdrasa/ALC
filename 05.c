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

    a = lerMatriz(n);
    b = lerVetor(n);

    u = criaMatriz(n);

    l = criaMatrizI(n);

    x = criaVetor(n);

    imprimeMatriz(a, n);

    imprimeVetor(b, n);

    lu(a, l, u, n);

    b2 = forwardSub(l, b, n);

    x = backSub(u, b2, n);

    imprimeMatriz(u, n);

    imprimeMatriz(l, n);

    imprimeVetor(b2, n);

    imprimeVetor(x, n);
	
	getchar();

    return 0;
}
