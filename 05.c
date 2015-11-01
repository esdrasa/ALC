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
    double **a, *b, **l, **u, *x, *x2, *b2, *r1, *r2, tol = 0.000001;
    int n;
    unsigned long int erro = 10000;

    printf("Digite o tamanho da matriz: ");
    scanf("%d%*c", &n);

    a = lerMatriz(n, n);
    b = lerVetor(n);

    u = criaMatriz(n, n);

    l = criaMatrizI(n);

    lu(a, l, u, n);
    
    b2 = criaVetor(n);
    forwardSub(l, b, b2, n);
    
    x = criaVetor(n);
    backSub(u, b2, x, n);
    
    printf("Matriz U:\n");
    imprimeMatriz(u, n, n);
    
    printf("Matriz L:\n");
    imprimeMatriz(l, n, n);
    
    printf("Solucao obtida com o metodo LU:\n");
    imprimeVetor(x, n);
    
    x2 = criaVetor(n);
    while(!gaussSeidel(a, b, x2, tol, erro, n))
    {
	if(tol != 1)
	    tol *= 10;
	else
	    tol += 0.1;
    }
    
    printf("Solucao obtida com o metodo de Gauss-Seidel.\nCom tolerancia %lf e condicao de parada %lu:\n", tol, erro);
    
    imprimeVetor(x2, n);
    
    r1 = criaVetor(n);
    residuo(a, x, b, r1, n);
    
    r2 = criaVetor(n);
    residuo(a, x2, b, r2, n);
    
    printf("Norma do residuo da solucao LU: %lf\n", normaDois(r1, n));
    
    printf("Norma do residuo da solocao de Gauss-Seidel: %lf\n", normaDois(r2, n));

    getchar();

    return 0;
}