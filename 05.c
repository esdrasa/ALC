#include <stdio.h>
#include "operacoes.h"
#include "metodos.h"
#include "saida.h"
#include "entrada.h"

int main()
{
    double **a, *b, **l, **u, *x, *x2, *b2, *r1, *r2, tol = 0.0001;
    int n;
    unsigned long int erro = 1000000;

    printf("Digite o tamanho da matriz: ");
    scanf("%d%*c", &n);

    a = lerMatriz(n, n);
    b = lerVetor(n);

    u = criaMatriz(n, n);

    l = criaMatrizI(n);
    
    /*
     * Dado um sitema A*x = b, teremos L*U*x = b
     * 
     * U*x = b2
     * 
     * L*b2 = b (b2 é resolvido com substituição para frente)
     * 
     * Após descobrir b2, se resolverá x:
     * U*x = b2 (aplica-se substituição para trás e x será resolvido)
    */
    
    //O método abaixo calcula as matrizes "l" e "u" à partir da matriz "a"
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
    
    liberaMatriz(a, n);
    liberaMatriz(l, n);
    liberaMatriz(u, n);
    liberaVetor(b);
    liberaVetor(x);
    liberaVetor(x2);
    liberaVetor(b2);
    liberaVetor(r1);
    liberaVetor(r2);

    return 0;
}