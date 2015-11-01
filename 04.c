#include <stdio.h>
#include "operacoes.h"
#include "metodos.h"
#include "saida.h"
#include "entrada.h"

int main()
{
    double *polinomio, *raiz, tol, p0, p1;
    int n, i, r = 0;
    unsigned long int max;
    
    printf("Digite o grau do polinomio: ");
    scanf("%d%*c", &n);
    
    polinomio = criaVetor(n + 1);
    
    printf("Digite os coeficientes do polinomio, em ordem: ");
    for(i = 0; i <= n; i++)
	scanf("%lf%*c", &polinomio[i]);
    
    printf("Digite uma tolerancia entre a raiz real e a raiz a ser calculada: ");
    scanf("%lf%*c", &tol);
    
    printf("Digite o numero maximo de iteracoes: ");
    scanf("%lu%*c", &max);
    
    raiz = criaVetor(1);
    
    if(newton(polinomio, raiz, tol, max, n))
    {
	printf("Existe pelo menos uma raiz no polinomio digitado:\n");
	printf("Raiz = %lf\n", *raiz);
    }
    else
    {
	printf("Nao foram encontradas raizes neste polinomio.\n");
    }
    
    p0 = ordenada(polinomio, 0, n);
    p1 = ordenada(polinomio, 1, n);
    
    printf("Os seguintes pontos pertencem ao polinomio:\n(%d, %lf)\n(%d, %lf)\n", 0, p0, 1, p1);
    
    getchar();
    
    return 0;
}