#include <stdio.h>
#include <float.h>
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
    
    //Um polinômio de n graus possui n+1 coeficientes
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
    
    //Calcula p0 no ponto (0, p0) do polinômio dado.
    p0 = ordenada(polinomio, 0, n);
    
    //Calcula p1 no ponto (1, p1) do polinômio dado.
    p1 = ordenada(polinomio, 1, n);
    
    printf("Os seguintes pontos pertencem ao polinomio:\n(%d, %lf)\n(%d, %lf)\n", 0, p0, 1, p1);
    
    getchar();
    
    return 0;
}