#include <stdio.h>
#include "operacoes.h"
#include "metodos.h"
#include "saida.h"
#include "entrada.h"

int main()
{
    double **a, *b, *x, *solucao, *r1, *r2;
    double tol;
    unsigned long int max;
    int n;
    
    printf("Digite o tamanho do sistema: ");
    scanf("%d%*c", &n);
    printf("Digite a tolerancia: ");
    scanf("%lf%*c", &tol);
    printf("Digite a quantida maxima de iteracoes: ");
    scanf("%lu%*c", &max);
    
    a = lerMatriz(n, n);
    b = lerVetor(n);
    solucao = lerVetorSolucao(n);
    
    x = criaVetor(n);
    
    gradienteConjugado(a, b, x, tol, max, n);
    
    printf("Solucao obtida com o metodo do Gradiente Cojudado:\n");
    imprimeVetor(x, n);
    
    r1 = criaVetor(n);
    r2 = criaVetor(n);
    
    residuo(a, x, b, r1, n);
    residuo(a, solucao, b, r2, n);
    
    printf("Norma do residuo com a solucao fornecida: %lf\n", normaDois(r2, n));
    printf("Norma do residuo com a solucao obtida pelo gradiente conjugado: %lf\n", normaDois(r1, n));
    
    getchar();
    
    return 0;
}