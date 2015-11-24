#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "operacoes.h"
#include "metodos.h"
#include "saida.h"
#include "entrada.h"

int main(){
    
    int n, nr;
    double **A, *b, *x, alfa, ls, li;
    
    srand((unsigned) time(NULL));
    
    printf("Digite o tamanho do sistema: ");
    
    scanf("%d%*c", &n);
    
    A = lerMatriz(n, n);
    
    b = lerVetor(n);
    
    x = criaVetor(n);
    
    printf("Digite a constante de variacao do limite superior: ");
    
    scanf("%lf%*c", &alfa);
    
    printf("Digite o limite superior: ");
    
    scanf("%lf%*c", &ls);
    
    printf("Digite o limite inferior: ");
    
    scanf("%lf%*c", &li);
    
    printf("Digite o numero de repeticoes para cada iteracao: ");
    
    scanf("%d%*c", &nr);
    
    simulatedAnnealing(A, b, x, alfa, ls, li, nr, n);
    
    printf("Vetor solucao encontrado:\n");
    
    imprimeVetor(x, n);
    
    getchar();
    
    return 0;
}